
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

include_guard()

function( ncrystal_ensure_valid_pkgname pkgname )
  string(MAKE_C_IDENTIFIER "${pkgname}" "tmp" )
  string(TOLOWER "${tmp}" "tmp" )
  string(FIND "${tmp}" "__" "badidx" )
  if ( NOT "${pkgname}" STREQUAL "${tmp}" OR NOT "${badidx}" STREQUAL "-1" )
    message( FATAL_ERROR "Not a valid name of NCrystal component:"
      " \"${pkgname}\" (must be lowercase, start with a-z and only contain"
      " \"a-z0-9_\" and no \"__\"")
  endif()
endfunction()

function( ncrystal_parse_dep_txt filename out_deplist )
  set( deplist "" )
  file(STRINGS "${filename}" "tmp" )
  foreach( "line" "${tmp}" )
    string(REGEX REPLACE "[ \t\r\n]" ";" "parts" "${line}" )
    foreach( "part" ${parts} )
      if ( part AND NOT "${part}" IN_LIST deplist )
        list(APPEND deplist "${part}" )
      endif()
    endforeach()
  endforeach()
  list(SORT deplist)
  set( "${out_deplist}" "${deplist}" PARENT_SCOPE )
endfunction()

function( ncrystal_detectsrc )
  #Code is in "packages", with:
  #  Headers: ncrystal_core/include/NCrystal/[internal/]<pgkname>
  #  Sources: ncrystal_core/src/<pkgname>
  #  Inter-package dependencies: ncrystal_core/src/<pkgname>/dep.txt

  set( "incdir" "${PROJECT_SOURCE_DIR}/include")

  #First look for the src directories:
  set( "pkglist" "" )
  file(
    GLOB "srcdir_candidates" LIST_DIRECTORIES true
    CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*"
  )
  foreach( "srcdir" ${srcdir_candidates} )
    if ( IS_DIRECTORY "${srcdir}" AND EXISTS "${srcdir}/dep.txt" )
      get_filename_component( "pkgname" "${srcdir}" NAME )
      ncrystal_ensure_valid_pkgname( "${pkgname}" )
      list( APPEND pkglist "${pkgname}" )
      ncrystal_parse_dep_txt( "${srcdir}/dep.txt" deplist )
      set( "ncpkg__${pkgname}__deplist_direct" "${deplist}" )
      file(
        GLOB "ncpkg__${pkgname}__srcfiles" LIST_DIRECTORIES false
        CONFIGURE_DEPENDS "${srcdir}/*.c"  "${srcdir}/*.cc"
      )

      file(
        GLOB "ncpkg__${pkgname}__privhdrfiles" LIST_DIRECTORIES false
        CONFIGURE_DEPENDS "${srcdir}/*.h"  "${srcdir}/*.hh"
      )

      set( tmp_api_type "absent" )
      set( pkg_incdir "" )

      if ( IS_DIRECTORY "${incdir}/NCrystal/${pkgname}")
        set( tmp_api_type "public" )
        set( pkg_incdir "${incdir}/NCrystal/${pkgname}" )
      endif()
      if ( IS_DIRECTORY "${incdir}/NCrystal/internal/${pkgname}")
        if ( NOT "${tmp_api_type}" STREQUAL "absent" )
          message(FATAL_ERROR "NCrystal component \"${pkgname}\" has"
            " both public and internal headers")
        endif()
        set( tmp_api_type "internal" )
        set( pkg_incdir "${incdir}/NCrystal/internal/${pkgname}" )
      endif()
      if(  "${tmp_api_type}" STREQUAL "absent" )
        set( "ncpkg__${pkgname}__headerfiles" "" )
        set( tmp_api_type "internal" )
      else()
        file(
          GLOB "tmp_hdrs" LIST_DIRECTORIES false
          CONFIGURE_DEPENDS "${pkg_incdir}/*.h" "${pkg_incdir}/*.hh"
        )
        if ( NOT tmp_hdrs )
          message(FATAL_ERROR "No header files found in ${pkg_incdir}")
        endif()
        set( "ncpkg__${pkgname}__headerfiles" "${tmp_hdrs}" )

      endif()
      set( "ncpkg__${pkgname}__apitype" "${tmp_api_type}" )
    endif()
  endforeach()

  #Now expand the dependencies transitively. We do do this in a big loop rather
  #than through recursive function calls, to avoid scope issues (PARENT_SCOPE in
  #a nested function call is not the scope of the current top-level function):
  set( pkg_pending "${pkglist}" )
  set( loopcount 1 )
  while(pkg_pending)
    #Simplistic guard against circular dependencies:
    math(EXPR loopcount "${loopcount}+1" )
    if ( loopcount GREATER "10000" )
      message(FATAL_ERROR "NCrystal components have circular dependencies!" )
    endif()

    #Check first item in pkg_pending:
    list(GET pkg_pending 0 pkgname)

    #Check whether all direct deps have already been processed:
    set( pkg_deplist "${ncpkg__${pkgname}__deplist_direct}" )
    if ( "${pkgname}" IN_LIST ncpkg__${pkgname}__deplist_direct )
      message(FATAL_ERROR "NCrystal component \"${pkgname}\" has a"
        " direct dependency on itself!")
    endif()

    foreach( dep ${ncpkg__${pkgname}__deplist_direct} )
      if ( "${dep}" IN_LIST pkg_pending )
        if ( "${pkgname}" IN_LIST ncpkg__${dep}__deplist_direct )
          message(FATAL_ERROR "NCrystal components depend on each other:"
            " \"${pkgname}\" and \"${dep}\".")
        endif()
        #Oups, we are not ready. Let us wait with ${pkgname} and move it to the
        #end of pkg_pending:
        set( pkg_deplist "<<notready>>" )
        list(REMOVE_ITEM pkg_pending "${pkgname}")
        list(APPEND pkg_pending "${pkgname}")
        break()#break out of foreach loop
      endif()
    endforeach()
    if ( NOT "${pkg_deplist}" STREQUAL "<<notready>>" )
      foreach( dep ${ncpkg__${pkgname}__deplist_direct} )
        foreach( dep2 ${ncpkg__${dep}__deplist} )
          if ( NOT "${dep2}" IN_LIST pkg_deplist )
            list( APPEND pkg_deplist "${dep2}" )
          endif()
        endforeach()
      endforeach()
      list(REMOVE_ITEM pkg_pending "${pkgname}")
      list( SORT pkg_deplist )
      set( ncpkg__${pkgname}__deplist "${pkg_deplist}" )
    endif()
  endwhile()

  #Final results:
  if ( OFF )
    foreach( "pkgname" ${pkglist} )
      message(STATUS "COMPONENT: ${pkgname}")
      message(STATUS "   abitype : ${ncpkg__${pkgname}__apitype}" )
      message(STATUS "   directdeplist : ${ncpkg__${pkgname}__deplist_direct}" )
      message(STATUS "   deplist : ${ncpkg__${pkgname}__deplist}" )
      message(STATUS "   srcfiles: ${ncpkg__${pkgname}__srcfiles}" )
      message(STATUS "   privhdrfiles: ${ncpkg__${pkgname}__privhdrfiles}" )
      message(STATUS "   headers : ${ncpkg__${pkgname}__headerfiles}" )
    endforeach()
  endif()

  set( "ncpkg_names" "${pkglist}" PARENT_SCOPE )
  foreach( "pkgname" ${pkglist} )
    foreach( varname
        "apitype" "deplist_direct" "deplist" "srcfiles" "privhdrfiles" "headerfiles" )
      set( tmp "ncpkg__${pkgname}__${varname}" )
      set( "${tmp}" "${${tmp}}" PARENT_SCOPE )
    endforeach()
  endforeach()

endfunction()

function( ncrystal_impl_setup_properties_for_pkg pkgname )
  if ( ncpkg__${pkgname}__headerfiles )
    foreach ( hf ${ncpkg__${pkgname}__headerfiles} )
      file(RELATIVE_PATH "hfrp" "${PROJECT_SOURCE_DIR}/include" "${hf}" )
      get_filename_component( hfrpdir "${hfrp}" DIRECTORY )
      ncinstall(FILES "${hf}" DESTINATION "${NCrystal_INCDIR}/${hfrpdir}" )
    endforeach()
  endif()

  set( "hdrfiles"
    ${ncpkg__${pkgname}__privhdrfiles}
    ${ncpkg__${pkgname}__headerfiles}
  )
  foreach( dep ${ncpkg__${pkgname}__deplist} )
    list( APPEND hdrfiles ${ncpkg__${dep}__headerfiles} )
  endforeach()
  set( "incdirs" "${ncapi_include_path};${PROJECT_SOURCE_DIR}/include" )
  if ( ncpkg__${pkgname}__privhdrfiles )
    list( APPEND incdirs "${PROJECT_SOURCE_DIR}/src/${pkgname}" )
  endif()

  set_source_files_properties(
    ${ncpkg__${pkgname}__srcfiles}
    PROPERTIES
    INCLUDE_DIRECTORIES "${incdirs}"
    LANGUAGE "CXX"
    OBJECT_DEPENDS "${hdrfiles}"
  )
endfunction()

function( ncrystal_process_core_srcfiles )
  set( "allsrcfiles" )
  foreach( pkgname ${ncpkg_names} )
    ncrystal_impl_setup_properties_for_pkg( ${pkgname} )
    list(APPEND allsrcfiles ${ncpkg__${pkgname}__srcfiles} )
  endforeach()
  set( "ncpkg_all_ncrystal_lib_srcfiles" "${allsrcfiles}" PARENT_SCOPE )
endfunction()
