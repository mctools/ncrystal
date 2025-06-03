
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

# Utilities for building up a testing infrastructure based on C++/C apps and
# libraries as well as python scripts. The API functions from this file are:
#
#   mctools_testutils_add_tests_pyscripts
#   mctools_testutils_add_tests_apps
#   mctools_testutils_add_test_libs (utility libs for the apps)
#   mctools_testutils_add_test_modules (shared libs for python/ctypes access)
#
# This file (and the associated mctools_testlauncher.py file) is intended to be
# eventually shared between at least NCrystal and MCPL projects.

include_guard()

include( mctools_utils )

set(
  "mctools_launcher_file"
  "${CMAKE_CURRENT_LIST_DIR}/mctools_testlauncher.py"
)
if ( NOT EXISTS "${mctools_launcher_file}" )
  message( FATAL_ERROR "Did not find expected ${mctools_launcher_file}" )
endif()

if ( DEFINED MCTOOLS_REQUIRE_ALL_TEST_DEPS )
  if ( MCTOOLS_REQUIRE_ALL_TEST_DEPS )
    set( MCTOOLS_REQUIRE_ALL_TEST_DEPS "ON" )
  else()
    set( MCTOOLS_REQUIRE_ALL_TEST_DEPS "OFF" )
  endif()
endif()

if ( NOT DEFINED MCTOOLS_TESTUTILS_PYTHON_EXECUTABLE
     AND NOT DEFINED Python3_EXECUTABLE )
  set( MCTOOLS_TESTUTILS_PYTHON_EXECUTABLE "auto" )
endif()

function( mctools_testutils_add_tests_pyscripts scriptsdir envmod ignorepattern )
  file(
    GLOB pyscriptlist
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS "${scriptsdir}/*.py"
  )
  list(
    APPEND envmod
    "PYTHONPYCACHEPREFIX=set:${PROJECT_BINARY_DIR}/tests/test_pycache"
    "MCTOOLS_TESTMODULES_LOCDIR=set:${PROJECT_BINARY_DIR}/tests/modlib_locations_$<CONFIG>"
  )
  set( pyexec "" )
  foreach(pyscript ${pyscriptlist})
    if ( NOT pyexec )
      mctools_testutils_internal_getpyexec( "pyexec" )
    endif()
    get_filename_component(bn "${pyscript}" NAME_WE)
    get_filename_component(psdir "${pyscript}" DIRECTORY)

    set(missingdeps "")
    mctools_testutils_internal_detectpydeps( "deplist" "${pyscript}" )
    if ( deplist )
      mctools_testutils_internal_missingpydeps( "missingdeps" "${deplist}" )
    endif()

    set( reflog "${psdir}/${bn}.log" )
    if ( EXISTS "${reflog}" )
      set( testname "py_rl_${bn}" )
    else()
      set( reflog "" )
      set( testname "py_${bn}" )
    endif()

    if ( ignorepattern AND "${bn}" MATCHES "${ignorepattern}" )
      message( STATUS "Skipping test ${testname}" )
      continue()
    endif()

    if ( missingdeps )
      string(REPLACE ";" " " md_pretty "${missingdeps}")
      if ( MCTOOLS_REQUIRE_ALL_TEST_DEPS )
        message( FATAL_ERROR "Test ${testname} has missing deps: ${md_pretty}")
      endif()
      message( WARNING "Skipping py test ${testname} due to missing deps: ${md_pretty}")
      continue()
    endif()

    mctools_testutils_internal_addtest( "${testname}" "${pyscript}" "${reflog}" )
    if ( envmod )
      set_property(
        TEST "${testname}"
        PROPERTY ENVIRONMENT_MODIFICATION "${envmod}"
      )
    endif()
  endforeach()
endfunction()

#NB: a project with both mcpl and ncrystal subdirs would perhaps get a name
#clash here? (of course in this particular case, that would be a feature)
set_property(GLOBAL PROPERTY mctools_testutils_internal_pydepspresent "")
set_property(GLOBAL PROPERTY mctools_testutils_internal_pydepsabsent "")

function( mctools_testutils_internal_haspydep resvar pydep )
  get_property(
    known_present
    GLOBAL PROPERTY mctools_testutils_internal_pydepspresent
  )
  if( "${pydep}" IN_LIST known_present )
    set( "${resvar}" "ON" PARENT_SCOPE )
    return()
  endif()
  get_property(
    known_absent
    GLOBAL PROPERTY mctools_testutils_internal_pydepsabsent
  )
  if( "${pydep}" IN_LIST known_absent )
    set( "${resvar}" "OFF" PARENT_SCOPE )
    return()
  endif()

  #Ok, first time encountering this particular dependency. We will detect its
  #presence by importing it in Python. As a special case, for matplotlib we
  #import matplotlib.pyplot and not just matplotlib. This has the advantage in
  #that it ensures the matplotlib font cache building is triggered, so we don't
  #get unexpected stderr output in the middle of our tests later.

  mctools_testutils_internal_getpyexec( "pyexec" )

  set( cmdtocheckflags "" )
  set( cmdtocheck "" )

  set( skip_status_msg OFF )

  if ( "x${pydep}" STREQUAL "xtoml" )
    #"tomllib" in python 3.11 and later (where it is ALWAYS present), otherwise
    #require "tomli".
    execute_process(
      COMMAND "${pyexec}" "-c"
      "import sys; print('tomllib' if sys.version_info[0:2]>=(3,11) else 'tomli')"
      RESULT_VARIABLE "res" OUTPUT_VARIABLE "pymodtoimport"
      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
    )
    if( NOT "x${res}" STREQUAL "x0" OR "x${pymodtoimport}" STREQUAL "x" )
      message(
        FATAL_ERROR "Problems invoking Python to"
        " determine tomllib/tomli name"
      )
    endif()
  elseif ( "x${pydep}" STREQUAL "xmatplotlib" )
    #One advantage of importing matplotlib.pyplot and not just matplotlib, is
    #that it triggers font cache building, so it won't clobber test output
    #later.
    set( pymodtoimport "matplotlib.pyplot" )
  elseif ( "x${pydep}" STREQUAL "xpyyaml" )
    set( pymodtoimport "yaml" )
  elseif ( "x${pydep}" STREQUAL "xcmake" )
    #cmake is of course present since we are running cmake...
    set( pymodtoimport "" )
    set( cmdtocheck "" )
  elseif ( "x${pydep}" STREQUAL "xruff" )
    set( pymodtoimport "" )
    set( cmdtocheck "ruff" )
    set( cmdtocheckflags "--version" )
  else()
    set( pymodtoimport "${pydep}" )
  endif()
  if ( pymodtoimport )
    if ( "${pymodtoimport}" STREQUAL "${pydep}" )
      set( "tmp" "${pydep}" )
    else()
      set( "tmp" "${pymodtoimport} (for ${pydep})" )
    endif()
    message(STATUS "Testing for presence of python module: ${tmp}")
    execute_process(
      COMMAND "${pyexec}" "-c" "import ${pymodtoimport}"
      RESULT_VARIABLE "res" OUTPUT_QUIET ERROR_QUIET
    )
  elseif( cmdtocheck )
    if ( "${cmdtocheck}" STREQUAL "${pydep}" )
      set( "tmp" "${pydep}" )
    else()
      set( "tmp" "${cmdtocheck} (for ${pydep})" )
    endif()
    message(STATUS "Testing for presence of command: ${tmp})")
    execute_process(
      COMMAND "${cmdtocheck}" ${cmdtocheckflags}
      RESULT_VARIABLE "res" OUTPUT_QUIET ERROR_QUIET
    )
  else()
    #Nothing to check (e.g. cmake):
    set( res 0 )
    set( skip_status_msg ON )
  endif()
  if( "x${res}" STREQUAL "x0" )
    set_property(
      GLOBAL PROPERTY mctools_testutils_internal_pydepspresent
      "${known_present};${pydep}"
    )
    set( "${resvar}" "ON" PARENT_SCOPE )
    if ( NOT skip_status_msg )
      message(STATUS "... PRESENT.")
    endif()
  else()
    set_property(
      GLOBAL PROPERTY mctools_testutils_internal_pydepsabsent
      "${known_absent};${pydep}"
    )
    set( "${resvar}" "OFF" PARENT_SCOPE )
    if ( NOT skip_status_msg )
      message(STATUS "... ABSENT.")
    endif()
  endif()
endfunction()

function( mctools_testutils_internal_missingpydeps resvar pydeps )
  set( missing "")
  foreach( pydep ${pydeps} )
    mctools_testutils_internal_haspydep( present "${pydep}" )
    if ( NOT present )
      list( APPEND missing "${pydep}" )
    endif()
  endforeach()
  set( "${resvar}" "${missing}" PARENT_SCOPE )
endfunction()

function(
    mctools_testutils_add_test_modules
    librootdir extra_link_libs extra_inc_dirs
  )
  mctools_testutils_internal_add_test_libs(
    "${librootdir}" "${extra_link_libs}" "${extra_inc_dirs}" "ON"
  )
endfunction()

function(
    mctools_testutils_add_test_libs
    librootdir extra_link_libs extra_inc_dirs
  )
  mctools_testutils_internal_add_test_libs(
    "${librootdir}" "${extra_link_libs}" "${extra_inc_dirs}" "OFF"
  )
endfunction()

function(
    mctools_testutils_internal_add_test_libs
    librootdir extra_link_libs extra_inc_dirs is_module
  )
  file(
    GLOB libdirs
    LIST_DIRECTORIES true
    CONFIGURE_DEPENDS "${librootdir}/lib_*"
  )
  foreach(libdir ${libdirs})
    get_filename_component(bn "${libdir}" NAME)
    string(SUBSTRING "${bn}" 4 -1 "bn")

    if ( is_module )
      set( libprefix "TestMod_" )
    else()
      set( libprefix "TestLib_" )
    endif()
    set( name "${libprefix}${bn}" )
    mctools_testutils_internal_getsrcfiles( srcfiles is_cxx "${libdir}" is_module )
    mctools_testutils_internal_detectlibdeps( "deplist" "${srcfiles}" )
    foreach( dep ${deplist} )
      if ( NOT "x${dep}" STREQUAL "x${bn}" )
        #Private because not headers (TODO: Add headers as explicit deps?):
        target_link_libraries( ${bn} PRIVATE "TestLib_${dep}" )
      endif()
    endforeach()
    if ( is_module )
      #Shared library to be loaded by python ctypes
      add_library( ${name} MODULE ${srcfiles} )
      if ( is_cxx )
        target_compile_features( ${name} PRIVATE cxx_std_11 )
      endif()
      set_target_properties(
        ${name} PROPERTIES
        CXX_VISIBILITY_PRESET "hidden"
        C_VISIBILITY_PRESET "hidden"
        VISIBILITY_INLINES_HIDDEN "ON"
      )
      file (GENERATE
        OUTPUT "${PROJECT_BINARY_DIR}/tests/modlib_locations_$<CONFIG>/module_loc_${name}.txt"
        CONTENT "$<TARGET_FILE:${name}>"
        TARGET ${name}
      )
    else()
      #Library to be linked by test applications:
      add_library( ${name} ${srcfiles} )
    endif()

    target_link_libraries( ${name} PRIVATE ${extra_link_libs} )
    target_include_directories( ${name} ${extra_inc_dirs} )
    mctools_apply_strict_comp_properties( ${name} )

    if ( EXISTS "${libdir}/include" )
      if ( is_module )
        message( FATAL_ERROR
          "Modules can not export headers and should therefore not"
          " contain include dirs. Offending directory is: ${libdir}/include")
      else()
        if ( NOT EXISTS "${libdir}/include/TestLib_${bn}" )
          message(
            FATAL_ERROR
            "${libdir}/include does not include TestLib_${bn} subdir"
          )
        endif()
        file( GLOB tmp LIST_DIRECTORIES true "${libdir}/include/*" )
        list(LENGTH tmp tmp)
        if ( NOT tmp EQUAL 1 )
          message(
            FATAL_ERROR
            "${libdir}/include can only contain an TestLib_${bn} subdir"
          )
        endif()
        target_include_directories( ${name} PUBLIC "${libdir}/include" )
      endif()
    endif()
  endforeach()
endfunction()

function(
    mctools_testutils_add_tests_apps
    approotdir extra_link_libs extra_inc_dirs envmod
  )
  set( testsbindir "${PROJECT_BINARY_DIR}/mctools_tests_bin" )
  file( MAKE_DIRECTORY "${testsbindir}" )
  file(
    GLOB appdirs
    LIST_DIRECTORIES true
    CONFIGURE_DEPENDS "${approotdir}/app_*"
  )
  foreach(appdir ${appdirs})
    get_filename_component(bn "${appdir}" NAME_WE)
    string(SUBSTRING "${bn}" 4 -1 "bn")
    mctools_testutils_internal_getsrcfiles( srcfiles is_cxx "${appdir}" "OFF" )
    add_executable( ${bn} ${srcfiles})
    if ( is_cxx )
      target_compile_features( ${bn} PRIVATE cxx_std_11 )
    endif()
    mctools_apply_strict_comp_properties( ${bn} )
    mctools_testutils_internal_detectlibdeps( "deplist" "${srcfiles}" "" )
    foreach( dep ${deplist} )
      target_link_libraries( ${bn} PRIVATE "TestLib_${dep}" )
    endforeach()
    target_link_libraries( ${bn} PRIVATE ${extra_link_libs} )
    target_include_directories( ${bn} ${extra_inc_dirs} )

    #Most apps won't need it, but we let an executable know its own path in the
    #build area, so we can unit test an executable's ability to find its own
    #path:
    target_compile_definitions(
      ${bn} PRIVATE "MCTOOLS_TESTAPP_FILE=\"$<TARGET_FILE:${bn}>\""
    )

    set( reflog "${appdir}/test.log" )
    if ( EXISTS "${reflog}" )
      set( testname "app_rl_${bn}" )
    else()
      set( testname "app_${bn}" )
      set( reflog "" )
    endif()
    mctools_testutils_internal_addtest(
      "${testname}"
      "$<TARGET_FILE:${bn}>"
      "${reflog}"
    )
    #A bit controversial perhaps, but a pain to maintain windows tests
    #otherwise:
    if ( WIN32 )
      list(
        APPEND envmod
        "PYTHONIOENCODING=set:UTF-8"
        "PYTHONLEGACYWINDOWSSTDIO=set:UTF-8"
      )
    endif()
    set_property(
      TEST "${testname}"
      PROPERTY ENVIRONMENT_MODIFICATION "${envmod}"
    )
  endforeach()
endfunction()

######################################
##### INTERNAL FUNCTIONS FOLLOWS #####
######################################

function( mctools_testutils_internal_getpyexec resvar )
  get_property(
    pyexec
    GLOBAL PROPERTY mctools_testutils_internal_pyexec
  )
  if ( pyexec )
    #Already found:
    set( "${resvar}" "${pyexec}" PARENT_SCOPE )
    return()
  endif()
  #Must figure out which one to use:
  if ( MCTOOLS_TESTUTILS_PYTHON_EXECUTABLE )
    if ( "x${MCTOOLS_TESTUTILS_PYTHON_EXECUTABLE}" STREQUAL "xauto" )
      find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter)
      set_property(
        GLOBAL PROPERTY
        mctools_testutils_internal_pyexec "${Python3_EXECUTABLE}"
      )
      set( "${resvar}" "${Python3_EXECUTABLE}" PARENT_SCOPE )
    else()
      set_property(
        GLOBAL PROPERTY mctools_testutils_internal_pyexec
        "${MCTOOLS_TESTUTILS_PYTHON_EXECUTABLE}"
      )
      set( "${resvar}" "${MCTOOLS_TESTUTILS_PYTHON_EXECUTABLE}" PARENT_SCOPE )
    endif()
  elseif( Python3_EXECUTABLE )
    set_property(
      GLOBAL PROPERTY
      mctools_testutils_internal_pyexec "${Python3_EXECUTABLE}"
    )
    set( "${resvar}" "${Python3_EXECUTABLE}" PARENT_SCOPE )
  else()
    message(
      FATAL_ERROR
      "Could not find Python3 executable (you can set "
      "MCTOOLS_TESTUTILS_PYTHON_EXECUTABLE or Python3_EXECUTABLE "
      "before calling functions in mctools_testutils)."
    )
  endif()
endfunction()

function( mctools_testutils_internal_detectlibdeps resvar files )
  #Determine dependencies on test libs by looking
  #for '#include #"TestLib_<DEPNAME>/'
  set(deplist "")
  foreach( file ${files} )
    file( STRINGS "${file}" inclines REGEX "^#include *\"TestLib_.*/.*\"" )
    foreach( entry ${inclines} )
      string(FIND "${entry}" "\"TestLib_" loc)
      math(EXPR loc "${loc}+9")#9 is length of '"TestLib_'
      string(SUBSTRING "${entry}" ${loc} -1 entry)
      string(FIND "${entry}" "/" loc)
      string(SUBSTRING "${entry}" 0 ${loc} entry)
      if( NOT "${entry}" IN_LIST deplist )
        list( APPEND deplist "${entry}" )
      endif()
    endforeach()
  endforeach()
  set( "${resvar}" "${deplist}" PARENT_SCOPE )
endfunction()

function( mctools_testutils_internal_detectpydeps resvar pyfile )
  #Determine dependencies on optional python modules by looking for lines like:
  #   "# NEEDS: numpy spglib"
  #Their presence will then be determined or not based on whether "import
  #<name>" works in python.
  set(deplist "")
  file( STRINGS "${pyfile}" needs_lines REGEX "^# NEEDS: " )
  set( comment_marker_seen OFF )
  foreach( entry ${needs_lines} )
    string(SUBSTRING "${entry}" 9 -1 entry)#9 is length of "# NEEDS: "

    #Support comment '#' char after the deplist:
    string(REPLACE "#" " # " entry "${entry}")#add space around comment marker
    string(REPLACE " " ";" entry "${entry}")#Turn into CMake list
    foreach(depname ${entry})
      if ( "x${depname}" STREQUAL "x#" )
        break()#break on comment marker
      endif()
      if (NOT "${depname}" IN_LIST deplist )
        list(APPEND deplist "${depname}")
      endif()
    endforeach()
  endforeach()
  set( "${resvar}" "${deplist}" PARENT_SCOPE )
endfunction()

function( mctools_testutils_internal_getsrcfiles
    resvar_srcfiles resvar_iscxx dir is_module  )
  file(
    GLOB srcfiles_cxx LIST_DIRECTORIES false CONFIGURE_DEPENDS
    "${dir}/*.cc"
  )
  file(
    GLOB srcfiles_c LIST_DIRECTORIES false CONFIGURE_DEPENDS
    "${dir}/*.c"
  )
  if ( NOT srcfiles_cxx AND NOT srcfiles_c )
    message(FATAL_ERROR "No source files found in ${dir}")
  endif()
  if ( srcfiles_cxx AND srcfiles_c )
    message(FATAL_ERROR "Both C and C++ source files found in ${dir}")
  endif()
  set_source_files_properties(
    ${srcfiles_c}
    PROPERTIES LANGUAGE C
  )
  if ( is_module )
    set_source_files_properties(
      ${srcfiles_cxx} ${srcfiles_c}
      PROPERTIES VISIBILITY_INLINES_HIDDEN "ON"
    )
    set_source_files_properties(
      ${srcfiles_cxx}
      PROPERTIES CXX_VISIBILITY_PRESET "hidden"
    )
    set_source_files_properties(
      ${srcfiles_c}
      PROPERTIES C_VISIBILITY_PRESET "hidden"
    )
  endif()
  set( res "" )
  list( APPEND res ${srcfiles_cxx} )
  list( APPEND res ${srcfiles_c} )
  set( ${resvar_srcfiles} "${res}" PARENT_SCOPE )
  if ( srcfiles_cxx )
    set( ${resvar_iscxx} "ON" PARENT_SCOPE )
  else()
    set( ${resvar_iscxx} "OFF" PARENT_SCOPE )
  endif()
endfunction()

function( mctools_testutils_internal_addtest name cmd_file reflog )
  #Make number of tests added available in global variable (need to use internal
  #cache var since there might be nested function calls, so PARENT_SCOPE is not
  #enough):
  if ( NOT DEFINED mctools_tests_count )
    set( tmpcount "1" )
  else()
    math( EXPR tmpcount "${mctools_tests_count}+1" )
  endif()
  set(
    mctools_tests_count "${tmpcount}"
    CACHE INTERNAL "number of tests added" FORCE
  )
  #Wrap all tests in our launcher script:
  mctools_testutils_internal_getpyexec( "pyexec" )
  add_test(
    NAME "${name}"
    COMMAND "${pyexec}" "${mctools_launcher_file}" "${cmd_file}" ${reflog}
  )
  set( wd "${PROJECT_BINARY_DIR}/tests/rundirs_$<CONFIG>/${name}" )
  if ( CMAKE_CONFIGURATION_TYPES )
    set( cfgtypelist "${CMAKE_CONFIGURATION_TYPES}" )
  else()
    set( cfgtypelist "${CMAKE_BUILD_TYPE}" )
  endif()
  foreach( cfgval ${cfgtypelist} )
    file( MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/tests/rundirs_${cfgval}/${name}" )
  endforeach()
  set_property( TEST "${name}" PROPERTY WORKING_DIRECTORY "${wd}" )

  #Set a timeout to prevent hanging jobs and test times getting out of hand
  #(exact values can be revisited). We allow longer jobs in Debug mode:
  set_property( TEST "${name}" PROPERTY TIMEOUT "$<IF:$<CONFIG:Debug>,480,240>" )
  #TODO: Support tests/costs.txt for adding cost properties to tests?
endfunction()
