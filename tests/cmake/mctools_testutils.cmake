
# Utilities for building up a testing infrastructure based on C++/C apps and
# libraries as well as python scripts. The API functions from this file are:
#
#   mctools_testutils_add_tests_pyscripts
#   mctools_testutils_add_test_libs
#   mctools_testutils_add_tests_apps
#
# This file (and the associated mctools_testlauncher.py file) is intended to be
# eventually shared between at least NCrystal and MCPL projects.

include_guard()

set(
  "mctools_launcher_file"
  "${CMAKE_CURRENT_LIST_DIR}/mctools_testlauncher.py"
)
if ( NOT EXISTS "${mctools_launcher_file}" )
  message( FATAL_ERROR "Did not find expected ${mctools_launcher_file}" )
endif()

if ( NOT DEFINED MCTOOLS_PYTHON_EXECUTABLE )
  set( MCTOOLS_PYTHON_EXECUTABLE "auto" )
endif()

function( mctools_testutils_add_tests_pyscripts scriptsdir pyenvmod )
  file(
    GLOB pyscriptlist
    LIST_DIRECTORIES false
    CONFIGURE_DEPENDS "${scriptsdir}/*.py"
  )
  set( pyexec "" )
  foreach(pyscript ${pyscriptlist})
    if ( NOT pyexec )
      mctools_testutils_internal_getpyexec( "pyexec" )
    endif()
    get_filename_component(bn "${pyscript}" NAME_WE)
    get_filename_component(psdir "${pyscript}" DIRECTORY)
    if ( EXISTS "${psdir}/${bn}.log" )
      set( testname "py_rl_${bn}" )
      mctools_testutils_internal_addreflogtest(
        "${testname}"
        "${pyscript}"
        "${psdir}/${bn}.log"
      )
    else()
      add_test( NAME "py_${bn}" COMMAND "${pyexec}" "${pyscript}" )
      mctools_testutils_internal_settestprops( "py_${bn}" )
    endif()
    if ( pyenvmod )
      set_property(
        TEST "${testname}"
        PROPERTY ENVIRONMENT_MODIFICATION "${pyenvmod}"
      )
    endif()
  endforeach()
endfunction()

function( mctools_testutils_add_test_libs librootdir extra_link_libs )
  file(
    GLOB libdirs
    LIST_DIRECTORIES true
    CONFIGURE_DEPENDS "${librootdir}/lib_*"
  )
  foreach(libdir ${libdirs})
    get_filename_component(bn "${libdir}" NAME)
    string(SUBSTRING "${bn}" 4 -1 "bn")
    set( name "TestLib_${bn}" )
    mctools_testutils_internal_getsrcfiles( srcfiles "${libdir}" )
    mctools_testutils_internal_detectlibdeps( "deplist" "${srcfiles}" )
    foreach( dep ${deplist} )
      if ( NOT "x${dep}" STREQUAL "x${bn}" )
        #Private because not headers (TODO: Add headers as explicit deps?):
        target_link_libraries( ${bn} PRIVATE "TestLib_${dep}" )
      endif()
    endforeach()
    add_library( ${name} ${srcfiles} )
    target_link_libraries( ${name} PRIVATE ${extra_link_libs} )
    if ( EXISTS "${libdir}/include" )
      if ( NOT EXISTS "${libdir}/include/TestLib_${bn}" )
        message(
          FATAL_ERROR
          "${libdir}/include does not include TestLib_${bn} subdir"
        )
      endif()
      file( GLOB tmp LIST_DIRECTORIES true "${libdir}/include/*"
      )
      list(LENGTH tmp tmp)
      if ( NOT tmp EQUAL 1 )
        message(
          FATAL_ERROR
          "${libdir}/include can only contain an TestLib_${bn} subdir"
        )
      endif()
      target_include_directories( ${name} PUBLIC "${libdir}/include" )
    endif()
  endforeach()
endfunction()

function( mctools_testutils_add_tests_apps approotdir )
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
    mctools_testutils_internal_getsrcfiles( srcfiles "${appdir}" )
    add_executable( ${bn} ${srcfiles})
    mctools_testutils_internal_detectlibdeps( "deplist" "${srcfiles}" "" )
    foreach( dep ${deplist} )
      target_link_libraries( ${bn} PRIVATE "TestLib_${dep}" )
    endforeach()
    target_link_libraries( ${bn} PRIVATE ${extra_link_libs} )
    set_target_properties(
      ${bn} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${testsbindir}/${bn}"
      RUNTIME_OUTPUT_NAME  "${bn}"
    )
    if ( EXISTS "${appdir}/test.log" )
      mctools_testutils_internal_addreflogtest(
        "app_rl_${bn}"
        "${testsbindir}/${bn}/${bn}"
        "${appdir}/test.log"
      )
    else()
      add_test( NAME "app_${bn}" COMMAND ${bn} )
      mctools_testutils_internal_settestprops( "app_${bn}" )
    endif()
  endforeach()
endfunction()

######################################
##### INTERNAL FUNCTIONS FOLLOWS #####
######################################

function( mctools_testutils_internal_getpyexec resvar )
  if ( MCTOOLS_PYTHON_EXECUTABLE )
    if ( "x${MCTOOLS_PYTHON_EXECUTABLE}" STREQUAL "xauto" )
      find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter)
      set( "MCTOOLS_PYTHON_EXECUTABLE" "${Python3_EXECUTABLE}" PARENT_SCOPE )
      set( "${resvar}" "${Python3_EXECUTABLE}" PARENT_SCOPE )
    else()
      set( "${resvar}" "${MCTOOLS_PYTHON_EXECUTABLE}" PARENT_SCOPE )
    endif()
  elseif( Python3_EXECUTABLE )
    set( "${resvar}" "${Python3_EXECUTABLE}" PARENT_SCOPE )
  else()
    message(
      FATAL_ERROR
      "MCTOOLS_PYTHON_EXECUTABLE or Python3_EXECUTABLE"
      "must be set before calling functions in mctools_testutils"
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

function( mctools_testutils_internal_getsrcfiles resvar_srcfiles dir )
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
  set_source_files_properties(${srcfiles_cxx} PROPERTIES LANGUAGE CXX)
  set_source_files_properties(${srcfiles_c} PROPERTIES LANGUAGE C)
  set( res "" )
  list( APPEND res ${srcfiles_cxx} )
  list( APPEND res ${srcfiles_c} )
  set( ${resvar_srcfiles} "${res}" PARENT_SCOPE )
endfunction()

function( mctools_testutils_internal_settestprops name )
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/rundirs/run_${name}")
  set_property(
    TEST "${name}"
    PROPERTY WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/rundirs/run_${name}"
  )
  #Default to 5 second timeout (at some point, we can increase this value for
  #some tests):
  set_property( TEST "${name}" PROPERTY TIMEOUT 5 )
endfunction()

function( mctools_testutils_internal_addreflogtest name cmd_file reflog )
  mctools_testutils_internal_getpyexec( "pyexec" )
  add_test(
    NAME "${name}"
    COMMAND "${pyexec}" "${mctools_launcher_file}" "${cmd_file}" "${reflog}"
  )
  mctools_testutils_internal_settestprops( "${name}" )
endfunction()