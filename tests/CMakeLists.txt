################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

#Tests require newer cmake (NB: Synchronize minimum value with that found in the
#check for NCRYSTAL_ENABLE_TESTING support in the top-level CMakeLists.txt
#file):
cmake_minimum_required(VERSION 3.28...3.30)

project( NCrystalTests LANGUAGES CXX C )

if ( NOT TARGET NCrystal::NCrystal )
  if ( NOT NCrystal_DIR )
    #Support pip-installed NCrystal:
    execute_process(
      COMMAND ncrystal-config --show cmakedir
      OUTPUT_VARIABLE NCrystal_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  endif()
  message(
    STATUS
    "Testing against existing NCrystal installation in ${NCrystal_DIR}"
  )
  find_package( NCrystal 3.9.7 REQUIRED )
  set( NCRYSTAL_ENABLE_PYTHON ON )
  set( NCTESTS_STANDALONE ON )
else()
  message(STATUS "Testing as part of the NCrystal build")
  set( NCTESTS_STANDALONE OFF )
endif()

enable_testing()

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake")
include( mctools_testutils )

if ( NOT NCRYSTAL_ENABLE_PYTHON )
  message( WARNING "Skipping python tests since NCRYSTAL_ENABLE_PYTHON=OFF")
else()
  set( nctest_pyenvmod "" )
  if ( NOT NCTESTS_STANDALONE )
    add_custom_target(NCrystalTestsCopyPyModules ALL
      COMMAND ${CMAKE_COMMAND} -E copy_directory
      "${NCrystal_SOURCE_DIR}/NCrystal"
      "${PROJECT_BINARY_DIR}/nctests_pypath/NCrystal"
      DEPENDS "${NCrystal_SOURCE_DIR}/NCrystal" )
    list(
      APPEND nctest_pyenvmod
      "PYTHONPATH=path_list_prepend:${PROJECT_BINARY_DIR}/nctests_pypath"
    )
    list(
      APPEND nctest_pyenvmod
      "NCRYSTAL_LIB=set:${NCrystal_BINARY_DIR}/${NCrystal_LIBNAME}"
    )
  endif()
  mctools_testutils_add_tests_pyscripts(
    "${CMAKE_CURRENT_LIST_DIR}/scripts" "${nctest_pyenvmod}"
  )
endif()
mctools_testutils_add_test_libs(
  "${CMAKE_CURRENT_LIST_DIR}" NCrystal::NCrystal
)
mctools_testutils_add_tests_apps(
  "${CMAKE_CURRENT_LIST_DIR}" NCrystal::NCrystal
)
