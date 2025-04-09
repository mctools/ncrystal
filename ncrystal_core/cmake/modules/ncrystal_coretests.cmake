
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


# Tests which can be run when only ncrystal_core is built, and which is even OK
# to use before an installation step (this is not the case with the tests in
# <reporoot>/tests).

include_guard()

function( ncrystal_core_setup_tests )
  set_source_files_properties(
    "${NCrystal_SOURCE_DIR}/app_test/main.cc"
    PROPERTIES LANGUAGE CXX
  )
  add_executable( "coretestapp" "${NCrystal_SOURCE_DIR}/app_test/main.cc" )
  target_link_libraries( "coretestapp" PRIVATE NCrystal::NCrystal )
  add_test( NAME "coretestapp" COMMAND "coretestapp" )
  set_property( TEST "coretestapp" PROPERTY TIMEOUT "$<IF:$<CONFIG:Debug>,20,7>" )
endfunction()

ncrystal_core_setup_tests()
