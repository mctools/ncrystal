
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

cmake_minimum_required(VERSION 3.20...3.31)

# Note: There is no find_package(NCrystal) in this file, since the point of this
# example is to show how one can access NCrystal features via a plugin mechanism
# at run-time.

project( ExampleProject VERSION 0.0.1 LANGUAGES CXX )

add_executable( "ncrystal-virtapi-example" "src/main.cc" "src/ncrystal_load.cc" )
if ( UNIX AND CMAKE_DL_LIBS )
  target_link_libraries( "ncrystal-virtapi-example" PRIVATE ${CMAKE_DL_LIBS} )
endif()

#Require C++11:
target_compile_features( "ncrystal-virtapi-example" PUBLIC cxx_std_11 )

install (TARGETS "ncrystal-virtapi-example" RUNTIME DESTINATION bin)
