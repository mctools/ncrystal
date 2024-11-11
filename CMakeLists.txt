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

cmake_minimum_required(VERSION 3.16...3.27)

#Fixme fix CMAKE_BUILD_TYPE here as well as in ncrystal_core subdir?

project( NCrystalMeta LANGUAGES C CXX )

add_subdirectory("${PROJECT_SOURCE_DIR}/ncrystal_core")
if ( NCRYSTAL_ENABLE_TESTING )
  enable_testing()
  add_subdirectory("${PROJECT_SOURCE_DIR}/tests")
endif()
