
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

import pathlib
tmp = pathlib.Path(__file__).resolve().absolute()
tmp = tmp.parent.parent.parent

genroot = tmp.joinpath('autogen')
reporoot = tmp.parent.parent

srcroot = reporoot.joinpath('ncrystal_core')
srcincroot = reporoot.joinpath('ncrystal_core/include')
testroot = reporoot.joinpath('tests')
pysrcroot = reporoot.joinpath('ncrystal_python')
exsrcroot = reporoot.joinpath('examples')
datadir = reporoot.joinpath('data')
ncg4srcroot = reporoot.joinpath('ncrystal_geant4')

del tmp
