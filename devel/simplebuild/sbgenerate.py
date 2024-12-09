
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

# Make the sbgen module (from ./pypath/sbgen) and the ncrystal_repo_tools module
# (from (../pypath/ncrystal_repo_tools) available:
import pathlib
import sys
d = pathlib.Path(__file__).resolve().absolute().parent
d1 = d / 'pypath'
d2 = d.parent / 'pypath'
assert d1.is_dir()
assert d2.is_dir()
sys.path.insert(0,str(d2))
sys.path.insert(0,str(d1))
#Now launch the main generator entry point:
import sbgen.main # noqa E402
sbgen.main.main()

