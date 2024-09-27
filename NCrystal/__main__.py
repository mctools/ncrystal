
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

"""

This __main__.py module is here to ensure that NCrystal commandline script
can be accessed even when not terminal scripts (entry points) are
installed. This is for instance done via:

$> python3 -mNCrystal nctool <args>
$> python3 -mNCrystal ncmat2cpp <args>
$> python3 -mNCrystal ncrystal_ncmat2cpp <args>

Note that the two last commands above run the same command line script.

"""
from .cliutils import _run_from_main_init
_run_from_main_init()
