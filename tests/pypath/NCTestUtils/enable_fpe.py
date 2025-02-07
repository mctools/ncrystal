
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
"""

   Importing this module enables FPE detection.

"""

__all__=[]

_can_catch = [False]
def is_catching_fpe():
    return _can_catch[0]

def _enable_fpe():
    from .loadlib import Lib
    lib = Lib('fpe')
    lib.nctest_catch_fpe()
    if lib.nctest_can_catch_fpe():
        _can_catch[0] = True
    return lib

__keepalive = _enable_fpe()
del _enable_fpe
