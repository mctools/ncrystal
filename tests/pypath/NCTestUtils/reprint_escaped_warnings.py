
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

Catch escaped warnings and print them without a stack tracke (since the stack
traces tend to make the test output irreprodicible.

"""

__all__=[]

import warnings
def printWarnings( message, category, *args, **kwargs ):
    import NCrystalDev._common as nc_common
    msg,cat = nc_common.WarningSpy._fmtwarning( message, category )
    print(f'CAUGHT ESCAPED WARNING ({cat}): {msg}')
warnings.showwarning = printWarnings
