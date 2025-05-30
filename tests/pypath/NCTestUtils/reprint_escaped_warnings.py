
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

Catch escaped warnings and print them without a stack trace (since stack traces
tend to make test output reproducible).

Additionally, use the NCrystalDev._common._add_warn_counts_to_msgs hook to
ensure that any NCrystal warnings gets a running warning count appended to the
message, to prevent the warning framework from silencing repeated warnings.

"""

__all__=[]

import warnings
import NCrystalDev._common as _c

_c._add_warn_counts_to_msgs[0] = True

def printWarnings( message, category, *args, **kwargs ):
    import NCrystalDev._common as nc_common
    msg,cat = nc_common.WarningSpy._fmtwarning( message, category )
    print(f'CAUGHT ESCAPED WARNING ({cat}): {msg}')

warnings.showwarning = printWarnings

del _c
del warnings
del printWarnings
