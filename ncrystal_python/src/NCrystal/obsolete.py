
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

Obsolete functions

"""

from ._common import ncgetenv_bool as _ncgetenv_bool

def decodecfg_packfact(cfgstr):
    """OBSOLETE FUNCTION (always returns 1.0 now)."""
    from . import _common as nc_common
    nc_common.warn('The decodecfg_packfact function is obsolete and now always'
                   ' returns 1.0. It will be removed in a future release.')
    return 1.0

def getFileContents(name):
    """OBSOLETE FUNCTION: Use createTextData(..).rawData instead."""
    from . import _common as nc_common
    nc_common.warn('The getFileContents(name) function is obsolete. Please use the'
                   ' createTextData function instead (specifically calling '
                   'createTextData(name).rawData will produce the same results getFileContents(name).')
    from .core import createTextData
    return createTextData(name).rawData

def clearInfoCaches():
    """Deprecated. Does the same as clearCaches()"""
    from . import _common as nc_common
    nc_common.warn('The clearInfoCaches function is deprecated. Please'
                   ' call the clearCaches() function instead.')
    from .core import clearCaches
    clearCaches()

def disableCaching():
    """Obsolete function. Instead call clearCaches() as needed."""
    raise RuntimeError('The disableCaching() function has been obsoleted and no longer works. Users can'
                       +' instead call the clearCaches() function if really needed to clear the caches.')

def enableCaching():
    """Obsolete function. Instead call clearCaches() as needed."""
    raise RuntimeError('The enableCaching function has been removed. Users can'
                       +' instead call the clearCaches function if really needed to clear the caches.')

if not _ncgetenv_bool('NOPYOBSOLETE'):
    _ = globals()
    for _k in [_k for _k in _.keys() if ( hasattr(_k,'startswith') and not _k.startswith('_') )]:
        del _[_k]
