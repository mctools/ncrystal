
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
"""

   Importing this module clips FP precision of NCrystal._common.ncpprint.

"""

def _clip_floats(obj):
    #Pass all floats in data structure (assumed loaded from JSON) through
    #'%.12g' to reduce FP fluctuations.
    if isinstance(obj, float):
        return float('%.12g'%obj)
    if isinstance(obj, dict):
        return {k: _clip_floats(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clip_floats(v) for v in obj]
    return obj

def doit():
    import NCrystalDev._common as _c
    _ncpprint_orig = _c.ncpprint
    def ncpprint_stablefp( obj, **kwargs ):
        _ncpprint_orig(_clip_floats(obj),**kwargs)
    _c.ncpprint = ncpprint_stablefp

doit()
del doit

