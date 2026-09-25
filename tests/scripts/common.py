#!/usr/bin/env python3

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

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev._common as nc_common
import os

def require(b):
    if not b:
        raise RuntimeError('check failed')

for v in [(0.25,'1/4'),
          (0.13137,'.13137'),
          (0.75,'3/4'),
          (-0.0,'0'),
          (0.9999999999999,'.9999999999999'),
          (0.999999999999999889,'.999999999999999889'),
          (0.9999999999999999999889,'1'),
          (0.00,'0')]:
    vfmt = nc_common.prettyFmtValue(v[0])
    print( f'nc_common.prettyFmtValue({v[0]}):',repr(vfmt))
    require( vfmt == v[1] )
    if '/' in vfmt:
        _ = vfmt.split('/')
        fmtval = int(_[0]) / int(_[1])
    else:
        fmtval = float(vfmt)
    require( abs(fmtval-v[0]) < 1e-6 )
    require( (fmtval==1.0) == (v[0]==1.0) )
    require( (fmtval==0.0) == (v[0]==0.0) )

def testcf( c ):
    print(f'format_chemform({c}):', repr(nc_common.format_chemform(c)) )

testcf( [('Al',0.99),('Cr',0.005),('B10',0.005)] )
testcf( [('Al',0.9),('Cr',0.1)])
testcf( [('Al',1/3),('Cr',2/3)])

#ncpprint (must go via nc_common.print and support do_sort):
_d = {'b':[1,2],'a':{'z':1,'y':2}}
for _ds in (False,True):
    with nc_common.capture_print_ctxmgr() as _cap:
        nc_common.ncpprint(_d,do_sort=_ds)
    print(f'ncpprint(do_sort={_ds}):',repr(''.join(_cap.data)))

#ncsetenv must use the same env var names as ncgetenv (and C++), also for
#namespaced builds without namespaced env vars (simulated by forcing the flag):
_orig_nsev = nc_common._cache_nsev[0]
try:
    for _flag in (True, False):
        nc_common._cache_nsev[0] = _flag
        nc_common.ncsetenv('TESTENVROUNDTRIP','17')
        _v = nc_common.ncgetenv('TESTENVROUNDTRIP')
        nc_common.ncsetenv('TESTENVROUNDTRIP',None)
        print(f'ncsetenv+ncgetenv (namespaced env vars={_flag}):',repr(_v))
        require( _v == '17' )
        require( nc_common.ncgetenv('TESTENVROUNDTRIP') is None )
        require( not any('TESTENVROUNDTRIP' in k for k in os.environ) )
finally:
    nc_common._cache_nsev[0] = _orig_nsev

#NCMAT header comment extraction (all-empty comment lines once looped forever):
from NCrystalDev._ncmatimpl import _extractInitialHeaderCommentsFromNCMATData

for _hdr in ('#\n#  a\n#   b\n#\n', '#\n#\n', ''):
    _res = _extractInitialHeaderCommentsFromNCMATData(
        'NCMAT v7\n' + _hdr + '@DENSITY\n  1 g_per_cm3\n' )
    print(f'Header comments of {_hdr!r}:', _res)
