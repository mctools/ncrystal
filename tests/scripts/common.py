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
