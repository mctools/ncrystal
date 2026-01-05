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

# NEEDS: numpy mpmath

import NCrystalDev as NC
from NCTestUtils.pwlindistmoments import PWLinDistMoments
import mpmath
mp = mpmath.mp
mp.dps = 100
mpf = mp.mpf
import NCTestUtils.enable_fpe # noqa F401

def test(cfgstr):
    info = NC.createInfo(cfgstr)
    for di in info.dyninfos:
        lbl = di.atomData.displayLabel()
        print()
        print('-'*80)
        print(f'Investigating "{cfgstr}" / {lbl}')
        print('-'*80)
        print()
        class fakemp:
            def mpf(self,x):
                return float(x)
            def fsum(self,alist):
                import math
                return math.fsum(alist)
            def sqrt(self,x):
                import math
                return math.sqrt(x)
        mp = fakemp()
        if not hasattr(di,'extract_Gn'):
            print('  Ignoring unsuitable dyninfo type.')
            continue
        def calcgn(n):
            return di.extract_Gn(n, expand_egrid=False, without_xsect=True)
        for n in [ 1, 2, 5, 8, 20, 50 ]:
            #TDOO: This is slow, since we redo the whole procedure for each n!
            #We should provide a way to get multiple Gns in one call to the C++
            #layer! This actually makes the test very slow for Debug builds!
            egrid,d = calcgn(n)
            PWLinDistMoments(mp,egrid,d).dump(7,fp_format='%.7g',title=f'G{n}')
            #di.plot_Gn(n,without_xsect=True)

def main():
    PWLinDistMoments.unit_test(mp)
    for filename in ('Au_sg225',
                     'Polyethylene_CH2',
                     'LiH_sg225_LithiumHydride',
                     'C_sg194_pyrolytic_graphite' ):
        for temp in (10,600):
            for vdoslux in (0,3):
                test(f'stdlib::{filename}.ncmat'
                     f';temp={temp:g}K'
                     ';comp=inelas'
                     f';vdoslux={vdoslux}')

if __name__ == '__main__':
    main()
