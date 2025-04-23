#!/usr/bin/env python3

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

# NEEDS: mpmath

from NCTestUtils.loadlib import Lib
import mpmath
mp = mpmath.mp
mp.dps = 100
mpf = mp.mpf
import NCTestUtils.enable_fpe # noqa F401

lib = Lib('extinction')
lib.dump()
assert hasattr(lib,'nctest_calcSabineA')
assert hasattr(lib,'nctest_calcSabineB')

def ref_calcSabineA(y):
    y = mpf(y)
    if y == mpf(0):
        return mpf(1)
    return mp.exp(-y)*mp.sinh(y)/y

def ref_calcSabineB(y):
    y = mpf(y)
    if y == mpf(0):
        return mpf(1)
    #Evaluate 1/y - exp(-y)/sinh(y) for y>=0
    return mpf(1)/y - mp.exp(-y)/mp.sinh(y)

def cmp( fct, reffct, x, fctname ):
    y = mpf(fct(x))
    print(f"   {fctname}({float(x):.14g}) = {float(y):.14g}")
    yr = mpf(reffct(x))
    if not mpmath.almosteq(y,yr, rel_eps=1e-15, abs_eps=1e-99):
        rd = y/yr - mpf(1)
        raise SystemExit(f'Function produced wrong output at argument {x}.'
                         f' Got {y} but expected {yr} (reldiff: {float(rd):e})')

def testAB():
    yvals = [ 0.0, 1e-99, 1e-20, 1e-3, 0.1,
              0.2-1e-3, 0.2, 0.2+1e-3,
              0.3-1e-3, 0.3, 0.3+1e-3,
              0.3, 0.7, 0.9, 20, 1e3, 1e5, 1e10, 1e50, 1e99, 1e198, 1e300 ]
    for y in yvals:
        cmp(lib.nctest_calcSabineA, ref_calcSabineA,y,"calcSabineA" )
    for y in yvals:
        cmp(lib.nctest_calcSabineB, ref_calcSabineB,y,"calcSabineB" )

def main():
    testAB()

if __name__=='__main__':
    main()
