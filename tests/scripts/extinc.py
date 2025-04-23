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
mp.dps = 200
mpf = mp.mpf
import NCTestUtils.enable_fpe # noqa F401

lib = Lib('extinction')
lib.dump()

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

def ref_calcSabineEb(x,y):
    x = mpf(x)
    y = mpf(y)
    A = ref_calcSabineA( y )
    B = ref_calcSabineB( y )
    return A / mp.sqrt( mpf(1) + B * x )

def ref_calcSabineEl(x,y):
    x = mpf(x)
    y = mpf(y)
    if x <= mpf(1):
        k = ( mpf(1)-(x/mpf(2))+((x**2)/mpf(4))
              -((x**3)*mpf(5)/mpf(48.))+((x**4)*mpf(7)/mpf(192)) )
    else:
        k = ( mp.sqrt( mpf(2)/(mp.pi*x))*
              ( mpf(1)
                - mpf(1)/(mpf(8)*x)
                - mpf(3)/(mpf(128)*(x**2))
                - mpf(15) / (mpf(1024)*(x**3)) )
             )
    return k * mp.exp( -y )

def cmp( fct, reffct, x, fctname, x2 = None ):
    def fmt(v):
        return f'{float(v):.14g}'
    argstr = fmt(x) if x2 is None else '%s, %s'%(fmt(x), fmt(x2))
    args = (x,) if x2 is None else (x,x2)
    y = mpf(fct(*args))
    print(f"   {fctname}({argstr}) = {fmt(y)}")
    yr = mpf(reffct(*args))
    if not mpmath.almosteq(y,yr, rel_eps=1e-15, abs_eps=1e-99):
        rd = y/yr - mpf(1)
        raise SystemExit(f'Function produced wrong output at argument ({argstr}).'
                         f' Got {fmt(y)} but expected {fmt(yr)} (reldiff: {float(rd):e})')

def testAB():
    yvals = [ 0.0, 1e-99, 1e-20, 1e-3, 0.1,
              0.2-1e-3, 0.2, 0.2+1e-3,
              0.3-1e-3, 0.3, 0.3+1e-3,
              0.3, 0.7, 0.9, 20, 1e3, 1e5, 1e10, 1e50, 1e99, 1e198, 1e300 ]
    xvals = [ 0.0, 1e-99, 0.1, 0.5, 0.999999, 1.0, 1.0000001, 10.0, 1e50,
              1e300 ]

    for y in yvals:
        cmp(lib.nctest_calcSabineA, ref_calcSabineA,y,"calcSabineA" )
    for y in yvals:
        cmp(lib.nctest_calcSabineB, ref_calcSabineB,y,"calcSabineB" )

    for y in yvals:
        for x in xvals:
            cmp(lib.nctest_calcSabineEb, ref_calcSabineEb,x,"calcSabineEb", x2 = y)

    for y in yvals:
        for x in xvals:
            cmp(lib.nctest_calcSabineEl, ref_calcSabineEl,x,"calcSabineEl", x2 = y)

def main():
    testAB()

if __name__=='__main__':
    main()
