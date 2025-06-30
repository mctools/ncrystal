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

def ref_calcSabineEb(x,y=0.0):
    x = mpf(x)
    y = mpf(y)
    A = ref_calcSabineA( y )
    B = ref_calcSabineB( y )
    return A / mp.sqrt( mpf(1) + B * x )

def ref_calcSabineElOriginal(x,y=0.0):
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

def ref_exact_El( x, y = 0 ):
    # Redoing Sabine's calculations leading to eq. 6.4.5.3 and 6.4.5.4, I (TK)
    # found the following exact form:
    #
    # El = exp(-y) * exp(-x) * ( I0(x)+I1(x) )
    #
    # Where I0(x) and I1(x) are the modified bessel functions of the first kind,
    # for nu=0 and nu=1 respectively. In mpmath these are given as besseli(0,x)
    # and besseli(1,x) (in sagemath the function name is bessel_I).
    x, y = mpf(x), mpf(y)
    return mp.exp( -y ) * mp.exp( -x ) * ( mp.besseli(0,x) + mp.besseli(1,x) )

def cmp( fct, reffct, x, fctname, x2 = None, rel_eps=1e-15, abs_eps=1e-99,
         fmtstr_res='%.14g' ):
    def fmt(v):
        return f'{float(v):.14g}'
    argstr = fmt(x) if x2 is None else '%s, %s'%(fmt(x), fmt(x2))
    args = (x,) if x2 is None else (x,x2)
    y = mpf(fct(*args))
    print(f"   {fctname}({argstr}) = {fmtstr_res}"%y)
    yr = mpf(reffct(*args))
    if not mpmath.almosteq(y,yr, rel_eps=rel_eps, abs_eps=abs_eps):
        rd = y/yr - mpf(1)
        raise SystemExit(f'Function produced wrong output at argument ({argstr}).'
                         f' Got {fmt(y)} but expected {fmt(yr)} (reldiff: {float(rd):e})')

def testAB():
    yvals = [ 0.0, 1e-99, 1e-20, 1e-3, 0.1,
              0.2-1e-3, 0.2, 0.2+1e-3,
              0.3-1e-3, 0.3, 0.3+1e-3,
              0.3, 0.7, 0.9, 20, 1e3, 1e5, 1e10, 1e50, 1e99, 1e198, 1e300 ]
    xvals = [ 0.0, 1e-99, 0.1, 0.5, 0.999999, 1.0, 1.0000001,
              1.3, 1.5, 1.57, 1.6, 2.0, 2.5, 4.0, 5.0, 6.0,
              6.39,6.4, 10.0, 1e50, 1e300 ]

    for y in yvals:
        cmp(lib.nctest_calcSabineA, ref_calcSabineA,y,"calcSabineA" )
    for y in yvals:
        cmp(lib.nctest_calcSabineB, ref_calcSabineB,y,"calcSabineB" )

    for y in yvals:
        for x in xvals:
            cmp(lib.nctest_calcSabineEb, ref_calcSabineEb,x,
                "calcSabineEb", x2 = y)

    for y in yvals:
        for x in xvals:
            cmp(lib.nctest_calcSabineElOriginal, ref_calcSabineElOriginal,x,
                "calcSabineElOriginal", x2 = y)

    for y in yvals:
        for x in xvals:
            cmp(lib.nctest_calcSabineEb_CachedAB, ref_calcSabineEb,x,
                "calcSabineEb_CachedAB", x2 = y)

    for x in xvals:
        cmp( lib.nctest_calcSabineEb_y0, ref_calcSabineEb,x,"calcSabineEb_y0" )

    for x in xvals:
        cmp( lib.nctest_calcSabineElOriginal_y0, ref_calcSabineElOriginal,x,
             "calcSabineElOriginal_y0" )

    for x in xvals:
        #we are comparing the exact formula with the approximation used in
        #NCrystal, so we reduce the rel_eps accordingly.
        eps = 1e-9
        if ( 1.06 < x < 2.35 ) or ( 5.65 < x < 11.04 ):
            #used expansions less precise in this area, this is a known feature
            eps = 1e-6
        cmp( lib.nctest_calcSabineEl_y0, ref_exact_El,x,
             "calcSabineEl_y0", rel_eps=eps )#, fmtstr_res='%.6g' )

def main():
    testAB()

if __name__=='__main__':
    main()
