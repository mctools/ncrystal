
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

#import NCrystalDev as NC
import math
import mpmath
mp = mpmath.mp
mp.dps = 50
mpf = mp.mpf
mp_one = mpf(1)

def bc_std_yp_of_x( theta_deg ):
    cos2th = math.cos( theta_deg * math.pi / 90 )
    A = 0.20 + 0.45 * cos2th
    B = 0.22 - 0.12 * ( 0.5 - cos2th )**2
    #print(f"USING A={A} B={B}")
    import numpy as np
    def yp( x ):
        return 1.0 / np.sqrt( 1.0 + 2.0 * x + A * x * x / ( 1 + B*x ) )
    return yp

class BCYPCalc:

    def __init__( self, theta_deg ):
        self._mp = mpmath.mp
        mindps = 30#fixme increase!!
        if not self._mp.dps or self._mp.dps<mindps:
            print(f"Increasing mp.dps to {mindps}...")
            self._mp.dps = mindps
        self._mpf = self._mp.mpf


        self._theta = mpf(theta_deg)
        self._sinth = mp.sin(self._theta * mp.pi/mpf(180) )
        self._sqrtsinth = mp.sqrt( self._sinth )
        self._sinth_to_3div2 = self._sinth * self._sqrtsinth

    @property
    def mp( self ):
        print(self._mp)
        return self._mp

    @property
    def theta_deg( self):
        return self._theta

    @property
    def sintheta( self):
        return self._sinth

    def phi( self, sigma_r ):
        #BC 1974 eq. 34 ( sr = sigma*r )
        sr = mpf(sigma_r)
        sinth = self._sinth
        assert sr >= 0.0
        if not sinth:
            return phi0( sr )
        if sinth == mp_one:
            return phipi( sr )
        u = sr * self._sqrtsinth
        return phi0( sr ) + self._sinth_to_3div2 * ( phipi(u) - phi0(u) )

    def integrand_fct_eq36( self, x ):
        x = mpf(x)
        phi = self.phi
        def integrand(eta):
            f = f_of_eta(eta)
            return f * phi( x * f )
        return integrand

    def calc_yp( self, x ):
        #print("...calc_yp( x=",x,")")
        k = mpf(6)/ (mpf(4)*mp.pi)
        i = self.integrand_fct_eq36( x )
        def calc( a, b ):
            return mp.quad( i, [ a, b ] )
        return k * (
            calc( mpf(0), mpf(1) )
            + calc( mpf(1), mpf(5) )
            + calc( mpf(5), mpf(10) )
            + calc( mpf(10), mpf(30) )
            + calc( mpf(30), mpf(200) )
            + calc( mpf(200), mpf(100000) )
        )

def phi0( sr ):
    #BC 1974 eq. 31 ( sr = sigma*r )
    sr = mpf(sr)
    assert sr >= 0
    if sr < 0.1:
        c0 = mpf("1")
        c1 = mpf("-3/2")
        c2 = mpf("8/5")
        c3 = mpf("-4/3")
        c4 = mpf("32/35")
        c5 = mpf("-8/15")
        c6 = mpf("256/945")
        c7 = mpf("-64/525")
        c8 = mpf("512/10395")
        c9 = mpf("-256/14175")
        c10 = mpf("4096/675675")
        c11 = mpf("-2048/1091475")
        c12 = mpf("16384/30405375")
        c13 = mpf("-2048/14189175")
        c14 = mpf("131072/3618239625")
        c15 = mpf("-16384/1915538625")
        c16 = mpf("131072/68746552875")
        c17 = mpf("-65536/162820783125")
        c18 = mpf("1048576/12993098493375")
        c19 = mpf("-524288/34029543673125")
        c20 = mpf("4194304/1494206326738125")
        return ( c0+sr*(c1+sr*(c2+sr*(c3+sr*(c4+sr*(c5+sr*(c6+sr*(
            c7+sr*(c8+sr*(c9+sr*(c10+sr*(c11+sr*(c12+sr*(c13+sr*(
                c14+sr*(c15+sr*(c16+sr*(c17+sr*(c18+sr*(c19+sr*c20)
                                                )))))))))))))))))) )
    foursr = mpf(4) * sr
    srsq = sr**2
    e = mp.exp( - foursr)
    #a = ( mpf(3)/(mpf(64)*srsq*sr) )
    a = mpf(3) / ( mpf(64)*(sr**3) )
    b = mpf(8) * srsq + foursr*e - ( mpf(1) - e )
    return a * b

def phipi( sr ):
    #BC 1974 eq. 32 ( sr = sigma*r )
    #ERRATA DISCOVERED THROUGH PAINFUL PROCESS: Factor of 1/2 has to be moved from
    #second term in the nominator to the third (logarithmic) term.
    sr = mpf(sr)
    assert sr >= 0
    if sr < 0.1:
        c0 = mpf("1")
        c1 = mpf("-3/2")
        c2 = mpf("12/5")
        c3 = mpf("-4")
        c4 = mpf("48/7")
        c5 = mpf("-12")
        c6 = mpf("64/3")
        c7 = mpf("-192/5")
        c8 = mpf("768/11")
        c9 = mpf("-128")
        c10 = mpf("3072/13")
        c11 = mpf("-3072/7")
        c12 = mpf("4096/5")
        c13 = mpf("-1536")
        c14 = mpf("49152/17")
        c15 = mpf("-16384/3")
        c16 = mpf("196608/19")
        c17 = mpf("-98304/5")
        c18 = mpf("262144/7")
        c19 = mpf("-786432/11")
        c20 = mpf("3145728/23")
        return ( c0+sr*(c1+sr*(c2+sr*(c3+sr*(c4+sr*(c5+sr*(c6+sr*(
            c7+sr*(c8+sr*(c9+sr*(c10+sr*(c11+sr*(c12+sr*(c13+sr*(
                c14+sr*(c15+sr*(c16+sr*(c17+sr*(c18+sr*(c19+sr*c20)
                                                )))))))))))))))))) )
    srsq = sr**2
    a = ( mpf(3)/(mpf(4)*(sr**3)) )
    #ORIGINAL: b = srsq - mpf('1/2')*sr + mp.log(mpf(1)+mpf(2)*sr)
    b = srsq - sr + mpf('1/2')*mp.log(mpf(1)+mpf(2)*sr) #FIXED
    return a * b

def _f_of_eta_taylor( eta ):
    #To be used for 0<=eta<=0.125 if using double precision:
    #constexpr double c0 = 1;
    #constexpr double c2 = -2./9;
    #constexpr double c4 = 1./45;
    #constexpr double c6 = -2./1575;
    #constexpr double c8 = 2./42525;
    #constexpr double c10 = -4./3274425;
    #constexpr double c12 = 1./42567525;
    c0 = mpf(1)
    c2 = mpf("-2/9")
    c4 = mpf("1/45")
    c6 = mpf("-2/1575")
    c8 = mpf("2/42525")
    c10 = mpf("-4/3274425")
    c12 = mpf("1/42567525")
    x = eta*eta
    return c0+x*(c2+x*(c4+x*(c6+x*(c8+x*(c10+x*c12)))))

def f_of_eta( eta, taylor = True ):
    # About f(eta) in equation 36. It is defined at the end of page 134 via: f(eta)
    # = sigma(eta1)/sigma(0), followed by definitions on top of page 135:
    #
    # sigma(0) = (3/4)*Q*beta
    # eta = pi*eta1*beta
    # beta = 2*r*sin2theta/lambda
    #
    # That leaves sigma(eta1) which I guess is given by eq. 29? So:
    #
    # f(eta) = ( eta^2-eta*sin(2*eta)+sin^2(eta) ) / eta^4
    eta = mpf(eta)
    if eta < 0:
        eta = -eta
    #Using lower threshold for taylor expansion due to mpmath:
    if not eta:
        taylor = True#avoid direct zero division
    if taylor and eta < mpf(0.125):
        return _f_of_eta_taylor( eta )

    #return mpf(1)/(mpf(1)+eta**2) #FIXME TESTING
    s = mp.sin(eta)
    s2 = mp.sin(2*eta)
    etasq = eta**2
    return  ( etasq - eta*s2 + s**2 ) / (etasq**2)
