#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

import NCrystal as NC
import mpmath
mp = mpmath.mp
mp.dps = 200
mpf = mp.mpf
import NCTestUtils.enable_fpe # noqa F401

class PiecewiseLinDist:

    """Utility class which takes a piecewise linear distribution function and
    (with mpmath), calculates standard moments. To be more easily interpretable,
    the highest order moment are scale independent standard moments, and the
    lowest order moments are provided in the following forms:

    n = 0 : The integral of the original distribution. The distribution will be
            normalised before calculating any higher order moment.
    n = 1 : The mean value of the distribution. Higher order moments will be
            central moments around this value.
    n = 2 : The square root of the variance (a.k.a. root-mean-square). All
            higher order moments will be normalised to a suitable power of the
            variance, to make them scale independent standard moments.
    n = 3 : skewness (scale independent)
    n = 4 : kurtosis (scale independent)

    """

    def __init__( self, mp, xgrid, y_values ):
        self.nbins = len(y_values)-1
        self.xvals, self.dx = self._unpack_grid( mp, xgrid, self.nbins )
        self._mp = mp
        mpf = mp.mpf
        self.yvals = [ mpf(e) for e in y_values ]
        #find integral and normalise:
        contribs = [ self.dx[i] * ( self.yvals[i] + self.yvals[i+1] )
                     for i in range(self.nbins) ]
        contribs.sort()
        self.integral = mp.fsum( contribs ) / 2
        normfact = 1.0 / self.integral
        for i in range(len(self.yvals)):
            self.yvals[i] *= normfact
        #Calculate dy, dy/dx and a=y1-x1*dy/dx in each bin:
        self.dy = [ ( self.yvals[i+1] - self.yvals[i] )
                    for i in range( self.nbins ) ]
        self.dydx = [ ( self.dy[i] / self.dx[i] )
                      for i in range(self.nbins) ]
        self._a = [ (self.yvals[i] - self.xvals[i] * self.dydx[i] )
                    for i in range(self.nbins) ]
        #find mean:
        xvals_pow2 = [ e**2 for e in self.xvals ]
        xvals_pow3 = [ e**3 for e in self.xvals ]
        contribs = []
        for i in range(self.nbins):
            k = self._a[i] / 2
            contribs.append( k * xvals_pow2[i+1] )
            contribs.append( - k * xvals_pow2[i] )
            k = self.dydx[i] / 3
            contribs.append( k * xvals_pow3[i+1] )
            contribs.append( - k * xvals_pow3[i] )
        contribs.sort()
        self.mean = mp.fsum( contribs )

        #offset with -mean for central moments:
        # u = x - mean
        self._u = [ x-self.mean for x in self.xvals ]
        # cache of (u**n)/n:
        self._cache_utondivn = {}
        self._cache_moments = {}

    def get_moment( self, n ):
        """Get the n'th moment of the distribution, for n>=0. Refer to the class
        description for details of what is returned.
        """

        assert isinstance(n,int)
        m = self._cache_moments.get(n)
        if m is not None:
            return m
        if n == 0:
            m = self.integral
        elif n == 1:
            m = self.mean
        else:
            m = self._calc_moment( n )
            if n == 2:
                m = m.sqrt() # variance -> rms
            else:
                #Standard moment:
                m /= self.rms**n
        self._cache_moments[n] = m
        return m

    @property
    def variance( self ):
        return self.get_moment(2) ** 2

    @property
    def rms( self ):
        return self.get_moment(2)

    @property
    def skewness( self ):
        return self.get_moment(3)

    @property
    def kurtosis( self ):
        return self.get_moment(4)

    def dump( self,
              n_max_moment = 4,
              *,
              title = None,
              fp_format = '%g' ):
        if title is not None:
            print(f'Distribution "{title}" with:')
        else:
            print("Distribution with:")
        names_and_vals = [ ('integral',self.integral),
                           ('mean',self.mean),
                           ('rms',self.rms) ]
        for n in range(3,n_max_moment+1):
            name = { 3 : 'skewness',
                     4 : 'kurtosis' }.get(n,f'{n}th moment (std)')
            names_and_vals.append( ( name, self.get_moment(n) ) )
        m = max( len(e[0]) for e in names_and_vals )
        for n,v in names_and_vals:
            print('  %s : %s'%( n.ljust(m), fp_format % v ) )

    @staticmethod
    def unit_test( mp ):
        rms = mp.mpf('1/3').sqrt()
        expected_moments = [ 6, 2, rms,
                             0, mp.mpf('1/5')/rms**4,
                             0, mp.mpf('1/7')/rms**6,
                             0 ]
        pw = PiecewiseLinDist( mp, [1.0,3.0],[3.0,3.0])
        assert all( mp.almosteq( pw.get_moment(i), m )
                    for i,m in enumerate(expected_moments) )
        pw = PiecewiseLinDist( mp, [1.0,3.0],[3.0,3.0,3.0,3.0,3.0])
        assert all( mp.almosteq( pw.get_moment(i), m )
                    for i,m in enumerate(expected_moments) )
        pw = PiecewiseLinDist( mp, [1.0,2.5,3.0],[3.0,3.0,3.0])
        assert all( mp.almosteq( pw.get_moment(i), m )
                    for i,m in enumerate(expected_moments) )

    def _unpack_grid( self, mp, xgrid, nbins ):
        assert nbins > 0
        mpf = mp.mpf
        if len(xgrid)==2 and nbins>1:
            #fixed grid:
            xmin, xmax = mpf(xgrid[0]), mpf(xgrid[1])
            dx = ( xmax - xmin ) / nbins
            xvals = [ (xmin + i * dx) for i in range(nbins) ]
            xvals.append( xmax )
            dxvals = [ dx for i in range(nbins) ]
        else:
            #flex grid:
            assert len(xgrid) == nbins + 1
            xvals = [ mpf(e) for e in xgrid ]
            dxvals = [ (xvals[i+1]-xvals[i]) for i in range(nbins) ]
        return xvals, dxvals

    def _get_utondivn( self, n ):
        assert isinstance(n,int) and n>=2
        #Get (u ** n) / n
        cv = self._cache_utondivn.get(n)
        if cv:
            return cv
        un = [ (e**n)/n for e in self._u ]
        self._cache_utondivn[n] = un
        return un

    def _calc_moment( self, n ):
        assert isinstance(n,int)
        unp1 = self._get_utondivn( n + 1 )
        unp2 = self._get_utondivn( n + 2 )
        contribs = []
        for i in range(self.nbins):
            k = self._a[i]
            contribs.append( k * unp1[i+1] )
            contribs.append( - k * unp1[i] )
            k = self.dydx[i]
            contribs.append( k * unp2[i+1] )
            contribs.append( - k * unp2[i] )
        contribs.sort()
        return self._mp.fsum( contribs )

def main():

    PiecewiseLinDist.unit_test(mp)

    cfgstr='stdlib::Au_sg225.ncmat;temp=20K;comp=inelas;vdoslux=3'
    info = NC.createInfo(cfgstr)
    di = info.findDynInfo('Au')
    def calcgn(n):
        return di.extract_Gn(n, expand_egrid=False, without_xsect=True)
    print()
    print('-'*80)
    print(f"Investigating {cfgstr}")
    print('-'*80)
    print()
    for n in [ 1, 2, 3, 5, 8, 14, 20, 50 ]:
        egrid,d = calcgn(n)
        PiecewiseLinDist(mp,egrid,d).dump(7,fp_format='%.7g',title=f'G{n}')
        #di.plot_Gn(n,without_xsect=True)




if __name__ == '__main__':
    main()
