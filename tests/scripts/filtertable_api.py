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

# NEEDS: numpy

# Test the Python API for filter tables (NCrystal.filter.NCrystalFilter), which
# uses the C API (ncrystal_filtertable).

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.filter import NCrystalFilter
from NCrystalDev.misc import evaluate_query as ncquery
from NCrystalDev.constants import wl2ekin
from NCrystalDev.exceptions import NCBadInput
from NCTestUtils.common import ensure_error
import bisect
import numpy as np

def lookup( wl, xs, w ):
    #Independent (point-by-point) implementation of the evaluation:
    if not w > wl[0]:
        return xs[0]
    n = len(wl)
    if w >= wl[-1]:
        slope = ( xs[-1] - xs[-2] ) / ( wl[-1] - wl[-2] )
        return xs[-1] if slope == 0.0 else max( 0.0, xs[-1] + ( w - wl[-1] ) * slope )
    i = bisect.bisect_right( wl, w ) - 1
    assert 0 <= i < n - 1 and wl[i+1] > wl[i]
    return xs[i] + ( w - wl[i] ) / ( wl[i+1] - wl[i] ) * ( xs[i+1] - xs[i] )

def test( cfg ):
    f = NCrystalFilter( cfg )
    assert f.cfgstr == cfg and f.options is None
    wl, macroxs = f.table
    #Same as the table of the JSON query, converted to 1/cm:
    r = ncquery( ['filtertable', cfg] )
    tol = r['tol']
    assert np.array_equal( wl, np.asarray(r['wl']) )
    np.testing.assert_allclose( macroxs, np.asarray(r['xs']) * r['numberdensity'],
                                rtol = 1e-15, atol = 0.0 )
    #Within the tolerance of the exact macroscopic cross section:
    rng = np.random.default_rng( 1234 )
    w = np.exp( rng.uniform( np.log(1e-6), np.log(wl[-1]), 20000 ) )
    exact = ( NC.createScatter(cfg).xsect(wl=w) + NC.createAbsorption(cfg).xsect(wl=w)
              ) * NC.createInfo(cfg).numberdensity
    tab = f.xsect( wl = w )
    floor = 1e-12 * r['numberdensity']#1e-12 barn per atom, in 1/cm
    assert np.all( np.abs( tab - exact ) <= tol * np.maximum( exact, floor ) )
    #Agrees with the point-by-point lookup, also at table points (including
    #pairs at discontinuities), at 0, when extrapolating, and for ekin=0:
    w = np.concatenate( [ w[:500], wl, [ -1.0, 0.0, wl[-1], 1.2*wl[-1], 1e6, np.inf ] ] )
    tab = f.xsect( wl = w )
    for v, t in zip( w, tab ):
        ref = lookup( wl, macroxs, float(v) )
        assert t == ref or abs( t - ref ) <= 1e-14 * max( abs(ref), floor ), (v, t, ref)
        assert f.xsect( wl = float(v) ) == t
    #Energies instead of wavelengths, and scalars:
    np.testing.assert_allclose( f.xsect( ekin = wl2ekin( w[:500] ) ), tab[:500],
                                rtol = 1e-12, atol = 1e-12 * floor )
    assert isinstance( f.xsect( wl = 1.0 ), float )
    assert isinstance( f.xsect( ekin = 0.025 ), float )
    assert f.xsect( ekin = 0.0 ) == tab[-1]
    print(f'{cfg}: npts={len(wl)} xsect(0)={f.xsect(wl=0.0):.6g}/cm'
          f' xsect(4 Aa)={f.xsect(wl=4.0):.6g}/cm'
          f' xsect(1000 Aa, extrapolated)={f.xsect(wl=1000.0):.6g}/cm: OK')

def main():
    test( 'stdlib::Al_sg225.ncmat' )
    test( 'stdlib::Be_sg194.ncmat;temp=80K' )
    test( 'stdlib::Polyethylene_CH2.ncmat' )
    test( 'stdlib::void.ncmat' )
    f = NCrystalFilter( 'stdlib::Al_sg225.ncmat', '' )
    assert f.options is None
    for arr in f.table:
        with ensure_error(ValueError,'assignment destination is read-only'):
            arr[0] = 1.0
    with ensure_error(NCBadInput,'Please provide exactly one of the "ekin" or'
                      ' "wl" parameters.'):
        f.xsect()
    with ensure_error(NCBadInput,'Please provide exactly one of the "ekin" or'
                      ' "wl" parameters.'):
        f.xsect( ekin = 1.0, wl = 1.0 )
    with ensure_error(NCBadInput,'ncrystal_filtertable: no options are supported'
                      ' (got "foo=1")'):
        NCrystalFilter( 'stdlib::Al_sg225.ncmat', 'foo=1' )
    with ensure_error(NCBadInput,'Filter table: only isotropic materials are'
                      ' supported (the cfg-string must not specify a crystal'
                      ' orientation)'):
        NCrystalFilter( 'stdlib::Al_sg225.ncmat;mos=0.3deg'
                        ';dir1=@crys_hkl:0,0,1@lab:0,0,1'
                        ';dir2=@crys_hkl:0,1,0@lab:0,1,0' )

if __name__ == '__main__':
    main()
