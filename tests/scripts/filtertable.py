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

# Test the "filtertable" query, which provides piecewise linear tables of total
# cross sections vs. wavelength.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query as ncquery
from NCrystalDev.exceptions import NCBadInput, NCCalcError
from NCTestUtils.common import ensure_error
import numpy as np

def exact_xs( cfg, wl ):
    return ( NC.createScatter(cfg).xsect(wl=wl)
             + NC.createAbsorption(cfg).xsect(wl=wl) )

def bragg_edges( cfg ):
    res = []
    def collect( info ):
        if info.isMultiPhase():
            for _,ph in info.phases:
                collect(ph)
        elif info.hasHKLInfo():
            res.extend( 2.0*h.d for h in info.hklObjects() )
    collect( NC.createInfo(cfg) )
    return np.unique(res)

def table_eval( wl_t, xs_t, wl ):
    #Linear interpolation, where a pair of identical wavelengths (a Bragg edge)
    #uses the second value (above the edge) for wl >= edge:
    i = np.clip( np.searchsorted(wl_t, wl, side='right') - 1, 0, len(wl_t)-2 )
    x0, x1, y0, y1 = wl_t[i], wl_t[i+1], xs_t[i], xs_t[i+1]
    t = np.where( x1 > x0, (wl-x0)/np.where(x1>x0,x1-x0,1.0), 0.0 )
    return y0 + t*(y1-y0)

def test( cfg, *opts ):
    r = ncquery( ['filtertable',cfg] + list(opts) )
    wl, xs = np.asarray(r['wl']), np.asarray(r['xs'])
    tol, wlmax = r['tol'], r['wlmax']
    assert len(wl) == r['npts'] == len(xs) and len(wl) >= 2
    assert wl[0] == 0.0 and wl[-1] == wlmax
    assert np.all( np.diff(wl) >= 0.0 )
    assert np.all( xs >= 0.0 )
    #Identical wavelengths only come in pairs, at discontinuities (for the
    #supported processes, these are the Bragg edges):
    idup = np.where( np.diff(wl) == 0.0 )[0]
    assert np.all( np.diff(idup) > 1 )
    edges = bragg_edges( cfg )
    for e in wl[idup]:
        assert np.abs(edges-e).min() <= 1e-9 * e
    #The first point is the limit for wavelength -> 0 (which might be ~0):
    assert np.isclose( xs[0], exact_xs( cfg, np.asarray([1e-14]) )[0],
                       rtol = 1e-6, atol = 1e-12 )
    #The other table points (except at edges, which have the values just
    #below or above the edge) are exact cross sections:
    if len(edges):
        j = np.clip( np.searchsorted( edges, wl ), 1, len(edges)-1 )
        dist = np.minimum( np.abs(wl-edges[j-1]), np.abs(wl-edges[j]) )
        iexact = np.where( dist > 1e-9 * wl )[0]
    else:
        iexact = np.arange(len(wl))
    iexact = iexact[ iexact > 0 ]
    xsx = exact_xs( cfg, wl[iexact] )
    assert np.allclose( xs[iexact], xsx, rtol = 1e-12, atol = 0.0 )
    #Interpolated values are within the tolerance, at random wavelengths
    #(log-uniform, and uniform at the shortest wavelengths) and close to the
    #Bragg edges:
    rng = np.random.default_rng( 12345 )
    wlt = np.concatenate( [ np.exp( rng.uniform( np.log(1e-9), np.log(wlmax), 20000 ) ),
                            rng.uniform( 0.0, 1e-5, 2000 ),
                            rng.uniform( 0.0, 1e-7, 2000 ),
                            rng.uniform( 0.0, 1e-9, 2000 ) ] )
    e = edges[ (edges > 1e-6) & (edges < wlmax/1.001) ]
    wlt = np.concatenate( [wlt] + [ e*f for f in (1-1e-8,1+1e-8,1-1e-5,1+1e-5) ] )
    xst = exact_xs( cfg, wlt )
    xsi = table_eval( wl, xs, wlt )
    #The tolerance is relative to max(xs,1e-12 barn):
    maxrelerr = ( np.abs( xsi-xst ) / np.maximum( xst, 1e-12 ) ).max()
    assert maxrelerr <= tol, f'max relative error {maxrelerr} > tol={tol}'
    print(f'{cfg} {" ".join(opts)}'.strip())
    print(f'   npts={r["npts"]} ndiscontinuities={r["ndiscontinuities"]}'
          f' ({len(idup)} in table) tol={tol:g} wlmax={wlmax:g} algo={r["algo"]}'
          f' maxrelerr<tol: OK')
    print(f'   processes: {" ".join(r["processes"])}')
    print(f'   first: ({wl[0]:.6g} Aa, {xs[0]:.6g} b)'
          f' last: ({wl[-1]:.6g} Aa, {xs[-1]:.6g} b)')

def main():
    test('stdlib::Al_sg225.ncmat')
    test('stdlib::Al_sg225.ncmat','tol=1e-2')
    test('stdlib::Al_sg225.ncmat','tol=1e-4')
    test('stdlib::Al_sg225.ncmat','algo=dp')
    test('stdlib::Al_sg225.ncmat','wlmax=10')
    test('stdlib::Al_sg225.ncmat','wlmax=4.7')
    test('stdlib::Al_sg225.ncmat;density=0.5x')
    test('stdlib::Al_sg225.ncmat;inelas=0')
    test('stdlib::Be_sg194.ncmat;temp=80K')
    test('stdlib::Polyethylene_CH2.ncmat')
    test('stdlib::Y2SiO5_sg15_YSO.ncmat')
    test('stdlib::B4C_sg166_BoronCarbide.ncmat')
    test('gasmix::air')
    test('phases<0.3*stdlib::Al_sg225.ncmat&0.7*stdlib::Cu_sg225.ncmat>')
    test('stdlib::void.ncmat')

    base = ['filtertable','stdlib::Al_sg225.ncmat']
    oriented = ('stdlib::Al_sg225.ncmat;mos=0.3deg;dir1=@crys_hkl:0,0,1@lab:0,0,1'
                ';dir2=@crys_hkl:0,1,0@lab:0,1,0')
    with ensure_error(NCBadInput,'Filter table: only isotropic materials are'
                      ' supported (the cfg-string must not specify a crystal'
                      ' orientation)'):
        ncquery(['filtertable',oriented])
    #Processes which have not been validated for filter tables give a warning:
    test('stdlib::Polyethylene_CH2.ncmat;ucnmode=remove')
    #Discontinuities which are not known in advance are detected (here by
    #disabling the known ones):
    with ensure_error(NCCalcError,'Filter table: the cross section has an'
                      ' unexpected discontinuity (or an extremely sharp'
                      ' feature) at wavelength 0.707627 Aa'):
        ncquery(base+['discontinuities=0'])
    with ensure_error(NCBadInput,'Filter table: invalid tolerance 0 (must be'
                      ' in the range (0,1))'):
        ncquery(base+['tol=0'])
    with ensure_error(NCBadInput,'Filter table: invalid wlmax 1e-09 (must be'
                      ' larger than 1e-08 Aa)'):
        ncquery(base+['wlmax=1e-9'])
    with ensure_error(NCBadInput,'Invalid filtertable query: ["filtertable",'
                      '"stdlib::Al_sg225.ncmat","wlmin=1"] (unknown option wlmin)'):
        ncquery(base+['wlmin=1'])
    with ensure_error(NCBadInput,'Invalid filtertable query: ["filtertable",'
                      '"stdlib::Al_sg225.ncmat","tol"] (options must have the'
                      ' form KEY=VALUE)'):
        ncquery(base+['tol'])

if __name__ == '__main__':
    main()
