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

# Tests related to the generation of lists of HKL planes.

import math
import re
import time

import NCrystalDev as NC
import NCTestUtils.enable_fpe  # noqa: F401
import numpy as np
from NCTestUtils.env import ncsetenv


# Without a space group, HKL planes are grouped into families by their
# (d-spacing, |F|^2) values. Check that this never splits the true symmetry
# families (i.e. those found when the space group is known).
def test_nosgfamilies():
    def families( ncmat, cfg ):
        info = NC.directLoad( ncmat, cfg + ';comp=bragg', doScatter = False,
                              doAbsorption = False ).info
        res = []
        for e in info.hklObjects():
            hkl = set( zip( e.h.tolist(), e.k.tolist(), e.l.tolist() ) )
            hkl |= { (-a,-b,-c) for a,b,c in hkl }
            assert len(hkl) == e.mult
            res.append( ( float(e.d), hkl ) )
        return res

    def check( fn, cfg ):
        txt = NC.createTextData( f'stdlib::{fn}' ).rawData
        nosg_txt = re.sub( r'@SPACEGROUP\s*\n\s*\d+\s*\n', '', txt )
        assert nosg_txt != txt
        sgfams = families( txt, cfg )
        nosgfams = families( nosg_txt, cfg )
        hkl2fam = {}
        for i, ( _, hkl ) in enumerate( nosgfams ):
            for e in hkl:
                hkl2fam[e] = i
        for d, hkl in sgfams:
            ids = { hkl2fam.get( e ) for e in hkl }
            assert len(ids) == 1, f'{fn}: symmetry family at d={d} split'
        print(f'{fn:<40} [{cfg}]: {len(sgfams)} families with space group,'
              f' {len(nosgfams)} without')

    #Materials where some families used to be split (depending on rounding of the
    #d-spacing):
    for fn, cfg in [ ( 'GaSe_sg194_GalliumSelenide.ncmat', 'dcutoff=0.4' ),
                     ( 'AlN_sg186_AluminumNitride.ncmat', 'dcutoff=0.15' ),
                     ( 'BeF2_sg152_Beryllium_Fluoride.ncmat', 'dcutoff=0.3' ),
                     ( 'BaO_sg225_BariumOxide.ncmat', 'dcutoff=0.35' ),
                     ( 'NaF_sg225_SodiumFlouride.ncmat', 'dcutoff=0.19' ),
                     ( 'Al_sg225.ncmat', 'dcutoff=0.3' ) ]:
        check( fn, cfg )

# Tests selection of a single HKL plane with the FILLHKL_SELECTHKL environment
# variable (for validation plots), with and without space group. The selected
# plane must get the values of its family in the full list of planes.
def test_select():
    def strip_sg( txt ):
        return re.sub( r'@SPACEGROUP\s*\n\s*\d+\s*\n', '', txt )

    def planes( txt, cfg ):
        NC.clearCaches()
        info = NC.directLoad( txt, cfg + ';comp=bragg', doScatter = False,
                              doAbsorption = False ).info
        return [ ( [ (int(a),int(b),int(c)) for a,b,c in zip(e.h,e.k,e.l) ],
                   e.mult, float(e.d), float(e.f2) ) for e in info.hklObjects() ]

    def select( txt, cfg, hkl ):
        ncsetenv( 'FILLHKL_SELECTHKL', hkl )
        try:
            return planes( txt, cfg )
        finally:
            ncsetenv( 'FILLHKL_SELECTHKL', None )

    def errmsg( e ):
        #Remove the build dependent env var prefix:
        return re.sub( r'NCRYSTAL[A-Z]*_FILLHKL', 'NCRYSTAL_FILLHKL', str(e) )

    def test( fn, cfg, hkls ):
        txt = NC.createTextData( f'stdlib::{fn}' ).rawData
        for label, t in ( ( 'SG', txt ), ( 'no SG', strip_sg( txt ) ) ):
            full = planes( t, cfg )
            for hkl in hkls:
                sel = tuple( int(x) for x in hkl.split(',') )
                msel = tuple( -x for x in sel )
                fam = [ p for p in full if sel in p[0] or msel in p[0] ]
                assert len(fam) <= 1
                try:
                    res = select( t, cfg, hkl )
                except NC.NCCalcError as e:
                    assert not fam
                    assert 'is not present' in str(e)
                    print(f'  {fn} [{label}] {hkl}: CalcError ({errmsg(e)[:40]}...)')
                    continue
                if not res:
                    #Plane outside the d-spacing range:
                    assert not fam
                    info = NC.directLoad( t, cfg + ';comp=bragg',
                                          doScatter = False,
                                          doAbsorption = False ).info
                    d = info.dspacingFromHKL( *sel )
                    assert not ( info.hklDLower() <= d <= info.hklDUpper() )
                    print(f'  {fn} [{label}] {hkl}: outside d-spacing range')
                    continue
                assert fam, 'selected plane not in full list'
                assert len(res) == 1
                hkllist, mult, d, f2 = res[0]
                assert hkllist == [ sel ] and mult == 2
                _, _, dref, f2ref = fam[0]
                if label == 'SG':
                    #Same code path and representative, so identical values:
                    assert d == dref and f2 == f2ref
                else:
                    #Family averages over (possibly) fewer members:
                    assert abs( d - dref ) < 1e-12 * dref
                    assert abs( f2 - f2ref ) < 1e-6 * f2ref
                print(f'  {fn} [{label}] {hkl}: d={d:.6g} F2={f2:.6g}')

    test( 'Al_sg225.ncmat', 'dcutoff=0.3',
          [ '1,1,1', '-1,-1,-1', '2,0,0', '0,0,-2', '0,-2,0', '-3,1,1',
            '1,-3,1', '1,1,0', '5,5,5', '10,10,10' ] )
    test( 'GaSe_sg194_GalliumSelenide.ncmat', 'dcutoff=0.4;dcutoffup=6',
          [ '1,0,0', '0,1,0', '-1,1,0', '1,-2,0', '2,-1,3', '0,0,2', '0,0,1',
            '0,0,4', '0,0,-4' ] )
    test( 'Al2O3_sg167_Corundum.ncmat', 'dcutoff=0.5',
          [ '1,0,4', '0,1,-4', '-1,0,-4', '1,1,3', '1,1,-3', '0,0,6' ] )

    #Invalid values:
    t = NC.createTextData( 'stdlib::Al_sg225.ncmat' ).rawData
    for bad in ( '1,1', '1,1,1,1', 'a,b,c', '1.5,1,1', '0,0,0', ' ' ):
        for label, tt in ( ( 'SG', t ), ( 'no SG', strip_sg( t ) ) ):
            try:
                select( tt, 'dcutoff=0.3', bad )
            except NC.NCBadInput as e:
                assert 'FILLHKL_SELECTHKL' in str(e)
                print(f'  invalid "{bad}" [{label}]: BadInput')
            else:
                raise RuntimeError(f'invalid value "{bad}" accepted')

# Tests that physics models work with a single plane selected via the
# FILLHKL_SELECTHKL environment variable (in case it is ever exposed as a cfg
# parameter): single crystals and layered crystals must scatter exactly on the
# selected (h,k,l) and (-h,-k,-l) planes, and powders must scale as d*|F|^2.
def test_selectphys():
    def strip_sg( txt ):
        return re.sub( r'@SPACEGROUP\s*\n\s*\d+\s*\n', '', txt )

    def load( txt, cfg, sel ):
        ncsetenv( 'FILLHKL_SELECTHKL', sel )
        try:
            NC.clearCaches()
            return NC.directLoad( txt, 'comp=bragg;' + cfg, doAbsorption = False )
        finally:
            ncsetenv( 'FILLHKL_SELECTHKL', None )

    def norm( v ):
        n = math.sqrt( sum( x*x for x in v ) )
        return tuple( x / n for x in v )

    def bragg_direction( g, d, wl ):
        """Incident direction satisfying the Bragg condition for the plane with
    normal g (in the lab frame) and d-spacing d."""
        g = norm( g )
        a = ( 0.0, 1.0, 0.0 ) if abs( g[0] ) > 0.9 else ( 1.0, 0.0, 0.0 )
        perp = norm( ( g[1]*a[2] - g[2]*a[1], g[2]*a[0] - g[0]*a[2],
                       g[0]*a[1] - g[1]*a[0] ) )
        s = wl / ( 2 * d )
        c = math.sqrt( 1 - s*s )
        return tuple( -s*a + c*b for a, b in zip( g, perp ) )

    def check_single_crystal( fn, cfg, gdir_lab, wl, selections, layered = False ):
        """For each selected plane, gdir_lab(hkl) gives its normal in the lab
    frame. Scattering at the exact Bragg geometry must have momentum transfer
    along the normal, with |q|=2pi/d. For single crystals, the planes of the
    full list must not add anything (they are not in the Bragg condition), but
    for layered crystals they might (due to rotation around the c-axis)."""
        txt = NC.createTextData( f'stdlib::{fn}' ).rawData
        full = NC.directLoad( txt, 'comp=bragg;' + cfg, doAbsorption = False )
        ekin = NC.wl2ekin( wl )
        k = 2 * math.pi / wl
        for sel, other in selections:
            hkl = tuple( int(x) for x in sel.split(',') )
            xs = []
            for label, t in ( ( 'SG', txt ), ( 'no SG', strip_sg( txt ) ) ):
                m = load( t, cfg, sel )
                assert m.info.nHKL() == 1
                e = next( m.info.hklObjects() )
                assert e.mult == 2 and ( int(e.h[0]), int(e.k[0]), int(e.l[0]) ) == hkl
                d = float( e.d )
                g = norm( gdir_lab( hkl ) )
                u = bragg_direction( g, d, wl )
                x = float( m.scatter.xsect( wl = wl, direction = u ) )
                assert x > 0.0
                xs.append( x )
                xfull = float( full.scatter.xsect( wl = wl, direction = u ) )
                if layered:
                    assert x <= xfull * ( 1 + 1e-9 )
                else:
                    assert abs( x - xfull ) < 1e-6 * x
                #No scattering at the Bragg geometry of another plane of the family:
                ho = tuple( int(v) for v in other.split(',') )
                uo = bragg_direction( norm( gdir_lab( ho ) ), d, wl )
                assert float( m.scatter.xsect( wl = wl, direction = uo ) ) == 0.0
                #Scattered neutrons have q along +-g, and |q|=2pi/d:
                ef, uf = m.scatter.sampleScatter( ekin, u, repeat = 100 )
                for ef_i, uf_i in zip( ef, zip( *uf ) ):
                    assert abs( ef_i - ekin ) < 1e-12 * ekin
                    q = tuple( k * ( a - b ) for a, b in zip( uf_i, u ) )
                    qmag = math.sqrt( sum( x*x for x in q ) )
                    assert abs( qmag - 2 * math.pi / d ) < 1e-6 * qmag
                    cosang = abs( sum( a*b for a, b in zip( q, g ) ) ) / qmag
                    assert cosang > math.cos( math.radians( 5.0 ) )
            assert abs( xs[0] - xs[1] ) < 1e-9 * xs[0]
            print(f'  {fn} {sel}: d={d:.5g}Aa, xsect at Bragg geometry {xs[0]:.4g} barn,'
                  f' zero for ({other}) (same without space group)')

    #Al (cubic), crystal axes aligned with the lab axes:
    check_single_crystal( 'Al_sg225.ncmat',
                          'mos=0.3deg;dir1=@crys_hkl:0,0,1@lab:0,0,1;'
                          'dir2=@crys_hkl:1,0,0@lab:1,0,0',
                          lambda hkl : hkl, 2.0,
                          [ ( '1,1,1', '1,1,-1' ), ( '0,0,2', '2,0,0' ),
                            ( '-3,1,1', '1,3,1' ) ] )

    #Pyrolytic graphite (layered crystal model), (0,0,l) planes along lab z:
    check_single_crystal( 'C_sg194_pyrolytic_graphite.ncmat',
                          'mos=2deg;lcaxis=0,0,1;dir1=@crys_hkl:0,0,1@lab:0,0,1;'
                          'dir2=@crys_hkl:1,0,0@lab:1,0,0',
                          lambda hkl : hkl, 3.0,
                          [ ( '0,0,2', '1,0,0' ), ( '0,0,-4', '1,0,0' ) ],
                          layered = True )

    #Powder: cross sections below the Bragg thresholds of two selected planes
    #must be proportional to d*|F|^2 (multiplicity 2 for both):
    al = NC.createTextData( 'stdlib::Al_sg225.ncmat' ).rawData
    for t in ( al, strip_sg( al ) ):
        res = []
        for sel in ( '1,1,1', '2,0,0', '-3,1,1' ):
            m = load( t, '', sel )
            e = next( m.info.hklObjects() )
            assert m.info.braggthreshold == 2 * float( e.d )
            assert float( m.scatter.xsect( wl = 2.01 * float( e.d ) ) ) == 0.0
            res.append( ( float( m.scatter.xsect( wl = 1.0 ) ),
                          float( e.d ) * float( e.f2 ) ) )
        for x, df2 in res[1:]:
            assert abs( ( x / res[0][0] ) / ( df2 / res[0][1] ) - 1 ) < 1e-9
    print('  Powder cross sections of selected planes scale as d*|F|^2')

# Tests the bookkeeping table used when generating HKL planes with space group
# symmetry: its memory limit (FILLHKL_MEMLIM env var), and that results are
# correct with large hkl indices (compared with results without space group).
def test_memlim():
    def strip_sg( txt ):
        return re.sub( r'@SPACEGROUP\s*\n\s*\d+\s*\n', '', txt )

    def load( txt, cfg ):
        NC.clearCaches()
        return NC.directLoad( txt, cfg + ';comp=bragg', doScatter = False,
                              doAbsorption = False ).info

    def families( info ):
        res = []
        for e in info.hklObjects():
            s = set( zip( e.h.tolist(), e.k.tolist(), e.l.tolist() ) )
            s |= { (-a,-b,-c) for a,b,c in s }
            assert len(s) == e.mult
            res.append( ( float(e.d), float(e.f2), s ) )
        return res

    def errmsg( e ):
        #Remove the build dependent env var prefix:
        return re.sub( r'NCRYSTAL[A-Z]*_FILLHKL', 'NCRYSTAL_FILLHKL', str(e) )

    def compare_sg_nosg( label, txt, cfg ):
        """Symmetry families must be exactly those found without space group,
    except that the latter might merge several of them."""
        sgf = families( load( txt, cfg ) )
        nosgf = families( load( strip_sg( txt ), cfg ) )
        h2i = {}
        for i, ( _, _, s ) in enumerate( nosgf ):
            for hkl in s:
                h2i[hkl] = i
        #F2 sums (weighted with multiplicity) must match in each nosg family:
        f2sum = [ 0.0 ] * len(nosgf)
        npts = 0
        maxhkl = 0
        for d, f2, s in sgf:
            ids = { h2i.get( hkl ) for hkl in s }
            assert len(ids) == 1 and None not in ids
            i = ids.pop()
            assert abs( nosgf[i][0] - d ) < 1e-6 * d
            f2sum[i] += f2 * len(s)
            npts += len(s)
            maxhkl = max( maxhkl, max( max(abs(x) for x in hkl) for hkl in s ) )
        assert npts == len(h2i)
        for ( _, f2, s ), f2s in zip( nosgf, f2sum ):
            assert abs( f2 * len(s) - f2s ) < 1e-6 * f2s
        print(f'{label}: {len(sgf)} symmetry families (max |h|,|k|,|l| = {maxhkl}),'
              f' {len(nosgf)} families without space group, all consistent')

    #Large c axis, so indices >128 (which used to require a 67MB table):
    large_cell = '''NCMAT v7
@CELL
  lengths 3 3 40
  angles 90 90 90
@SPACEGROUP
  123
@ATOMPOSITIONS
  Al 0 0 0
  O 1/2 1/2 1/2
@DYNINFO
  element Al
  fraction 1/2
  type vdosdebye
  debye_temp 300
@DYNINFO
  element O
  fraction 1/2
  type vdosdebye
  debye_temp 300
'''
    compare_sg_nosg( 'Large tetragonal cell', large_cell, 'dcutoff=0.25' )
    for fn, cfg in [ ( 'BO3H3_sg2_BoricAcid.ncmat', 'dcutoff=0.5' ),
                     ( 'GaSe_sg194_GalliumSelenide.ncmat', 'dcutoff=0.3' ),
                     ( 'Al2O3_sg167_Corundum.ncmat', 'dcutoff=0.3' ) ]:
        compare_sg_nosg( fn, NC.createTextData( f'stdlib::{fn}' ).rawData, cfg )

    #Memory limit (default 20MB), which must trigger before any time is spent:
    yag = NC.createTextData( 'stdlib::Y3Al5O12_sg230_YAG.ncmat' ).rawData
    al = NC.createTextData( 'stdlib::Al_sg225.ncmat' ).rawData
    def expect_error( errtype, txt, cfg, memlim = None ):
        ncsetenv( 'FILLHKL_MEMLIM', memlim )
        t0 = time.time()
        try:
            load( txt, cfg ).nHKL()
        except errtype as e:
            assert time.time() - t0 < 10.0
            assert 'NCRYSTAL_FILLHKL_MEMLIM' in errmsg(e)
            print(f'  {cfg} (memlim={memlim}): {errtype.__name__}:'
                  f' {errmsg(e)[:95]}...')
        else:
            raise RuntimeError('expected error did not happen')
        finally:
            ncsetenv( 'FILLHKL_MEMLIM', None )

    expect_error( NC.NCCalcError, yag, 'dcutoff=0.02' )
    expect_error( NC.NCCalcError, al, 'dcutoff=0.1', memlim = '0.01' )
    for bad in ( '0', '-1', 'abc' ):
        expect_error( NC.NCBadInput, al, 'dcutoff=0.1', memlim = bad )

    #Works when the limit is raised, and does not affect crystals without space
    #group:
    ncsetenv( 'FILLHKL_MEMLIM', '1' )
    print('  Al (memlim=1):', load( al, 'dcutoff=0.1' ).nHKL(), 'families')
    ncsetenv( 'FILLHKL_MEMLIM', '0.01' )
    print('  Al without space group (memlim=0.01):',
          load( strip_sg( al ), 'dcutoff=0.1' ).nHKL(), 'families')
    ncsetenv( 'FILLHKL_MEMLIM', None )

# Checks |F|^2 of all planes against an independent calculation, for materials
# with strongly damped light atoms at small d-spacings. This includes weak
# planes near fsquarecut, where contributions from heavily damped atoms used
# to be skipped (giving errors of several percent).
def test_weakf2():
    def check( fn, cfg ):
        info = NC.createInfo( f'stdlib::{fn};{cfg}' )
        hkl, d, f2 = [], [], []
        for e in info.hklObjects():
            #All planes of the family, with the family values:
            for h, k, l in zip( e.h, e.k, e.l ):  # noqa: E741
                hkl.append( ( h, k, l ) )
                d.append( float( e.d ) )
                f2.append( float( e.f2 ) )
        hkl = np.asarray( hkl, dtype = float )
        d, f2 = np.asarray( d ), np.asarray( f2 )
        q2 = ( 2 * math.pi / d )**2
        F = np.zeros( len(d), dtype = complex )
        for ai in info.atominfos:
            b = ai.atomData.coherentScatLen()
            pos = np.asarray( ai.positions, dtype = float )
            phases = np.exp( 2j * math.pi * ( hkl @ pos.T ) ).sum( axis = 1 )
            F += b * np.exp( -0.5 * ai.msd * q2 ) * phases
        f2ref = np.abs( F )**2
        relerr = np.abs( f2 - f2ref ) / f2ref
        nweak = int( ( f2 < 1e-4 ).sum() )
        assert relerr.max() < 1e-9, ( fn, float( relerr.max() ) )
        print(f'{fn:<36} {cfg}: {len(d)} planes ({nweak} with F2<1e-4 barn)'
              f' agree with independent F2 calculation')

    check( 'BO3H3_sg2_BoricAcid.ncmat', 'dcutoff=0.25' )
    check( 'MgH2_sg136_MagnesiumHydride.ncmat', 'dcutoff=0.2' )
    check( 'SiO2-beta_sg180_BetaQuartz.ncmat', 'dcutoff=0.22' )
    check( 'CaH2_sg62_CalciumHydride.ncmat', 'dcutoff=0.2' )

def main():
    test_nosgfamilies()
    test_select()
    test_selectphys()
    test_memlim()
    test_weakf2()

if __name__ == '__main__':
    main()
