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

# NEEDS: numpy spglib

# Validates the handling of unit cells (volume, density, d-spacings, structure
# factors, symmetry and plane directions) for many synthetic cells of all
# shapes, against independent metric-tensor calculations and against spglib:
# the same crystal in different cell settings must give identical physics.
# Also cross-checks all crystalline materials in the standard data library.
#
# Background: github issue #383 (general cells had a wrong real-space lattice
# matrix, giving wrong volumes and d-spacings).

import cmath
import math
import random

import NCrystalDev as NC
import NCTestUtils.enable_fpe  # noqa: F401
import numpy as np
import spglib
from NCTestUtils.env import ncsetenv

#Keep planes with tiny F2, so plane lists can be compared exactly (the
#variable is read when planes are calculated):
ncsetenv('FILLHKL_IGNOREFSQCUT','1')

_msd = { 'Al': 0.009, 'O': 0.013, 'Fe': 0.005, 'V': 0.007, 'Ni': 0.006 }

def metric_tensor( a, b, c, alpha, beta, gamma ):
    ca, cb, cg = ( math.cos(math.radians(x)) for x in (alpha,beta,gamma) )
    return np.array( [ [ a*a, a*b*cg, a*c*cb ],
                       [ a*b*cg, b*b, b*c*ca ],
                       [ a*c*cb, b*c*ca, c*c ] ] )

def cellpars_of( si ):
    return tuple( si[k] for k in ('a','b','c','alpha','beta','gamma') )

def cellpars_from_lattice( lattice ):
    #Lattice vectors are the rows (spglib convention). Results are rounded, to
    #snap e.g. gamma=120.00000000000001 to 120:
    va, vb, vc = ( np.asarray(v,dtype=float) for v in lattice )
    def ang( u, v ):
        cosang = float(u@v) / math.sqrt( float((u@u)*(v@v)) )
        return round( math.degrees( math.acos( cosang ) ), 10 )
    def length( v ):
        return round( float(np.linalg.norm(v)), 12 )
    return ( length(va), length(vb), length(vc),
             ang(vb,vc), ang(va,vc), ang(va,vb) )

def lattice_from_cellpars( a, b, c, alpha, beta, gamma ):
    ca, cb, cg = ( math.cos(math.radians(x)) for x in (alpha,beta,gamma) )
    sg = math.sin(math.radians(gamma))
    cy = (ca-cb*cg)/sg
    return np.array( [ [ a, 0.0, 0.0 ],
                       [ b*cg, b*sg, 0.0 ],
                       [ c*cb, c*cy, c*math.sqrt(1.0-cb*cb-cy*cy) ] ] )

def compose( cellpars, atoms, spacegroup = None ):
    #Atoms are (element,x,y,z). NCMAT data is only verified with spglib by the
    #NCMATComposer when a spacegroup is provided.
    a, b, c, alpha, beta, gamma = cellpars
    comp = NC.NCMATComposer()
    comp.set_cellsg( a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma,
                     spacegroup = spacegroup )
    comp.set_atompos( atoms )
    for e in sorted({ e for e,*_ in atoms }):
        comp.set_dyninfo_msd( e, _msd[e], temperature = 293.15 )
    return comp.create_ncmat(
        verify_crystal_structure = spacegroup is not None )

def load_info( ncmat, dcut ):
    return NC.directLoad( ncmat, f'dcutoff={dcut}',
                          doScatter = False, doAbsorption = False ).info

def cohelas_xs_per_atom( ncmat, dcut, wls ):
    sc = NC.directLoad( ncmat, f'dcutoff={dcut};comp=coh_elas',
                        doInfo = False, doAbsorption = False ).scatter
    return [ sc.xsect( wl = wl ) for wl in wls ]

def all_hkl( G, dcut ):
    #All (h,k,l) with d>=dcut (one of each Friedel pair), with d values.
    Gi = np.linalg.inv( G )
    #|h| <= |a|/d, since a.G_hkl = 2pi*h:
    hmax = [ math.floor( math.sqrt(G[i,i]) / dcut ) + 1
             for i in range(3) ]
    for h in range(-hmax[0],hmax[0]+1):
        for k in range(-hmax[1],hmax[1]+1):
            for ll in range(-hmax[2],hmax[2]+1):
                if (h,k,ll) <= (0,0,0):
                    continue
                v = np.array( (h,k,ll), dtype = float )
                d = 1.0 / math.sqrt( v@Gi@v )
                if d >= dcut:
                    yield (h,k,ll), d

def ref_planes( info, dcut ):
    #Independent calculation of d-spacings and structure factors (in barn,
    #including Debye-Waller factors) of all planes with d>=dcut.
    G = metric_tensor( *cellpars_of( info.structure_info ) )
    atoms = [ ( ai.atomData.coherentScatLen(), ai.msd, p )
              for ai in info.atominfos for p in ai.positions ]
    out = {}
    for hkl, d in all_hkl( G, dcut ):
        q2 = ( 2*math.pi/d )**2
        F = sum( b * math.exp(-0.5*msd*q2)
                 * cmath.exp( 2j*math.pi*float(np.dot(hkl,p)) )
                 for b, msd, p in atoms )
        out[hkl] = ( d, abs(F)**2 )
    return out

def safe_dcutoff( cellpars, nplanes ):
    #A dcutoff giving roughly nplanes planes, placed in a gap in the list of
    #d-spacings (so rounding at the boundary can never affect plane counts).
    G = metric_tensor( *cellpars )
    V = math.sqrt( np.linalg.det(G) )
    dc = ( 4*math.pi*V/(6*nplanes) )**(1/3.)
    ds = sorted( d for _,d in all_hkl( G, 0.8*dc ) )
    for d1, d2 in zip( ds, ds[1:] ):
        if d1 >= dc and d2 > d1*(1+1e-6):
            return 0.5*(d1+d2)
    raise RuntimeError('no gap found')

def check_cell( name, info, dcut ):
    #Checks volume, density, and every individual plane against independent
    #calculations. Returns number of planes.
    si = info.structure_info
    V = math.sqrt( np.linalg.det( metric_tensor( *cellpars_of(si) ) ) )
    assert abs( si['volume']/V - 1.0 ) < 1e-12, (name,si['volume'],V)
    natoms = sum( len(ai.positions) for ai in info.atominfos )
    assert natoms == si['n_atoms']
    mass = sum( ai.atomData.averageMassAMU()*len(ai.positions)
                for ai in info.atominfos )
    amu_per_aa3_in_gcm3 = 1.66053906660
    assert abs( info.density/(mass*amu_per_aa3_in_gcm3/V) - 1.0 ) < 1e-6
    assert abs( info.numberdensity*V/natoms - 1.0 ) < 1e-12
    ref = ref_planes( info, dcut )
    f2max = max( f2 for _,f2 in ref.values() )
    symeqv = info.hklIsSymEqvGroup()
    def f2eq( x, y ):
        return abs( x - y ) < 1e-9*max(1e-3*f2max,y)
    def merge_compatible( x, y ):
        #NCrystal's criterion for grouping planes into one family, when not
        #using symmetry (relative tolerance 1e-6):
        return abs( x - y ) < 1.01e-6*( x + y )
    seen = set()
    for e in info.hklObjects():
        assert 2*len(e.h) == e.mult
        nexact = 0
        for hkl in zip( e.h, e.k, e.l ):
            hkl = tuple( int(x) for x in hkl )
            key = hkl if hkl > (0,0,0) else tuple( -x for x in hkl )
            assert key in ref and key not in seen, (name,hkl)
            seen.add( key )
            d_ref, f2_ref = ref[key]
            #Each individual plane, via the C-API:
            assert abs( info.dspacingFromHKL(*hkl)/d_ref - 1.0 ) < 1e-11
            if abs( e.d/d_ref - 1.0 ) < 1e-11 and f2eq( e.f2, f2_ref ):
                nexact += 1
            elif symeqv:
                raise RuntimeError(f'{name}: wrong d or F2 for {hkl}: {e.d},'
                                   f' {e.f2} (expected {d_ref}, {f2_ref})')
            else:
                assert merge_compatible( e.d, d_ref ), (name,hkl)
                assert merge_compatible( e.f2, f2_ref ), (name,hkl)
        assert nexact >= 1, (name,str(e))
    #Only planes with vanishing F2 (systematic absences) may be absent:
    missing = [ k for k in set(ref)-seen if ref[k][1] > 1e-12*f2max ]
    assert not missing, (name,sorted(missing)[:5])
    return len(seen)

def check_symmetry_groups( name, info, dataset ):
    #Groups of symmetry-equivalent planes must be exactly the orbits of the
    #spglib point group operations (plus Friedel pairs).
    assert info.hklIsSymEqvGroup()
    rots = [ np.asarray( r, dtype = int ) for r in dataset.rotations ]
    for e in info.hklObjects():
        group = set()
        for hkl in zip( e.h, e.k, e.l ):
            hkl = tuple( int(x) for x in hkl )
            group.update( ( hkl, tuple( -x for x in hkl ) ) )
        orbit = set()
        for r in rots:
            hr = tuple( int(x) for x in np.array( e.hkl_label ) @ r )
            orbit.update( ( hr, tuple( -x for x in hr ) ) )
        assert group == orbit, (name,e.hkl_label)

def check_orientation( name, info, ncmat, dcut ):
    #Plane directions: a single crystal orientation defined via two hkl
    #directions must be accepted when the lab directions are exactly the
    #metric-tensor angle apart, and rejected when off by just 1e-6 rad.
    Gi = np.linalg.inv( metric_tensor( *cellpars_of(info.structure_info) ) )
    def invd( h ):
        return math.sqrt( h@Gi@h )
    hkls = sorted( ( np.array(e.hkl_label) for e in info.hklObjects() ),
                   key = invd )
    h1 = hkls[0]
    h2 = next( h for h in hkls[1:] if np.linalg.norm(np.cross(h1,h)) > 0 )
    phi = math.acos( (h1@Gi@h2) / ( invd(h1)*invd(h2) ) )
    def load( phi_lab ):
        NC.directLoad( ncmat,
                       f'dcutoff={dcut};mos=1deg;dirtol=1e-9;'
                       f'dir1=@crys_hkl:{h1[0]},{h1[1]},{h1[2]}@lab:0,0,1;'
                       f'dir2=@crys_hkl:{h2[0]},{h2[1]},{h2[2]}@lab:0,'
                       f'{math.sin(phi_lab):.17g},{math.cos(phi_lab):.17g}',
                       doInfo = False, doAbsorption = False )
    load( phi )
    try:
        load( phi + 1e-6 )
    except NC.NCBadInput:
        pass
    else:
        raise RuntimeError(f'{name}: wrong orientation accepted')

def powder_spectrum( info, dcut ):
    #Invariant (for a given crystal) list of (d, sum(mult*F2)/natoms^2).
    n2 = info.structure_info['n_atoms']**2
    out = []
    for d, w in sorted( ( e.d, e.mult*e.f2 ) for e in info.hklObjects() ):
        if out and d < out[-1][0]*(1+3e-6):#above hkl merge tolerance
            out[-1][1] += w
        else:
            out.append( [d, w] )
    return [ (d,w/n2) for d, w in out ]

def cmp_spectra( name, s1, s2 ):
    #Absent planes (e.g. in centred cells) only appear in some settings:
    wmax = max( w for d,w in s1 )
    s1, s2 = ( [ e for e in s if e[1] > 1e-9*wmax ] for s in (s1,s2) )
    assert len(s1) == len(s2), (name,len(s1),len(s2))
    for (d1,w1),(d2,w2) in zip( s1, s2 ):
        assert abs( d1/d2 - 1.0 ) < 3e-6, (name,d1,d2)
        assert abs( w1 - w2 ) < 1e-5*wmax, (name,d1,w1,w2)

def spglib_cell( cellpars, atoms ):
    elems = sorted( { e for e,*_ in atoms } )
    return ( ( lattice_from_cellpars( *cellpars ),
               [ (x,y,z) for e,x,y,z in atoms ],
               [ elems.index(e)+1 for e,*_ in atoms ] ),
             { i+1 : e for i, e in enumerate(elems) } )

def other_settings( cellpars, atoms ):
    #The same crystal in other settings: spglib standardised and primitive
    #cells, and the Niggli reduced cell.
    cell, num2elem = spglib_cell( cellpars, atoms )
    ds = spglib.get_symmetry_dataset( cell, symprec = 1e-5 )
    pl, pp, pn = spglib.standardize_cell( cell, to_primitive = True,
                                          symprec = 1e-5 )
    nl = spglib.niggli_reduce( cell[0] )
    npos = [ np.array(p) @ cell[0] @ np.linalg.inv(nl) for p in cell[1] ]
    for vname, vl, vp, vn, sg in ( ( 'std', ds.std_lattice, ds.std_positions,
                                     ds.std_types, ds.number ),
                                   ( 'prim', pl, pp, pn, None ),
                                   ( 'niggli', nl, npos, cell[2], None ) ):
        vatoms = [ ( num2elem[int(n)], )+tuple( float(x)%1.0 for x in p )
                   for p, n in zip( vp, vn ) ]
        vds = spglib.get_symmetry_dataset( ( vl, vp, vn ), symprec = 1e-5 )
        yield vname, cellpars_from_lattice( vl ), vatoms, sg, vds

def check_crystal( name, cellpars, atoms, nplanes = 300 ):
    ncmat = compose( cellpars, atoms )
    dcut = safe_dcutoff( cellpars, nplanes )
    info = load_info( ncmat, dcut )
    nchecked = check_cell( name, info, dcut )
    check_orientation( name, info, ncmat, dcut )
    spec = powder_spectrum( info, dcut )
    vpa = info.structure_info['volume'] / info.structure_info['n_atoms']
    dmax = max( e.d for e in info.hklObjects() )
    wls = [ 2*dcut + f*2*(dmax-dcut) for f in (0.01,0.2,0.5,0.8,0.99) ]
    xs = cohelas_xs_per_atom( ncmat, dcut, wls )
    sgno = spglib.get_symmetry_dataset( spglib_cell( cellpars, atoms )[0],
                                        symprec = 1e-5 ).number
    for vname, vcell, vatoms, sg, vds in other_settings( cellpars, atoms ):
        vn = f'{name}/{vname}'
        vncmat = compose( vcell, vatoms, spacegroup = sg )
        vinfo = load_info( vncmat, dcut )
        nchecked += check_cell( vn, vinfo, dcut )
        vsi = vinfo.structure_info
        assert abs( vsi['volume']/vsi['n_atoms']/vpa - 1.0 ) < 1e-9, vn
        assert abs( vinfo.density/info.density - 1.0 ) < 1e-9, vn
        cmp_spectra( vn, spec, powder_spectrum( vinfo, dcut ) )
        for x1, x2 in zip( xs, cohelas_xs_per_atom( vncmat, dcut, wls ) ):
            assert abs( x2/x1 - 1.0 ) < 1e-5, (vn,x1,x2)
        if sg is not None:
            assert vds.number == sg, vn
            check_symmetry_groups( vn, vinfo, vds )
    return sgno, nchecked

def valid_angles( alpha, beta, gamma ):
    ca, cb, cg = ( math.cos(math.radians(x)) for x in (alpha,beta,gamma) )
    #Normalised squared volume (1 for orthogonal axes), avoid flat cells:
    return 1.0 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg > 0.2**2

def gen_cells( rng ):
    def L():
        return round( rng.uniform(2.5,9.0), 5 )
    def A( lo = 60.0, hi = 120.0 ):
        return round( rng.uniform(lo,hi), 4 )
    def angles( lo = 60.0, hi = 120.0, gamma = None ):
        while True:
            res = ( A(lo,hi), A(lo,hi), gamma or A(lo,hi) )
            if valid_angles( *res ):
                return res
    a = L()
    fcc_p, bcc_p = 4.04932/math.sqrt(2), 2.8665*math.sqrt(3)/2
    bcc_ang = math.degrees( math.acos(-1/3) )
    yield 'triclinic', ( L(), L(), L() ) + angles()
    yield 'triclinic', ( L(), L(), L() ) + angles()
    yield 'triclinic-acute', ( L(), L(), L() ) + angles(62,89)
    yield 'triclinic-obtuse', ( L(), L(), L() ) + angles(91,118)
    yield 'monoclinic-a', ( L(), L(), L(), A(95,120), 90, 90 )
    yield 'monoclinic-b', ( L(), L(), L(), 90, A(95,120), 90 )
    yield 'monoclinic-c', ( L(), L(), L(), 90, 90, A(95,118) )
    yield 'orthorhombic', ( L(), L(), L(), 90, 90, 90 )
    yield 'tetragonal', ( a, a, L(), 90, 90, 90 )
    yield 'cubic', ( a, a, a, 90, 90, 90 )
    yield 'hexagonal-120', ( a, a, L(), 90, 90, 120 )
    yield 'hexagonal-60', ( a, a, L(), 90, 90, 60 )
    x = A(45,105)
    yield 'rhombohedral', ( a, a, a, x, x, x )
    yield 'fcc-primitive', ( fcc_p, fcc_p, fcc_p, 60, 60, 60 )
    yield 'bcc-primitive', ( bcc_p, bcc_p, bcc_p, bcc_ang, bcc_ang, bcc_ang )
    yield 'gamma120-alphabeta-random', ( L(), L(), L() ) + angles(gamma=120)
    yield 'alpha-beta-90-gamma-random', ( L(), L(), L(), 90, 90, A() )
    yield 'alpha-gamma-90-beta-random', ( L(), L(), L(), 90, A(), 90 )
    #Angles at or just outside the "exactly 90 or 120" tolerance used in the
    #C++ code (1e-14 rad):
    yield 'near-90-inside-tol', ( L(), L(), L(), 90+1e-13, 90-1e-13, 90 )
    yield 'near-90-outside-tol', ( L(), L(), L(), 90+1e-9, 90-1e-9, 90 )
    yield 'near-120-outside-tol', ( a, a, L(), 90, 90, 120+1e-9 )

def gen_atoms( rng, single ):
    if single:
        return [ ( rng.choice(sorted(_msd)), 0.0, 0.0, 0.0 ) ]
    return [ ( rng.choice(sorted(_msd)), round(rng.random(),4),
               round(rng.random(),4), round(rng.random(),4) )
             for i in range(rng.randint(2,4)) ]

def test_synthetic():
    print('Synthetic cells (crystal, spglib space group, planes checked):')
    rng = random.Random( 123456 )
    ntot = 0
    for name, cellpars in gen_cells( rng ):
        for single in ( True, False ):
            atoms = gen_atoms( rng, single )
            fullname = f'{name}-{len(atoms)}atom{"" if single else "s"}'
            sgno, n = check_crystal( fullname, cellpars, atoms )
            ntot += n
            print(f'  {fullname:38} SG-{sgno:<3} {n}')
    print(f'Synthetic cells OK ({ntot} planes checked)')

def test_stdlib():
    n = 0
    for fe in NC.browseFiles( factory = 'stdlib' ):
        info = NC.createInfo( f'stdlib::{fe.name}' )
        if not info.isSinglePhase() or not info.hasStructureInfo():
            continue
        cellpars = cellpars_of( info.structure_info )
        dcut = safe_dcutoff( cellpars, 150 )
        info = load_info( NC.createTextData(f'stdlib::{fe.name}').rawData,
                          dcut )
        check_cell( fe.name, info, dcut )
        atoms = [ ( ai.atomData.displayLabel(), )+tuple(p)
                  for ai in info.atominfos for p in ai.positions ]
        ds = spglib.get_symmetry_dataset( spglib_cell( cellpars, atoms )[0],
                                          symprec = 1e-5 )
        assert ds.number == info.structure_info['spacegroup'], fe.name
        check_symmetry_groups( fe.name, info, ds )
        n += 1
    print(f'Standard library crystals OK ({n} files)')

def main():
    test_synthetic()
    test_stdlib()

if __name__ == '__main__':
    main()
