
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

"""Utilities for generating synthetic crystal structures and CIF data (both
valid CIF data in various styles, and invalid CIF data which should be
rejected), for testing CIF-related functionality. Requires numpy and spglib.
"""

import math
from fractions import Fraction

import numpy as np
import spglib

_elements = ( 'Al', 'O', 'Fe', 'Ni', 'V', 'Cu', 'Si', 'C' )

class Structure:
    """Crystal structure in the spglib standard setting of a space group (or
    P1), with atoms given as (element,(x,y,z)) tuples for all atoms in the
    cell."""

    def __init__( self, cellpars, atoms, spacegroup = 1, hall_number = None ):
        self.cellpars = tuple( float(e) for e in cellpars )
        self.atoms = [ ( e, tuple( float(x)%1.0 for x in p ) )
                       for e, p in atoms ]
        self.spacegroup = spacegroup
        self.hall_number = hall_number or default_hall_number( spacegroup )

    def lattice( self ):
        return lattice_from_cellpars( *self.cellpars )

    def spglib_cell( self ):
        elems = self.elements()
        return ( self.lattice(), [ p for _,p in self.atoms ],
                 [ elems.index(e)+1 for e,_ in self.atoms ] )

    def elements( self ):
        return sorted( { e for e,_ in self.atoms } )

    def symops( self ):
        return symops_from_hall( self.hall_number )

    def unique_sites( self ):
        """One representative for each set of symmetry-equivalent atoms."""
        ds = spglib.get_symmetry_dataset( self.spglib_cell(), symprec = 1e-5,
                                          hall_number = self.hall_number )
        assert ds.number == self.spacegroup
        return [ self.atoms[i] for i in sorted(set(ds.equivalent_atoms)) ]

    def with_origin_shift( self, shift ):
        """Same crystal, with all coordinates shifted (i.e. a non-standard
        origin, which to_cif also takes into account for the symmetry
        operations)."""
        s = Structure( self.cellpars,
                       [ (e, tuple( x+d for x,d in zip(p,shift) ) )
                         for e,p in self.atoms ],
                       self.spacegroup, self.hall_number )
        s.origin_shift = np.asarray( shift, dtype = float )
        s.unshifted = self
        return s

def lattice_from_cellpars( a, b, c, alpha, beta, gamma ):
    ca, cb, cg = ( math.cos(math.radians(x)) for x in (alpha,beta,gamma) )
    sg = math.sin(math.radians(gamma))
    cy = (ca-cb*cg)/sg
    return np.array( [ [ a, 0.0, 0.0 ],
                       [ b*cg, b*sg, 0.0 ],
                       [ c*cb, c*cy, c*math.sqrt(1.0-cb*cb-cy*cy) ] ] )

def cellpars_from_lattice( lattice ):
    va, vb, vc = ( np.asarray(v,dtype=float) for v in lattice )
    def ang( u, v ):
        return math.degrees( math.acos( float(u@v)
                                        / math.sqrt(float((u@u)*(v@v))) ) )
    return ( float(np.linalg.norm(va)), float(np.linalg.norm(vb)),
             float(np.linalg.norm(vc)), ang(vb,vc), ang(va,vc), ang(va,vb) )

_hall_cache = {}
def default_hall_number( spacegroup ):
    """spglib's default (first) Hall number of a space group."""
    if not _hall_cache:
        for h in range( 1, 531 ):
            _hall_cache.setdefault( spglib.get_spacegroup_type(h).number, h )
    return _hall_cache[spacegroup]

def symops_from_hall( hall_number ):
    """List of (R,t) symmetry operations (including centring)."""
    d = spglib.get_symmetry_from_database( hall_number )
    return [ ( np.asarray(r,dtype=int), np.asarray(t,dtype=float) )
             for r, t in zip( d['rotations'], d['translations'] ) ]

def expand( symops, sites ):
    """Expand (element,position) sites with symmetry operations."""
    R = np.array( [ r for r,_ in symops ], dtype = float )
    T = np.array( [ t for _,t in symops ], dtype = float )
    out, seen = [], set()
    for e, p in sites:
        imgs = ( R @ np.asarray( p, dtype = float ) + T ) % 1.0
        for q in imgs:
            key = tuple( int(k) for k in np.round( q*1e6 ) % 1000000 )
            if key not in seen:
                seen.add( key )
                out.append( ( e, tuple( float(x) for x in q ) ) )
    return out

def min_interatomic_distance( lattice, atoms ):
    """Smallest distance between two atoms (also across cell boundaries)."""
    L = np.asarray( lattice, dtype = float )
    pos = np.array( [ p for _,p in atoms ], dtype = float )
    diff = pos[:,None,:] - pos[None,:,:]
    diff -= np.round( diff )#nearest image (fine for non-extreme cells)
    cells = np.array( list( np.ndindex(3,3,3) ), dtype = float ) - 1.0
    d = np.linalg.norm( ( diff[:,:,None,:] + cells[None,None,:,:] ) @ L,
                        axis = -1 )
    d[ np.arange(len(pos)), np.arange(len(pos)), 13 ] = np.inf#self (n=0)
    return float( d.min() )

def random_cellpars( rng, spacegroup ):
    """Random lattice parameters, compatible with the (spglib standard
    setting of the) space group."""
    def L():
        return round( rng.uniform( 3.0, 7.0 ), 4 )
    def A():
        return round( rng.uniform( 95.0, 115.0 ), 3 )
    a = L()
    if spacegroup <= 2:
        return ( L(), L(), L(), round(rng.uniform(80,100),3),
                 round(rng.uniform(80,100),3), A() )
    if spacegroup <= 15:
        return ( L(), L(), L(), 90.0, A(), 90.0 )
    if spacegroup <= 74:
        return ( L(), L(), L(), 90.0, 90.0, 90.0 )
    if spacegroup <= 142:
        return ( a, a, L(), 90.0, 90.0, 90.0 )
    if spacegroup <= 194:
        return ( a, a, L(), 90.0, 90.0, 120.0 )
    return ( a, a, a, 90.0, 90.0, 90.0 )

_special_points = [ (0,0,0), (0.5,0.5,0.5), (0,0,0.5), (0.5,0,0),
                    (0.25,0.25,0.25) ]
_rhombohedral_sgs = ( 146, 148, 155, 160, 161, 166, 167 )

_special_points_cache = {}
def special_points( spacegroup ):
    """Candidate special positions. All coordinates are exact with 4
    decimals, except (1/3,2/3,z) points which are only used when they are
    special positions (so that rounded coordinates in CIF data can be
    restored exactly by symmetrisation)."""
    if spacegroup in _special_points_cache:
        return _special_points_cache[spacegroup]
    res = list( _special_points )
    if 143 <= spacegroup <= 194 and spacegroup not in _rhombohedral_sgs:
        ops = symops_from_hall( default_hall_number( spacegroup ) )
        def mult( p ):
            return len( expand( ops, [ ('X',p) ] ) )
        #x and y must be fixed: any in-plane displacement increases the
        #multiplicity (otherwise e.g. (x,2x,0) with x=1/3 is a free value):
        eps = 1e-3
        for p in ( (1/3,2/3,0.25), (1/3,2/3,0) ):
            if all( mult( (p[0]+eps*dx, p[1]+eps*dy, p[2]) ) > mult(p)
                    for dx, dy in ( (1,0), (0,1), (1,1), (1,2), (2,1),
                                    (1,-1) ) ):
                res.append( p )
    _special_points_cache[spacegroup] = res
    return res

def random_structure( rng, spacegroup, nsites = 2, max_atoms = 64,
                      min_dist = 1.2, max_tries = 1000 ):
    """Random crystal structure with the given space group (as verified with
    spglib), with nsites symmetry-independent sites (of different elements),
    each either at a special or a general position."""
    hall = default_hall_number( spacegroup )
    ops = symops_from_hall( hall )
    for itry in range( max_tries ):
        if itry == max_tries//2:
            nsites = 1#fall back, in case of high multiplicities
        cellpars = random_cellpars( rng, spacegroup )
        L = lattice_from_cellpars( *cellpars )
        elems = rng.sample( _elements, nsites )
        sites = []
        for e in elems:
            r = rng.random()
            x, y, z = ( round( rng.random(), 4 ) for i in range(3) )
            if r < 0.35:
                p = rng.choice( special_points( spacegroup ) )
            elif r < 0.75:
                #Partially special positions (one or two free parameters):
                p = rng.choice( [ (x,0,0), (x,x,x), (0,y,z), (x,x,z) ] )
            else:
                p = ( x, y, z )
            sites.append( ( e, p ) )
        atoms = expand( ops, sites )
        if len(atoms) > max_atoms:
            continue
        dmin = min_interatomic_distance( L, atoms )
        if dmin < 0.5 * min_dist:
            continue#(nearly) coinciding atoms, do not scale the cell too much
        if dmin < min_dist:
            #Scale up cell (the lengths are rounded to 4 decimals):
            f = 1.05 * min_dist / dmin
            cellpars = tuple( round( x*f, 4 ) for x in cellpars[:3] ) \
                + cellpars[3:]
        s = Structure( cellpars, atoms, spacegroup, hall )
        #Must have the correct symmetry with both strict and the loose
        #tolerances (the latter is used by e.g. cif2ncmat):
        if all( getattr( spglib.get_symmetry_dataset( s.spglib_cell(),
                                                      symprec = sp ),
                         'number', None ) == spacegroup
                for sp in ( 1e-5, 1e-2 ) ):
            return s
    raise RuntimeError(f'could not generate structure in SG-{spacegroup}')

def primitive_cell( structure ):
    """Same crystal described in a primitive cell (listed as P1)."""
    lat, pos, num = spglib.find_primitive( structure.spglib_cell(),
                                           symprec = 1e-5 )
    elems = structure.elements()
    return Structure( cellpars_from_lattice( lat ),
                      [ ( elems[n-1], p ) for p, n in zip( pos, num ) ] )

def transformed_basis( structure, M ):
    """Same crystal (listed as P1) in the basis given by the rows of the
    unimodular integer matrix M acting on the lattice vectors."""
    M = np.asarray( M, dtype = float )
    assert abs( abs(np.linalg.det(M)) - 1.0 ) < 1e-9
    L2 = M @ structure.lattice()
    Minv = np.linalg.inv( M )
    return Structure( cellpars_from_lattice( L2 ),
                      [ ( e, tuple( np.asarray(p) @ Minv ) )
                        for e, p in structure.atoms ] )

def _fmt_coord( x, round_digits ):
    x = float(x) % 1.0
    if round_digits is not None:
        return f'{x:.{round_digits}f}'
    return f'{x:.12f}'

def _fmt_symop( R, t ):
    out = []
    for i in range(3):
        s = ''
        for j, v in enumerate( 'xyz' ):
            c = int(R[i][j])
            if c:
                s += ( '+' if c > 0 and s else '' ) + ( '-' if c < 0 else '' )
                s += ( str(abs(c)) if abs(c) != 1 else '' ) + v
        f = Fraction( float(t[i]) % 1.0 ).limit_denominator( 48 )
        assert abs( float(f) - float(t[i]) % 1.0 ) < 1e-9
        if f:
            s += f'+{f.numerator}/{f.denominator}'
        out.append( s )
    return ','.join( out )

def to_cif( structure, *, style = 'asym', round_digits = None,
            type_symbols = True, uiso = None, name = 'synthetic' ):
    """Create CIF data. Style 'asym' lists unique sites and all symmetry
    operations of the space group, style 'p1' lists all atoms in P1."""
    if style == 'p1':
        ops, sites, sgno, hm = [ ( np.eye(3,dtype=int), np.zeros(3) ) ], \
            structure.atoms, 1, 'P 1'
    else:
        assert style == 'asym'
        base = getattr( structure, 'unshifted', structure )
        shift = getattr( structure, 'origin_shift', np.zeros(3) )
        sites = [ ( e, tuple( np.asarray(p)+shift ) )
                  for e,p in base.unique_sites() ]
        #x' = x+s => ops become R.x' + (t + s - R.s):
        ops = [ ( R, t + shift - R @ shift ) for R,t in base.symops() ]
        sgno = structure.spacegroup
        sgtype = spglib.get_spacegroup_type( structure.hall_number )
        hm = sgtype.international_full
        if sgtype.choice:
            hm += ':' + sgtype.choice#origin choice or axes (e.g. "H")
    a, b, c, alpha, beta, gamma = structure.cellpars
    lines = [ f'data_{name}' ]
    for k, v in ( ('length_a',a), ('length_b',b), ('length_c',c),
                  ('angle_alpha',alpha), ('angle_beta',beta),
                  ('angle_gamma',gamma) ):
        lines.append( f'_cell_{k} {v:.12g}' )
    lines += [ f'_space_group_IT_number {sgno}',
               f"_space_group_name_H-M_alt '{hm}'",
               'loop_', '_space_group_symop_operation_xyz' ]
    lines += [ f"'{_fmt_symop(R,t)}'" for R, t in ops ]
    lines += [ 'loop_', '_atom_site_label' ]
    if type_symbols:
        lines.append( '_atom_site_type_symbol' )
    lines += [ '_atom_site_fract_x', '_atom_site_fract_y',
               '_atom_site_fract_z' ]
    if uiso is not None:
        lines.append( '_atom_site_U_iso_or_equiv' )
    for i, ( e, p ) in enumerate( sites ):
        cols = [ f'{e}{i+1}' ] + ( [ e ] if type_symbols else [] )
        cols += [ _fmt_coord( x, round_digits ) for x in p ]
        if uiso is not None:
            cols.append( f'{uiso[e]:.6g}' )
        lines.append( ' '.join( cols ) )
    return '\n'.join( lines ) + '\n'

#Invalid CIF data (should be rejected):

def invalid_cifs( structure ):
    """Dictionary of named, invalid, CIF data based on the structure."""
    good = to_cif( structure )
    out = {}
    out['missing-cell-length'] = '\n'.join( ll for ll in good.splitlines()
                                            if not ll.startswith('_cell_length_c') )
    i = good.index( 'loop_\n_atom_site_label' )
    out['no-atoms'] = good[:i]
    #Unknown element (both as label and symbol, so it can not be guessed):
    e0 = structure.unique_sites()[0][0]
    out['unknown-element'] = good.replace( f'\n{e0}1 {e0} ', '\nXx1 Xx ' )
    assert out['unknown-element'] != good
    out['impossible-cell-angles'] = '\n'.join(
        ( ll.split()[0] + ' 150' if ll.startswith('_cell_angle') else ll )
        for ll in good.splitlines() ) + '\n'
    #Two different elements related by an inversion centre which is declared,
    #but then coinciding with each other:
    p = ( 0.13, 0.27, 0.31 )
    pm = tuple( (-x)%1.0 for x in p )
    out['overlapping-atoms'] = to_cif( Structure(
        ( 4.0, 4.5, 5.0, 90, 90, 90 ), [ ('Al',p), ('O',pm) ] ),
        style = 'p1' ).replace( "'x,y,z'", "'x,y,z'\n'-x,-y,-z'" ).replace(
            '_space_group_IT_number 1', '_space_group_IT_number 2' ).replace(
                "'P 1'", "'P -1'" )
    return out
