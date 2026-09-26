
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

"""Tests of ncrystal_cif2ncmat with synthetic CIF data from cifgen: valid CIF
data in many styles must reproduce the original crystal, and invalid CIF data
must be rejected."""

import pathlib
import random

import NCrystalDev as NC
import NCrystalDev.cli as nc_cli

from . import cifgen
from .common import work_in_tmpdir
from .env import ncsetenv

_uiso_temp = 200.0

def convert( cif, extra_args = () ):
    """Run cif2ncmat on the CIF data, returning the NCMAT data."""
    pathlib.Path('in.cif').write_text( cif )
    nc_cli.run( 'cif2ncmat', 'in.cif', '-o', 'out.ncmat', '-q', *extra_args )
    return pathlib.Path('out.ncmat').read_text()

def reference_ncmat( structure, uiso = None ):
    """NCMAT data of the exact structure (no CIF or spglib involved)."""
    c = NC.NCMATComposer()
    a, b, cc, alpha, beta, gamma = structure.cellpars
    c.set_cellsg( a=a, b=b, c=cc, alpha=alpha, beta=beta, gamma=gamma )
    c.set_atompos( [ (e,*p) for e,p in structure.atoms ] )
    for e in structure.elements():
        if uiso:
            c.set_dyninfo_msd( e, uiso[e], temperature = _uiso_temp )
        else:
            c.set_dyninfo_debyetemp( e, 300.0 )
    return c.create_ncmat( verify_crystal_structure = False )

_wavelengths = ( 1.3, 2.0, 3.0, 4.5, 7.0 )

def physics( ncmat ):
    """Space group, volume/atom, density and coherent elastic cross sections
    (per atom) at a few wavelengths."""
    m = NC.directLoad( ncmat, 'dcutoff=0.6;comp=coh_elas',
                       doAbsorption = False )
    si = m.info.structure_info
    return ( si['spacegroup'], si['volume']/si['n_atoms'], m.info.density,
             [ float(x) for x in m.scatter.xsect( wl = _wavelengths ) ] )

def check_same_physics( name, ref, res, spacegroup ):
    sg, vpa, dens, xs = res
    _, vpa_ref, dens_ref, xs_ref = ref
    assert sg == spacegroup, (name,sg,spacegroup)
    assert abs( vpa/vpa_ref - 1.0 ) < 1e-9, (name,vpa,vpa_ref)
    assert abs( dens/dens_ref - 1.0 ) < 1e-9, (name,dens,dens_ref)
    xsmax = max( xs_ref )
    for x, xr in zip( xs, xs_ref ):
        assert abs( x - xr ) < 1e-6*xsmax, (name,xs,xs_ref)

def valid_variants( s, all_variants = True ):
    """Named variants of valid CIF data for structure s, as (cif, extra
    cif2ncmat args, uiso values)."""
    dt = ( '--debyetemp', '300' )
    v = { 'asym-unit+symops': ( cifgen.to_cif( s ), dt, None ) }
    if not all_variants:
        return v
    uiso = { e : round( 0.005 + 0.002*i, 4 )
             for i, e in enumerate( s.elements() ) }
    v['all-atoms-P1'] = ( cifgen.to_cif( s, style = 'p1' ), dt, None )
    v['primitive-cell-P1'] = ( cifgen.to_cif( cifgen.primitive_cell(s),
                                              style = 'p1' ), dt, None )
    v['other-basis-P1'] = ( cifgen.to_cif( cifgen.transformed_basis(
        s, [[1,0,0],[1,1,0],[0,1,1]] ), style = 'p1' ), dt, None )
    v['rounded-coords'] = ( cifgen.to_cif( s, round_digits = 4 ), dt, None )
    v['origin-shift'] = ( cifgen.to_cif( s.with_origin_shift(
        (0.125,0.25,0.375) ) ), dt, None )
    #Elements must be deduced from labels (like "Fe1"):
    v['no-type-symbols'] = ( cifgen.to_cif( s, type_symbols = False ),
                             dt, None )
    v['uiso'] = ( cifgen.to_cif( s, uiso = uiso ),
                  ( '--uisotemp', str(_uiso_temp) ), uiso )
    return v

def tolerated_variants( s ):
    """Inconsistent but recoverable CIF data (the explicitly listed symmetry
    operations must be used, and are expected to be completed if needed)."""
    dt = ( '--debyetemp', '300' )
    cif = cifgen.to_cif( s )
    ops = [ ll for ll in cif.splitlines() if ll.startswith("'") ]
    v = {}
    if len(ops) > 2:#(with only 2, the result is just P1 and not recoverable)
        v['incomplete-symops'] = ( cif.replace( ops[-1]+'\n', '' ), dt )
    v['wrong-HM-symbol'] = ( cif.replace(
        cif.split("_space_group_name_H-M_alt ")[1].split('\n')[0],
        "'P 1'" if s.spacegroup != 1 else "'P -1'" ), dt )
    return v

def check_structure( s, *, all_variants, tolerated ):
    ref = physics( reference_ncmat( s ) )
    variants = valid_variants( s, all_variants )
    for name, ( cif, args, uiso ) in variants.items():
        ncmat = convert( cif, args )
        check_same_physics( name, physics( reference_ncmat( s, uiso ) )
                            if uiso else ref, physics( ncmat ), s.spacegroup )
    ntol = 0
    if tolerated:
        for name, ( cif, args ) in tolerated_variants( s ).items():
            check_same_physics( name, ref, physics( convert(cif,args) ),
                                s.spacegroup )
            ntol += 1
    return len(variants), ntol

def check_invalid( s ):
    for name, cif in cifgen.invalid_cifs( s ).items():
        try:
            convert( cif, ( '--debyetemp', '300' ) )
        except NC.NCBadInput as e:
            print(f'  invalid CIF "{name}" rejected: {str(e)[:60]}')
        else:
            raise RuntimeError(f'invalid CIF "{name}" was not rejected')

def run( spacegroups, *, seed, full_for = None, invalid_for = () ):
    """Check a random structure in each space group. For those in full_for
    (default: all), all valid and tolerated CIF variants are checked,
    otherwise just the asymmetric unit + symops style."""
    ncsetenv( 'CIF2NCMAT_UNITTEST_NOPLOT', '1' )
    ncsetenv( 'ONLINEDB_FORBID_NETWORK', '1' )
    rng = random.Random( seed )
    nconv = 0
    with work_in_tmpdir():
        for sg in spacegroups:
            try:
                s = cifgen.random_structure( rng, sg, max_atoms = 100 )
            except RuntimeError:
                #Some cubic space groups require many atoms:
                s = cifgen.random_structure( rng, sg, max_atoms = 200 )
            full = full_for is None or sg in full_for
            nv, nt = check_structure( s, all_variants = full,
                                      tolerated = full )
            nconv += nv + nt
            print(f'SG-{sg:<3} ({len(s.atoms)} atoms): {nv} valid'
                  f' and {nt} tolerated CIF variants OK')
            if sg in invalid_for:
                check_invalid( s )
    print(f'All OK ({nconv} conversions)')
