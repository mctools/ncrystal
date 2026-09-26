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

# NEEDS: spglib
import NCrystalDev as NC
from NCTestUtils.common import ensure_error


def main():
    print('\n\n  ================> Al example 1\n\n')

    c_Al = NC.NCMATComposer()
    c_Al.set_cellsg_cubic( 4.05 )
    c_Al.set_atompos( [ ('Al',0,0,0),
                        ('Al',0,1/2,1/2),
                        ('Al',1/2,0,1/2),
                       ('Al',1/2,1/2,0)])
    c_Al.allow_fallback_dyninfo()
    with ensure_error(NC.NCBadInput,
                      "Must provide a space group number (or invoke"
                      " .refine_crystal_structure()) before it is"
                      " possible to verify a crystal structure"):
        c_Al.create_ncmat()

    c_Al.refine_crystal_structure()

    print( c_Al() )
    c_Al.load().dump()

    print('\n\n  ================> Al example 2\n\n')
    c_Al2 = NC.NCMATComposer()
    c_Al2.set_cellsg_cubic( 4.05, spacegroup=225 )
    c_Al2.set_atompos( [('Al',0,0,0),('Al',0,1/2,1/2),('Al',1/2,0,1/2),('Al',1/2,1/2,0)])
    c_Al2.allow_fallback_dyninfo()
    c_Al2.set_composition('Al','0.99 Al 0.01 Cr')
    print(c_Al2())

    print('\n\n  ================> Al example 3\n\n')

    c_Al3 = NC.NCMATComposer()
    c_Al3.set_cellsg_cubic( 4.05 )
    c_Al3.set_atompos( [('tight_atom',0,0,0),('loose_atom',0,1/2,1/2),('loose_atom',1/2,0,1/2),('loose_atom',1/2,1/2,0)])
    c_Al3.set_dyninfo_msd('tight_atom',msd=0.005, temperature=200)
    c_Al3.set_dyninfo_msd('loose_atom',msd=0.02, temperature=200)
    c_Al3.set_composition('tight_atom','Al')
    c_Al3.set_composition('loose_atom','Al')
    c_Al3.refine_crystal_structure()#Detect spacegroup
    print(c_Al3())

    print('\n\n  ================> Al SANS example\n\n')

    c=NC.NCMATComposer('Al_sg225.ncmat') #<--- Can init from cfg-string.
    c.add_secondary_phase(0.01,'void.ncmat')
    c.add_hard_sphere_sans_model(50)
    print(c())
    c.load().dump()


def test_verify_origin_and_triclinic():
    #verify_crystal_structure must not reject valid structures just because
    #spglib standardises with a different origin (arbitrary for P1 and polar
    #space groups) or to a different (reduced) triclinic basis:
    print('\n\n  ================> verify_crystal_structure\n\n')
    def verify( cell, atoms, sg ):
        c = NC.NCMATComposer()
        a, b, cc, al, be, ga = cell
        c.set_cellsg( a=a, b=b, c=cc, alpha=al, beta=be, gamma=ga,
                      spacegroup = sg )
        c.set_atompos( atoms )
        c.allow_fallback_dyninfo()
        c.verify_crystal_structure( quiet = True )
    tricl = (5.1,6.3,7.2,81,97,103)#mixed acute/obtuse angles
    p1 = [('Al',0.13,0.71,0.29),('O',0.62,0.05,0.93),('Fe',0.41,0.38,0.57)]
    verify( tricl, p1, 1 )
    pm1 = [('Al',0.1,0.2,0.3),('Al',0.9,0.8,0.7),
           ('O',0.35,0.05,0.6),('O',0.65,0.95,0.4)]
    verify( tricl, pm1, 2 )
    p = 4.05/2**0.5 #cell where spglib keeps the basis but shifts the origin:
    verify( (p,p,p,60,60,60), [('Ni',0.0424,0.8626,0.7802),
                               ('Al',0.0126,0.2354,0.7916)], 1 )
    fcc = [('Al',0,0,0),('Al',0,.5,.5),('Al',.5,0,.5),('Al',.5,.5,0)]
    verify( (4.05,)*3+(90,)*3,
            [ (e,(x+.123)%1,(y+.456)%1,(z+.789)%1) for e,x,y,z in fcc ], 225 )
    print('Valid structures verified OK')
    #Wrong space groups must still be rejected:
    for cell, atoms, sg in ( ( tricl, p1, 2 ),
                             ( (4.05,)*3+(90,)*3, fcc, 229 ) ):
        with ensure_error(NC.NCBadInput):
            verify( cell, atoms, sg )

def test_refine_corrections():
    #refine_crystal_structure must only report corrections when spglib
    #actually moves atoms, not for mere basis changes (e.g. triclinic
    #reduction, primitive to centred cells) or origin shifts. Anisotropic
    #properties can only be kept if the basis is unchanged.
    print('\n\n  ================> refine corrections\n\n')
    from NCrystalDev._ncmatimpl import _spglib_refine_cell
    def refine( cell, atoms ):
        c = NC.NCMATComposer()
        a, b, cc, al, be, ga = cell
        c.set_cellsg( a=a, b=b, c=cc, alpha=al, beta=be, gamma=ga )
        c.set_atompos( atoms )
        d = _spglib_refine_cell( c.as_spglib_cell()[0] )
        return d['sgno'], d['warnings'], d['can_keep_anisotropic_properties']
    p = 4.05/2**0.5
    fcc = [('Al',0,0,0),('Al',0,.5,.5),('Al',.5,0,.5),('Al',.5,.5,0)]
    for descr, cell, atoms, expect in (
            ( 'P1, triclinic basis change', (5.1,6.3,7.2,81,97,103),
              [('Al',0.13,0.71,0.29),('O',0.62,0.05,0.93),
               ('Fe',0.41,0.38,0.57)], (1,[],False) ),
            ( 'P1, repeated reductions', (p,p,p,60,60,60),
              [('Al',0.24,0.54,0.38),('O',0.62,0.06,0.94),
               ('Fe',0.40,0.83,0.11)], (1,[],False) ),
            ( 'Imm2, from primitive cell', (p,p,p,60,60,60),
              [('Al',.3,.1,.7),('O',.9,.5,.2)], (44,[],False) ),
            ( 'fcc, shifted origin', (4.05,)*3+(90,)*3,
              [ (e,(x+.1)%1,(y+.2)%1,(z+.3)%1) for e,x,y,z in fcc ],
              (225,[],True) ),
            ( 'fcc, displaced atom', (4.05,)*3+(90,)*3,
              [('Al',0.001,0,0)]+fcc[1:],
              (225,[('Structure received minor corrections by spglib'
                     ' (at the 0.1% level)')],True) ) ):
        res = refine( cell, atoms )
        print(f'{descr}: SG-{res[0]} warnings={res[1]} keep_aniso={res[2]}')
        assert res == expect, (descr,res)

if __name__ == '__main__':
    main()
    test_verify_origin_and_triclinic()
    test_refine_corrections()
