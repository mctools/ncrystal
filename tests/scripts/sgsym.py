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

# NEEDS: numpy spglib gemmi

# Tests of the sgsym component (space groups and their symmetries),
# mostly by comparing with spglib and gemmi.

from fractions import Fraction

import gemmi
import NCTestUtils.enable_fpe  # noqa: F401
import numpy as np
import spglib
from NCTestUtils.loadlib import Lib

lib = Lib('testsgsym')

# Verifies NCrystal's hardwired table of the 530 space group settings (Hall
# numbers, space group numbers, setting choice codes and Hall symbols) against
# spglib and gemmi, including that the Hall symbols give the same symmetry
# operations as those in spglib's database.
def test_table():
    def spglib_type( hn ):
        t = spglib.get_spacegroup_type( hn )
        get = ( lambda k : getattr( t, k ) ) if hasattr( t, 'number' ) else t.get
        return get('number'), get('choice'), ' '.join( get('hall_symbol').split() )

    def ops_set( rots, trans ):
        return frozenset( ( tuple( int(x) for x in np.asarray(r).flatten() ),
                            tuple( Fraction( round( float(x)*24 ) % 24, 24 )
                                   for x in t ) )
                          for r, t in zip( rots, trans ) )

    def spglib_ops( hn ):
        d = spglib.get_symmetry_from_database( hn )
        return ops_set( d['rotations'], d['translations'] )

    def gemmi_ops( hall ):
        ops = gemmi.symops_from_hall( hall )
        return ops_set( [ np.array(o.rot)//gemmi.Op.DEN for o in ops ],
                        [ np.array(o.tran)/gemmi.Op.DEN for o in ops ] )

    gemmi_itb = list( gemmi.spacegroup_table_itb() )
    assert len( gemmi_itb ) == 530
    first_hn = {}
    for hn in range( 1, 531 ):
        number = lib.nctest_sgsym_number( hn )
        choice = lib.nctest_sgsym_choice( hn )
        hall = lib.nctest_sgsym_hallsymbol( hn )
        #spglib:
        assert ( number, choice, hall ) == spglib_type( hn ), ( hn, number, choice,
                                                                hall )
        assert gemmi_ops( hall ) == spglib_ops( hn ), hn
        #gemmi:
        g = gemmi_itb[hn-1]
        gchoice = ( g.ext + g.qualifier ).strip('\x00 ')
        assert ( g.number, gchoice, ' '.join( g.hall.split() ) ) == ( number, choice,
                                                                     hall ), hn
        #String round trips:
        s = lib.nctest_sgsym_tostring( hn )
        assert s == ( f'{number}:{choice}' if choice else str(number) )
        assert lib.nctest_sgsym_parse( s ) == hn
        first_hn.setdefault( number, hn )
        assert lib.nctest_sgsym_parse( str(number) ) == first_hn[number]
    assert sorted( first_hn ) == list( range( 1, 231 ) )
    assert lib.nctest_sgsym_parse( '227:3' ) == 0
    print('All 530 space group settings consistent with spglib and gemmi')

def main():
    test_table()

if __name__ == '__main__':
    main()
