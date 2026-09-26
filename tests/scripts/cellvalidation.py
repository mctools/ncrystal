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

# Tests validation of unit cells and atom positions in NCMAT data (invalid input
# must be rejected rather than silently modified or accepted).

import NCrystalDev as NC
import NCTestUtils.enable_fpe  # noqa: F401
from NCTestUtils.common import ensure_error


def ncmat( lengths, angles, spacegroup, atoms ):
    elems = sorted( { e for e,*_ in atoms } )
    out = [ 'NCMAT v7', '@CELL', f'  lengths {lengths}', f'  angles {angles}' ]
    if spacegroup:
        out += [ '@SPACEGROUP', f'  {spacegroup}' ]
    out += [ '@ATOMPOSITIONS' ] + [ f'  {e} {pos}' for e, pos in atoms ]
    for e in elems:
        n = sum( 1 for ee,_ in atoms if ee == e )
        out += [ '@DYNINFO', f'  element {e}', f'  fraction {n}/{len(atoms)}',
                 '  type vdosdebye', '  debye_temp 300' ]
    return '\n'.join( out ) + '\n'

def load( *args ):
    return NC.directLoad( ncmat( *args ), doScatter = False,
                          doAbsorption = False ).info

def test_hexagonal_angles():
    #Hexagonal/trigonal space groups require gamma=120 (other values were
    #once silently changed to 120 when below 120):
    si = load( '3 3 5', '90 90 120', 191, [('Al','0 0 0')] ).structure_info
    print(f'SG-191 with gamma=120: loaded with gamma={si["gamma"]:g}')
    for gamma in ( 60, 100, 130 ):
        with ensure_error( NC.NCBadInput,
                           'Spacegroup (191) requires alpha=beta=90'
                           ' and gamma=120' ):
            load( '3 3 5', f'90 90 {gamma}', 191, [('Al','0 0 0')] )

def main():
    test_hexagonal_angles()

if __name__ == '__main__':
    main()
