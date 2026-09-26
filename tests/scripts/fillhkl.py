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

# Tests related to the generation of lists of HKL planes.

import re

import NCrystalDev as NC
import NCTestUtils.enable_fpe  # noqa: F401


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

def main():
    test_nosgfamilies()

if __name__ == '__main__':
    main()
