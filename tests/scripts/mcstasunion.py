#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
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

import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.common import fix_ncrystal_version_printouts
import NCrystalDev.cli as nc_cli
import shlex
import pathlib

_some_ncmat_data="""NCMAT v5
# Here is a comment
@DENSITY
  1 g_per_cm3#And another one
@DYNINFO
  element H
  fraction 1
  type vdosdebye
  debye_temp 300
"""

def test_cli( *args ):
    print(f"============= CLI >>{shlex.join(args)}<< ====================")
    nc_cli.run('mcstasunion',*args)
    print("===========================================")

def main():
    from NCrystalDev._common import print
    fix_ncrystal_version_printouts()
    test_cli('--help')
    test_cli('--split','MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5')
    test_cli('--split','MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=10.0')
    test_cli('--split','MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5;comp=bragg')
    test_cli('--split','MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5;inelas=0')
    test_cli('MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5')
    test_cli('MyMat','phases<0.65*Al_sg225.ncmat;dcutoff=0.4&0.35'
             '*MgO_sg225_Periclase.ncmat;dcutoff=0.8>;temp=100K;vdoslux=0',
             '--split')

    test_cli('MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5',
             '-o','myblafile.c')

    print("Contents of myblafile.c:",'='*40)
    print(pathlib.Path('myblafile.c').read_text())

    import NCrystalDev.mcstasutils as ncmu

    print("TEST PyAPI1",'='*60)
    o = ncmu.cfgstr_2_union_instrument_code(
        cfgstr='Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5;inelas=0',
        name='BlaBla',
        split_by_physics = True
    )
    print(o)
    print("TEST PyAPI2",'='*60)
    o = ncmu.cfgstr_2_union_instrument_code(
        cfgstr='Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5;inelas=0',
        name='BlaBla',
        split_by_physics = False
    )
    print(o)

if __name__ == '__main__':
    main()
