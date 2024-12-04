#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

#FIXME: Update NEEDS

# NEEDS: numpy

import NCTestUtils.enable_fpe
import NCrystalDev as NC
import NCrystalDev.cliutils as nc_cliutils

from NCTestUtils.common import ( ensure_error,
                                 work_in_tmpdir )

import pathlib
import shlex
import os

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

#FIXME: Also test python API, once we create one...

def test_cli( *args ):
    print(f"============= CLI >>{shlex.join(args)}<< ====================")
    nc_cliutils.run('mcstasunion',*args)
    print("===========================================")

def main():
    test_cli('--help')
    test_cli('--split','MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5')
    test_cli('--split','MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=10.0')
    test_cli('--split','MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5;comp=bragg')
    test_cli('--split','MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5;inelas=0')
    test_cli('MyMat','Al_sg225.ncmat;temp=20K;vdoslux=0;dcutoff=0.5')
    test_cli('MyMat','phases<0.65*Al_sg225.ncmat;dcutoff=0.4&0.35'
             '*MgO_sg225_Periclase.ncmat;dcutoff=0.8>;temp=100K;vdoslux=0',
             '--split')

    #FIXME also test the output when writing to files!!

if __name__ == '__main__':
    main()

