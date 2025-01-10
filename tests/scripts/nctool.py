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

# NEEDS: numpy matplotlib

#gives problems with matplotlib: import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev as NC
import NCrystalDev.cliutils as nc_cliutils
#from NCrystalDev.misc import AnyTextData
from NCTestUtils.env import ncsetenv

from NCTestUtils.common import ( ensure_error,
                                 work_in_tmpdir )

import pathlib
import shlex

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

def test_cli( args, *, nstart = 30, nend = 20 ):
    print(f"============= CLI >>{shlex.join(args)}<< ====================")
    nc_cliutils.run('nctool',*args)
    print("===========================================")

def main():
    test_cli(['--help'])

    test_cli(['SrF2_sg225_StrontiumFluoride.ncmat','-d'])
    with ensure_error(NC.NCFileNotFound,
                      'Could not find data: "some_nonexistant_file.ncmat"'):
        test_cli(['some_nonexistant_file.ncmat','-d'])

    NC.registerInMemoryFileData('dummyAl.ncmat',_some_ncmat_data)
    test_cli(['dummyAl.ncmat','-d'])
    test_cli(['dummyAl.ncmat','--cfg'])
    test_cli(['dummyAl.ncmat;temp=300K','--cfg'])

    test_cli(['--test'])

    test_cli(['--doc'])
    test_cli(['--doc','--doc'])
    test_cli(['--extract','dummyAl.ncmat'])
    with ensure_error(NC.NCFileNotFound,
                      'Could not find data: "dummyAl.ncmat;temp=300K"'):
        test_cli(['--extract','dummyAl.ncmat;temp=300K'])
    import argparse
    with ensure_error(argparse.ArgumentError,
                      'argument --extract: expected one argument'):
        test_cli(['dummyAl.ncmat','--extract'])

    with work_in_tmpdir():
        print("Testing work_in_tmpdir()")
        test_cli(['dummyAl.ncmat;temp=300K','--cfg'])
        test_cli(['--extract','dummyAl.ncmat'])
        test_cli(['stdlib::Al_sg225.ncmat;temp=300K','--cfg'])
        with ensure_error(NC.NCFileNotFound,
                          'Could not find data: "missing.ncmat"'):
            test_cli(['missing.ncmat;temp=300K','--cfg'])
        pathlib.Path('./missing.ncmat').write_text(_some_ncmat_data)
        test_cli(['missing.ncmat;temp=300K','--cfg'])
    with ensure_error(NC.NCFileNotFound,
                      'Could not find data: "missing.ncmat"'):
        test_cli(['missing.ncmat;temp=300K','--cfg'])

    ncsetenv('TOOL_UNITTESTS','1')
    ncsetenv('DPI','75')
    test_cli(['--bench','Al_sg225.ncmat'])
    test_cli(['Al_sg225.ncmat;dcutoff=0.5;vdoslux=2'])
    test_cli(['Al_sg225.ncmat;dcutoff=0.5;vdoslux=2','-a','--energy','-x','1e-3:1e2'])
    with ensure_error(argparse.ArgumentError,
                      'argument --mc: expected 2 arguments'):
        test_cli(['--mc'])
    test_cli(['--mc','-h'])
    test_cli(['--mc','-h'])

    test_cli(['--mc','2Aa','2cm','Al_sg225.ncmat;temp=250K'])
    test_cli(['--mc','constant;wl=2','sphere;r=0.02',
              'Al_sg225.ncmat;vdoslux=1;comp=inelas;inelas=0'])

    for a in ['coh_elas','bragg','incoh_elas','sans','elastic','inelastic']:
        test_cli(['--cfg','Al_sg225.ncmat;dcutoff=0.5;vdoslux=2',f'--{a}'])

    test_cli(['Al_sg225.ncmat;vdoslux=1;temp=1000K',
              'Al_sg225.ncmat;vdoslux=1;temp=10K',
              '--logy'])
    test_cli(['Al_sg225.ncmat;vdoslux=1;temp=1000K',
              'Al_sg225.ncmat;vdoslux=1;temp=10K',
              '--common','dcutoff=1.0Aa','--liny'])

    #This doesn't work, since we have issue #202:
    #assert not pathlib.Path("./ncrystal.pdf").is_file()
    #    test_cli(['Al_sg225.ncmat;vdoslux=1;temp=1000K','--pdf'])
    #    assert pathlib.Path("./ncrystal.pdf").is_file()

    with ensure_error(RuntimeError,'Currently we do not support switching back'
                      ' and forth between pdf plotting in the same process (see'
                      ' also https://github.com/mctools/ncrystal/issues/202)'):
        test_cli(['Al_sg225.ncmat;vdoslux=1;temp=1000K','--pdf'])

    test_cli(['--cfg',
              'phases<0.65*Al_sg225.ncmat&0.35'
              '*MgO_sg225_Periclase.ncmat>;temp=100K'])

    test_cli(['--phases',
              'phases<0.65*Al_sg225.ncmat&0.35'
              '*MgO_sg225_Periclase.ncmat>;temp=100K'])

    test_cli(['Al_sg225.ncmat;temp=10 K;vdoslux=1',
              'Al_sg225.ncmat;vdoslux=1;temp=10K',
              '--common','dcutoff=1.0Aa','--liny'])

if __name__ == '__main__':
    main()
