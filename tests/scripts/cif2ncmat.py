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

# NEEDS: gemmi spglib ase

from argparse import ArgumentError

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
from NCTestUtils.env import ncsetenv
from NCTestUtils.common import ( print_text_file_with_snipping,
                                 ensure_error,
                                 work_in_tmpdir )
import pathlib
import contextlib
import shlex

def create_fake_onlinedb_cache_dir():
    #create fake NCRYSTAL_ONLINEDB_CACHEDIR and populate with 3 entries, so we are
    #able to exercise codid:: and mpid:: capabilities in tests.
    d = pathlib.Path('fake_onlinedb_cachedir_testcif2ncmat')
    ncsetenv('ONLINEDB_CACHEDIR',str(d.resolve().absolute()))
    if d.is_dir():
        return
    d.mkdir(parents=True)
    from NCTestUtils.dirs import test_data_dir
    for f in ['cod_1000257.cif','cod_9005777.cif','mp_55.cif']:
        (d/f).write_text((test_data_dir/f'fake_{f}').read_text())

#TODO: Also test python API if not done elsewhere.

def test_cli( args, *, nstart = 100, nend = 20, outfile = None, in_tmp_dir = True ):
    import NCrystalDev.cli as nc_cli
    if isinstance(args,str):
        args = shlex.split(args)
    hr=f"============= CLI >>{shlex.join(args)}<< ===================="
    print(hr)
    ctx = work_in_tmpdir if in_tmp_dir else contextlib.nullcontext
    with ctx():
        nc_cli.run('cif2ncmat',*args)
        if outfile not in ('stdout',None):
            content = pathlib.Path(outfile).read_text()
            print_text_file_with_snipping( content,
                                           nstart=nstart,
                                           nend=nend,
                                           prefix='OUT>')
    print('='*len(hr))

def main():
    ncsetenv('CIF2NCMAT_UNITTEST_NOPLOT','1')
    ncsetenv('ONLINEDB_FORBID_NETWORK','1')
    create_fake_onlinedb_cache_dir()

    test_cli(['--help'])
    test_cli('codid::9005777 --showcif')
    test_cli('codid::9005777 --nomerge',outfile='autogen_CaO3Si_sg2_cod9005777.ncmat')
    test_cli('codid::9005777',outfile='autogen_CaO3Si_sg2_cod9005777.ncmat')
    test_cli('codid::9005777 --uisotemp=50',outfile='autogen_CaO3Si_sg2_cod9005777.ncmat')
    test_cli('codid::9005777 --debyetemp=50',outfile='autogen_CaO3Si_sg2_cod9005777.ncmat')
    test_cli('codid::1000257',outfile='autogen_F7FeNa2Ni_sg74_cod1000257.ncmat')
    test_cli('codid::1000257 --uisotemp=50',outfile='autogen_F7FeNa2Ni_sg74_cod1000257.ncmat')
    test_cli('codid::1000257 --debyetemp=150',outfile='autogen_F7FeNa2Ni_sg74_cod1000257.ncmat')
    test_cli('codid::1000257 --debyetemp=150 --uisotemp=50',outfile='autogen_F7FeNa2Ni_sg74_cod1000257.ncmat')
    test_cli('mpid::55',outfile='autogen_Sn_sg139_mp55.ncmat')
    test_cli('mpid::55 --noval',outfile='autogen_Sn_sg139_mp55.ncmat')
    test_cli('mpid::55 --dynamics stdlib::Sn_sg141.ncmat',outfile='autogen_Sn_sg139_mp55.ncmat')
    test_cli('codid::9005777 --dynamics stdlib::CaSiO3_sg2_Wollastonite.ncmat',outfile='autogen_CaO3Si_sg2_cod9005777.ncmat')
    test_cli('codid::9005777 --showcif')
    test_cli('QE_pw_Al.out --via-ase',outfile='autogen_Al_sg225.ncmat')
    test_cli(['QE_pw_Al.out','--via-ase','--remap','Al is 0.99 Al 0.01 Cr'],
             outfile='autogen_Al99Cr_sg225.ncmat')

    test_cli('dummy_D.cif',outfile='autogen_C3D2_sg221.ncmat')
    test_cli(['dummy_D.cif','--remap','D is H'],outfile='autogen_C3H2_sg221.ncmat')

    with ensure_error(ArgumentError,'invalid --remap syntax in "D = H"'):
        test_cli(['dummy_D.cif','--remap','D = H'],outfile='autogen_C3H2_sg221.ncmat')

    with ensure_error(ArgumentError,'Conflicting options'):
        test_cli('--valplot stdlib::Al_sg225.ncmat --showcif')
    with ensure_error(ArgumentError,'Conflicting options'):
        test_cli('--valplot stdlib::Al_sg225.ncmat --via-ase')

    test_cli(['dummy_D.cif','--remap','D is 0.6 H 0.4 D'],outfile='autogen_C15D4H6_sg221.ncmat')
    test_cli('codid::9005777 -ostdout --via-ase')
    test_cli('mpid::55 --dynamics stdlib::Sn_sg141.ncmat -obla.ncmat',outfile='bla.ncmat')
    test_cli(['codid::9005777','--atomdata',"Si 28.09u::4.1491fm 0.004b 0.171b"],outfile='autogen_CaO3Si_sg2_cod9005777.ncmat')
    test_cli(['codid::9005777','--atomdata',"Si 28.09u 0fm 0b 0b"],outfile='autogen_CaO3Si_sg2_cod9005777.ncmat')

    with ensure_error(ArgumentError,'invalid --atomdata syntax in "Si 28.09 0fm 0b 0b"'):
        test_cli(['codid::9005777','--atomdata',"Si 28.09 0fm 0b 0b"])

    test_cli(['codid::9005777','--atomdata',"Si 28.09u 0fm 0b 0b",
              '--atomdata',"Ca 128.09u 0fm 0b 0b"],
             outfile='autogen_CaO3Si_sg2_cod9005777.ncmat')

    test_cli('--valplot stdlib::Al_sg225.ncmat')
    test_cli('--valplot stdlib::Al_sg225.ncmat -opdf')
    test_cli('--valplot stdlib::Al_sg225.ncmat -obla.pdf')
    test_cli('--valplot stdlib::Al_sg225.ncmat stdlib::Be_sg194.ncmat')
    test_cli('--valplot stdlib::Al_sg225.ncmat stdlib::Be_sg194.ncmat -opdf')
    test_cli('--valplot stdlib::Al_sg225.ncmat stdlib::Be_sg194.ncmat -obla.pdf')
    test_cli('Wyckoff1963_Sulphor_0011255.cif',outfile='autogen_S_sg70.ncmat')
    test_cli(['Wyckoff1963_Sulphor_0011255.cif','--spacegroup','F d d d:1'],outfile='autogen_S_sg70.ncmat')
    test_cli(['Wyckoff1963_Sulphor_0011255.cif','--spacegroup','F d d d:2'],outfile='autogen_S_sg70.ncmat')

    test_cli('mp_55_unrefined_conventional.cif',outfile='autogen_Sn_sg139.ncmat')
    test_cli('mp_55_unrefined_primitive.cif',outfile='autogen_Sn_sg139.ncmat')

if __name__ == '__main__':
    main()
