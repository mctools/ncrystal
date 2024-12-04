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

# NEEDS: numpy

import NCTestUtils.enable_fpe
import NCrystalDev as NC
import NCrystalDev.cliutils as nc_cliutils
from NCrystalDev.ncmat2cpp import ncmat2cpp
from NCrystalDev.misc import AnyTextData

from NCTestUtils.common import ( print_text_file_with_snipping,
                                 ensure_error,
                                 work_in_tmpdir,
                                 fmt_args_as_str )

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

def test_pyapi( *filelist, nstart = 30, nend = 20, **kwargs ):
    args_str = fmt_args_as_str( *filelist, **kwargs )
    print(f"============= PyAPI >>{args_str}<< ====================")
    cpp = ncmat2cpp(*filelist,**kwargs )
    print_text_file_with_snipping(cpp,nstart=nstart,nend=nend,prefix='CPP>')
    print("===========================================")

def test_cli( args, outfile, *, nstart = 30, nend = 20 ):
    print(f"============= CLI >>{shlex.join(args)}<< ====================")
    nc_cliutils.run('ncmat2cpp',*args)
    if outfile not in ('stdout',None):
        cpp = pathlib.Path(outfile).read_text()
        print_text_file_with_snipping(cpp,nstart=nstart,nend=nend,prefix='CPP>')
    print("===========================================")

def main():
    test_pyapi('Al_sg225.ncmat')
    test_pyapi('Al_sg225.ncmat',compact=True)

    with ensure_error(NC.NCBadInput,
                      'Can not accept unnamed text data when converting to C++'
                      ' code, since a string key will be needed in the code.'):
        test_pyapi(_some_ncmat_data)
    NC.registerInMemoryFileData('dummyAl.ncmat',_some_ncmat_data)
    with ensure_error(NC.NCBadInput,
                      'Can not accept unnamed text data when converting to C++'
                      ' code, since a string key will be needed in the code.'):
        test_pyapi(NC.createTextData('virtual::dummyAl.ncmat').rawData)

    test_pyapi(NC.createTextData('virtual::dummyAl.ncmat'))
    test_pyapi( AnyTextData(_some_ncmat_data, name = 'some_name.ncmat' ) )

    test_pyapi('LiquidHeavyWaterD2O_T293.6K.ncmat')
    test_pyapi('Al_sg225.ncmat',
               'stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat',
               nstart=60,
               compact=True)
    test_pyapi('stdlib::Al_sg225.ncmat')

    with ensure_error(NC.NCBadInput,'No files or text data provided'):
        test_pyapi()#empty file list

    test_cli(['--help'],None)
    test_cli(['SrF2_sg225_StrontiumFluoride.ncmat','-o','stdout'],'stdout')

    import argparse
    with ensure_error(argparse.ArgumentError,
                      'the following arguments are required: FILE'):
        test_cli(['-o','stdout'],'stdout')

    with work_in_tmpdir():
        print("Testing work_in_tmpdir()")
    print("work_in_tmpdir() gave no errors (which is good)")

    with work_in_tmpdir():
        pathlib.Path("somefile.ncmat").write_text(_some_ncmat_data)
        assert pathlib.Path("somefile.ncmat").is_file()
        assert pathlib.Path("./somefile.ncmat").is_file()
        test_cli(['somefile.ncmat','-o','bla.cc'],'bla.cc')
        test_cli(['somefile.ncmat','-o','stdout'],None)
        test_cli(['Al_sg225.ncmat','-o','foo.cc'],'foo.cc')
        test_cli(['Al_sg225.ncmat','-o','foo.cc'],'foo.cc')#OK to overwrite
        with ensure_error(NC.NCFileNotFound,'No such file: "./Al_sg225.ncmat"'):
            test_cli(['./Al_sg225.ncmat','-o','foo2.cc'],'foo2.cc')
        test_cli(['stdlib::Al_sg225.ncmat','-o','foo3.cc'],'foo3.cc')
        test_cli(['gasmix::Ar/1gcm3/20K','-o','foo4.cc'],'foo4.cc')
        test_cli(['gasmix::Ar/1gcm3/20K','-o','foo4.cc','--compact'],'foo4.cc')

    test_pyapi('dummyAl.ncmat',
               cppfunctionname = 'somefctname',
               compact=True,
               width=30,
               validate=True,
               regfctname='Some::Name::Space::some_reg_fct',
               extra_includes=['SomeHeader.hh']
               )
    test_pyapi('dummyAl.ncmat',
               cppfunctionname = 'somefctname',
               compact=True,
               width=30,
               validate=True,
               regfctname='Some::Name::Space::some_reg_fct',
               extra_includes=['SomeHeader.hh'],
               quiet = True
               )

    #FIXME: It would be useful with a test that we can actually compile the
    #output (but... we do use the functionality for embedding). Perhaps we need
    #a special/dedicated CTest?

if __name__ == '__main__':
    main()
