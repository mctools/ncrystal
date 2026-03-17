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

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev.cli as nc_cli
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev._common as nc_common
from NCTestUtils.common import ensure_error
from argparse import ArgumentError
import shlex
import re
import pathlib
import sys
import io
import json

_re_find_values = re.compile(r'(-)?\d+(\.\d+)?([eE][+-]?\d+)?')

class FilterHistDumps:
    def __p( self, *args, **kwargs ):
        patterns = ('tallied ray count:',
                    'integral :','mean :','rms :',
                    'minfilled :','maxfilled :',
                    'underflow :','overflow :')
        args=list(args)
        for i in range(len(args)):
            a = args[i]
            newa = []
            for e in a.splitlines(keepends=True):
                e_norm = ' '.join(e.split())
                if any(n in e_norm for n in patterns):
                    newa.append(_re_find_values.sub('<HIDDEN>',e))
                else:
                    newa.append(e)
            args[i] = ''.join(newa)
        return self.__op( *newa, **kwargs )
    def __enter__( self ):
        self.__op = nc_common.set_ncrystal_print_fct( self.__p )
    def __exit__( self, *a, **kw ):
        nc_common.set_ncrystal_print_fct( self.__op )

class FilterCatchStdout:
    def __enter__(self):
        self.__buf = io.StringIO()
        self.__orig = sys.stdout
        sys.stdout = self.__buf
        return self
    def __exit__(self, *a, **kw):
        sys.stdout = self.__orig
    def getvalue(self):
        return self.__buf.getvalue()

def test( *args, test_catch_stdout = False ):
    print(f"============= CLI >>{shlex.join(args)}"
          "<< ====================")
    retval = None
    if test_catch_stdout:
        with FilterCatchStdout() as buf:
            nc_cli.run('minimc',*args)
        retval = buf.getvalue()
    elif '-d' in args or '--dump' in args:
        with FilterHistDumps():
            nc_cli.run('minimc',*args)
    else:
        nc_cli.run('minimc',*args)
    print("===========================================")
    return retval

def main():
    with ensure_error(ArgumentError,
                      'Missing arguments. Use -h, --help or --full-help'
                      ' for information about proper usage.'):
        test()

    test('--hel')
    test('--full')
    test('-e','help')
    test('--srccfg=help')
    test('--geomcfg','help')
    test('solid::GdO3/1gcm3','-d')
    test('solid::GdO3/1gcm3','2Aa on 2cm','-d')
    test('solid::GdO3/1gcm3','2Aa on 2cm','-e','ignoremiss=1','-d')
    test('solid::GdO3/1gcm3','--decode')
    test('solid::GdO3/1gcm3','-g',"sphere;r=0.0102578",
         '-s',
         "  \t circular;  ;r=0.0102578;ekin=0.0252617;"
         "  z=-0.0102578 ; n=100.0e4",
         '-e', '', '--decode')

    test('solid::GdO3/1gcm3','-o','bla.json.gz')
    with ensure_error(ArgumentError,
                      'Do not specify CFGSTR when using --input'):
        test('solid::GdO3/1gcm3','-i','bla.json.gz','--dump')
    test('-i','bla.json.gz','--dump')
    test('-i','bla.json.gz','--decode')

    test('solid::GdO3/1gcm3','1e3 times','-o','bla_qemu.json.gz',
         '-e','nthreads=1;tallybins=e:50:0:0.5,q:50:0:30','-t','q,e,mu')
    test('-i','bla_qemu.json.gz','--decode',)
    test('-i','bla_qemu.json.gz','--decode','--quiet')
    with ensure_error(ArgumentError,
                      'Incompatible options: --decode and --dump'):
        test('-i','bla_qemu.json.gz','--decode','--dump')
    with ensure_error(ArgumentError,
                      'Incompatible options: --decode and --plot'):
        test('-i','bla_qemu.json.gz','--decode','--plot')
    test('-i','bla_qemu.json.gz','--dump')
    test('-i','bla_qemu.json.gz','--dump','--tally=e,mu')
    test('-i','bla_qemu.json.gz','--dump','--tally=e,l')
    buf = test('-i','bla_qemu.json.gz','-o','stdout',test_catch_stdout=True)
    json.loads(buf)#<-- check for valid json
    from NCrystalDev.minimc import MMCResults
    res = MMCResults(buf)
    with FilterHistDumps():
        res.dump()
    test('-i','bla_qemu.json.gz','-o','bla_qemu_2.json.gz')

    assert ( pathlib.Path('bla_qemu.json.gz').read_bytes()
             == pathlib.Path('bla_qemu_2.json.gz').read_bytes() )

    test('solid::GdO3/1gcm3','1kT on 1mfp','--decode')
    test( '-e','tally=q;nthreads=97',
          '--decode',
          'solid::GdO3/1gcm3',
          '1kT on 1mfp',
         )
    test('solid::GdO3/1gcm3','1kT on 1mfp','--decode','--quiet')

if __name__ == '__main__':
    main()
