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
import NCrystalDev._common as nc_common
from NCTestUtils.common import ensure_error
from argparse import ArgumentError
import shlex
import re

_re_find_values = re.compile(r'\d+(\.\d+)?([eE][+-]?\d+)?')

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

def test( *args ):
    print(f"============= CLI >>{shlex.join(args)}"
          "<< ====================")
    if '-d' in args or '--dump' in args:
        with FilterHistDumps():
            nc_cli.run('minimc',*args)
    else:
        nc_cli.run('minimc',*args)
    print("===========================================")

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

    test('solid::GdO3/1gcm3','1kT on 1mfp','--decode')
    test( '-e','tally=q;nthreads=97',
          '--decode',
          'solid::GdO3/1gcm3',
          '1kT on 1mfp',
         )
    #fixme: make it possible to test two files for compatibility


if __name__ == '__main__':
    main()
