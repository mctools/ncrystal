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

# NEEDS: numpy

import sys
if '--plot' not in sys.argv[1:]:
    import NCTestUtils.enable_fpe # noqa F401

import NCrystalDev.cli as nc_cli
import NCrystalDev._common as nc_common
from NCrystalDev.exceptions import NCBadInput

import NCTestUtils.reprint_escaped_warnings # noqa F401
from NCTestUtils.env import ncsetenv
from NCTestUtils.common import ensure_error

from argparse import ArgumentError
import shlex
import re
import pathlib
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

def main(do_plot):
    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')


    with ensure_error(ArgumentError,
                      'Missing arguments. Use -h, --help or --full-help'
                      ' for information about proper usage.'):
        test()

    test('--hel')
    test('--full')
    test('-e','help')
    test('--srccfg=help')
    with ensure_error(ArgumentError,
                      'Use -h/--help for instructions'
                      ' (or --full-help for even more)'):
        test('help')
    test('--geomcfg','help')
    test('','-t','help')
    test('-t','help')
    with ensure_error(ArgumentError,
                      'Missing file: missing.json'):
        test('solid::GdO3/1gcm3','-i','missing.json')
    test('solid::GdO3/1gcm3','-d')
    test('solid::GdO3/1e-20gcm3','-d')
    test('solid::GdO3/1gcm3','help')
    test('solid::GdO3/1gcm3','2Aa on 2cm','-d')
    test('solid::GdO3/1gcm3;comp=','2Aa','-d')
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
    with ensure_error(ArgumentError,
                      'Incompatible options: --doc and --input'):
        test('-i','bla.json.gz','--doc=engine')
    with ensure_error(ArgumentError,
                      'Do not specify --enginecfg when using --input (use'
                      ' --tally instead if you are trying to filter tallies).'):
        test('-i','bla.json.gz','-e','')
    with ensure_error(ArgumentError,
                      'Do not specify --geomcfg or --srccfg when using --input'):
        test('-i','bla.json.gz','-s','')
    with ensure_error(ArgumentError,
                      'Do not specify --geomcfg or --srccfg when using --input'):
        test('-i','bla.json.gz','-g','')
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

    with ensure_error(ArgumentError,
                      'Output file aready exists: bla_qemu.json.gz'):
        test('solid::GdO3/1gcm3','1e3 times','-o','bla_qemu.json.gz')
    pathlib.Path('./bla_dummydir').mkdir()
    with ensure_error(ArgumentError,
                      'Output file is a directory: ./bla_dummydir'):
        test('solid::GdO3/1gcm3','1e3 times','-o','./bla_dummydir')
    with ensure_error(ArgumentError,
                      'Output directory not found: bla_dummydir2'):
        test('solid::GdO3/1gcm3','1e3 times','-o','./bla_dummydir2/bla.json')
    with ensure_error(ArgumentError,
                      'Invalid filename (must end with .json or .json.gz): bla_123.txt'):
        test('solid::GdO3/1gcm3','1e3 times','-o','./bla_123.txt')

    with ensure_error(ArgumentError,
                      '--output=stdout is incompatible with --decode'):
        test('solid::GdO3/1gcm3','1kT on 1mfp','--decode','--output=stdout')
    with ensure_error(ArgumentError,
                      '--output=stdout is incompatible with --dump'):
        test('solid::GdO3/1gcm3','1kT on 1mfp','--dump','--output=stdout')
    with ensure_error(ArgumentError,
                      '--output=stdout is incompatible with --doc'):
        test('solid::GdO3/1gcm3','1kT on 1mfp','--doc','engine',
             '--output=stdout')
    with ensure_error(ArgumentError,
                      '--output=stdout is incompatible with --plot'):
        test('solid::GdO3/1gcm3','1kT on 1mfp','--plot','--output=stdout')
    with ensure_error(ArgumentError,
                      'Incompatible options: --decode and --doc'):
        test('solid::GdO3/1gcm3','1kT on 1mfp','--doc','engine','--decode')
    #The magic "help" strings do not trigger the usual checking of other
    #arguments being consistent, they are merely a handy shortcut which should
    #always work:
    test('solid::GdO3/1gcm3','help','--decode')
    test('solid::GdO3/1gcm3','-t','help')

    with ensure_error(ArgumentError,'Missing CFGSTR argument'):
        test('-g','sphere;r=1','-s','constant;wl=1.8')
    with ensure_error(ArgumentError,
                      'Do not supply CFGSTR argument with --doc'):
        test('void.ncmat','--doc=engine')
    with ensure_error(ArgumentError,
                      'Do not supply --srccfg or --geomcfg if also supplying'
                      ' a SCENARIO string'):
        test('void.ncmat','2Aa','-s','constant;wl=1.8')
    with ensure_error(ArgumentError,
                      'Do not supply --srccfg or --geomcfg if also supplying'
                      ' a SCENARIO string'):
        test('void.ncmat','2Aa','-g','sphere;r=1')

    with ensure_error(ArgumentError,
                      'Options --srccfg and --geomcfg must always be supplied'
                      ' together (--enginecfg is optional and will default to'
                      ' an empty string).'):
        test('void.ncmat','-g','sphere;r=1')
    with ensure_error(ArgumentError,
                      'Option --doc and --output can not be used together.'):
        test('--doc=engine','-o','foobla_qemu.json.gz')
    with ensure_error(ArgumentError,
                      'Unsupported tally flag: foo'):
        test('void.ncmat','2Aa','-t','q,q,e,   foo, mu')

    with ensure_error(ArgumentError,
                      'Can not use --tally when the --enginecfg'
                      ' also contains "tally=..."'):
        test('void.ncmat','2Aa','-t','q','-e','   tally =   ;ddsfsdf')

    with ensure_error(NCBadInput,
                      'Invalid parameter for chosen engine: "sdftally"'):
        test('void.ncmat','2Aa','-t','q','-e','   sdftally =   ')

    test('void.ncmat','2Aa','-o','bla_hello.json')
    with ensure_error(ArgumentError,
                      'Output file aready exists: bla_hello.json'):
        test('void.ncmat','2Aa','-o','bla_hello.json')
    test('void.ncmat','2Aa','-o','bla_hello.json','--force')


    test('void.ncmat','2Aa')#should auto plot
    test('void.ncmat','2Aa','-e','tally=')#should auto plot NO tallies
    test('void.ncmat','2Aa','-e','tallybreakdown=0','--plot')

    test('void.ncmat','2Aa','-t','q,mu,e','-o','bla17.json.gz')
    test('-i','bla17.json.gz','--plot')
    test('-i','bla17.json.gz','-t','mu,q','--plot')#skips e


if __name__ == '__main__':
    main(do_plot = '--plot' in sys.argv[1:])
