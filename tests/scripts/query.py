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
import NCTestUtils.stabilise_ncpprint # noqa F401
import NCrystalDev.cli as nc_cli
from NCrystalDev.exceptions import NCBadInput
from NCTestUtils.env import ncsetenv
from NCTestUtils.common import ensure_error
from argparse import ArgumentError
import shlex
import pathlib
from NCrystalDev.misc import evaluate_query as ncquery

def safefmt(s):
    return str(s).encode('unicode_escape').decode('utf-8')

def test(*query):
    print()

    print(f">>> Sending query: {[safefmt(e) for e in query]!r}")
    r = ncquery(query)
    import pprint
    pprint.pprint(r)

def test_cli( *args ):
    print(f"============= CLI >>{shlex.join(safefmt(a) for a in args)}"
          "<< ====================")
    nc_cli.run('query',*args)
    print("===========================================")

def main():
    from NCrystalDev._common import print
    ncsetenv('FIX_QUERY_VERSION_FOR_TESTS','1')

    with ensure_error(NCBadInput,
                      'Missing or empty query'):
        test()
    with ensure_error(NCBadInput,
                      'Invalid JSON query key: "foo"'):
        test('foo')

    with ensure_error(NCBadInput,
                      'Invalid JSON query key: ""'):
        test('')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc",""]'):
        test('mmc','')

    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: "-bla"'):
        test('mmc','-bla')
    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: "-bla"'):
        test('mmc','-bla')
    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: " -bla"'):
        test('mmc',' -bla')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","bla-bla"]'):
        test('mmc','bla-bla')

    with ensure_error(NCBadInput,
                      "Invalid query (not all entries are strings):"
                      " ('foo', 117)"):
        test('foo',117)

    for a in ('bar'+'\x07','\x07','\x07\x07','\x07'+'bar'):
        with ensure_error(NCBadInput,
                          'Can not use character 0x07 (ASCII BEL) in JSON'
                          ' query strings'):
            test('foo',a)

    test_cli('--help')
    test_cli('version')
    test_cli('version','-j')
    fn = 'tmp_test_file.json'
    if pathlib.Path(fn).exists():
        pathlib.Path(fn).unlink()
    test_cli('version','-o',fn)
    print(f"Contents of {fn}:",'='*40)
    print(pathlib.Path(fn).read_text())
    with ensure_error(RuntimeError,
                      'File already exists (use --force to overwrite):'
                      ' tmp_test_file.json'):
        test_cli('list','-o',fn)

    test_cli('list','-o',fn,'--force')
    print(f"Contents of {fn}:",'='*40)
    print(pathlib.Path(fn).read_text())

    with ensure_error(NCBadInput,
                      'Invalid JSON query key: "foo"'):
        test_cli('foo')

    test_cli('list')
    test_cli('mmc','list')
    test_cli('util','list')
    test_cli('util','wl2ekin','1.8')
    test_cli('util','ekin2wl','0.025')
    test_cli('util','eqrefl','149','1,2,3','1,1,0')
    test_cli('util','eqrefl','150','1,2,3','1,1,0')

    with ensure_error(ArgumentError,
                      'the following arguments are required: STR'):
        test_cli()
    with ensure_error(NCBadInput,
                      'Invalid JSON query key: "foo"'):
        test_cli('foo')
    with ensure_error(NCBadInput,
                      'Invalid JSON query key: ""'):
        test_cli('')

    with ensure_error(NCBadInput,
                      'Invalid JSON query key: ""'):
        test_cli('')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc",""]'):
        test_cli('mmc','')

    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: " -bla"'):
        test_cli('mmc',' -bla')

    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","bla-bla"]'):
        test_cli('mmc','bla-bla')

    for a in ('bar'+'\x07','\x07','\x07\x07','\x07'+'bar'):
        with ensure_error(NCBadInput,
                          'Can not use character 0x07 (ASCII BEL) in JSON'
                          ' query strings'):
            test_cli('foo',a)

    test_cli('sab','refeval','0.025',
             'stdlib::Al2O3_sg167_Corundum.ncmat','3','Al')
    test_huge_arrays()

def test_huge_arrays():
    #In the mode where huge arrays are transferred separately from the JSON
    #data, arrays which are never collected must be discarded when caches are
    #cleared, and when the mode is enabled the next time.
    import json

    from NCrystalDev._chooks import _get_raw_cfcts
    from NCrystalDev.core import clearCaches
    rawfct = _get_raw_cfcts()
    query = [ 'sab', 'sglcell', '@40.0', '@158.28456094789',
              '@160.60010716438572', '@-2.2701714000978403e-14',
              '@-2.120250108329995e-17', '@0.02833179339754119',
              '@0.027532855700519738', '@0.027532855700519426',
              '@0.026748257342144074', '3', '12462' ]
    def find_keys( data ):
        if isinstance( data, dict ):
            data = list( data.values() )
        if isinstance( data, list ):
            return [ k for e in data for k in find_keys(e) ]
        is_key = ( isinstance( data, str )
                   and data.startswith('__ncrystal__dblarray::') )
        return [ data ] if is_key else []
    def query_without_collecting():
        rawfct['enablejsonarr'](True)
        res = rawfct['jsonquery'](query)
        rawfct['enablejsonarr'](False)
        keys = find_keys( json.loads( res ) )
        assert len(keys) >= 2
        return keys
    print()
    keys = query_without_collecting()
    assert len( rawfct['getjsonarr']( keys[0] ) ) > 0
    print("Collected huge array")
    clearCaches()
    def errmsg( key ):
        return 'Invalid key for JSON Huge Vector: "{}"'.format(
            key.split('::')[1] )
    with ensure_error(NCBadInput,errmsg(keys[1])):
        rawfct['getjsonarr']( keys[1] )
    keys = query_without_collecting()
    rawfct['enablejsonarr'](True)
    rawfct['enablejsonarr'](False)
    with ensure_error(NCBadInput,errmsg(keys[0])):
        rawfct['getjsonarr']( keys[0] )
    res = ncquery( query, huge_arrays = True )
    assert not find_keys( res )
    print("Normal query with huge arrays worked afterwards")


if __name__ == '__main__':
    main()
