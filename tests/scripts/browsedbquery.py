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

# Test the ["util","browsedb",...] and ["util","browsefactories"] JSON
# queries.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query
from NCTestUtils.common import ensure_error
import pprint

_crystal = """NCMAT v7
#
#   A small Al crystal.
#
@CELL
 cubic 4.04958
@SPACEGROUP
 225
@ATOMPOSITIONS
 Al 0 1/2 1/2
 Al 0 0 0
 Al 1/2 1/2 0
 Al 1/2 0 1/2
@DEBYETEMPERATURE
 Al 400
"""

_gas = """NCMAT v7
#
#
@STATEOFMATTER
  gas
@DENSITY
  0.5 kg_per_m3
@DYNINFO
  element He
  fraction 1
  type freegas
"""

_broken = """NCMAT v7
# Broken file
@DENSITY
  -1 g_per_cm3
"""

def q( *args ):
    return evaluate_query( ['util'] + list(args) )

def rounded( x ):
    #Round floats, for reproducibility across platforms:
    if isinstance( x, float ):
        return float( '%.10g'%x )
    if isinstance( x, list ):
        return [ rounded(e) for e in x ]
    if isinstance( x, dict ):
        return dict( (k,rounded(v)) for k,v in x.items() )
    return x

def show( title, data ):
    print(f'==> {title}:')
    pprint.pp( rounded(data), width = 78 )

def run_scenario( name, code, envval = None ):
    #Run code in fresh process (the NCRYSTAL_FACTORY_THREADS env var is
    #only read once per process), printing the factorythreads state:
    import subprocess
    import sys

    from NCTestUtils.env import ncsetenv
    full = ( 'import NCrystalDev as NC\n'
             'from NCrystalDev.misc import evaluate_query as q\n'
             'def ft():\n'
             '    d = q(["util","factorythreads"])\n'
             '    return d["nthreads"], d["user_configured"]\n'
             'def load(n):\n'
             '    q(["util","browsedb","stdlib","nthreads=%i"%n])\n'
             + code )
    ncsetenv('FACTORY_THREADS',envval)
    try:
        rv = subprocess.run( [ sys.executable, '-c', full ],
                             capture_output = True, text = True,
                             check = True )
    finally:
        ncsetenv('FACTORY_THREADS',None)
    return rv.stdout.strip()

def test_factorythreads():
    ft = q('factorythreads')
    assert set(ft) == set(['threads_available','nthreads','user_configured'])
    avail = ft['threads_available']
    def n( nthreads ):
        return str( nthreads if avail else 1 )
    #Nothing configured: loading with temporary threads leaves no trace, but
    #an explicit disabling is respected:
    res = run_scenario('A','print(ft()); load(4); print(ft());'
                       ' NC.enableFactoryThreads(1); load(4); print(ft())')
    assert res.split() == ['(1,','False)','(1,','False)','(1,','True)']
    #Env var is respected (also 0, which means disabled):
    res = run_scenario('B','print(ft()); load(8); print(ft())','3')
    assert res.split() == [ f'({n(3)},', 'True)' ]*2
    res = run_scenario('C','print(ft()); load(8); print(ft())','0')
    assert res.split() == [ '(1,', 'True)' ]*2
    #Explicit calls before the env var is read take precedence:
    res = run_scenario('D','NC.enableFactoryThreads(2); print(ft())','3')
    assert res.split() == [ f'({n(2)},', 'True)' ]
    print('Factory thread queries OK')

def main():
    NC.removeAllDataSources()
    NC.enableStandardDataLibrary()
    NC.registerInMemoryFileData('crystal.ncmat',_crystal)
    NC.registerInMemoryFileData('gas.ncmat',_gas)
    NC.registerInMemoryFileData('broken.ncmat',_broken)
    NC.registerInMemoryFileData('notncmat.laz','whatever')
    #Hides stdlib::Al_sg225.ncmat:
    NC.registerInMemoryFileData('Al_sg225.ncmat',_gas)

    assert 'browsedb' in q('list') and 'browsefactories' in q('list')

    #Overviews (contents otherwise depend on build and environment):
    counts = q('browsedb')
    assert counts['virtual'] == 5 and counts['stdlib'] > 100
    facts = q('browsefactories')
    assert set(facts) == set(['textdata','info','scatter','absorption'])
    assert 'virtual' in facts['textdata'] and 'stdlib' in facts['textdata']
    assert 'stdncmat' in facts['info']
    print('Overview queries OK')

    full = q('browsedb','virtual')
    show('browsedb virtual',full)
    cheap = q('browsedb','virtual','cheap')
    show('browsedb virtual cheap',cheap)
    #Cheap mode is the same, minus info and load errors:
    for e_full, e_cheap in zip(full,cheap):
        e = dict( (k,v) for k,v in e_full.items()
                  if k not in ('info','error') )
        assert e == e_cheap

    al = [ e for e in q('browsedb','stdlib') if e['name']=='Al_sg225.ncmat' ]
    assert len(al) == 1
    al = al[0]
    al['source'] = '<stdlib-location>'
    al['comments'] = al['comments'][0:3]
    show('stdlib::Al_sg225.ncmat (source replaced, first comments only)',al)

    #Chunks concatenate to the full list:
    for fact in ('virtual','stdlib'):
        ref = q('browsedb',fact,'cheap')
        for n in (1,2,3,7,200):
            chunks = [ q('browsedb',fact,str(i),str(n),'cheap')
                       for i in range(n) ]
            assert [ e for c in chunks for e in c ] == ref
            assert max(len(c) for c in chunks) - min(len(c) for c in chunks) <= 1
    assert q('browsedb','virtual','1','2') == full[2:]
    print('Chunked queries OK')

    #Loading with factory threads gives identical results:
    ref = q('browsedb','stdlib')
    NC.enableFactoryThreads(4)
    try:
        assert q('browsedb','stdlib') == ref
    finally:
        NC.enableFactoryThreads(1)
    print('Loading with factory threads OK')

    test_factorythreads()

    def bad( msg, *args ):
        with ensure_error(NC.NCBadInput,msg):
            q(*args)
    bad('Unknown TextData factory in browsedb query: "nonexistent"',
        'browsedb','nonexistent')
    badchunk = ('Invalid chunk specification in browsedb query'
                ' (must have 0<=I<N): ')
    bad(badchunk+'I=2, N=2','browsedb','virtual','2','2')
    bad(badchunk+'I=0, N=0','browsedb','virtual','0','0')
    bad(('Invalid browsedb query (usage: ["util","browsedb",FACTNAME,'
         '(I,N,)("cheap",)("nthreads=N")])'),'browsedb','virtual','1')
    bad('Invalid nthreads in browsedb query: "x"',
        'browsedb','virtual','nthreads=x')
    bad('Invalid chunk index I in browsedb query: "a"',
        'browsedb','virtual','a','2')
    bad(('Invalid util query: ["util","browsefactories","virtual"] (no'
         ' arguments should come after: ["util","browsefactories"])'),
        'browsefactories','virtual')

if __name__ == '__main__':
    main()
