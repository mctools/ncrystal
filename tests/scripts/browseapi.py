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

# Test the NCrystal.browse Python API (query_data, DataBrowser, DataEntry,
# PhysicsProps).

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.stabilise_ncpprint # noqa F401
import NCrystalDev as NC
import NCrystalDev.browse as nb
from NCrystalDev.misc import evaluate_query
from NCTestUtils.common import ensure_error

_crystal = """NCMAT v7
#
#   A small Al crystal.
#
#   Mentions Togo.
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
# ------------------
#  Helium gas.
#  Second line.
#
#  Other paragraph.
@STATEOFMATTER
  gas
@DENSITY
  0.5 kg_per_m3
@DYNINFO
  element He
  fraction 1
  type freegas
"""

def nthreads():
    return evaluate_query(['util','factorythreads'])['nthreads']

def show( title, b ):
    print(f'{title}: {b!r} {b.names()}')

def lazlau( fmt ):
    #Small .laz/.lau data (high dcutoff):
    from NCrystalDev.mcstasutils import cfgstr_2_hkl
    return '\n'.join( cfgstr_2_hkl( cfgstr = ( 'stdlib::Al_sg225.ncmat;'
                                               'dcutoff=1.0' ),
                                    tgtformat = fmt, verbose = False,
                                    fp_format = '%.8g' ) ) + '\n'

def test_query_data():
    import NCrystalDev._common as nc_common
    d = nb.query_data('virtual',physics=False)
    assert 'info' not in d[0] and d[0]['name'] == 'Al_sg225.ncmat'
    d = nb.query_data('virtual')
    assert [ e['name'] for e in d ] == ['Al_sg225.ncmat','broken.ncmat',
                                        'crystal.ncmat','gas.ncmat']
    print('query_data("virtual")[2]:')
    nc_common.ncpprint( d[2] )#floats clipped by stabilise_ncpprint
    import json
    assert json.loads( nb.query_data('virtual',as_json=True) ) == d
    progress = []
    nb.query_data('stdlib',progress=lambda a,b : progress.append((a,b)))
    ntot = nb.list_factories()['stdlib']
    assert progress[0] == (0,ntot) and progress[-1] == (ntot,ntot)
    assert len(progress) > 2
    assert [ a for a,b in progress ] == sorted( a for a,b in progress )
    print('query_data OK')

def test_entries():
    b = nb.DataBrowser( factory = 'virtual' )
    print('==> Entries (not loaded):')
    for e in b:
        print(f'  {e} name={e.name!r} factory={e.factory!r}'
              f' source={e.source!r} priority={e.priority!r}'
              f' hidden={e.hidden} datatype={e.datatype!r} path={e.path!r}')
        print(f'     description={e.description!r}')
        print(f'     comments={e.comments!r}')
    print('==> Entries (loaded):')
    for e in b.with_physics():
        print(f'  {e.display_name}: {e.props}' if e.props
              else f'  {e.display_name}: ERROR {e.error!r}')
    lb = b.with_physics()
    p = lb.match('crystal')[0].props
    assert p.sg == 225 and p.as_dict()['natoms'] == 4
    assert p.composition == ((13,((0,1.0),)),)
    assert p.crystalsystem == 'cubic' and p.a == p.b == p.c == 4.04958
    assert p.alpha == 90.0 and abs( p.volume - 4.04958**3 ) < 1e-9
    assert set(p.debyetemps) == {'Al'} and p.debyetemps['Al'] == 400.0
    assert set(p.msds) == {'Al'} and p.customsections == frozenset()
    assert abs( p.braggthreshold - 2*4.04958/3**0.5 ) < 1e-9
    assert abs( p.mass - 26.9815 ) < 1e-3 and p.incohxs > 0.0
    with ensure_error(AttributeError,'PhysicsProps has no attribute "foo"'):
        _ = p.foo
    g = lb.match('gas.ncmat')[0].props
    assert g.a is None and g.debyetemps is None and g.crystalsystem is None
    assert g.braggthreshold is None
    d = lb.match('crystal')[0].as_dict()
    print('as_dict keys:',list(d))
    assert d['props']['dyninfo'] == ['vdosdebye'] and d['error'] is None
    assert d['description'] == 'A small Al crystal.' and d['path'] is None
    d = lb.match('broken')[0].as_dict()
    assert d['props'] is None and d['error']
    e = b.match('crystal')[0]
    print('matching_lines:',e.matching_lines(['togo','A SMALL']))
    print('matching_lines (regex):',e.matching_lines('to+go$',regex=True))
    m = e.load(';temp=100K')
    assert m.info.getTemperature() == 100.0 and m.info.hasStructureInfo()
    assert e.textdata().rawData == _crystal
    print('Entries OK')

def test_selection():
    b = nb.DataBrowser()
    assert b.factories()[0] == 'virtual'#highest priority first
    al = b.match('Al_sg225')
    show('Al_sg225',al)
    assert [ e.hidden for e in al ] == [False,True]
    assert al[1].fullkey == 'stdlib::Al_sg225.ncmat'
    v = b.from_factory('virtual')
    show('virtual',v)
    assert len(v) == 4 and bool(v) and not v.match('nonexistent')
    show('slice',v[1:3])
    assert isinstance( v[0], nb.DataEntry )
    show('*.NCMAT',v.match('*.NCMAT'))
    show('stdlib::Be*',b.match('stdlib::Be*'))
    show('regex',b.from_factory('stdlib').match('^b.*_sg1[0-9]{2}',
                                                regex=True))
    show('two patterns',v.match('crystal','gas'))
    show('search togo',v.search('togo'))
    show('search togo small',v.search('togo','small'))
    show('search regex',v.search('^ *mentions',regex=True))
    show('mixed search',v.search('togo').search('^ *mentions',regex=True))
    show('where He',v.where("'He' in elements"))
    show('where several',v.where('crystal','sg==225','natoms>3'))
    show('where None',v.where('sg > 200'))#None comparisons are false
    show('where fct',v.where(lambda p : p.state == 'gas'))
    show('where debyetemps',v.where('max(debyetemps.values()) > 300'))
    show('filter',v.filter(lambda e : e.name.startswith('g')))
    show('sorted absxs',v.sorted('absxs'))
    show('sorted sg reverse',v.sorted('sg',reverse=True))
    show('sorted name reverse',v.sorted('name',reverse=True))
    show('sorted fct',v.sorted(lambda e : len(e.name)))
    show('chain',v.match('*.ncmat').search('gas').where('absxs < 1')
         .sorted('name',reverse=True))
    #Selections never modify the original:
    assert len(v) == 4 and v.names()[0] == 'Al_sg225.ncmat'
    #Constructor arguments:
    show('constructor',nb.DataBrowser('gas','crystal',factory='virtual',
                                      search='small',where='crystal'))
    assert nb.DataBrowser(factory='virtual',physics=True)[0].props
    print('Selection OK')

def test_output():
    import NCrystalDev._common as nc_common
    v = nb.DataBrowser(factory='virtual')
    print('==> dump():')
    v.dump()
    print('==> dump(comments=True) of search:')
    v.search('togo').dump(comments=True)
    print('==> dump(props=True,linewidth=60):')
    v.match('crystal').dump(props=True,linewidth=60)
    print('==> format_listing(highlight=...,truncate=False):')
    print(v.search('gas').format_listing(highlight=lambda s : f'[{s}]',
                                         truncate=False),end='')
    for fmt in ('text','csv','html'):
        #NB: CSV has full precision, so no derived values like density:
        cols = ( 'formula,sg,absxs,dyninfo,description' if fmt == 'csv'
                 else 'formula,sg,density,dyninfo,description' )
        print(f'==> table(fmt={fmt!r}):')
        print(v.sorted('density').table(cols,fmt=fmt),end='')
    print('==> table(fmt="json"):')
    import json
    nc_common.ncpprint(json.loads(v.table(['formula','density'],fmt='json')))
    assert v.to_csv('sg') == v.table('sg',fmt='csv')
    assert v.to_html('sg') == v.table('sg',fmt='html')
    dicts = v.to_dicts()
    assert json.loads(v.to_json()) == dicts and dicts[2]['props']['sg'] == 225
    print('==> info():')
    v.match('crystal').dump_info()
    print('Suggestions:',v.suggestions('crystl'),v.suggestions('qwertyzzz'))
    assert repr(nb.DataBrowser(factory='virtual').match('crystal')) == (
        'DataBrowser(1 entry from 1 factory)' )
    print('Output OK')

def test_errors():
    def bad( msg, fct ):
        with ensure_error(NC.NCBadInput,msg):
            fct()
    v = nb.DataBrowser(factory='virtual')
    bad('Unknown TextData factory: "nonexistent"',
        lambda : nb.DataBrowser(factory='nonexistent'))
    bad('Unknown TextData factory: "nonexistent"',
        lambda : nb.query_data('nonexistent'))
    bad('Invalid regular expression "(": missing ), unterminated subpattern'
        ' at position 0',lambda : v.search('(',regex=True))
    bad('Unknown name "foo" in where expression "foo > 1"',
        lambda : v.where('foo > 1'))
    bad('Invalid where expression "elements.__class__" (private attributes'
        ' are not allowed)',lambda : v.where('elements.__class__'))
    bad('Invalid where expression "absxs >": invalid syntax',
        lambda : v.where('absxs >'))
    bad('Error evaluating where expression "absxs/0 > 1": division by'
        ' zero',lambda : v.where('absxs/0 > 1'))
    bad('Error evaluating where expression "elements.foo": \'frozenset\''
        ' object has no attribute \'foo\'',lambda : v.where('elements.foo'))
    sortable = ', '.join( n for n,d in nb.physics_props_doc()
                          if n not in ('debyetemps','msds') )
    bad('Invalid sort key "foo" (must be "name" or one of: '+sortable+')',
        lambda : v.sorted('foo'))
    bad('Invalid column "foo" (must be "description" or one of: '
        + ', '.join( n for n,d in nb.physics_props_doc() ) + ')',
        lambda : v.table('sg,foo'))
    bad('Invalid table format: "xml" (must be "text", "csv", "json", or'
        ' "html")',lambda : v.table('sg',fmt='xml'))
    nophys = nb.DataBrowser(factory='virtual',physics=False)
    bad('Physics properties are not available, since the DataBrowser was'
        ' created with physics=False',lambda : nophys.where('crystal'))
    print('Errors OK')

def test_lazy_loading():
    #Physics is loaded on demand, once per factory, shared by derived
    #browsers:
    progress = []
    b = nb.DataBrowser( progress = lambda a,b : progress.append(b) )
    assert not progress
    s1 = b.from_factory('virtual').where('crystal')
    assert set(progress) == {4}#only the virtual factory loaded
    s2 = b.from_factory('virtual').sorted('density')
    s3 = b.match('virtual::gas*').table('formula')
    assert set(progress) == {4} and s1 and s2 and s3
    b.match('stdlib::Al_sg225').where('crystal')
    assert set(progress) == {4,nb.list_factories()['stdlib']}
    #Physics is also loaded on demand when accessed via an entry:
    progress.clear()
    b = nb.DataBrowser( progress = lambda a,b : progress.append(b) )
    e = b.match('virtual::crystal*')[0]
    assert not progress
    assert e.props.sg == 225 and e.error is None and progress
    progress.clear()
    assert b.match('virtual::broken*')[0].error and not progress#cached
    e = nb.DataBrowser( factory = 'virtual', physics = False )[0]
    assert e.props is None and e.error is None
    print('Lazy loading OK')

def test_lazlau():
    #Non-NCMAT data, both in-memory and on-disk (in current directory):
    import pathlib

    from NCTestUtils.common import work_in_tmpdir
    NC.registerInMemoryFileData('mem.laz',lazlau('laz'))
    NC.registerInMemoryFileData('mem.lau',lazlau('lau'))
    NC.enableRelativePaths(True)
    with work_in_tmpdir():
        pathlib.Path('disk.laz').write_text(lazlau('laz'))
        pathlib.Path('disk.lau').write_text(lazlau('lau'))
        b = nb.DataBrowser()
        for e in ( list( b.from_factory('virtual').match('*.la?')
                         .with_physics() )
                   + list( b.from_factory('relpath').with_physics() ) ):
            path = ( e.path and pathlib.Path(e.path).name )#abs path varies
            print(f'  {e.display_name}: datatype={e.datatype}'
                  f' comments={e.comments} descr={e.description!r}'
                  f' path(basename)={path} error={e.error}')
            p = e.props
            print(f'     sg={p.sg} crystalsystem={p.crystalsystem}'
                  f' formula={p.formula} a={p.a:g}'
                  f' braggthreshold={p.braggthreshold:.6g}'
                  f' dyninfo={sorted(p.dyninfo)}')
            assert e.load().info.hasStructureInfo()
    NC.enableRelativePaths(False)
    print('LAZ/LAU data OK')

def test_atomdb():
    import json

    import NCrystalDev._common as nc_common
    #The static database is only queried once (and cached):
    nqueries = [0]
    orig_q = nb._q
    def counting_q( *args ):
        nqueries[0] += ( args == ('atomdb',) )
        return orig_q( *args )
    nb._q = counting_q
    try:
        d = nb.query_atomdb()
        assert json.loads( nb.query_atomdb(as_json=True) ) == d
        nb.AtomDBBrowser('He')
        nb.AtomDBBrowser().where('absxs > 1')
        assert nqueries[0] <= 1
    finally:
        nb._q = orig_q
    #Returned data is a copy (modifying it does not affect the cache):
    d[0]['mass'] = -1.0
    assert nb.query_atomdb()[0]['mass'] > 0.0
    print('atomdb fields:',[ n for n,_ in nb.atomdb_fields_doc() ])
    a = nb.AtomDBBrowser()
    assert len(a) == len(d) and repr(a) == f'AtomDBBrowser({len(d)} entries)'
    print('He + b10 (case-insensitive):',nb.AtomDBBrowser('He','b10').names())
    print('glob:',a.match('Li*').names())
    print('where:',a.where('absxs > 1000 and natural')
          .sorted('absxs',reverse=True).names())
    print('where fct:',a.where(lambda e : e['z'] == 1).names())
    print('sorted fct:',a.match('H').sorted(lambda e : -e['a']).names())
    print('slice:',a.match('H')[1:3].names(),'item:',a.match('Al')[0]['label'])
    sel = nb.AtomDBBrowser('He','B10')
    print('==> text table:')
    sel.dump()
    print('==> text table (columns):')
    print(sel.table('mass,absxs'),end='')
    print('==> csv:')
    print(sel.to_csv('a,cohsl,incohxs,absxs'),end='')#no mass: FP
    print('==> html:')
    print(sel.to_html('natural,absxs'),end='')
    print('==> json:')
    nc_common.ncpprint( json.loads( sel.to_json('element,cohxs') ) )
    assert sel.to_dicts()[1]['label'] == 'He3'
    def bad( msg, fct ):
        with ensure_error(NC.NCBadInput,msg):
            fct()
    names = ', '.join( n for n,_ in nb.atomdb_fields_doc() )
    bad('Unknown name "foo" in where expression "foo > 1"',
        lambda : a.where('foo > 1'))
    bad('Invalid sort key "foo" (must be one of: '+names+')',
        lambda : a.sorted('foo'))
    bad('Invalid column "label" (must be one of: '
        + names.replace('label, ','') + ')', lambda : a.table('label'))
    bad('Invalid table format: "xml" (must be "text", "csv", "json", or'
        ' "html")',lambda : a.table(fmt='xml'))
    print('AtomDB OK')

def test_threads():
    #Factory threads are only changed temporarily (NB: must be last, since
    #it configures factory threads):
    assert nthreads() == 1
    nb.DataBrowser(factory='stdlib',nthreads=4,physics=True)
    assert nthreads() == 1
    NC.enableFactoryThreads(2)
    nb.DataBrowser(factory='stdlib',nthreads=4,physics=True)
    assert nthreads() == ( 2 if evaluate_query(['util','factorythreads'])
                           ['threads_available'] else 1 )
    NC.enableFactoryThreads(1)
    print('Factory thread handling OK')

def main():
    NC.removeAllDataSources()
    NC.enableStandardDataLibrary()
    NC.registerInMemoryFileData('crystal.ncmat',_crystal)
    NC.registerInMemoryFileData('gas.ncmat',_gas)
    NC.registerInMemoryFileData('broken.ncmat',
                                'NCMAT v7\n@DENSITY\n -1 g_per_cm3\n')
    #Hides stdlib::Al_sg225.ncmat:
    NC.registerInMemoryFileData('Al_sg225.ncmat',_gas)

    facts = nb.list_factories()
    assert facts['virtual'] == 4 and facts['stdlib'] > 100
    assert 'stdncmat' in nb.list_all_factories()['info']
    print('Properties:',[ n for n,d in nb.physics_props_doc() ])
    test_query_data()
    test_entries()
    test_selection()
    test_output()
    test_errors()
    test_lazy_loading()
    test_lazlau()
    test_atomdb()
    test_threads()

if __name__ == '__main__':
    main()
