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

# Test the NCrystal.browse Python API.

import NCTestUtils.enable_fpe # noqa F401
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

def names( entries ):
    return [ e.display_name for e in entries ]

def nthreads():
    return evaluate_query(['util','factorythreads'])['nthreads']

def lazlau( fmt ):
    #Small .laz/.lau data (high dcutoff):
    from NCrystalDev.mcstasutils import cfgstr_2_hkl
    return '\n'.join( cfgstr_2_hkl( cfgstr = ( 'stdlib::Al_sg225.ncmat;'
                                               'dcutoff=1.0' ),
                                    tgtformat = fmt, verbose = False,
                                    fp_format = '%.8g' ) ) + '\n'

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
        for e in ( nb.find('*.la?',factory='virtual',load=True)
                   + nb.browse('relpath',load=True) ):
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

def main():
    NC.removeAllDataSources()
    NC.enableStandardDataLibrary()
    NC.registerInMemoryFileData('crystal.ncmat',_crystal)
    NC.registerInMemoryFileData('gas.ncmat',_gas)
    NC.registerInMemoryFileData('broken.ncmat','NCMAT v7\n@DENSITY\n -1 g_per_cm3\n')
    #Hides stdlib::Al_sg225.ncmat:
    NC.registerInMemoryFileData('Al_sg225.ncmat',_gas)

    facts = nb.list_factories()
    assert facts['virtual'] == 4 and facts['stdlib'] > 100
    assert 'stdncmat' in nb.list_all_factories()['info']
    print('Properties:',[ n for n,d in nb.physics_props_doc() ])

    print('==> Cheap browsing of virtual factory:')
    for e in nb.browse('virtual'):
        assert e.props is None and e.error is None
        print(f'  {e} name={e.name!r} factory={e.factory!r}'
              f' source={e.source!r} priority={e.priority!r}'
              f' hidden={e.hidden} datatype={e.datatype!r}')
        print(f'     description={e.description!r}')
        print(f'     comments={e.comments!r}')

    print('==> Loaded:')
    for e in nb.browse('virtual',load=True):
        print(f'  {e.display_name}: {e.props}' if e.props
              else f'  {e.display_name}: ERROR {e.error!r}')
    p = [ e for e in nb.browse('virtual',load=True)
          if e.name=='crystal.ncmat' ][0].props
    assert p.sg == 225 and p.as_dict()['natoms'] == 4
    assert p.composition == ((13,((0,1.0),)),)
    with ensure_error(AttributeError,'PhysicsProps has no attribute "foo"'):
        p.foo

    all_entries = nb.browse(load=True)
    al = [ e for e in all_entries if e.name == 'Al_sg225.ncmat' ]
    print('Al_sg225.ncmat entries:',names(al),[e.hidden for e in al])
    assert al[1].props.dyninfo == frozenset(['vdos'])
    #Sorted by priority first:
    assert all_entries[0].factory == 'virtual'

    def find( *a, **kw ):
        res = names( nb.find( *a, **kw ) )
        kwstr = dict( (k,( '<function>' if callable(v) else v ))
                      for k,v in kw.items() )
        print(f'find{a if a else ""}{kwstr if kwstr else ""}:',res)
        return res
    find('crystal')
    find('*.NCMAT',factory='virtual')
    find('stdlib::Be*')
    find('^b.*_sg1[0-9]{2}',regex=True,factory='stdlib')
    find(search='togo',factory='virtual')
    find(search=['togo','small'],factory='virtual')
    find(search=['^ *mentions'],regex=True,factory='virtual')
    find(where="'He' in elements",factory='virtual')
    find(where=['crystal','sg==225','natoms>3'],factory='virtual')
    find(where='sg > 200',factory='virtual')#None comparisons are false
    find(where=lambda p : p.state == 'gas',factory='virtual')
    find('gas',where='absxs < 1',factory='virtual')

    #Sorting and dicts:
    loaded = nb.browse('virtual',load=True)
    print('sorted by absxs:',names(nb.sort_entries(loaded,'absxs')))
    print('sorted by sg (reverse):',
          names(nb.sort_entries(loaded,'sg',reverse=True)))
    print('sorted by name (reverse):',
          names(nb.sort_entries(loaded,'name',reverse=True)))
    sortable = ', '.join( n for n,d in nb.physics_props_doc()
                          if n not in ('debyetemps','msds') )
    with ensure_error(NC.NCBadInput,'Invalid sort key "foo" (must be "name"'
                      ' or one of: ' + sortable + ')'):
        nb.sort_entries(loaded,'foo')
    p = [ e for e in loaded if e.name == 'crystal.ncmat' ][0].props
    assert p.crystalsystem == 'cubic' and p.a == p.b == p.c == 4.04958
    assert p.alpha == 90.0 and abs( p.volume - 4.04958**3 ) < 1e-9
    assert set(p.debyetemps) == set(['Al']) and p.debyetemps['Al'] == 400.0
    assert set(p.msds) == set(['Al']) and p.customsections == frozenset()
    assert abs( p.braggthreshold - 2*4.04958/3**0.5 ) < 1e-9
    assert abs( p.mass - 26.9815 ) < 1e-3 and p.incohxs > 0.0
    g = [ e for e in loaded if e.name == 'gas.ncmat' ][0].props
    assert g.a is None and g.debyetemps is None and g.crystalsystem is None
    assert g.braggthreshold is None
    print('Newer properties OK')
    d = [ e for e in loaded if e.name == 'crystal.ncmat' ][0].as_dict()
    print('as_dict keys:',list(d))
    assert d['path'] is None and all( e.path is None for e in loaded )
    assert d['props']['dyninfo'] == ['vdosdebye'] and d['error'] is None
    assert d['description'] == 'A small Al crystal.'
    d = [ e for e in loaded if e.name == 'broken.ncmat' ][0].as_dict()
    assert d['props'] is None and d['error']

    e = nb.find('crystal')[0]
    print('matching_lines:',e.matching_lines(['togo','A SMALL']))
    print('matching_lines (regex):',e.matching_lines('to+go$',regex=True))
    m = e.load(';temp=100K')
    assert m.info.getTemperature() == 100.0 and m.info.hasStructureInfo()
    assert e.textdata().rawData == _crystal

    progress = []
    nb.browse('stdlib',load=True,progress=lambda a,b : progress.append((a,b)))
    ntot = facts['stdlib']
    assert progress[0] == (0,ntot) and progress[-1] == (ntot,ntot)
    assert len(progress) > 2
    assert [ a for a,b in progress ] == sorted( a for a,b in progress )
    print('Progress reporting OK')

    test_lazlau()

    #Factory threads are only changed temporarily:
    assert nthreads() == 1
    nb.browse('stdlib',load=True,nthreads=4)
    assert nthreads() == 1
    NC.enableFactoryThreads(2)
    nb.browse('stdlib',load=True,nthreads=4)
    assert nthreads() == ( 2 if evaluate_query(['util','factorythreads'])
                           ['threads_available'] else 1 )
    NC.enableFactoryThreads(1)
    print('Factory thread handling OK')

    def bad( msg, *a, **kw ):
        with ensure_error(NC.NCBadInput,msg):
            nb.find( *a, **kw )
    bad('Unknown TextData factory: "nonexistent"',factory='nonexistent')
    bad('Invalid regular expression "(": missing ), unterminated subpattern'
        ' at position 0',search='(',regex=True)
    bad('Unknown name "foo" in where expression "foo > 1"',where='foo > 1')
    bad('Invalid where expression "elements.__class__" (private attributes'
        ' are not allowed)',where='elements.__class__')
    bad('Invalid where expression "absxs >": invalid syntax',where='absxs >')
    bad('Error evaluating where expression "absxs/0 > 1": float division by'
        ' zero',where='absxs/0 > 1',factory='virtual')
    bad('Where conditions require entries loaded with load=True',
        where='absxs > 1',factory='virtual',load=False)

if __name__ == '__main__':
    main()
