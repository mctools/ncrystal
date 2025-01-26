
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

from .dirs import reporoot

toml_files =  dict( monolith = 'pyproject.toml',
                    core = 'ncrystal_core/pyproject.toml',
                    coreempty = 'ncrystal_core/empty_pypkg/pyproject.toml',
                    py = 'ncrystal_python/pyproject.toml',
                    plugmgr = 'ncrystal_pypluginmgr/pyproject.toml',
                    meta = 'ncrystal_metapkg/pyproject.toml' )

_data_cache = {}
def load_data( key ):
    if key not in _data_cache:
        _data_cache[key] = _actual_load_data( toml_files[key] )
    return _data_cache[key]

def _actual_load_data(subpath):
    from .toml import parse_toml, can_parse_toml

    if not can_parse_toml():
        #Acceptable simplification to simply abort, since already we mostly have
        #python >=3.11 in dev envs, and CI certainly will catch it with python
        #>=3.11 somewhere.
        print("Silent abort: tomli not present and python < 3.11")
        raise SystemExit(0)

    d = parse_toml(reporoot / subpath)
    #sneak in the source location:
    assert '__file__' not in d
    d['__srcloc__'] = subpath
    print(f"Loaded: {describe(d)}")
    return d

def describe( data ):
    return '<root>/%s'%data['__srcloc__']

def cmp_common_entries( keypath, dict1, dict2, allow_diff = [] ):
    d1, d2 = dict1, dict2
    for k in keypath.split('.'):
        d1, d2 = d1[k], d2[k]

    keys = set(d1.keys()).intersection(set(d2.keys())) - set(allow_diff)
    ok = True
    for k in keys:
        if d1[k] != d2[k]:
            print()
            print(f'Inconsistency found in key "{keypath}.{k}:"')
            print()
            print(f'   {describe(dict1)} has value: {repr(d1[k])}')
            print()
            print(f'   {describe(dict2)} has value: {repr(d2[k])}')
            print()
            ok = False
    if not ok:
        raise SystemExit(1)

def check_metadata():
    #Check that there is consistency between the three pyproject.toml files,
    #except where expected.

    data_monolith = load_data( 'monolith' )
    data_core = load_data( 'core' )
    data_coreempty = load_data( 'coreempty' )
    data_py = load_data( 'py' )
    data_meta = load_data( 'meta' )
    data_plugmgr = load_data( 'plugmgr' )

    #Check that there are no unexpected sections:
    toplvlkeys_notool = set(['build-system','project','__srcloc__'])
    toplvlkeys = toplvlkeys_notool.union( set(['tool']) )
    assert toplvlkeys == set( data_monolith.keys() )
    assert toplvlkeys == set( data_core.keys() )
    assert toplvlkeys_notool == set( data_coreempty.keys() )
    assert toplvlkeys == set( data_py.keys() )
    assert toplvlkeys_notool == set( data_meta.keys() )
    assert toplvlkeys == set( data_plugmgr.keys() )

    #Check 'project' section
    version = data_core['project']['version']
    assert data_coreempty['project']['version']==version
    assert data_monolith['project']['version']==version
    assert data_meta['project']['version']==version
    projkeys_monolith = set( data_monolith['project'].keys() )
    projkeys_core = set( data_core['project'].keys() )
    projkeys_coreempty = set( data_coreempty['project'].keys() )
    projkeys_py = set( data_py['project'].keys() )
    projkeys_meta = set( data_meta['project'].keys() )
    projkeys_plugmgr = set( data_plugmgr['project'].keys() )
    assert 'dependencies' in projkeys_monolith
    assert 'dependencies' not in projkeys_core
    assert 'dependencies' not in projkeys_coreempty
    assert 'dependencies' in projkeys_py
    assert 'dependencies' in projkeys_meta
    assert 'dependencies' in projkeys_plugmgr
    assert data_plugmgr['project']['dependencies']==[]

    assert ( set(data_meta['project']['dependencies'])
             == set([f'ncrystal-core=={version}',
                     f'ncrystal-python=={version}']) )
    assert ( projkeys_monolith - projkeys_core ) == set(['dependencies'])
    assert ( projkeys_core - projkeys_monolith ) == set([])
    assert ( projkeys_core - projkeys_meta ) == set(['scripts'])
    assert ( projkeys_meta - projkeys_core ) == set(['dependencies',
                                                     'optional-dependencies'])
    assert ( projkeys_py - projkeys_monolith ) == set(['dynamic'])
    assert ( projkeys_monolith - projkeys_py ) == set(['version'])
    assert ( projkeys_core - projkeys_meta ) == set(['scripts'])
    assert ( projkeys_coreempty - projkeys_core ) == set([])
    assert ( projkeys_core - projkeys_coreempty ) == set(['scripts'])

    cmp_common_entries( 'project', data_coreempty, data_core )
    cmp_common_entries( 'build-system', data_coreempty, data_py )
    assert ( data_coreempty['build-system'] == data_py['build-system'] )

    cmp_common_entries( 'project', data_coreempty, data_core,
                        allow_diff = [] )

    cmp_common_entries( 'project', data_meta, data_core,
                        allow_diff = ['name'] )
    cmp_common_entries( 'project', data_py, data_plugmgr,
                        allow_diff = ['name','description',
                                      'dependencies','scripts'] )
    cmp_common_entries( 'project', data_monolith, data_core,
                        allow_diff = ['name','scripts','description'] )
    cmp_common_entries( 'project', data_monolith, data_py,
                        allow_diff = ['name','scripts','description'] )
    #NB: projects.scripts checked elsewhere in more detail!

    #Check 'build-system' section
    cmp_common_entries( 'build-system', data_monolith, data_core,
                        allow_diff = [] )
    cmp_common_entries( 'build-system', data_py, data_meta,
                        allow_diff = [] )
    cmp_common_entries( 'build-system', data_py, data_plugmgr,
                        allow_diff = [] )

    #Check 'tool' section
    assert set( data_core['tool'].keys() ) == set(['scikit-build'])
    assert set( data_monolith['tool'].keys() ) == set(['scikit-build'])
    cmp_common_entries( 'tool.scikit-build', data_monolith, data_core,
                        allow_diff = ['sdist','wheel'] )
    cmp_common_entries( 'tool.scikit-build.sdist', data_monolith, data_core,
                        allow_diff = ['include'] )
    cmp_common_entries( 'tool.scikit-build.wheel', data_monolith, data_core,
                        allow_diff = ['packages'] )


def check_all_project_scripts():
    data_monolith = load_data( 'monolith' )
    data_core = load_data( 'core' )
    data_py = load_data( 'py' )
    data_meta = load_data( 'meta' )
    data_plugmgr = load_data( 'plugmgr' )

    extra_mono = [
        ('ncrystal-config',
         "_ncrystal_core_monolithic.info:_ncrystal_config_cli_wrapper"),
    ]
    extra_core = [
        ('ncrystal-config',
         "_ncrystal_core.info:_ncrystal_config_cli_wrapper"),
    ]
    extra_plugmgr = [
        ('ncrystal-pluginmanager',
         "_ncrystal_pypluginmgr.cli:main"),
    ]
    ok = True
    if not _check_project_scripts_impl( data_py,
                                        cli_scripts = True,
                                        extra=[] ):
        ok = False
    if not _check_project_scripts_impl( data_monolith,
                                        cli_scripts = True,
                                        extra = extra_mono ):
        ok = False
    if not _check_project_scripts_impl( data_core,
                                        cli_scripts = False,
                                        extra = extra_core ):
        ok = False
    if not _check_project_scripts_impl( data_meta,
                                        cli_scripts = False,
                                        extra = [] ):
        ok = False
    if not _check_project_scripts_impl( data_plugmgr,
                                        cli_scripts = False,
                                        extra = extra_plugmgr ):
        ok = False
    if not ok:
        raise SystemExit('Failures detected in project.scripts')


def _check_project_scripts_impl( data, *, cli_scripts, extra ):
    #We already check project.scripts consistency between the three
    #pyproject.toml files in check_metadata(), so here we simply have to verify
    #that ncrystal_python/project.toml has exactly 1 entry points for each of
    #the actual NCrystal/_cli_*.py modules.
    from_toml = data['project'].get('scripts',{})
    expected = dict( e for e in extra )
    if cli_scripts:
        for f in (reporoot
                  / 'ncrystal_python/NCrystal').glob('_cli_*.py'):
            bn = f.name[:-3]
            n = bn[5:]
            n = f'ncrystal_{n}' if n!='nctool' else n
            expected[n] = f'NCrystal.{bn}:main'
    if expected == from_toml:
        return True#all ok
    descr = describe( data )
    for k,v in expected.items():
        fmt = f'{k} = "{v}"'
        v_toml = from_toml.get(k)
        if v != v_toml:
            print(f"ERROR: Missing line in {descr}:\n\n   {fmt}\n")
            if v_toml is not None:
                print(f'Got instead:\n\n   {k} = "{v_toml}"\n')
    for k,v in from_toml.items():
        if k not in expected:
            print(f'ERROR: Unexpected line in {descr}:\n\n   {k}="{v}"\n')
    return False

def main():
    check_all_project_scripts()
    check_metadata()
    print("All OK!")

if __name__=='__main__':
    main()
