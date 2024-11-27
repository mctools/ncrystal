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

import sys
import pathlib
sys.path.insert(0,str(pathlib.Path(__file__).parent/'../pypath'))
import common

if not common.can_parse_toml():
    #Acceptable simplification to simply abort, since already we mostly have
    #python >=3.11 in dev envs, and CI certainly will catch it with python
    #>=3.11 somewhere.
    print("Silent abort: tomli not present and python < 3.11")
    raise SystemExit(0)

def load_data(subpath):
    d = common.parse_toml(common.reporoot / subpath)
    #sneak in the source location:
    assert '__file__' not in d
    d['__srcloc__'] = subpath
    print(f"Loaded: {describe(d)}")
    return d

def describe( data ):
    return '<root>/%s'%data['__srcloc__']

data_monolith = load_data( 'pyproject.toml' )
data_core = load_data( 'ncrystal_core/pyproject.toml' )
data_py = load_data( 'ncrystal_python/pyproject.toml' )

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

    #Check that there are no unexpected sections:
    toplvlkeys = set(['build-system','project','tool','__srcloc__'])
    assert toplvlkeys == set( data_monolith.keys() )
    assert toplvlkeys == set( data_core.keys() )
    assert toplvlkeys == set( data_py.keys() )

    #Check 'project' section
    assert data_monolith['project']['version']==data_core['project']['version']
    projkeys_monolith = set( data_monolith['project'].keys() )
    projkeys_core = set( data_core['project'].keys() )
    projkeys_py = set( data_py['project'].keys() )
    assert 'dependencies' in projkeys_monolith
    assert 'dependencies' not in projkeys_core
    assert 'dependencies' in projkeys_py
    assert ( projkeys_monolith - projkeys_core ) == set(['dependencies'])
    assert ( projkeys_core - projkeys_monolith ) == set([])
    assert ( projkeys_py - projkeys_monolith ) == set(['dynamic'])
    assert ( projkeys_monolith - projkeys_py ) == set(['version'])
    cmp_common_entries( 'project', data_monolith, data_core,
                        allow_diff = ['name','scripts','description'] )
    cmp_common_entries( 'project', data_monolith, data_py,
                        allow_diff = ['name','scripts','description'] )
    #NB: projects.scripts checked elsewhere in more detail!

    #Check 'build-system' section
    cmp_common_entries( 'build-system', data_monolith, data_core,
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
    extra_mono = [
        ('ncrystal-config',
         "_ncrystal_core_monolithic.info:_ncrystal_config_cli_wrapper"),
    ]
    extra_core = [
        ('ncrystal-config',
         "_ncrystal_core.info:_ncrystal_config_cli_wrapper"),
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
    if not ok:
        raise SystemExit('Failures detected in project.scripts')


def _check_project_scripts_impl( data, *, cli_scripts, extra ):
    #We already check project.scripts consistency between the three
    #pyproject.toml files in check_metadata(), so here we simply have to verify
    #that ncrystal_python/project.toml has exactly 1 entry points for each of
    #the actual NCrystal/_cli_*.py modules.
    from_toml = data['project']['scripts']
    expected = dict( e for e in extra )
    if cli_scripts:
        for f in (common.reporoot
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
