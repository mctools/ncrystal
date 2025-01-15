
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

#[project.optional-dependencies]
#composer = [ 'spglib' ]
#cif = [ 'ase<3.24.0', 'gemmi', 'spglib' ]
#plot = [ 'matplotlib' ]
#all = [ 'ase<3.24.0', 'gemmi', 'spglib', 'matplotlib' ]
#devel = [ 'simple-build-system', 'ruff', 'cppcheck', 'mpmath', 'tomli', 'pybind11' ]

def get_base_deps():
    from .dirs import reporoot
    from .toml import parse_toml
    t1 = parse_toml( reporoot.joinpath('ncrystal_python/pyproject.toml') )
    t2 = parse_toml( reporoot.joinpath('ncrystal_python/pyproject.toml') )
    return set( t1['project']['dependencies'] + t2['project']['dependencies'])

def get_master_source():
    from .dirs import reporoot
    from .toml import parse_toml
    frel = 'ncrystal_metapkg/pyproject.toml'
    t = parse_toml( reporoot.joinpath(frel) )
    assert len(t['project']['dependencies'])==2#just ncrystal-* pkgs
    t = t['project']['optional-dependencies']

    #Make sure the same dependency is listed with same version reqs everywhere:
    seen = {}
    for k,v in t.items():
        for e in v:
            bn = e.split('<')[0].split('=')[0].split('>')[0]
            assert bn.replace('-','_').isidentifier(), f"basename issues: {bn}"
            if bn in seen:
                assert e == seen[bn], f"same dep ({bn}) seen with multiple reqs"
            else:
                seen[bn] = e

    #Verify that 'all' contains all the rest (except what is in 'devel'), and
    #that nothing in devel is in 'all':
    assert 'devel' in t
    assert 'all' in t
    assert set(t['devel']).isdisjoint(set(t['all']))
    deps_expected_all = set()
    for k,v in t.items():
        if k in ('devel','all'):
            continue
        deps_expected_all.update(set(v))
    assert set(t['all']) == deps_expected_all
    print(f"  Verified internal consistency of {frel}")
    return t

def produce_pipreqfile( deplist ):
    return '\n'.join(sorted(set(deplist)))+'\n'

def check_pipreqfiles(src):
    #FIXME: Missing global deps like numpy, which we should dig out of ncrystal_core+ncrystal_python. Also, add all_plus_devel with both combined?
    base_deps = get_base_deps()
    from .dirs import reporoot
    files = {}
    for k,v in src.items():
        frel = f'devel/reqs/requirement_{k}.txt'
        assert not frel in files
        files[frel] = produce_pipreqfile( set(v + list(base_deps)) )
    def print_expected_content(frel):
        print('-------- Expected content ------------')
        print(files[frel],end='')
        print('--------------------------------------')
    found = set()
    for f in reporoot.joinpath('devel','reqs').glob('requirement_*.txt'):
        frel = str(f.relative_to(reporoot))
        if not frel in files:
            raise SystemExit(f'ERROR: Excess file: {frel}')
        if not files[frel]==f.read_text():
            print(f'ERROR: Unexpected content in file: {frel}')
            print_expected_content(frel)
            raise SystemExit(1)
        found.add(frel)
    for frel in files:
        if not frel in found:
            print(f'ERROR: Missing file: {frel}')
            print_expected_content(frel)
            raise SystemExit(1)
        print(f'  Verified {frel}')

def main():
    src = get_master_source()
    #Fixme: actually produce these files (+conda files) in devel/envs:
    check_pipreqfiles(src)

if __name__=='__main__':
    main()
