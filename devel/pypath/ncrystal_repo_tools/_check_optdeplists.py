
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

def get_master_source():
    from .dirs import reporoot
    from .toml import parse_toml
    frel = 'ncrystal_metapkg/pyproject.toml'
    t = parse_toml( reporoot.joinpath(frel) )
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

def main():
    src = get_master_source()
    #Fixme: actually produce these files (+conda files) in devel/envs:
    for k,v in src.items():
        print(k)
        print(produce_pipreqfile(v))

if __name__=='__main__':
    main()
