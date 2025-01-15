
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

#Check that the optional dependencies in ncrystal_metapkg/pyproject.toml are
#sensible, and that the devel/reqs/requirements_*.txt and devel/reqs/conda_*.yml
#files are always kept in synch with what is specified in
#ncrystal_metapkg/pyproject.toml.

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
    return deplist, '\n'.join(sorted(set(deplist)))+'\n'

def _check_files( files, pattern  ):
    from .dirs import reporoot
    def print_expected_content(frel):
        print('-------- Expected content ------------')
        print(files[frel][1],end='')
        print('--------------------------------------')
    found = set()
    for f in reporoot.joinpath('devel','reqs').glob(pattern):
        frel = str(f.relative_to(reporoot)).replace('\\','/')
        if frel not in files:
            raise SystemExit(f'ERROR: Excess file: {frel}')
        if not files[frel][1]==f.read_text():
            print(f'ERROR: Unexpected content in file: {frel}')
            print_expected_content(frel)
            raise SystemExit(1)
        found.add(frel)
    for frel in files:
        if frel not in found:
            print(f'ERROR: Missing file: {frel}')
            print_expected_content(frel)
            raise SystemExit(1)
        print(f'  Verified {frel}')

def check_pipreqfiles(src):
    all_deps_combined = set()
    base_deps = get_base_deps()
    files = {}
    for k,v in src.items():
        frel = f'devel/reqs/requirement_{k}.txt'
        assert frel not in files
        deps = set(v + list(base_deps))
        all_deps_combined.update(deps)
        files[frel] = produce_pipreqfile( deps )
    frel_all_combined = 'devel/reqs/requirement_all_combined.txt'
    assert frel_all_combined not in files
    files[frel_all_combined] = produce_pipreqfile( all_deps_combined )
    _check_files(files,'requirement_*.txt')
    return files

def pip_2_conda_dep( depname ):
    #NB: I am not sure that all version specifications are the same in
    #pyproject.toml, requirements.txt and conda.yml. If not, we might be able to
    #salvage something here.
    return depname

def produce_conda_env( name, deplist ):
    res=f"""name: ncrystal_{name}
channels:
  - nodefaults
  - conda-forge
dependencies:
  - python
  - compilers
  - cmake
  - numpy
  - pip
"""
    for n in sorted(deplist):
        res += '  - %s\n'%pip_2_conda_dep(n)
    return res

def check_condaenvfiles(pipfiles):
    files = {}
    for frel, (deplist,_) in pipfiles.items():
        name = frel.split('/')[-1].split('.')[0].split('_',1)[1]
        frel = f'devel/reqs/conda_{name}.yml'
        files[frel] = ( deplist, produce_conda_env( name, deplist ) )
    _check_files(files,'conda_*.yml')
    return files

def main():
    src = get_master_source()
    files = check_pipreqfiles(src)
    check_condaenvfiles(files)

if __name__=='__main__':
    main()
