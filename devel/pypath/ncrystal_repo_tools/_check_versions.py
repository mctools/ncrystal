
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

def get_versionfile_version( path ):
    assert path.name=='VERSION'
    return path.read_text().strip()

def get_py_version( path ):
    import ast
    t = ( path ).read_text()
    for line in t.splitlines():
        if line.startswith('__version__'):
            var,val = line.split('=')
            assert var.rstrip()=='__version__'
            return ast.literal_eval(val.strip())

def get_toml_version( path ):
    from .toml import parse_toml
    assert path.name == 'pyproject.toml'
    return parse_toml( path )['project']['version']

def check_ncrystal_metapkg_version( path ):
    from .toml import parse_toml
    assert path.name == 'pyproject.toml'
    t = parse_toml( path )
    d = t['project']['dependencies']
    v = t['project']['version']
    for pn in 'ncrystal-core','ncrystal-python':
        e = f'{pn}=={v}'
        if e not in d:
            raise SystemExit(f'Missing pinned dependency on {e} in {path}')
    return v

def get_cmake_version( path ):
    assert path.name == 'CMakeLists.txt'
    lv = [ line for line in path.read_text().splitlines()
           if line.startswith('project(' ) ]
    assert len(lv)==1
    assert lv[0].startswith('project( NCrystal VERSION ')
    return lv[0][len('project( NCrystal VERSION '):].split()[0]

def check_versions():
    from .dirs import reporoot

    versions = [
        ( get_versionfile_version, 'VERSION' ),
        ( get_py_version, 'ncrystal_python/NCrystal/__init__.py' ),
        ( get_toml_version, 'pyproject.toml' ),
        ( get_toml_version, 'ncrystal_core/pyproject.toml' ),
        ( check_ncrystal_metapkg_version, 'ncrystal_metapkg/pyproject.toml' ),
        ( get_cmake_version, 'ncrystal_core/CMakeLists.txt' ),
    ]

    versions_found = []
    subpathmaxlen = max(len(subpath) for fct, subpath in versions)

    print("Extracted versions:")
    for fct, subpath in versions:
        v = fct( reporoot / subpath )
        print(f"   {subpath.ljust(subpathmaxlen)} : {v}")
        versions_found.append(v)

    if len(set(versions_found)) == 1:
        print("All OK!")
    else:
        raise SystemExit('ERROR: Inconsistencies detected')

def main():
    check_versions()

if __name__=='__main__':
    main()
