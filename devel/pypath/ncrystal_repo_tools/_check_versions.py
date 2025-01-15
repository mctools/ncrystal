
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

def get_testcfglog_version( path ):
    #Look for single line like "NCrystal::getVersion() = 3009080"
    assert path.name == 'test.log'
    assert path.parent.name == 'app_cfg'
    lv = [ line for line in path.read_text().splitlines()
           if line.startswith('NCrystal::getVersion() = ') ]
    assert len(lv)==1
    v = int( lv[0][len('NCrystal::getVersion() = '):].strip() )
    vmajor = v // 1000000
    vminor = (v%1000000)//1000
    vpatch = v%1000
    return f'{vmajor}.{vminor}.{vpatch}'

def get_generictestlog_version( path ):
    #Look for strings like "by NCrystal v3.9.80"
    assert path.name.endswith('.log')
    versions = set()
    for line in path.read_text().splitlines():
        if 'by NCrystal v' not in line:
            continue
        versions.add(line.split('by NCrystal v',1)[1].split()[0])
    if len(versions)>1:
        raise SystemExit(f'Multiple versions found in {path}: {versions}')
    if not versions:
        raise SystemExit(f'No versions found in {path}')
    return versions.pop()

def get_testtoollog_version( path ):
    #Look for single string like "with NCrystal (v3.9.80)"
    assert path.name == 'nctool.log'
    lv = [ line for line in path.read_text().splitlines()
           if 'with NCrystal (v' in line ]
    assert len(lv)==1
    return lv[0].split('with NCrystal (v',1)[1].split(')',1)[0]

def get_cmake_version( path ):
    assert path.name == 'CMakeLists.txt'
    lv = [ line for line in path.read_text().splitlines()
           if line.startswith('project(' ) ]
    assert len(lv)==1
    assert lv[0].startswith('project( NCrystal VERSION ')
    return lv[0][len('project( NCrystal VERSION '):].split()[0]

#VERSION:3.9.80

def check_versions():
    from .dirs import reporoot

    versions = [
        ( get_versionfile_version, 'VERSION' ),
        ( get_py_version, 'ncrystal_python/NCrystal/__init__.py' ),
        ( get_toml_version, 'pyproject.toml' ),
        ( get_toml_version, 'ncrystal_core/pyproject.toml' ),
        ( get_cmake_version, 'ncrystal_core/CMakeLists.txt' ),
        ( get_testcfglog_version, 'tests/src/app_cfg/test.log' ),
        ( get_testtoollog_version, 'tests/scripts/nctool.log' ),
        ( get_generictestlog_version, 'tests/scripts/mcstasunion.log' ),
        ( get_generictestlog_version, 'tests/scripts/ncmat2hkl.log' ),
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

    #TODO: ncrystal_core/CMakeLists.txt
    #TODO: ncrystal_core/include/NCrystal/NCVersion.hh
    #TODO: ncrystal_core/include/NCrystal/ncrystal.h
    #TODO: Git describe! (if .git present)

def main():
    check_versions()

if __name__=='__main__':
    main()
