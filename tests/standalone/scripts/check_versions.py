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

def get_py_version( path ):
    import ast
    t = ( path ).read_text()
    for line in t.splitlines():
        if line.startswith('__version__'):
            var,val = line.split('=')
            assert var.rstrip()=='__version__'
            return ast.literal_eval(val.strip())

def get_toml_version( path ):
    assert path.name == 'pyproject.toml'
    return common.parse_toml( path )['project']['version']

def check_versions():
    root = common.reporoot
    versions = [
        ( get_py_version, 'ncrystal_python/NCrystal/__init__.py' ),
        ( get_toml_version, 'pyproject.toml' ),
        ( get_toml_version, 'ncrystal_core/pyproject.toml' ),
    ]

    versions_found = []
    subpathmaxlen = max(len(subpath) for fct, subpath in versions)

    print("Extracted versions:")
    for fct, subpath in versions:
        v = fct( root / subpath )
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

if __name__=='__main__':
    check_versions()

