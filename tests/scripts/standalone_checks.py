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

# NEEDS: toml

def run( script_file ):
    import sys
    import subprocess
    rv = subprocess.run( [ sys.executable or 'python3',
                           '-BI',
                           script_file.resolve().absolute() ] )
    return rv.returncode

def run_all():
    import pathlib
    testdir = pathlib.Path(__file__).resolve().absolute().parent.parent
    scripts =sorted(testdir.joinpath('standalone','scripts').glob('*.py'))
    failures = []
    print('='*80)
    for sf in scripts:
        print(f"==> Running standalone script: {sf.name}\n")
        if run( sf ) != 0:
            print(f"ERROR: Standalone script {sf.name} failed")
            failures.append(sf.name)
        print('='*80)

    if failures:
        print()
        for f in failures:
            print(f"   -> Ended in failure: {f}")
        print()
        raise SystemExit("Standalone test script failure")

if __name__ == '__main__':
    run_all()
