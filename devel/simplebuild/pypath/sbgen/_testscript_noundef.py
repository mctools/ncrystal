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

import os
import platform
import pathlib
import subprocess

if platform.system() == 'Darwin':
    print('No ldd on OSX -- aborting.')
    raise SystemExit

libdir = pathlib.Path(os.environ['SBLD_LIB_DIR'])
libs = sorted( libdir.glob('libPKG__*') )
assert len(libs)>10
ok = True
for lib in libs:
    print(f"Checking {lib.name}")
    rv = subprocess.run( ['ldd','-r',lib],
                         capture_output = True,
                         check = True )
    assert not rv.stderr
    for line in rv.stdout.decode().splitlines():
        if line.startswith('undefined symbol:'):
            ok = False
            print('    ',line)
if not ok:
    raise SystemExit("Errors detected")
