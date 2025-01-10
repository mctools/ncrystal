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

import pathlib
import subprocess
from simplebuild.cfg import dirs

rv = subprocess.run(['sb_nccmd_config','-s'],
                    check=True,
                    capture_output=True)
assert not rv.stderr

vars = {}
varnames_ordered = []
for iline, line in enumerate(rv.stdout.decode().splitlines()):
    line = line.strip()
    if iline<2 or not line:
        continue
    print(line)
    tmp = line.split(' : ',1)
    varname, varval = [tmp[0],None] if len(tmp)==1 else tmp
    assert varname not in vars
    vars[varname] = varval
    varnames_ordered.append( varname )

if varnames_ordered != sorted(varnames_ordered):
    raise SystemExit('variable name list not sorted properly in'
                     ' ncrystal-config main.cc nccfg_show_item_list()')

sbinst = dirs.installprefix
reporoot = pathlib.Path(__file__).resolve().parent.parent.parent.parent.parent

assert sbinst.joinpath('bin').samefile( vars['bindir'] )
assert sbinst.joinpath('lib').samefile( vars['libdir'] )
assert sbinst.joinpath('lib').samefile( vars['shlibdir'] )
assert sbinst.joinpath('data','NCData').samefile( vars['datadir'] )
assert reporoot.joinpath('ncrystal_core','include').samefile( vars['includedir'] )
assert vars['libname'].startswith('libPKG__NCCInterface.')
assert pathlib.Path( vars['libpath'] ).is_file()
assert pathlib.Path( vars['shlibpath'] ).is_file()
