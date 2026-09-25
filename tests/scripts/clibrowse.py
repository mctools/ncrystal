#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

# Test the "ncrystal browse" command-line tool.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import NCrystalDev.cli as nc_cli
from NCrystalDev._common import capture_print_ctxmgr
from NCTestUtils.common import ensure_error
import re
import shlex

_ncmat_withcomments = """NCMAT v7
#
#   My test material (for testing only).
#   Second line of description which is long enough that it will have to be
#   truncated in the listing.
#
#   Some more comments mentioning Togo and VDOS.
#
@DENSITY
  1 g_per_cm3
@DYNINFO
  element H
  fraction 1
  type vdosdebye
  debye_temp 300
"""

_ncmat_nocomments = """NCMAT v7
@DENSITY
  2 g_per_cm3
@DYNINFO
  element O
  fraction 1
  type vdosdebye
  debye_temp 400
"""

def run( *args, show = True ):
    print(f"============= CLI >>browse {shlex.join(args)}<< =============")
    with capture_print_ctxmgr() as cap:
        nc_cli.run('browse',*args)
    out = ''.join(cap.data)
    #Location of stdlib depends on installation:
    out = re.sub(r'from "stdlib" \(.*, priority=',
                 'from "stdlib" (<stdlib-location>, priority=', out)
    if show:
        print(out,end='')
    return out

def main():
    NC.removeAllDataSources()
    NC.enableStandardDataLibrary()
    NC.registerInMemoryFileData('mytestmat.ncmat',_ncmat_withcomments)
    NC.registerInMemoryFileData('othermat.ncmat',_ncmat_nocomments)
    NC.registerInMemoryFileData('notncmat.laz','whatever')
    #Hide stdlib file by in-memory file with same name:
    NC.registerInMemoryFileData('Al_sg225.ncmat',_ncmat_nocomments)

    run('--help')
    run('-f','virtual')
    run('Al_sg225')
    run('stdlib::Be*')
    run('--names','-f','stdlib','*sg229*')
    run('--names','Al_sg225')
    run('-s','togo','-f','virtual')
    run('-s','TOGO','-s','vdos','-f','virtual')
    run('-s','togo','-s','nonexistentword')
    run('-s','togo','-f','stdlib','Al_sg225')
    run('-c','-f','virtual')
    run('-c','stdlib::Al_sg225')
    run('-x','mytestmat.ncmat')
    out = run('--plugins',show=False)
    assert 'plugins loaded' in out
    with ensure_error(NC.NCFileNotFound,
                      'Could not find data: "nonexistent.ncmat"'):
        run('-x','nonexistent.ncmat')
    import argparse
    with ensure_error(argparse.ArgumentError,
                      '--extract and --plugins can not be combined'
                      ' with other options.'):
        run('-x','mytestmat.ncmat','-f','virtual')
    with ensure_error(argparse.ArgumentError,
                      'Do not specify both --names and --comments.'):
        run('--names','-c')

if __name__ == '__main__':
    main()
