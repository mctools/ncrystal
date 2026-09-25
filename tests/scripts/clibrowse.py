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
import os
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

_ncmat_crystal = """NCMAT v7
# A small crystal.
@CELL
 cubic 4.04958
@SPACEGROUP
 225
@ATOMPOSITIONS
 Al 0 1/2 1/2
 Al 0 0 0
 Al 1/2 1/2 0
 Al 1/2 0 1/2
@DEBYETEMPERATURE
 Al 400
"""

def run( *args, show = True ):
    print(f"============= CLI >>browse {shlex.join(args)}<< =============")
    with capture_print_ctxmgr() as cap:
        nc_cli.run('browse',*args)
    out = ''.join(cap.data)
    #Location of stdlib depends on installation:
    out = re.sub(r'from "stdlib" \(.*, priority=',
                 'from "stdlib" (<stdlib-location>, priority=', out)
    out = out.replace('\x1b','<ESC>')#make color codes visible in log
    if show:
        print(out,end='')
    return out

_color_envvars = ('NO_COLOR','FORCE_COLOR','CLICOLOR_FORCE',
                  'GREP_COLORS','TERM')

def setenv( **kw ):
    for k,v in kw.items():
        if v is None:
            os.environ.pop(k,None)
        else:
            os.environ[k] = v

def test_colors():
    #NB: Output is captured, so "auto" mode means no colors unless forced.
    run('-s','togo','-s','my','-f','virtual')
    run('-s','togo','-s','my','-f','virtual','--color=always')
    run('-s','togo','-f','virtual','-c','--colour=yes')
    run('-s','togo','-f','virtual','--color=never')
    run('-s','togo','--color','mytestmat')#plain --color means auto
    setenv( FORCE_COLOR = '1' )
    run('-s','togo','-f','virtual')
    setenv( NO_COLOR = '1' )
    run('-s','togo','-f','virtual')
    run('-s','togo','-f','virtual','--color=always')
    setenv( NO_COLOR = None, FORCE_COLOR = None,
            GREP_COLORS = 'sl=1:ms=01;32:ln=35' )
    run('-s','togo','-f','virtual','--color=always')
    setenv( GREP_COLORS = None )
    #Regex highlighting, with overlapping matches merged:
    run('-E','-s','tog|nothing','-s','ogo','-f','virtual','--color=always')
    import argparse
    with ensure_error(argparse.ArgumentError,
                      "argument --color/--colour: invalid choice: 'blue'"
                      " (choose from 'always', 'auto', 'force', 'if-tty',"
                      " 'never', 'no', 'none', 'tty', 'yes')"):
        run('--color=blue')

def test_physics():
    run('--props','-f','virtual')
    run('--props','-f','stdlib','LiquidHeavyWater')
    run('-w','"H" in elements','-f','virtual')
    run('-w','crystal and sg==225 and natoms==4','-f','virtual')
    run('-w','sg > 200','-f','virtual')#None comparisons are false
    run('-w','"vdosdebye" in dyninfo','-w','density<1.5','-f','virtual')
    run('-w','absxs>0.2 and len(atoms)==1 and state=="solid"','-f','virtual')
    run('-w','formula in ("D2O","H2O")','-f','stdlib','Liquid')
    run('-w','"scatknl" in dyninfo and "O" in elements','-f','stdlib',
        'Liquid','--names')
    run('-w','elements <= {"Al","O"} and nphases==1','-f','virtual')
    import argparse
    def bad_where( expr, errmsg ):
        with ensure_error(argparse.ArgumentError,errmsg):
            run('-w',expr)
    bad_where('foo > 1', ('Unknown name "foo" in --where expression'
                          ' "foo > 1" (see --help for available properties)'))
    bad_where('open("x")', ('Unknown name "open" in --where expression'
                            ' "open(\"x\")" (see --help for available'
                            ' properties)'))
    bad_where('elements.__class__', ('Invalid --where expression'
                                     ' "elements.__class__" (private'
                                     ' attributes are not allowed)'))
    bad_where('absxs >', 'Invalid --where expression "absxs >": invalid syntax')
    with ensure_error(NC.NCBadInput,'Error evaluating --where expression'
                      ' "absxs/0 > 1": float division by zero'):
        run('-w','absxs/0 > 1','-f','virtual')

def main():
    #Colors in output must not depend on the environment of the test:
    orig_env = dict( (k,os.environ.get(k)) for k in _color_envvars )
    setenv( **dict( (k,None) for k in _color_envvars ) )
    try:
        main_impl()
        test_colors()
    finally:
        setenv( **orig_env )

def main_impl():
    NC.removeAllDataSources()
    NC.enableStandardDataLibrary()
    NC.registerInMemoryFileData('mytestmat.ncmat',_ncmat_withcomments)
    NC.registerInMemoryFileData('othermat.ncmat',_ncmat_nocomments)
    NC.registerInMemoryFileData('notncmat.laz','whatever')
    NC.registerInMemoryFileData('mycrystal.ncmat',_ncmat_crystal)
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
    #Literal vs. regex search and patterns:
    run('-s','togo|nomatch','-f','virtual')
    run('-E','-s','togo|nomatch','-f','virtual')
    run('-E','-s','test.*only','-s','^ *some','-f','virtual')
    run('-E','--names','^my.*mat','othermat\\.ncmat$','^other$')
    run('--names','my.*mat')
    run('-c','-f','virtual')
    test_physics()
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
                      'Do not specify --names together with --comments'
                      ' or --props.'):
        run('--names','-c')
    with ensure_error(argparse.ArgumentError,
                      'Invalid regular expression "(": missing ),'
                      ' unterminated subpattern at position 0'):
        run('-E','-s','(')

if __name__ == '__main__':
    main()
