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

def run( *args, show = True, sanitize_paths = False ):
    print(f"============= CLI >>browse {shlex.join(args)}<< =============")
    with capture_print_ctxmgr() as cap:
        nc_cli.run('browse',*args)
    out = ''.join(cap.data)
    #Location of stdlib depends on installation:
    out = re.sub(r'from "stdlib" \(.*, priority=',
                 'from "stdlib" (<stdlib-location>, priority=', out)
    out = re.sub(r'from "relpath" \(.*, priority=',
                 'from "relpath" (<current-dir>, priority=', out)
    if sanitize_paths:
        out = re.sub(r'(Source        : ).*( \(factory "(stdlib|relpath)")',
                     r'\1<location-dependent>\2', out)
        out = re.sub(r'(On-disk path  : ).*', r'\1<location-dependent>', out)
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

def test_lazlau():
    #Non-NCMAT data, both in-memory and on-disk (in current directory):
    import pathlib

    from NCrystalDev.mcstasutils import cfgstr_2_hkl
    from NCTestUtils.common import work_in_tmpdir
    def lazlau( fmt ):
        return '\n'.join( cfgstr_2_hkl( cfgstr = ( 'stdlib::Al_sg225.ncmat;'
                                                   'dcutoff=1.0' ),
                                        tgtformat = fmt, verbose = False,
                                        fp_format = '%.8g' ) ) + '\n'
    NC.registerInMemoryFileData('mem.laz',lazlau('laz'))
    NC.registerInMemoryFileData('mem.lau',lazlau('lau'))
    NC.enableRelativePaths(True)
    with work_in_tmpdir():
        pathlib.Path('disk.laz').write_text(lazlau('laz'))
        pathlib.Path('disk.lau').write_text(lazlau('lau'))
        run('*.la?')
        run('*.la?','--columns','formula,sg,a,braggthreshold,dyninfo')
        run('-s','disk')
        run('-f','relpath','--info','disk.laz',sanitize_paths=True)
        run('-f','virtual','-w','sg==225','--count')
    NC.enableRelativePaths(False)

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
    #Newer properties:
    run('-w','crystalsystem=="cubic" and braggthreshold > 4','-f','virtual')
    run('-w','max(debyetemps.values()) > 300','-f','virtual')#None is false
    run('-f','virtual','--columns','a,volume,debyetemps,msds,mass,cohxs')
    #HTML tables and atomdb:
    run('-f','virtual','--sort','sg','--columns','formula','--html')
    run('--atomdb','He','b10')
    run('--atomdb','-w','absxs > 1000 and natural','--sort','absxs',
        '--reverse','--columns','mass,absxs')
    run('--atomdb','Li*','--names')
    run('--atomdb','Li','--csv','--columns','a,cohsl,absxs')#no mass: FP
    run('--atomdb','Li6','--html','--columns','absxs')
    run('--atomdb','H','--json','--columns','a')
    run('--atomdb','xyz')
    run('--atomdb','-w','absxs > 1e9','--count')
    #Info view:
    run('-f','virtual','--info','mycrystal','notncmat')
    run('--info','stdlib::Al_sg225',sanitize_paths=True)#hidden entry
    run('-f','virtual','--info','mycrystl')
    #Counting, paths and suggestions:
    run('-f','virtual','--count')
    run('-f','virtual','-w','crystal','--count')
    run('-f','virtual','--path')#in-memory files have no path
    out = run('-f','stdlib','Al_sg225','--path',show=False).strip()
    #NB: Location depends on installation (might even be embedded):
    assert out == '' or out.replace('\\','/').endswith('/Al_sg225.ncmat')
    run('-f','virtual','mycrystl')
    run('-f','virtual','mycrystl','--columns','sg')
    run('-f','virtual','mycrystal','-w','absxs > 100')#no suggestions
    run('-f','virtual','qwertyzzz')
    #Tables, sorting and JSON:
    run('-f','virtual','--columns','formula,sg,density,dyninfo,description')
    run('-f','virtual','--sort','density','--reverse')
    run('-f','virtual','--sort','sg')#unavailable values last
    run('-f','virtual','--sort','name','--reverse')
    import json
    out = run('-f','virtual','mycrystal','--json',show=False)
    import NCTestUtils.stabilise_ncpprint # noqa F401
    import NCrystalDev._common as nc_common
    nc_common.ncpprint( json.loads(out) )#FP precision clipped
    out = run('-f','virtual','--sort','density','--columns','formula,sg,'
              'dyninfo,description','--json',show=False)
    nc_common.ncpprint( json.loads(out) )
    #CSV (full precision, so check values rather than printing them):
    out = run('-f','virtual','--sort','sg','--columns','formula,density,'
              'debyetemps,elements,description','--csv',show=False)
    import csv
    import io
    rows = list( csv.reader( io.StringIO(out) ) )
    print('CSV header:',rows[0])
    for r in rows[1:]:
        dens = float(r[2]) if r[2] else None
        print('CSV row:',r[0],r[1],None if dens is None else '%.6g'%dens,
              r[3].split(':')[0] if r[3] else None,r[4],repr(r[5]))
    #No truncation:
    run('-f','virtual','mytestmat','--no-truncate')
    run('-f','virtual','mytestmat','--columns','description','--no-truncate')
    run('-f','virtual','mytestmat','-s','truncated','--no-truncate')
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
    def bad_args( errmsg, *args ):
        with ensure_error(argparse.ArgumentError,errmsg):
            run(*args)
    from NCrystalDev.browse import physics_props_doc
    propnames = ', '.join( n for n,d in physics_props_doc() )
    propnames_sortable = ', '.join( n for n,d in physics_props_doc()
                                    if n not in ('debyetemps','msds') )
    bad_args('Invalid column "foo" (must be "description" or one of: '
             + propnames + ')', '--columns','sg,foo')
    bad_args('Invalid sort key "foo" (must be "name" or one of: '
             + propnames_sortable + ')', '--sort','foo')
    bad_args('--reverse requires --sort.','--reverse')
    bad_args('Do not specify both --names and --count.','--names','--count')
    bad_args('--search can not be used together with --atomdb.',
             '--atomdb','-s','x')
    bad_args('--html requires --columns or --sort.','--html')
    bad_args(('Invalid column "foo" (must be one of: z, a, element, natural,'
              ' mass, cohsl, cohxs, incohxs, scatxs, absxs)'),
             '--atomdb','--columns','foo')
    bad_args(('Unknown name "sg" in --where expression "sg > 1" (see --help'
              ' for available properties)'),'--atomdb','-w','sg > 1')
    bad_args('Do not specify both --info and --json.','--info','--json')
    bad_args(('Do not specify --info together with --comments, --props,'
              ' --columns, or --sort.'),'--info','--props')
    bad_args(('Do not specify --path together with --comments, --props,'
              ' --columns, or --sort.'),'--path','--sort','name')
    bad_args('Invalid sort key "debyetemps" (must be "name" or one of: '
             + propnames_sortable + ')', '--sort','debyetemps')
    with ensure_error(NC.NCBadInput,'Error evaluating --where expression'
                      ' "elements.foo": \'frozenset\' object has no'
                      ' attribute \'foo\''):
        run('-w','elements.foo','-f','virtual')
    bad_args('Do not specify --props together with --columns or --sort.',
             '--columns','sg','--props')
    bad_args(('Do not specify --json together with --names, --comments,'
              ' or --props.'),'--json','--props')
    bad_args('--csv requires --columns or --sort.','--csv')
    bad_args('Do not specify both --csv and --json.','--csv','--json',
             '--sort','sg')
    bad_args('Do not specify both --names and --json.','--json','--names')
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
    test_lazlau()
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
