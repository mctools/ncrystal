
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


def short_description():
    return 'Move files between C++ components'

def search_files( name2comp, filename ):
    bn = filename if filename.startswith('NC') else 'NC'+filename
    if bn.endswith('.hh'):
        bn = bn[:-3]
    if not bn.isidentifier():
        return
    bn_hh = bn + '.hh'
    bn_icc = bn + '.icc'
    bn_cc = bn + '.cc'
    comp = None
    def name_in_filelist( name, filelist ):
        for f in filelist:
            if name == f.name:
                return f
    for c in name2comp.values():
        if ( name_in_filelist( bn_hh, c.hdrfiles)
             or name_in_filelist( bn_cc, c.srcfiles) ):
            comp = c
            break
    if not comp:
        return
    info = {}
    name_in_filelist( bn_hh, comp.hdrfiles)
    info['hh'] = name_in_filelist( bn_hh, comp.hdrfiles)
    info['icc'] = name_in_filelist( bn_icc, comp.hdrfiles)
    info['local_hh'] = name_in_filelist( bn_hh, comp.local_hdrs)
    info['local_icc'] = name_in_filelist( bn_icc, comp.local_hdrs)
    info['cc'] = name_in_filelist( bn_cc, comp.srcfiles)
    info['comp'] = comp
    info['bn'] = bn
    return info

def determine_moves( info ):
    from .dirs import ( coreroot_include,
                        coreroot_src )

    tgtcomp_name = info['tgtcomp_name']
    tgtdir_hdr = coreroot_include.joinpath( f'NCrystal/internal/{tgtcomp_name}'
                                            if info['tgtcomp_is_internal'] else
                                            f'NCrystal/{tgtcomp_name}' )
    tgtdir_src = coreroot_src.joinpath( tgtcomp_name )

    moves = []
    f_hdrdir = [ info['hh'], info['icc'] ]
    f_srcdir = [ info['local_hh'], info['local_icc'], info['cc'] ]
    for f in f_hdrdir:
        if f:
            moves.append( ( f, tgtdir_hdr / f.name ) )
    for f in f_srcdir:
        if f:
            moves.append( ( f, tgtdir_src / f.name ) )
    return moves

def determine_replace_map( info ):
    if not info['hh']:
        return []
    from .dirs import coreroot_include
    tgtcomp_name = info['tgtcomp_name']
    nc = ( f'NCrystal/internal/{tgtcomp_name}'
           if info['tgtcomp_is_internal'] else
           f'NCrystal/{tgtcomp_name}' )
    fmtinc = 'include "%s"'
    rm = []
    for f in ( info['hh'], info['icc'] ):
        if not f:
            continue
        rm.append( ( fmtinc % f.relative_to(coreroot_include),
                     fmtinc%((f'{nc}/%s')%f.name) ) )
    return rm

def main( parser ):

    parser.init( """Move existing C++ files to different components.""" )

    parser.add_argument(
        'FILENAME',
        help="""File to move. Provide either like NCSomeFile.hh, NCSomeFile or
        SomeFile. Will move all associated files (.hh/.cc/.icc) and also update
        include statements and dependencies in all components."""
    )

    parser.add_argument(
        'TARGETCOMP',
        help="""TargetFile to move. Provide either like NCSomeFile.hh, NCSomeFile or
        SomeFile. Will move all associated files (.hh/.cc/.icc) and also update
        include statements and dependencies in all components."""
    )
    parser.add_argument(
        '--new', choices = ('internal','public'), metavar='TYPE',
        help="""Allow TARGETCOMP to specify an new component which is not yet
        present. TYPE must be either 'internal' or 'public'."""
    )
    parser.add_argument(
        '--dryrun', action = 'store_true',
        help='Show changes but do not modify anything'
    )

    args = parser.parse_args()
    setattr(args,'comp',args.TARGETCOMP)

    from .core_components import ( load_components,
                                   is_valid_component_name )
    from .dirs import coreroot_src

    if not args.comp or not is_valid_component_name(args.comp):
        parser.error('Invalid component name')

    name2comp = load_components()

    files_to_touch = []
    if args.comp not in name2comp:
        if not args.new:
            parser.error(f'Unknown component "{args.comp}" (supply '
                         '--new TYPE to create).')
        tgtcomp_is_internal = args.new == 'internal'
        files_to_touch.append( coreroot_src.joinpath(args.comp,'dep.txt') )
    else:
        tgtcomp_is_internal = name2comp[args.comp].is_internal

    info = search_files(name2comp, args.FILENAME)
    if not info:
        parser.error( 'Could not find files to move based'
                      f' on pattern "{args.FILENAME}"' )
    info['tgtcomp_name'] = args.comp
    info['tgtcomp_is_internal'] = tgtcomp_is_internal

    replacements = determine_replace_map( info )
    moves = determine_moves( info )

    if args.dryrun:
        if moves:
            print("Would move:")
        for k,v in moves:
            print()
            print(f"       {k}")
            print(f"    to {v}")
        if files_to_touch:
            print()
            print("Would create empty files:")
        for f in files_to_touch:
            print()
            print(f"       {f}")
        if replacements:
            print()
            print("Would replace:")
        for k,v in replacements:
            print()
            print(f"       {repr(k)}")
            print(f"  with {repr(v)}")
        print()
        return

    for k,v in moves:
        v.parent.mkdir(parents = True, exist_ok = True )
        v.write_text( k.read_text() )
        k.unlink()

    for f in files_to_touch:
        assert not f.exists()
        f.parent.mkdir(parents = True, exist_ok = True )
        f.write_text('')

    from ._cliutils_common import do_replace
    for k,v in replacements:
        do_replace( k, v, types = ['cpp'] )

    from .core_components import fix_deps
    fix_deps()
