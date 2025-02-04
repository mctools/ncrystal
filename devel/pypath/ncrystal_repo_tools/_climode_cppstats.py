
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
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
    return 'Provide information and statistics about particular C++ component'

print_comp_stats_path_mode_options = ('reltorepo','abs','rel')

def print_comp_stats( comp,
                      path_mode = 'reltorepo',
                      withcore = True,
                      showincs_to_comp = None ):
    #from .core_components import file_calc_sloc_count

    assert path_mode in print_comp_stats_path_mode_options
    if path_mode == 'reltorepo':
        from .dirs import reporoot
        path_fmt_root = reporoot
    elif path_mode == 'rel':
        from pathlib import Path as _
        path_fmt_root = _().absolute()
    else:
        assert path_mode == 'abs'
        path_fmt_root = None
    if path_fmt_root:
        import os.path

    def fmtpath( p ):
        if not p:
            return '<none>'
        return os.path.relpath(p,path_fmt_root) if path_fmt_root else str(p)

    def fmtpath_name( p ):
        if p and path_mode == 'reltorepo':
            return p.name
        return fmtpath(p)

    from .extract_includes import get_nccomp_include_statements as getncinc
    ignore_list = set([comp.name])
    if not withcore:
        ignore_list.add('core')
    def getinc_comps( f ):
        return set( c for i,c in getncinc(f, ignore_list = ignore_list) )
    def fmtinc_comps( f ):
        c = ' '.join(getinc_comps( f ))
        return ' (includes %s)'%c if c else ''


    print(f'Component "{comp.name}"')
    print()
    print('  Part of public API : %s'%('no' if comp.is_internal else 'yes'))
    print('  Include dir        : %s'%fmtpath(comp.hdrdir))
    print('  Source dir         : %s'%fmtpath(comp.srcdir))
    print('  Dependency file    : %s'%fmtpath(comp.depfile))
    print()
    def show_extra_incs( f ):
        if not showincs_to_comp:
            return
        for i,c in getncinc(f, ignore_list = ignore_list):
            if showincs_to_comp=='ALL' or c == showincs_to_comp:
                print(f'      VIA #include "{i}"')

    if comp.hdrfiles:
        print("Header files:")
        for f in comp.hdrfiles:
            print('   %s%s'%(fmtpath_name(f),fmtinc_comps(f)))
            show_extra_incs(f)
    if comp.srcfiles:
        print("Source files:")
        import itertools
        for f in itertools.chain(comp.local_hdrs,comp.srcfiles):
            print('   %s%s'%(fmtpath_name(f),fmtinc_comps(f)))
            show_extra_incs(f)

def main( parser ):
    parser.init( """Provide information and statistics about particular C++
    component""" )

    parser.add_argument(
        'COMP',
        help="""Code component name (e.g. "core", "utils", ...) that you want
        to investigate."""
    )
    _pmodedef = 'reltorepo'
    assert _pmodedef in print_comp_stats_path_mode_options
    parser.add_argument(
        '-p','--pathmode', choices = print_comp_stats_path_mode_options,
        default = _pmodedef,
        help=f"""How to display paths (default is {_pmodedef})."""
    )
    parser.add_argument( '--withcore', action='store_true', help="""Without
    this options, include statements for the core components will not be
    shown""" )
    parser.add_argument(
        '--showincs',type=str,
        help="""Show full include statements for provided component (use name
        ALL for all)"""
    )

    args = parser.parse_args()

    from .core_components import load_components
    #from .extract_includes import get_include_staments_from_file as getinc
    #from .dirs import coreroot, coreroot_include

    name2comp = load_components()
    if args.COMP not in name2comp:
        parser.error(f'Invalid component name: "{args.COMP}"')

    print_comp_stats( name2comp[ args.COMP ],
                      args.pathmode,
                      withcore = args.withcore,
                      showincs_to_comp = args.showincs )
    return
