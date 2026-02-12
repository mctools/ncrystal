
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

group_defs = [
    ('main','Most important modules'),
    ('conv','Modules related to material data conversions'),
    ('misc','Miscellaneous modules')
]

_cache_md = {}
def get_module_metadata( mode ):
    md=_cache_md.get(mode)
    if md:
        return md
    from ._cliimpl import _resolve_cmd_and_import_climod
    climod, _ = _resolve_cmd_and_import_climod( mode, [] )
    md = climod.climod_metadata()
    assert len(md)==3
    assert isinstance(md['descr'],str)
    assert isinstance(md['displayorder'],int)
    assert isinstance(md['displaygroup'],str)
    md['descr'] = ' '.join(md['descr'].split())
    assert md['descr'].endswith('.')
    _cache_md[mode] = md
    return md

def get_mode_list():
    from .cli import cli_tool_list
    return cli_tool_list( canonical_names = False )

def progname():
    return 'ncrystal'

def usage():
    print(generate_usage(return_list=False))

def collect_usage_data( mode_list = None ):
    if mode_list is None:
        mode_list = get_mode_list()

    group_keys = [k for k,v in group_defs]
    group_key2descr = dict(group_defs)

    group_key2data = dict( (k,[]) for k in group_keys )
    for mode in mode_list:
        md = get_module_metadata( mode )
        descr = md['descr']
        grp = md['displaygroup']
        order = md['displayorder']
        assert grp in group_key2descr, ( "invalid group (must be one"
                                         " of: %s)"%(' '.join(group_keys)) )
        group_key2data[grp].append( (order,mode,descr) )
    out = []
    for grp in group_keys:
        items=[ (mode,descr)
                for order,mode,descr in sorted(group_key2data[grp]) ]
        out.append( (grp,items) )
    return out

def generate_usage(*,return_list = True):
    import textwrap
    theprogname = progname()
    ml = get_mode_list()
    example_mode = 'minimc'
    assert example_mode in ml
    modelist_indent = '     '
    max_len_mode = max(len(mode) for mode in ml)
    width=80
    descrwidth = width - max_len_mode -3 - len(modelist_indent)
    group_key2descr = dict(group_defs)

    out = []
    out.append('Usage:')
    out.append('')
    out.append(f'  $> {theprogname} MODULE <module options>')
    out.append('')
    out.append('Here MODULE is the name of a specific command-line'
               ' module, like "nctool",')
    out.append('"minimc", etc.')
    out.append('')
    out.append('Use the -h or --help flag to get information about the'
               ' options available for a')
    out.append('given module. For example:')
    out.append('')
    out.append(f'  $> {theprogname} {example_mode} --help')
    out.append('')
    out.append('The full list of available modules are:')

    for grp, items in collect_usage_data(ml):
        assert len(items)>0,f"empty usage group: {grp}"
        grpdescr = group_key2descr[grp]
        out.append('')
        out.append(f'  {grpdescr}:')
        out.append('  '+'='*(len(grpdescr)+1))
        out.append('')
        for mode,descr in items:
            ld = textwrap.wrap(descr, width=descrwidth)
            out.append(f'{modelist_indent}{mode.ljust(max_len_mode)} : {ld[0]}')
            for e in ld[1:]:
                out.append(' '*(width-descrwidth)+e)

    return out if return_list else ('\n'.join(out)+'\n')

def main( argv = None ):
    if argv is None:
        import sys
        argv = sys.argv
    assert argv
    if len(argv)<2 or argv[1] in ('-h','--h','--he','--hel','--help'):
        usage()
        return
    if len(argv)==2 and argv[1] in ('-l','--l','--li','--lis','--list'):
        print( ' '.join(get_mode_list()) )
        return

    mode = argv[1]
    if mode not in get_mode_list():
        raise SystemExit(f'ERROR: Invalid module "{mode}". Run without'
                         ' arguments to list available modules.')

    from ._cliimpl import _resolve_cmd_and_import_climod
    climod, argv = _resolve_cmd_and_import_climod( mode, argv[2:] )
    argv[0] = f'{progname()} {mode}'
    climod.main( argv = argv )
