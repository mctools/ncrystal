
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

"""Command-line tool for browsing available data files and plugins."""

from ._cliimpl import cli_entry_point, create_ArgumentParser, print


def climod_metadata():
    return dict(
        displaygroup = 'main',
        displayorder = 15,
        descr = ( "Browse and search available data files (e.g. the"
                  " standard library of NCMAT files) and plugins." )
    )

def parseArgs( progname, arglist, return_parser = False ):
    import textwrap
    descr = textwrap.dedent("""
    Browse and search the data files available to NCrystal (e.g. files in the
    standard data library, in the current directory, or in-memory files),
    and the loaded plugins.

    By default a list of available files is printed, grouped by the source
    delivering them, and with a short description extracted from the header
    comments of any NCMAT files. The list can be narrowed down by providing
    one or more PATTERNs (case-insensitive substrings of file names, or glob
    patterns like "Al*.ncmat"), by requiring certain words to be present in
    the file names or NCMAT header comments (--search), or by only showing
    files from a given source (--factory). Files can also be selected based
    on their physics content with --where (see below), which requires all
    candidate files to be loaded.

    Files marked "(hidden)" are shadowed by files with the same name from a
    higher priority source, and are listed with their full name (like
    "stdlib::Al_sg225.ncmat"), which can be used to select them explicitly.
    """).strip()
    epilog = textwrap.dedent("""
    examples:
      %(prog)s                     # list all files
      %(prog)s Al                  # files with "al" in the name
      %(prog)s "*_sg225*"          # files matching a glob pattern
      %(prog)s -s vdos -s togo     # search names and header comments
      %(prog)s -E -s "boron|b4c"   # search with regular expression
      %(prog)s -c Al_sg225.ncmat   # show header comments of a file
      %(prog)s -x Al_sg225.ncmat   # show full content of a file
      %(prog)s --info Al_sg225.ncmat # show all info about a file
      %(prog)s --plugins           # list loaded plugins
      %(prog)s -w "'B' in elements and absxs > 100"
      %(prog)s -w "'vdos' in dyninfo" -w "crystal and sg == 225"
      %(prog)s --props Al_sg225.ncmat  # show the properties of a file
      %(prog)s -w "'B' in elements" --sort absxs --reverse
      %(prog)s -f stdlib --columns formula,sg,density
    """).strip()
    epilog += ( '\n\nphysics properties available in --where expressions'
                ' (and --props):\n' + _propdocs_str() + '\n\n' )
    epilog += textwrap.fill(
        '--where expressions are Python expressions using the properties'
        ' above, comparison and boolean operators, set/string/number'
        ' literals, and the functions: %s. Expressions failing due to'
        ' unavailable (None) values are considered false.'%(
            ', '.join(sorted(_where_fct_names()))), width = 79 )
    import argparse
    parser = create_ArgumentParser( prog = progname,
                                    description = descr,
                                    epilog = epilog,
                                    formatter_class
                                    = argparse.RawDescriptionHelpFormatter )
    parser.add_argument('pattern', type=str, nargs='*', metavar='PATTERN',
                        help='Only show files whose names match the pattern.')
    parser.add_argument('-s','--search', type=str, action='append',
                        default=[], metavar='WORD',
                        help=('Only show files which contain WORD in their'
                              ' name or NCMAT header comments (case'
                              '-insensitive). Can be specified multiple'
                              ' times, in which case all words must be'
                              ' present. Matching comment lines are shown.'))
    parser.add_argument('-E','--regex', action='store_true',
                        help=('Interpret search WORDs and PATTERNs as'
                              ' (case-insensitive) Python regular'
                              ' expressions, which match anywhere in the text,'
                              ' with ^ and $ matching at line boundaries'
                              ' (e.g. -E -s "boron|b4c").'))
    parser.add_argument('-f','--factory', type=str, default=None,
                        metavar='NAME',
                        help=('Only show files delivered by the named'
                              ' factory (e.g. "stdlib" or "virtual").'))
    parser.add_argument('-w','--where', type=str, action='append',
                        default=[], metavar='EXPR',
                        help=('Only show files for which the Python expression'
                              ' EXPR is true, based on physics properties of'
                              ' the loaded material (see below). Can be'
                              ' specified multiple times, in which case all'
                              ' expressions must be true.'))
    parser.add_argument('--props', action='store_true',
                        help=('Show the physics properties of each file (i.e.'
                              ' the values available in --where expressions).'
                              ))
    parser.add_argument('--columns', type=str, default=None, metavar='KEYS',
                        help=('Show a table with the given comma-separated'
                              ' physics properties (and "description")'
                              ' of each file, instead of the usual listing.'))
    parser.add_argument('--sort', type=str, default=None, metavar='KEY',
                        help=('Show a table sorted by KEY, which is "name"'
                              ' or a physics property (which is then also'
                              ' shown). Files without a value are listed'
                              ' last.'))
    parser.add_argument('--reverse', action='store_true',
                        help='Reverse the order of --sort.')
    parser.add_argument('--json', action='store_true',
                        help=('Output all information about the selected'
                              ' files (including physics properties) as'
                              ' JSON.'))
    parser.add_argument('-c','--comments', action='store_true',
                        help='Show full NCMAT header comments of the files.')
    parser.add_argument('--info', action='store_true',
                        help=('Show all available information about each'
                              ' selected file, including physics properties,'
                              ' header comments, and usage examples.'))
    parser.add_argument('--count', action='store_true',
                        help='Only print the number of selected files.')
    parser.add_argument('--path', action='store_true',
                        help=('Only print the on-disk paths of the selected'
                              ' files, one per line (files which are not'
                              ' on disk, e.g. in-memory files, are'
                              ' skipped).'))
    parser.add_argument('--names', action='store_true',
                        help=('Only print the names of the files, one per'
                              ' line (useful for scripting).'))
    parser.add_argument('-x','--extract', type=str, default=None,
                        metavar='DATANAME',
                        help=('Print the full content of DATANAME (e.g. a file'
                              ' name), using the same lookup mechanism as for'
                              ' data in NCrystal cfg-strings. This can'
                              ' therefore also be used to inspect in-memory'
                              ' (or on-demand created) data.'))
    parser.add_argument('--plugins', action='store_true',
                        help='List the currently loaded plugins.')
    parser.add_argument('--color','--colour', type=str,
                        default='auto', metavar='WHEN',
                        choices=sorted(_color_choices),
                        help=('Whether to highlight --search hits with colors'
                              ' (like grep): "always", "never", or "auto"'
                              ' (the default, only when printing to a'
                              ' terminal). In "auto" mode the NO_COLOR,'
                              ' FORCE_COLOR, CLICOLOR_FORCE, and TERM(=dumb)'
                              ' env vars are respected, and the color can be'
                              ' changed via GREP_COLORS (e.g. "mt=01;32").'
                              ' Like for grep, WHEN must be given as'
                              ' --color=WHEN, and a plain --color means'
                              ' --color=auto.'))
    if return_parser:
        return parser
    #Like grep, a plain --color (without "=WHEN") means --color=auto:
    arglist = [ ( '--color=auto' if a in ('--color','--colour') else a )
                for a in arglist ]
    args = parser.parse_args( arglist )
    nmodes = sum( bool(e) for e in ( args.extract, args.plugins ) )
    if nmodes > 1:
        parser.error('Do not specify both --extract and --plugins.')
    if nmodes and ( args.pattern or args.search or args.factory
                    or args.comments or args.names or args.where
                    or args.props or args.columns or args.sort
                    or args.json or args.count or args.path
                    or args.info ):
        parser.error('--extract and --plugins can not be combined'
                     ' with other options.')
    outmodes = [ o for o,v in [ ('--info',args.info),
                                ('--names',args.names),
                                ('--count',args.count),
                                ('--path',args.path),
                                ('--json',args.json) ] if v ]
    if len(outmodes) > 1:
        parser.error(f'Do not specify both {outmodes[0]} and'
                     f' {outmodes[1]}.')
    if args.names and ( args.comments or args.props ):
        parser.error('Do not specify --names together with --comments'
                     ' or --props.')
    if ( args.count or args.path or args.info ) and (
            args.comments or args.props or args.columns or args.sort ):
        parser.error(f'Do not specify {outmodes[0]} together with'
                     ' --comments, --props, --columns, or --sort.')
    table = bool( args.columns or args.sort )
    others = [ o for o,v in [ ('--names',args.names),
                              ('--comments',args.comments),
                              ('--props',args.props) ] if v ]
    if table and others:
        parser.error(f'Do not specify {others[0]} together with --columns'
                     ' or --sort.')
    if args.json and ( table or others ):
        parser.error('Do not specify --json together with --columns,'
                     ' --sort, --names, --comments, or --props.')
    if args.reverse and not args.sort:
        parser.error('--reverse requires --sort.')
    from .browse import physics_props_doc
    propnames = [ n for n,d in physics_props_doc() ]
    args.columns = [ c.strip() for c in ( args.columns or '' ).split(',')
                     if c.strip() ]
    for c in args.columns:
        if c not in propnames + ['description']:
            parser.error(f'Invalid column "{c}" (must be "description" or'
                         f' one of: {", ".join(propnames)})')
    from .browse import _dict_props
    sortnames = [ n for n in propnames if n not in _dict_props ]
    if args.sort and args.sort not in sortnames + ['name']:
        parser.error(f'Invalid sort key "{args.sort}" (must be "name" or'
                     f' one of: {", ".join(sortnames)})')
    if args.sort and args.sort != 'name' and args.sort not in args.columns:
        args.columns.append( args.sort )
    args.table = table
    #Validate patterns, search words and --where expressions up front, so
    #problems are reported as usage errors:
    from . import browse as nb
    from .exceptions import NCBadInput
    try:
        nb.filter_entries( [], patterns = args.pattern, search = args.search,
                           regex = args.regex, where = args.where )
    except NCBadInput as e:
        parser.error( _cli_msg( str(e) ) )
    args.search_re = nb._compile_words( args.search, args.regex )
    return args

def _where_fct_names():
    from .browse import _where_funcs
    return list( _where_funcs )

def _cli_msg( msg ):
    #Adapt error messages from the NCrystal.browse module to CLI usage:
    msg = msg.replace('where expression','--where expression')
    if msg.startswith('Unknown name '):
        msg += ' (see --help for available properties)'
    return msg

def _propdocs_str():
    import textwrap

    from .browse import physics_props_doc
    pd = physics_props_doc()
    w = max( len(n) for n,d in pd )
    out = []
    for n,d in pd:
        ll = textwrap.wrap( d, width = 76 - w - 5 )
        out.append( f'  {n.ljust(w)} : {ll[0]}' )
        out += [ ' '*(w+5)+e for e in ll[1:] ]
    return '\n'.join(out)

class _Progress:
    #Progress indicator on stderr, only shown on a terminal and only if the
    #work takes more than a second.
    def __init__( self, what ):
        import sys
        import time
        self.__t0 = time.time()
        self.__what = what
        self.__shown = False
        self.__enabled = ( hasattr(sys.stderr,'isatty')
                           and sys.stderr.isatty() )
    def update( self, n, ntotal ):
        if not self.__enabled:
            return
        import sys
        import time
        t = time.time()
        if not self.__shown and t - self.__t0 < 1.0:
            return
        if self.__shown and t - self.__tlast < 0.1:
            return#limit update rate
        self.__shown, self.__tlast = True, t
        sys.stderr.write(f'\r{self.__what}: {n}/{ntotal}')
        sys.stderr.flush()
    def done( self ):
        if self.__shown:
            import sys
            sys.stderr.write('\r\x1b[K')
            sys.stderr.flush()

#Same options and synonyms as GNU grep:
_color_choices = { 'always' : True, 'yes' : True, 'force' : True,
                   'never' : False, 'no' : False, 'none' : False,
                   'auto' : None, 'tty' : None, 'if-tty' : None }

def _use_color( when ):
    #Decide whether to use ANSI color codes. NB: These env vars are general
    #conventions (not NCrystal specific), so are never namespaced.
    v = _color_choices[when]
    if v is not None:
        return v
    import os
    import sys
    env = os.environ
    if env.get('NO_COLOR'):
        return False
    if any( env.get(k,'0') not in ('','0')
            for k in ('FORCE_COLOR','CLICOLOR_FORCE') ):
        return True
    if env.get('TERM') == 'dumb':
        return False
    from ._common import _builtin_print, get_ncrystal_print_fct
    if get_ncrystal_print_fct() is not _builtin_print():
        return False#output redirected (e.g. captured via cli.run)
    if not ( hasattr(sys.stdout,'isatty') and sys.stdout.isatty() ):
        return False
    if sys.platform == 'win32' and not ( 'WT_SESSION' in env
                                         or 'TERM' in env ):
        return False#classic Windows consoles might not handle ANSI codes
    return True

def _match_sgr():
    #SGR code for matches, default and GREP_COLORS handling like GNU grep.
    import os
    sgr = '01;31'
    for part in os.environ.get('GREP_COLORS','').split(':'):
        k,_,v = part.partition('=')
        if k in ('mt','ms') and v and all( c.isdigit() or c==';' for c in v ):
            sgr = v
    return sgr

class _Highlighter:
    #Wraps all (non-empty) matches of the compiled regexes in color codes.
    def __init__( self, regexes, enabled ):
        self.__res = list( regexes ) if enabled else []
        self.__start = f'\x1b[{_match_sgr()}m\x1b[K'
        self.__end = '\x1b[m\x1b[K'
    def __call__( self, s ):
        if not self.__res or not s:
            return s
        spans = sorted( m.span() for r in self.__res for m in r.finditer(s)
                        if m.end() > m.start() )
        merged = []
        for a,b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max( merged[-1][1], b )
            else:
                merged.append( [a,b] )
        out, pos = [], 0
        for a,b in merged:
            out += [ s[pos:a], self.__start, s[a:b], self.__end ]
            pos = b
        out.append( s[pos:] )
        return ''.join(out)

def create_argparser_for_sphinx( progname ):
    return parseArgs( progname, [], return_parser = True )

def _linewidth():
    #Terminal width if interactive, otherwise fixed (reproducible output):
    import shutil
    import sys
    if sys.stdout.isatty():
        return max( 80, shutil.get_terminal_size().columns )
    return 80

def _truncate( s, n ):
    return s if len(s) <= n else s[:max(0,n-3)].rstrip() + '...'

def _collect( args ):
    from . import browse as nb
    from .exceptions import NCBadInput
    load = bool( args.where or args.props or args.json or args.info
                 or [ c for c in args.columns if c != 'description' ]
                 or ( args.sort and args.sort != 'name' ) )
    progress = _Progress( 'Loading materials' ) if load else None
    entries = nb.browse( args.factory, load = load,
                         progress = progress.update if progress else None )
    if progress:
        progress.done()
    args.all_entries = entries#for suggestions
    try:
        return nb.filter_entries( entries, patterns = args.pattern,
                                  search = args.search, regex = args.regex,
                                  where = args.where )
    except NCBadInput as e:
        raise NCBadInput( _cli_msg( str(e) ) ) from e

def _print_props( entry, linewidth, indent = 8 ):
    pre = ' '*indent
    if entry.props is None:
        print(_truncate(f'{pre}[could not load: {entry.error}]',
                        linewidth))
        return
    import textwrap

    from .browse import _fmt_prop
    parts = [ f'{k}={_fmt_prop(v)}' for k,v in entry.props.as_dict().items() ]
    for line in textwrap.wrap( '  '.join(parts), width = linewidth - indent,
                               break_long_words = False,
                               break_on_hyphens = False ):
        print(f'{pre}{line}')

def _table_value( entry, col ):
    if col == 'description':
        return entry.description
    if entry.props is None:
        return '-'
    v = getattr( entry.props, col )
    if v is None:
        return '-'
    if isinstance( v, frozenset ):
        return ','.join( sorted(v) ) if v else '-'
    if isinstance( v, dict ):
        return ','.join( f'{k}:{x:g}' for k,x in v.items() ) if v else '-'
    if isinstance( v, float ):
        return '%g'%v
    return str(v)

def _print_table( entries, args ):
    #NB: Rows are never truncated (no data should be lost), except for the
    #free-text description column which is shortened to fit if possible.
    if not entries:
        _print_no_matches( args )
        return
    cols = args.columns
    header = [ 'NAME' ] + [ c.upper() for c in cols ]
    rows = [ [ e.display_name ] + [ _table_value(e,c) for c in cols ]
             for e in entries ]
    def width( k ):
        return max( len(r[k]) for r in rows + [header] )
    if 'description' in cols:
        kd = 1 + cols.index('description')
        other = sum( width(k) + 2 for k in range(len(header)) if k != kd )
        room = max( 20, _linewidth() - other )
        for r in rows:
            r[kd] = _truncate( r[kd], room )
    widths = [ width(k) for k in range(len(header)) ]
    def fmt( r ):
        return '  '.join( ( v.ljust(w) if k==0 or cols[k-1]=='description'
                            else v.rjust(w) )
                          for k,(v,w) in enumerate(zip(r,widths)) ).rstrip()
    print( fmt(header) )
    for r in rows:
        print( fmt(r) )

def _print_info( entries, args ):
    if not entries:
        _print_no_matches( args )
        return
    import textwrap
    linewidth = _linewidth()
    for n, e in enumerate( entries ):
        if n:
            print()
        name = e.display_name
        print(f'==> {name}')
        srcdescr = ( f'{e.source} (factory "{e.factory}",'
                     f' priority {e.priority})' )
        fields = [ ('Description', e.description or '-'),
                   ('Full key', e.fullkey),
                   ('Source', srcdescr),
                   ('Data type', e.datatype or '-'),
                   ('On-disk path', e.path or '-') ]
        if e.hidden:
            fields.append( ('Hidden', ( 'yes, by a higher priority entry'
                                        ' with the same name' ) ) )
        for k, v in fields:
            ll = textwrap.wrap( v, width = linewidth - 20,
                                break_long_words = False,
                                break_on_hyphens = False ) or ['']
            print(f'    {k:<14}: {ll[0]}')
            for x in ll[1:]:
                print(f'{"":20}{x}')
        print('  Physics properties:')
        _print_props( e, linewidth, indent = 4 )
        comments = _strip_empty( e.comments )
        if comments:
            print('  Header comments:')
            for ll in comments:
                print(f'    # {ll}'.rstrip())
        print('  Usage examples:')
        pycmd = ( f'python3 -c \'import NCrystal as NC;'
                  f' NC.load("{name}").dump()\'' )
        for d, c in [ ( 'Plot cross sections', f'nctool "{name}"' ),
                      ( 'Show full content',
                        f'ncrystal browse -x "{name}"' ),
                      ( 'Load in Python', pycmd ) ]:
            print(f'    # {d}:')
            print(f'    {c}')

def _print_no_matches( args ):
    print('No matching files found.')
    if not args.regex and any( c in w for w in args.search + args.pattern
                               for c in '|^$+()\\{}' ):
        print('Note: Search WORDs are matched literally. Use -E/--regex'
              ' for regular expressions (e.g. -E -s "boron|b4c").')
    sugg = _suggestions( args )
    if sugg:
        print('Did you mean: %s?'%( ', '.join(sugg) ))

def _suggestions( args ):
    #Suggest similar names if (plain) name patterns alone matched nothing:
    from .browse import _is_glob, filter_entries
    patterns = [ p for p in args.pattern if not _is_glob(p) ]
    if args.regex or not patterns or patterns != args.pattern:
        return []
    entries = getattr( args, 'all_entries', [] )
    if filter_entries( entries, patterns = patterns ):
        return []
    import difflib
    def stem( n ):
        return n.rsplit('.',1)[0].lower() if '.' in n else n.lower()
    #Compare with full names (without extension), and with their parts
    #(e.g. "diamond" in "C_sg227_Diamond.ncmat"):
    cands = {}
    for e in entries:
        st = stem(e.name)
        cands.setdefault( st, e.name )
        for part in st.split('_'):
            if len(part) >= 4:
                cands.setdefault( part, e.name )
    res = []
    for p in patterns:
        for m in difflib.get_close_matches( stem(p), list(cands), n = 3,
                                            cutoff = 0.7 ):
            if cands[m] not in res:
                res.append( cands[m] )
    return res

def _strip_empty( lines ):
    lines = list( lines or [] )
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return lines

def _print_listing( items, args ):
    linewidth = _linewidth()
    hl = _Highlighter( args.search_re, _use_color( args.color ) )
    groups = []
    for i in items:
        key = ( i.factory, i.source, i.priority )
        if not groups or groups[-1][0] != key:
            groups.append( ( key, [] ) )
        groups[-1][1].append( i )
    for (factname, source, priority), group in groups:
        n = len(group)
        src = f' ({source}, priority={priority})' if source else (
            f' (priority={priority})' )
        print(f'==> {n} file{"" if n==1 else "s"} from "{factname}"{src}:')
        namew = min( 40, max( len(i.display_name) for i in group ) )
        for i in group:
            name = i.display_name
            extra = ''
            if i.hidden:
                extra = ' (hidden)'
            descr = i.description
            #NB: Truncate before highlighting, so color codes do not count:
            padding = ' '*( max(0,namew-len(name)) )
            if descr and not args.comments:
                room = linewidth - 4 - max(namew,len(name)) - 2 - len(extra)
                line = ( f'    {hl(name)}{padding}  '
                         f'{hl(_truncate(descr,room))}{extra}' )
            else:
                line = f'    {hl(name)}{extra}'
            print(line.rstrip())
            matching = ( i.matching_lines( args.search, regex = args.regex )
                         if args.search else [] )
            for ll in matching:
                if not args.comments:
                    print('        | '+hl(_truncate(ll,linewidth-10)))
            if args.props:
                _print_props( i, linewidth )
            if args.comments:
                comments = _strip_empty( i.comments )
                for ll in comments:
                    print(f'        # {hl(ll)}'.rstrip())
                if comments:
                    print()
    if not items:
        _print_no_matches( args )

@cli_entry_point
def main( progname, arglist ):
    args = parseArgs( progname, arglist )
    if args.extract:
        from .core import createTextData
        print( createTextData( args.extract ).rawData, end='' )
        return
    if args.plugins:
        from .plugins import browsePlugins
        browsePlugins( dump = True )
        return
    items = _collect( args )
    if args.sort:
        from .browse import sort_entries
        items = sort_entries( items, args.sort, reverse = args.reverse )
    if args.count:
        print( len(items) )
        return
    if args.info:
        _print_info( items, args )
        return
    if args.path:
        for i in items:
            if i.path:
                print( i.path )
        return
    if args.json:
        import json
        print( json.dumps( [ i.as_dict() for i in items ], indent = 1 ) )
        return
    if args.table:
        _print_table( items, args )
        return
    if args.names:
        for i in items:
            print( i.display_name )
        return
    _print_listing( items, args )
