
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

"""Command-line tool for browsing available data and plugins."""

from ._cliimpl import cli_entry_point, create_ArgumentParser, print


def climod_metadata():
    return dict(
        displaygroup = 'main',
        displayorder = 15,
        descr = ( "Browse and search available data (e.g. the"
                  " standard library of NCMAT files) and plugins." )
    )

def parseArgs( progname, arglist, return_parser = False ):
    import textwrap
    descr = textwrap.dedent("""
    Browse and search the data available to NCrystal (e.g. files in the
    standard data library or in the current directory, in-memory data, or
    data created on-demand like "solid::B4C/2.52gcm3"), and the loaded
    plugins.

    By default a list of available data entries is printed, grouped by the
    source delivering them, and with a short description extracted from the
    header comments of any NCMAT data. The list can be narrowed down by
    providing one or more PATTERNs (case-insensitive substrings of names, or
    glob patterns like "Al*.ncmat"), by requiring certain words to be present
    in the names or NCMAT header comments (--search), or by only showing
    entries from a given source (--factory). Entries can also be selected
    based on their physics content with --where (see below), which requires
    all candidate entries to be loaded.

    Entries marked "(hidden)" are shadowed by entries with the same name from
    a higher priority source, and are listed with their full name (like
    "stdlib::Al_sg225.ncmat"), which can be used to select them explicitly.
    """).strip()
    epilog = textwrap.dedent("""
    examples:
      %(prog)s                     # list all data
      %(prog)s Al                  # data with "al" in the name
      %(prog)s "*_sg225*"          # data matching a glob pattern
      %(prog)s -s vdos -s togo     # search names and header comments
      %(prog)s -E -s "boron|b4c"   # search with regular expression
      %(prog)s -c Al_sg225.ncmat   # show header comments
      %(prog)s -x Al_sg225.ncmat   # show full content
      %(prog)s --info Al_sg225.ncmat # show all info
      %(prog)s --plugins           # list loaded plugins
      %(prog)s -w "'B' in elements and absxs > 100"
      %(prog)s -w "'vdos' in dyninfo" -w "crystal and sg == 225"
      %(prog)s --props Al_sg225.ncmat  # show the physics properties
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
                        help='Only show data whose names match the pattern.')
    parser.add_argument('-s','--search', type=str, action='append',
                        default=[], metavar='WORD',
                        help=('Only show data which contain WORD in their'
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
                        help=('Only show data delivered by the named'
                              ' factory (e.g. "stdlib" or "virtual").'))
    parser.add_argument('-w','--where', type=str, action='append',
                        default=[], metavar='EXPR',
                        help=('Only show data for which the Python expression'
                              ' EXPR is true, based on physics properties of'
                              ' the loaded material (see below). Can be'
                              ' specified multiple times, in which case all'
                              ' expressions must be true.'))
    parser.add_argument('--props', action='store_true',
                        help=('Show the physics properties of each entry (i.e.'
                              ' the values available in --where expressions).'
                              ))
    parser.add_argument('--columns', type=str, default=None, metavar='KEYS',
                        help=('Show a table with the given comma-separated'
                              ' physics properties (and "description")'
                              ' of each entry, instead of the usual listing.'))
    parser.add_argument('--sort', type=str, default=None, metavar='KEY',
                        help=('Show a table sorted by KEY, which is "name"'
                              ' or a physics property (which is then also'
                              ' shown). Entries without a value are listed'
                              ' last.'))
    parser.add_argument('--reverse', action='store_true',
                        help='Reverse the order of --sort.')
    parser.add_argument('--json', action='store_true',
                        help=('Output all information about the selected'
                              ' data (including physics properties) as'
                              ' JSON. With --columns or --sort, only the'
                              ' name and the table columns are output.'))
    parser.add_argument('--csv', action='store_true',
                        help=('Output the table of --columns or --sort in'
                              ' CSV format (with full numerical'
                              ' precision).'))
    parser.add_argument('--no-truncate', action='store_true',
                        help=('Never shorten long descriptions or other'
                              ' text to fit the line width.'))
    parser.add_argument('-c','--comments', action='store_true',
                        help='Show full NCMAT header comments.')
    parser.add_argument('--info', action='store_true',
                        help=('Show all available information about each'
                              ' selected entry, including physics properties,'
                              ' header comments, and usage examples.'))
    parser.add_argument('--count', action='store_true',
                        help='Only print the number of selected entries.')
    parser.add_argument('--path', action='store_true',
                        help=('Only print the on-disk paths of the selected'
                              ' files, one per line (entries which are not'
                              ' on disk, e.g. in-memory data, are'
                              ' skipped).'))
    parser.add_argument('--names', action='store_true',
                        help=('Only print the names of the entries, one per'
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
                    or args.info or args.csv ):
        parser.error('--extract and --plugins can not be combined'
                     ' with other options.')
    outmodes = [ o for o,v in [ ('--csv',args.csv),
                                ('--info',args.info),
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
    if args.json and others:
        parser.error('Do not specify --json together with --names,'
                     ' --comments, or --props.')
    if args.csv and not table:
        parser.error('--csv requires --columns or --sort.')
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
        if args.regex:
            nb._compile_words( args.pattern, True )
        nb._compile_words( args.search, args.regex )
        for w in args.where:
            nb._WhereExpr( w )
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

def _browser( args ):
    #Create DataBrowser with selection according to the arguments:
    from . import browse as nb
    from .exceptions import NCBadInput
    progress = _Progress( 'Loading materials' )
    b = nb.DataBrowser( factory = args.factory, progress = progress.update )
    args.all_browser = b#for suggestions
    try:
        sel = b.match( *args.pattern, regex = args.regex )
        sel = sel.search( *args.search, regex = args.regex )
        sel = sel.where( *args.where )
        if args.sort:
            sel = sel.sorted( args.sort, reverse = args.reverse )
    except NCBadInput as e:
        raise NCBadInput( _cli_msg( str(e) ) ) from e
    finally:
        progress.done()
    return sel

def _print_no_matches( args ):
    print('No matching data found.')
    if not args.regex and any( c in w for w in args.search + args.pattern
                               for c in '|^$+()\\{}' ):
        print('Note: Search WORDs are matched literally. Use -E/--regex'
              ' for regular expressions (e.g. -E -s "boron|b4c").')
    sugg = _suggestions( args )
    if sugg:
        print('Did you mean: %s?'%( ', '.join(sugg) ))

def _suggestions( args ):
    #Suggest similar names if (plain) name patterns alone matched nothing:
    from .browse import _is_glob
    patterns = [ p for p in args.pattern if not _is_glob(p) ]
    if args.regex or not patterns or patterns != args.pattern:
        return []
    b = args.all_browser
    if b.match( *patterns ):
        return []
    res = []
    for p in patterns:
        res += [ n for n in b.suggestions(p) if n not in res ]
    return res

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
    sel = _browser( args )
    truncate = not args.no_truncate
    if args.count:
        print( len(sel) )
        return
    if args.path:
        for e in sel:
            if e.path:
                print( e.path )
        return
    if args.names:
        for n in sel.names():
            print( n )
        return
    if args.json and not args.table:
        print( sel.to_json(), end = '' )
        return
    if not sel:
        _print_no_matches( args )
        return
    if args.info:
        print( sel.info( linewidth = _linewidth() ), end = '' )
        return
    if args.table:
        fmt = 'csv' if args.csv else ( 'json' if args.json else 'text' )
        print( sel.table( args.columns, fmt = fmt, truncate = truncate,
                          linewidth = _linewidth() ), end = '' )
        return
    hl = _Highlighter( args.search_re, _use_color( args.color ) )
    print( sel.format_listing( comments = args.comments,
                               props = args.props, truncate = truncate,
                               linewidth = _linewidth(), highlight = hl ),
           end = '' )
