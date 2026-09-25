
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
      %(prog)s --plugins           # list loaded plugins
      %(prog)s -w "'B' in elements and absxs > 100"
      %(prog)s -w "'vdos' in dyninfo" -w "crystal and sg == 225"
      %(prog)s --props Al_sg225.ncmat  # show the properties of a file
    """).strip()
    epilog += ( '\n\nphysics properties available in --where expressions'
                ' (and --props):\n' + _propdocs_str() + '\n\n' )
    epilog += textwrap.fill(
        '--where expressions are Python expressions using the properties'
        ' above, comparison and boolean operators, set/string/number'
        ' literals, and the functions: %s. Comparisons with unavailable'
        ' (None) values are considered false.'%(
            ', '.join(sorted(_where_funcs))), width = 79 )
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
    parser.add_argument('-c','--comments', action='store_true',
                        help='Show full NCMAT header comments of the files.')
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
                    or args.props ):
        parser.error('--extract and --plugins can not be combined'
                     ' with other options.')
    if args.names and ( args.comments or args.props ):
        parser.error('Do not specify --names together with --comments'
                     ' or --props.')
    try:
        args.where_code = [ _compile_where( w ) for w in args.where ]
    except _WhereError as e:
        parser.error(str(e))
    import re
    def compile_re( s ):
        try:
            #MULTILINE: ^ and $ match at line boundaries, like for grep:
            return re.compile( s if args.regex else re.escape(s),
                               re.IGNORECASE | re.MULTILINE )
        except re.error as e:
            parser.error(f'Invalid regular expression "{s}": {e}')
    args.search_re = [ compile_re(w) for w in args.search ]
    args.pattern_re = ( [ compile_re(p) for p in args.pattern ]
                        if args.regex else None )
    return args

#Physics properties of materials, available in --where expressions:
_propdocs = [
    ('elements', 'set of element names, e.g. {"Al","O"}'),
    ('atoms', 'set of atom labels, including isotopes, e.g. {"D","O"}'),
    ('nelements', 'number of different elements'),
    ('formula', 'chemical formula, e.g. "Al2O3"'),
    ('absxs', 'absorption cross section per atom at 2200m/s [barn]'),
    ('scatxs', 'free scattering cross section per atom [barn]'),
    ('density', 'density [g/cm3]'),
    ('numdens', 'number density [atoms/Aa3]'),
    ('temp', 'temperature [K]'),
    ('state', 'state of matter: "solid", "liquid", "gas", or "unknown"'),
    ('crystal', 'True if crystalline (for multiphase: any phase)'),
    ('sg', 'space group number (None if not available)'),
    ('natoms', 'number of atoms in unit cell (None if not available)'),
    ('dyninfo', ('set of dynamic info types present, among "vdos",'
                 ' "vdosdebye", "scatknl" (a full scattering kernel),'
                 ' "freegas", and "sterile"')),
    ('nphases', 'number of phases (1 for single-phase materials)'),
]

def _propdocs_str():
    import textwrap
    w = max( len(n) for n,d in _propdocs )
    out = []
    for n,d in _propdocs:
        ll = textwrap.wrap( d, width = 76 - w - 5 )
        out.append( f'  {n.ljust(w)} : {ll[0]}' )
        out += [ ' '*(w+5)+e for e in ll[1:] ]
    return '\n'.join(out)

_where_funcs = dict( len = len, min = min, max = max, any = any, all = all,
                     abs = abs, round = round, set = set, sorted = sorted )

class _WhereError(Exception):
    pass

def _compile_where( expr ):
    #Compile --where expression, only allowing known names and no private
    #attributes (as a basic safety measure and to catch typos early).
    import ast
    try:
        tree = ast.parse( expr, mode = 'eval' )
    except SyntaxError as e:
        raise _WhereError(f'Invalid --where expression "{expr}": {e.msg}')
    allowed = set( n for n,d in _propdocs ) | set( _where_funcs )
    for node in ast.walk( tree ):
        if isinstance( node, ast.Name ) and node.id not in allowed:
            raise _WhereError(f'Unknown name "{node.id}" in --where'
                              f' expression "{expr}" (see --help for'
                              ' available properties)')
        if isinstance( node, ast.Attribute ) and node.attr.startswith('_'):
            raise _WhereError(f'Invalid --where expression "{expr}"'
                              ' (private attributes are not allowed)')
    return expr, compile( tree, '<where>', 'eval' )

def _eval_where( where_code, props ):
    ns = dict( _where_funcs )
    ns['__builtins__'] = {}
    for expr, code in where_code:
        try:
            ok = eval( code, ns, dict(props) )
        except TypeError:
            ok = False#e.g. comparison with None
        except Exception as e:
            from .exceptions import NCBadInput
            raise NCBadInput(f'Error evaluating --where expression'
                             f' "{expr}": {e}') from e
        if not ok:
            return False
    return True

def _physics_props( info ):
    #Physics properties of loaded Info object (see _propdocs).
    from ._common import format_chemform
    from .atomdata import elementZToName
    phases = ( [ ph for fr,ph in info.phases ] if info.isMultiPhase()
               else [ info ] )
    leaves = []
    def add_leaves( i ):
        if not i.isMultiPhase():
            leaves.append( i )
        for fr,ph in ( i.phases if i.isMultiPhase() else [] ):
            add_leaves( ph )
    add_leaves( info )
    fc = info.getFlattenedComposition()
    elements = frozenset( elementZToName(Z) for Z,_ in fc )
    atomfracs = {}
    for Z, isotopes in fc:
        en = elementZToName(Z)
        for A, fr in isotopes:
            lbl = en if A == 0 else { (1,2):'D', (1,3):'T' }.get( (Z,A),
                                                                  f'{en}{A}' )
            atomfracs[lbl] = atomfracs.get(lbl,0.0) + fr
    ditypes = { 'DI_VDOS' : 'vdos', 'DI_VDOSDebye' : 'vdosdebye',
                'DI_ScatKnlDirect' : 'scatknl', 'DI_FreeGas' : 'freegas',
                'DI_Sterile' : 'sterile' }
    dyninfo = frozenset( ditypes.get( type(di).__name__, 'other' )
                         for i in leaves for di in i.dyninfos )
    struct = ( info.getStructureInfo() if ( not info.isMultiPhase()
                                            and info.hasStructureInfo() )
               else None )
    return dict(
        elements = elements,
        atoms = frozenset( atomfracs ),
        nelements = len( elements ),
        formula = format_chemform( sorted( atomfracs.items() ) ),
        absxs = info.getXSectAbsorption(),
        scatxs = info.getXSectFree(),
        density = info.getDensity(),
        numdens = info.getNumberDensity(),
        temp = info.getTemperature() if info.hasTemperature() else None,
        state = info.stateOfMatter().name.lower(),
        crystal = any( i.isCrystalline() for i in leaves ),
        sg = struct['spacegroup'] if struct else None,
        natoms = struct['n_atoms'] if struct else None,
        dyninfo = dyninfo,
        nphases = max( 1, len( phases ) ),
    )

def _fmt_prop( v ):
    if v is None:
        return 'None'
    if isinstance( v, frozenset ):
        return '{' + ','.join( repr(e) for e in sorted(v) ) + '}'
    if isinstance( v, float ):
        return '%g'%v
    return repr(v)

class _Progress:
    #Progress indicator on stderr, only shown on a terminal and only if the
    #work takes more than a second.
    def __init__( self, total, what ):
        import sys
        import time
        self.__t0 = time.time()
        self.__total, self.__what = total, what
        self.__shown = False
        self.__enabled = ( hasattr(sys.stderr,'isatty')
                           and sys.stderr.isatty() )
    def update( self, n ):
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
        sys.stderr.write(f'\r{self.__what}: {n}/{self.__total}')
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

def _is_glob( pattern ):
    return any( c in pattern for c in '*?[' )

def _name_matches( entry, pattern, pattern_re = None ):
    #Patterns with '::' are matched against the full key (e.g. stdlib::Al..):
    s = entry.fullKey if '::' in pattern else entry.name
    if pattern_re is not None:
        return bool( pattern_re.search( s ) )
    import fnmatch
    p, s = pattern.lower(), s.lower()
    return fnmatch.fnmatchcase( s, p ) if _is_glob( p ) else ( p in s )

def _header_comments( entry ):
    #Header comments of NCMAT data as list of lines (None if not NCMAT).
    if not entry.name.lower().endswith('.ncmat'):
        return None
    from ._ncmatimpl import _extractInitialHeaderCommentsFromNCMATData as f
    from .core import createTextData
    try:
        data = createTextData( entry.fullKey ).rawData
        #Only the part before the first section is needed:
        idx = data.find('\n@')
        return f( data if idx == -1 else data[:idx+1] )
    except Exception: # noqa BLE001
        #Unreadable data should not prevent browsing other files:
        return None

def _short_descr( comments ):
    #First paragraph of the header comments, as a single line.
    words = []
    for line in ( comments or [] ):
        if not line.strip():
            if words:
                break
            continue
        #Skip ascii-art rulers like "----" or "=====":
        if not any( c.isalnum() for c in line ):
            if words:
                break
            continue
        words += line.split()
    return ' '.join(words)

def _truncate( s, n ):
    return s if len(s) <= n else s[:max(0,n-3)].rstrip() + '...'

class _Item:
    def __init__( self, entry, hidden ):
        self.entry = entry
        self.hidden = hidden
        self.__comments = False
    @property
    def comments( self ):
        if self.__comments is False:
            self.__comments = _header_comments( self.entry )
        return self.__comments
    @property
    def display_name( self ):
        e = self.entry
        return ( e.fullKey if ( e.priority == 'OnlyOnExplicitRequest'
                                or self.hidden ) else e.name )

def _collect( args ):
    from .datasrc import browseFiles
    items, seen = [], set()
    for e in browseFiles():
        #NB: browseFiles returns entries sorted by priority, so entries with
        #names seen previously are hidden by higher priority entries.
        items.append( _Item( e, hidden = ( e.name in seen ) ) )
        seen.add( e.name )
    if args.factory is not None:
        items = [ i for i in items if i.entry.factName == args.factory ]
    if args.pattern:
        prs = args.pattern_re or [ None ]*len(args.pattern)
        items = [ i for i in items
                  if any( _name_matches( i.entry, p, pr )
                          for p, pr in zip( args.pattern, prs ) ) ]
    regexes = args.search_re
    for i in items:
        i.matching_lines = []
    if regexes:
        selected = []
        for i in items:
            lines = i.comments or []
            text = '\n'.join( [ i.entry.name ] + lines )
            if all( r.search( text ) for r in regexes ):
                i.matching_lines = [ ll.strip() for ll in lines
                                     if any( r.search(ll) for r in regexes ) ]
                selected.append( i )
        items = selected
    if args.where or args.props:
        _load_props( items )
    if args.where:
        items = [ i for i in items if i.props is not None
                  and _eval_where( args.where_code, i.props ) ]
    return items

def _load_props( items ):
    from ._msg import _suppress_msgs_ctx
    from .core import createInfo
    progress = _Progress( len(items), 'Loading materials' )
    #Loading might emit warnings (e.g. about @CUSTOM_ sections) which are not
    #relevant here:
    with _suppress_msgs_ctx():
        for n, i in enumerate( items ):
            progress.update( n )
            i.props, i.load_error = None, None
            try:
                i.props = _physics_props( createInfo( i.entry.fullKey ) )
            except Exception as e: # noqa BLE001
                #Files which can not be loaded are simply not selected:
                i.load_error = str(e) or e.__class__.__name__
    progress.done()

def _print_props( item, linewidth ):
    if item.props is None:
        print(_truncate(f'        [could not load: {item.load_error}]',
                        linewidth))
        return
    import textwrap
    parts = [ f'{k}={_fmt_prop(v)}' for k,v in item.props.items() ]
    for line in textwrap.wrap( '  '.join(parts), width = linewidth - 8,
                               break_long_words = False,
                               break_on_hyphens = False ):
        print(f'        {line}')

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
        e = i.entry
        key = ( e.factName, e.source, e.priority )
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
            descr = _short_descr( i.comments )
            #NB: Truncate before highlighting, so color codes do not count:
            padding = ' '*( max(0,namew-len(name)) )
            if descr and not args.comments:
                room = linewidth - 4 - max(namew,len(name)) - 2 - len(extra)
                line = ( f'    {hl(name)}{padding}  '
                         f'{hl(_truncate(descr,room))}{extra}' )
            else:
                line = f'    {hl(name)}{extra}'
            print(line.rstrip())
            for ll in getattr(i,'matching_lines',[]):
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
        print('No matching files found.')
        if not args.regex and any( c in w for w in args.search + args.pattern
                                   for c in '|^$+()\\{}' ):
            print('Note: Search WORDs are matched literally. Use -E/--regex'
                  ' for regular expressions (e.g. -E -s "boron|b4c").')

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
    if args.names:
        for i in items:
            print( i.display_name )
        return
    _print_listing( items, args )
