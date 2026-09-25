
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
    files from a given source (--factory).

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
      %(prog)s -c Al_sg225.ncmat   # show header comments of a file
      %(prog)s -x Al_sg225.ncmat   # show full content of a file
      %(prog)s --plugins           # list loaded plugins
    """).strip()
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
    parser.add_argument('-f','--factory', type=str, default=None,
                        metavar='NAME',
                        help=('Only show files delivered by the named'
                              ' factory (e.g. "stdlib" or "virtual").'))
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
    if return_parser:
        return parser
    args = parser.parse_args( arglist )
    nmodes = sum( bool(e) for e in ( args.extract, args.plugins ) )
    if nmodes > 1:
        parser.error('Do not specify both --extract and --plugins.')
    if nmodes and ( args.pattern or args.search or args.factory
                    or args.comments or args.names ):
        parser.error('--extract and --plugins can not be combined'
                     ' with other options.')
    if args.names and args.comments:
        parser.error('Do not specify both --names and --comments.')
    return args

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

def _name_matches( entry, pattern ):
    import fnmatch
    p = pattern.lower()
    #Patterns with '::' are matched against the full key (e.g. stdlib::Al..):
    s = ( entry.fullKey if '::' in p else entry.name ).lower()
    return fnmatch.fnmatchcase( s, p ) if _is_glob( p ) else ( p in s )

def _header_comments( entry ):
    #Header comments of NCMAT data as list of lines (None if not NCMAT).
    if not entry.name.lower().endswith('.ncmat'):
        return None
    from ._ncmatimpl import _extractInitialHeaderCommentsFromNCMATData as f
    from .core import createTextData
    try:
        return f( createTextData( entry.fullKey ).rawData )
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
        items = [ i for i in items
                  if any( _name_matches( i.entry, p ) for p in args.pattern ) ]
    words = [ w.lower() for w in args.search ]
    for i in items:
        i.matching_lines = []
    if words:
        selected = []
        for i in items:
            lines = i.comments or []
            text = '\n'.join( [ i.entry.name ] + lines ).lower()
            if all( w in text for w in words ):
                i.matching_lines = [ ll.strip() for ll in lines
                                     if any( w in ll.lower() for w in words ) ]
                selected.append( i )
        items = selected
    return items

def _strip_empty( lines ):
    lines = list( lines or [] )
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return lines

def _print_listing( items, args ):
    linewidth = _linewidth()
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
            if descr and not args.comments:
                room = linewidth - 4 - max(namew,len(name)) - 2 - len(extra)
                line = ( f'    {name.ljust(namew)}  '
                         f'{_truncate(descr,room)}{extra}' )
            else:
                line = f'    {name}{extra}'
            print(line.rstrip())
            for ll in getattr(i,'matching_lines',[]):
                if not args.comments:
                    print(_truncate(f'        | {ll}',linewidth))
            if args.comments:
                comments = _strip_empty( i.comments )
                for ll in comments:
                    print(f'        # {ll}'.rstrip())
                if comments:
                    print()
    if not items:
        print('No matching files found.')

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
