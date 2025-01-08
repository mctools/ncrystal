
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

"""Module implementing the closely related grep, find, and replace modes."""

import fnmatch

def main_grep( parser ):
    parser.init( 'Grep though files for specified PATTERNS.' )
    grepfindreplace_addargs( parser, is_grep = True )
    args = parser.parse_args()
    pos,neg = prepare_needles( args.PATTERN)
    files = iter_repo_files( args.types, args.pathfilter )
    if args.list:
        for f in files:
            if iter_nonempty(grep(f,pos,neg)):
                print(f)
    else:
        for f in files:
            for line in grep(f,pos,neg):
                print('%s:%s'%(f,line))

def main_grepl( parser ):
    parser.init(
        """Grep though files for specified PATTERNS and list files with hits."""
    )
    grepfindreplace_addargs( parser, is_grep = True )
    args = parser.parse_args()
    pos,neg = prepare_needles( args.PATTERN )
    for f in iter_repo_files( args.types, args.pathfilter ):
        if iter_nonempty(grep(f,pos,neg)):
            print(f)

def main_find( parser ):
    parser.init(
        """List repo files."""
    )
    grepfindreplace_addargs( parser )
    args = parser.parse_args()
    for f in iter_repo_files( args.types, args.pathfilter ):
        print(f)

def do_replace( fromstr, tostr,
                greppatterns = None,
                types = None,
                pathfilter = None ):
    assert not isinstance( types, str ), "give list of strings not a string"
    if fromstr == tostr:
        return
    pos,neg = prepare_needles( greppatterns or [] )
    for f in iter_repo_files( types, pathfilter ):
        replacement_lines = set()
        for line in grep(f,pos,neg):
            if fromstr in line:
                replacement_lines.add(line)
        if not replacement_lines:
            continue
        #At least one hit in the file!
        content = f.read_text().splitlines()#fixme: deal with binary files
                                            #(unicodedecodeerror).
        n = 0
        for i in range(len(content)):
            if content[i] in replacement_lines:
                n += 1
                content[i] = content[i].replace(fromstr,tostr)
        content.append('')#for join to add final newline (not strictly kosher
                          #for a generic replace, but we want to enforce these
                          #trailing newlines anyway.
        print(f"Replacing {n} lines in {f}")
        f.write_text('\n'.join(content))

def main_replace( parser ):
    parser.init(
        """Perform replacement (case sensitive,
        no wild cards or regexp syntax)"""
    )
    grepfindreplace_addargs( parser, is_replace = True )
    args = parser.parse_args()
    do_replace( args.FROMSTR, args.TOSTR,
                args.greppatterns,
                args.types,
                args.pathfilter )

def grepfindreplace_addargs( parser, *, is_grep = False, is_replace = False ):
    assert not ( is_grep and is_replace )
    if is_grep:
        parser.add_argument( 'PATTERN', nargs='+', help="""Grep pattern""")
    if is_replace:
        parser.add_argument( 'FROMSTR', help="""String to be substituted""")
        parser.add_argument( 'TOSTR', help="""Value of substitution""")
    parser.add_argument( '-t','--types', nargs='+',
                         help="""Patterns to select certain files. Either using
                         wildcards on the file paths, or special keys like 'py',
                         'cpp', etc. Provide a pattern called "HELP" in
                         order to get more information.""")
    parser.add_argument( '-p','--pathfilter', nargs='+', metavar=('PATH'),
                         help="""Use to operate only on certain subpaths of
                         the repository (no wildcards allowed, but can be
                         negated with !).""")
    if is_grep:
        parser.add_argument( '-l','--list', action='store_true',
                             help="""Simply list files with hits.""")
    if is_replace:
        parser.add_argument( '--greppatterns', nargs='+',metavar='PATTERN',
                             help="""If these grep patterns are specified,
                             replacements will only be performed in lines
                             matching these patterns.""")

def mode_replace( files, fromstr, tostr, needles ):
    if fromstr == tostr:
        return
    pos,neg = prepare_needles( needles)
    for f in files:
        replacement_lines = set()
        for line in grep(f,pos,neg):
            if fromstr in line:
                replacement_lines.add(line)
        if replacement_lines:
            content = f.read_text().splitlines()
            n = 0
            for i in range(len(content)):
                if content[i] in replacement_lines:
                    n += 1
                    content[i] = content[i].replace(fromstr,tostr)
            content.append('')#for join to add final newline
            print(f"Replacing {n} lines in {f}")
            f.write_text('\n'.join(content))

class StrMatch:
    def __init__( self, pattern ):
        self._case_insensitive = False
        if pattern.endswith('//i'):
            self._case_insensitive = True
            pattern = pattern[:-3].lower()
        self._pattern = pattern
        self._has_wildcards = any( c in pattern for c in '*?' )
        if self._has_wildcards:
            if not self._pattern.startswith('*'):
                self._pattern = '*'+self._pattern
            if not self._pattern.endswith('*'):
                self._pattern += '*'
    def do_match( self, s ):
        if self._case_insensitive:
            s = s.lower()
        return ( fnmatch.fnmatchcase(s,self._pattern)
                 if self._has_wildcards else
                 ( self._pattern in s ) )

def grep( f, pos_needles, neg_needles ):
    with f.open('rt') as fh:
        for line in fh:
            if line.endswith('\n'):
                line = line[:-1]
                if ( ( not pos_needles
                       or any( n.do_match(line) for n in pos_needles ) )
                     and not any( n.do_match(line) for n in neg_needles ) ):
                    yield line

def iter_nonempty( iterable ):
    return next(iterable, None) is not None

def prepare_needles( needles ):
    pos,neg=[],[]
    for n in (needles or []):
        if n.startswith('!'):
            neg.append(StrMatch(n[1:]))
        else:
            pos.append(StrMatch(n))
    return pos,neg

def iter_repo_files( types = None, pathfilter = None ):
    #types is the patterns one can provide to srciter.all_files_iter, and the
    #pathfilter is an additional selection on the paths.

    patterns = types or []
    if 'HELP' in patterns:
        from .srciter import special_patterns_db
        import shlex
        print('The following special types are predefined for convenience, and')
        print('also provide an idea of the syntax otherwise supported:')
        print()
        def fmt(k):
            #shell quote, but include at least quotes for pedagogical reasons
            f = shlex.quote(k)
            return "'%s'"%k if f==k else f
        n = max(len(fmt(k)) for k in special_patterns_db )
        for k,v in special_patterns_db.items():
            print("   %s expands to: %s"%( fmt(k).ljust(n),
                                           ' '.join(fmt(e) for e in v) ) )
        print()
        print("Additionally, add '!' in front of any pattern to exclude files")
        print('matching such patterns')
        raise SystemExit(0)

    from .srciter import all_files_iter
    files = all_files_iter( *patterns )
    if not pathfilter:
        return files

    import pathlib
    def filter_files(file_list):
        from .dirs import reporoot
        from .util import path_is_relative_to

        pos,neg = [], []
        for pf in pathfilter:
            if pf.startswith('!'):
                tgt = neg
                pf = pf[1:]
            else:
                tgt = pos
            pf = pathlib.Path(pf).resolve().absolute()
            if not path_is_relative_to(pf,reporoot):
                raise SystemExit(f'Path outside repo: {pf}')
            if not pf.exists():
                raise SystemExit(f'Missing file: {pf}')
            tgt.append( str(pf.relative_to(reporoot)) )

        for f in file_list:
            k = str(f.relative_to(reporoot))
            ok = True
            for n in neg:
                if k.startswith(n):
                    ok = False
                    break
            if not ok:
                continue
            if not pos:
                yield f
            for p in pos:
                if k.startswith(p):
                    yield f
                    break
    return filter_files(files)
