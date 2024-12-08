from .srciter import all_files_iter
from .dirs import reporoot
import fnmatch

def parse_args():
    import argparse
    import textwrap
    def wrap(t):
        return textwrap.fill(' '.join(t.split()),width=77)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Tool which can be used by NCrystal developers in various
        way, including for checking and editing files.""",
    )
    parser.add_argument( "PATTERN", nargs='*',
                         help=wrap("""Patterns to select certain files. Either
                         using wildcards on the file paths, or special keys like
                         'py', 'cpp', etc."""))
    parser.add_argument( "-c", "--check", nargs='*',metavar='CHECK',
                         help=wrap("""Run various tests. Defaults to running
                         all. Provide the name 'list' to simply see a list of
                         available checks."""))
    parser.add_argument( "-f","--find", action='store_true',
                         help=wrap("""List all files matching the provided
                         PATTERNs."""))
    parser.add_argument( "--grep", nargs='+',metavar='STR',
                         help=wrap("""Grep according to patterns"""))
    parser.add_argument( "--grepl", nargs='+',metavar='STR',
                         help=wrap("""Grep and just list files with any hits"""))
    parser.add_argument( "--replace", nargs=2,metavar=('FROMSTR','TOSTR'),
                         help=wrap("""Perform replacement (case sensitive, no
                         wild cards or regexp syntax)"""))
    parser.add_argument( "--linefilter", nargs='+', metavar=('STR'),
                         help=wrap("""Use with --replace to select only certain
                         lines (of course lines must already have the FROMSTR in
                         it, this is an additional restriction."""))
    parser.add_argument( "--pathfilter", nargs='+', metavar=('STR'),
                         help=wrap("""Use to operate only on certain subpaths of
                         the repository (no wildcards allowed, but can be
                         negated with !)."""))
    args = parser.parse_args()
    nmodes = sum((1 if e else 0) for e in ( (args.check is not None),
                                           args.find,
                                           args.grepl,
                                           args.grep,
                                           args.replace))
    if nmodes == 0:
        parser.error('Please specify flags to select '
                     'operation (run with -h to learn more)')
    if ( nmodes>1
         or ( (args.check is not None) and args.PATTERN )
         or args.linefilter and not args.replace() ):
        parser.error('inconsistent arguments')
    return args

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

def mode_grepl( files, needles ):
    pos,neg = prepare_needles( needles)
    for f in files:
        if iter_nonempty(grep(f,pos,neg)):
            print(f)

def mode_grep( files, needles ):
    pos,neg = prepare_needles( needles)
    for f in files:
        for line in grep(f,pos,neg):
            print('%s:%s'%(f,line))

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

def main():
    args = parse_args()
    files = all_files_iter( *args.PATTERN)
    if args.pathfilter:
        import pathlib
        def filter_files(file_list):
            pos,neg = [], []
            for pf in args.pathfilter:
                if pf.startswith('!'):
                    tgt = neg
                    pf = pf[1:]
                else:
                    tgt = pos
                pf = pathlib.Path(pf).resolve().absolute()
                if not pf.is_relative_to(reporoot):
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
        files = filter_files(files)

    if args.grepl:
        mode_grepl( files, args.grepl )
    elif args.grep:
        mode_grep( files, args.grep )
    elif args.replace:
        mode_replace( files,
                      args.replace[0],
                      args.replace[1],
                      args.linefilter )
    elif args.find:
        #List selected:
        for f in files:
            print(f)
    elif (args.check is not None):
        from . import check_runner
        if not args.check:
            check_runner.run_all_checks()
        elif 'list' in args.check:
            for t in check_runner.get_available_checks_list():
                print(t)
        else:
            for c in sorted(set(args.check)):
                check_runner.run_check(c)
    else:
        raise RuntimeError('logic error: should have been caught in parser')

if __name__ == '__main__':
    main()
