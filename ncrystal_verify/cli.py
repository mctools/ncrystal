
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

from contextlib import contextmanager as _ctxmgr

@_ctxmgr
def work_in_tmpdir():
    """Context manager for working in a temporary directory (automatically
    created+cleaned) and then switching back"""
    import os
    import tempfile
    the_cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            os.chdir(tmpdir)
            yield
        finally:
            os.chdir(the_cwd)#Important to leave tmpdir *before* deletion, to
                             #avoid PermissionError on Windows.

def prepend_to_path_var( env, varname, val ):
    import platform
    val = str(val)
    varsep = ':' if platform.system()!='Windows' else ';'
    if varname in env:
        env[varname] = '%s%s%s'%(val,varsep,env[varname])
    else:
        env[varname] = val


def run_test( script ):
    name = script.stem
    print(" .. running test %s"%name)
    pypath = script.parent.parent.joinpath('pypath')
    assert pypath.is_dir()
    import subprocess
    import sys
    import os
    env = os.environ.copy()
    prepend_to_path_var( env, 'PYTHONPATH', pypath )
    #Needed for windows, leaving them on all the time for now:
    env['PYTHONIOENCODING']='UTF-8'
    env['PYTHONLEGACYWINDOWSSTDIO']='UTF-8'
    with work_in_tmpdir():
        rv = subprocess.run( [sys.executable, str(script)],
                             env = env,
                             capture_output = True,
                              )
        if rv.returncode != 0:
            print(" .. failed !")
            print_lines_with_snipping(rv.stdout,prefix='stdout: ')
            print_lines_with_snipping(rv.stderr,prefix='stderr: ')
            return False
        else:
            logfile = script.parent.joinpath('%s.log'%name)
            if logfile.exists():
                diff = calc_diff_output( logfile.read_bytes(),
                                         rv.stdout,
                                         'REFERENCE',
                                         'ACTUAL' )
                if diff:
                    print('Log file differs!')
                    print(diff)
                    return False
            print(" .. success")
    return True #success

def as_lines( s ):
    return ( s.decode('utf8','backslashreplace')
             if isinstance(s,bytes)
             else s ).splitlines(keepends=True)

def print_lines_with_snipping( b, prefix ):
    lines = as_lines(b)
    n1,n2 = 20,80
    if len(lines)>n1+n2+20:
        lines = ( lines[0:n1]
                  + ['<SNIPPED %i lines>'%(len(lines)-n1-n2)]
                  +lines[-n2:] )
    print(prefix + prefix.join(lines))

def calc_diff_output( a, b, beforetxt, aftertxt ):
    if a == b:
        return ''
    import difflib
    ud = difflib.unified_diff( as_lines(a), as_lines(b),
                               fromfile=beforetxt, tofile=aftertxt )
    return ''.join(ud)


def extract_needs_statement( pyfile ):
    with pyfile.open('rt') as fh:
        for line in fh:
            if line.startswith('# NEEDS: '):
                return set(line[len('# NEEDS: '):].split())
    return None

def prepare_needs( scripts ):
    import importlib
    script2dep = {}
    alldeps = set()
    for script in scripts:
        name = script.stem
        assert name not in script2dep
        needs = extract_needs_statement( script )
        if needs:
            script2dep[name] = needs
            alldeps.update( needs )
    absent = set()
    for dep in sorted(alldeps):
        modtoimport = dep
        if dep=='matplotlib':
            modtoimport = 'matplotlib.pyplot'
        if dep=='toml':
            import sys
            if sys.version_info[0:2] >= (3,11):
                #tomllib always present in 3.11+
                continue
            else:
                #use optional module tomli in 3.10 and earlier:
                modtoimport = 'tomli'
        try:
            importlib.import_module(modtoimport)
        except ModuleNotFoundError:
            absent.add( dep )
    return dict( absent = absent,
                 script2dep = script2dep )


def parse_args():
    from argparse import ArgumentParser, RawTextHelpFormatter
    import textwrap
    def wrap(t,w=59):
        return textwrap.fill( ' '.join(t.split()), width=w )

    descr = """Validate NCrystal installation by launching several tests."""
    parser = ArgumentParser( description=wrap(descr,79),
                             formatter_class = RawTextHelpFormatter )
    parser.add_argument('-m','--allow-missing', action='append',
                        metavar='DEP',
                        help=wrap(
                            """Allow to skip tests which do not have the listed
                            optional dependencies installed (the special name
                            "all" means all of them). Any test missing a
                            dependency not listed with this flag will be marked
                            as a failure. Specify this flag multiple times or
                            use comma-separation to list more than one""" )
                        )

    parser.add_argument('-f','--filter', action='append',
                        metavar='PATTERN', dest='filters',
                        help=wrap(
                            """Select tests based on fnmatch'ing againt their
                            names, and only run these tests (default is to run
                            all). Specify this flag multiple times or use
                            comma-separation to list more than one pattern""" )
                        )

    args = parser.parse_args()

    def flatten_to_set( list_of_strs ):
        m = set()
        for a in list_of_strs or []:
            m.update(set(e.strip() for e in a.split(',')))
        return m

    args.allow_missing = flatten_to_set( args.allow_missing )
    args.filters = flatten_to_set( args.filters )

    return args

def main():
    args = parse_args()
    filterfct = None
    if args.filters:
        import fnmatch
        def filterfct( name ):
            return any( fnmatch.fnmatch( name, f )
                        for f in args.filters )


    import pathlib
    datadir = pathlib.Path(__file__).parent.joinpath('data')
    scripts = sorted( datadir.joinpath('scripts').glob('*.py'),
                      key = lambda p : (p.name, p) )
    depinfo = prepare_needs( scripts )
    jobs = []

    for script in scripts:
        name = script.stem
        needs = depinfo['script2dep'].get(name,set())
        missing = needs.intersection( depinfo['absent'] )
        jobs.append( (name, script, missing ) )

    failures = 0
    successes = 0
    dep_skipped_ok = 0
    dep_skipped_bad = 0
    hidden = 0
    ntot = 0
    for name, script, missing in jobs:
        ntot += 1
        if filterfct and not filterfct(name):
            hidden += 1
            continue
        if missing:
            ok = ( 'all' in args.allow_missing
                   or not (missing-args.allow_missing) )
            if ok:
                dep_skipped_ok += 1
                print(" .. skipping test %s (missing: %s)"%(name,' '.join(missing)))
            else:
                dep_skipped_bad += 1
                print(" .. can not run test %s (missing: %s)"%(name,' '.join(missing)))
        else:
            if not run_test(script):
                failures += 1
            else:
                successes += 1
    assert len(jobs) == ntot

    print()
    print('Test summary:')
    print()
    stats = [
        ('successful tests', successes ),
        ('failed tests', failures ),
        ('tests w/o deps (ok)', dep_skipped_ok),
        ('tests w/o deps (not ok)', dep_skipped_bad),
        ('tests hidden by filters', hidden ),
    ]
    stats.append( ( 'TOTAL', ntot) )

    w1 = max( len(s) for s,n in stats )
    w2 = max( len(str(n)) for s,n in stats )
    for s,n in stats:
        print('   %s:  %s'%(s.rjust(w1),str(n).rjust(w2)))
    print()

    assert ntot == ( successes
                     + failures
                     + dep_skipped_ok
                     + dep_skipped_bad
                     + hidden )
    assert sum( n for s,n in stats ) == 2*ntot

    nbad = failures + dep_skipped_bad
    if nbad:
        print(f'ERROR: Problems with {nbad} tests.')
        raise SystemExit(1)
    if not successes:
        print('ERROR: Requires at least a '
              'single test to run succesfully.')
        raise SystemExit(1)
    print('All selected tests ran successfully.')
