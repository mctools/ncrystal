
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
    import pathlib
    script_path = pathlib.Path(script)
    name = script_path.stem
    print(" .. running test %s"%name)
    pypath = script_path.parent.parent.joinpath('pypath')
    assert pypath.is_dir()
    import subprocess
    import sys
    import os
    env = os.environ.copy()
    prepend_to_path_var( env, 'PYTHONPATH', pypath )
    with work_in_tmpdir():
        rv = subprocess.run( [sys.executable, str(script)],
                             env = env,
                             capture_output = True,
                              )
        if rv.returncode != 0:
            print(" .. failed !")
            print('stdout:')
            print(rv.stdout.decode())
            print('stderr:')
            print(rv.stderr.decode())
            raise SystemExit(1)
        else:
            logfile = script_path.parent.joinpath('%s.log'%name)
            if logfile.exists():
                diff = calc_diff_output( logfile.read_bytes(),
                                         rv.stdout,
                                         'REFERENCE',
                                         'ACTUAL' )
                if diff:
                    print('Log file differs!')
                    print(diff)
                    raise SystemExit(1)
            print(" .. success")


def calc_diff_output( a, b, beforetxt, aftertxt ):
    if a == b:
        return ''
    def as_lines( s ):
        return ( s.decode('utf8','backslashreplace')
                 if isinstance(s,bytes)
                 else s ).splitlines(keepends=True)

    import difflib
    ud = difflib.unified_diff( as_lines(a), as_lines(b),
                               fromfile=beforetxt, tofile=aftertxt )
    return ''.join(ud)

def main():
    import pathlib
    datadir = pathlib.Path(__file__).parent.joinpath('data')
    scripts = sorted( str(p) for p in datadir.joinpath('scripts').glob('*.py') )
    for s in scripts:
        run_test(s)
