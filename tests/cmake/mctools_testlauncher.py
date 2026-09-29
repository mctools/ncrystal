
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

# Utility script needed by mctools_testutils.cmake for launching tests and
# comparing with reference output.

import pathlib
import platform
import shlex
import shutil
import subprocess
import sys

is_windows = (platform.system() == 'Windows')
ENCODING = sys.stdout.encoding

def run( app_file, reflogfile = None ):
    wd=pathlib.Path('./tmprundir')
    shutil.rmtree(wd,ignore_errors=True)
    wd.mkdir()
    cmd = [str(app_file)]
    if app_file.name.endswith('.py'):
        #-u so a crashing script does not lose its final (pipe-buffered)
        #stdout, which is essential for debugging from CI logs:
        cmd = [sys.executable,'-u'] + cmd
    print("MCTools TestLauncher running command:")
    for e in cmd:
        print(f"  {shlex.quote(e)}")
    print()
    sys.stdout.flush()
    sys.stderr.flush()
    #Stream the child's stdout live (line by line) while also collecting
    #it for the reference-log comparison below: the previous
    #buffer-then-print approach (capture_output=True) meant that a child
    #killed externally (e.g. by a ctest timeout after hanging) took all
    #its already-produced diagnostic output with it. TextIOWrapper gives
    #the same encoding/universal-newline treatment as the text-mode
    #subprocess.run did, and the stderr pipe is drained from a thread so
    #neither pipe can fill up and block the child:
    import io
    import threading
    p = subprocess.Popen( cmd, cwd = wd,
                          stdout = subprocess.PIPE,
                          stderr = subprocess.PIPE )
    stderr_parts = []
    def _drain_stderr():
        with io.TextIOWrapper( p.stderr, encoding=ENCODING,
                               errors='backslashreplace' ) as f:
            stderr_parts.extend( f )
    t = threading.Thread( target = _drain_stderr, daemon = True )
    t.start()
    stdout_parts = []
    with io.TextIOWrapper( p.stdout, encoding=ENCODING,
                           errors='backslashreplace' ) as f:
        for line in f:
            stdout_parts.append( line )
            sys.stdout.write( line )
            sys.stdout.flush()
    returncode = p.wait()
    t.join()
    sys.stdout.flush()
    sys.stderr.flush()
    print("MCTools TestLauncher done running command.")
    r_stdout = ''.join( stdout_parts )
    r_stderr = ''.join( stderr_parts )
    assert isinstance(r_stdout,str)
    assert isinstance(r_stderr,str)
    if r_stderr:
        for line in r_stderr.splitlines():
            sys.stderr.write('stderr> '+line+'\n')
        sys.stderr.flush()
        raise SystemExit('Error: Process emitted output on'
                         ' stderr (not supported yet with ref logs)')
    output_raw = r_stdout
    newout = pathlib.Path('./output.log').absolute()
    newout.unlink(missing_ok=True)
    assert not newout.exists()
    newout.write_bytes(output_raw.encode(ENCODING,errors='backslashreplace'))

    if returncode == 3221225781 and is_windows:
        raise SystemExit('Error: Command ended with exit'
                         f' code {returncode} (usually'
                         ' indicates "DLL not found")')
    if returncode != 0:
        raise SystemExit(f'Error: Command ended with exit code {returncode}')
    if reflogfile is None:
        return #Done!
    refoutput = reflogfile.read_text(encoding='utf-8')
    if output_raw == refoutput:
        sys.stdout.flush()
        print("Reference log-files are exact match!")
        sys.stdout.flush()
        return
    #output = output_raw.decode('utf-8').splitlines()
    #refoutput = refoutput.decode('utf-8').splitlines()
    output = output_raw.splitlines()
    refoutput = refoutput.splitlines()
    if output == refoutput:
        sys.stdout.flush()
        print("Reference log-files match!")
        return
    if len(output)==len(refoutput):
        for i,(o,r) in enumerate(zip(output,refoutput)):
            if o!=r:
                print(f"L{i+1} - {r}")
                print(f"L{i+1} + {o}")
    def qp( p ):
        return shlex.quote(str(p.absolute()))

    def explicit_unicode_char(c):
        #32 is space, <32 are control chars, 127 is DEL.
        return c if 32<=ord(c)<=126 else rf'\u{{{hex(ord(c))[2:]}}}'

    def explicit_unicode_str(s):
        return ''.join( explicit_unicode_char(c) for c in s)

    import difflib
    for line in difflib.unified_diff(refoutput,
                                     output,
                                     fromfile='BEFORE',
                                     tofile='AFTER',
                                     lineterm=''):
        print(f'DIFF> {explicit_unicode_str(line)}')

    raise SystemExit(f"""
ERROR: Output does not match that of the reference log.
New output is at: {newout.absolute()}
Reference output is at: {reflogfile.absolute()}
Unix commands to diff and update:

    colordiff -y {qp(reflogfile)} {qp(newout)} | less -r
    diff {qp(reflogfile)} {qp(newout)}
    cp {qp(newout)} {qp(reflogfile)}

""")

def main( ):
    assert len(sys.argv) in (2,3)
    app_file = pathlib.Path(sys.argv[1])
    if not app_file.is_file():
        raise SystemExit(f'File to run not found: {app_file}')
    reflogfile = None
    if len(sys.argv)==3:
        reflogfile = pathlib.Path(sys.argv[2])
        if not reflogfile.is_file():
            raise SystemExit(f'Reference log file not found: {reflogfile}')
    run(app_file, reflogfile)
    sys.stdout.flush()
    sys.stderr.flush()

if __name__=='__main__':
    main()
