#!/usr/bin/env python3

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

import sys
import pathlib
import shutil
import shlex
import platform
import subprocess

is_windows = (platform.system() == 'Windows')
ENCODING = sys.stdout.encoding

def run( app_file, reflogfile = None ):
    wd=pathlib.Path('./tmprundir')
    shutil.rmtree(wd,ignore_errors=True)
    wd.mkdir()
    cmd = [str(app_file)]
    if app_file.name.endswith('.py'):
        cmd = [sys.executable] + cmd
    print("MCTools TestLauncher running command:")
    for e in cmd:
        print(f"  {shlex.quote(e)}")
    print()
    sys.stdout.flush()

    sys.stdout.flush()
    sys.stderr.flush()
    r = subprocess.run( cmd,
                        encoding=ENCODING,
                        errors='backslashreplace',
                        capture_output = True,
                        cwd = wd )
    sys.stdout.flush()
    sys.stderr.flush()
    print("MCTools TestLauncher done running command.")
    r_stdout = ( r.stdout.decode(ENCODING,errors='backslashreplace')
                 if isinstance(r.stdout,bytes)
                 else ( r.stdout or '' ) )
    r_stderr = ( r.stderr.decode(ENCODING,errors='backslashreplace')
                 if isinstance(r.stderr,bytes)
                 else ( r.stderr or '' ) )
    assert isinstance(r_stdout,str)
    assert isinstance(r_stderr,str)
    sys.stdout.write(r_stdout)
    if r_stderr:
        #Todo support this, by merging them (probably needs subprocess.Popen not
        #subprocess.run):
        for line in r_stderr.splitlines():
            #sys.stderr.buffer.write(b'stderr> '+line+'\n'.encode())
            sys.stderr.write('stderr> '+line+'\n')
        sys.stderr.flush()
        raise SystemExit('Error: Process emitted output on'
                         ' stderr (not supported yet with ref logs)')
    output_raw = r_stdout
    newout = pathlib.Path('./output.log').absolute()
    newout.unlink(missing_ok=True)
    assert not newout.exists()
    newout.write_bytes(output_raw.encode(ENCODING,errors='backslashreplace'))

    if r.returncode == 3221225781:
        if is_windows:
            raise SystemExit('Error: Command ended with exit'
                             f' code {r.returncode} (usually'
                             ' indicates "DLL not found")')
    if r.returncode != 0:
        raise SystemExit(f'Error: Command ended with exit code {r.returncode}')
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
        return c if 32<=ord(c)<=126 else r'\u{%s}'%(hex(ord(c))[2:])

    def explicit_unicode_str(s):
        return ''.join( explicit_unicode_char(c) for c in s)

    import difflib
    for line in difflib.unified_diff(refoutput,
                                     output,
                                     fromfile='BEFORE',
                                     tofile='AFTER',
                                     lineterm=''):
        print('DIFF> %s'%explicit_unicode_str(line))

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
