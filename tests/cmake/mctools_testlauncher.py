#!/usr/bin/env python3

# Utility script needed by mctools_testutils.cmake for launching tests and
# comparing with reference output.

import pathlib
import shutil
import shlex

def run( app_file, reflogfile = None ):
    wd=pathlib.Path('./wd')
    shutil.rmtree(wd,ignore_errors=True)
    wd.mkdir()
    import subprocess
    import sys
    cmd = [str(app_file)]
    if app_file.name.endswith('.py'):
        cmd = [sys.executable] + cmd
    print("MCTools TestLauncher running command:")
    for e in cmd:
        print(f"  {shlex.quote(e)}")
    print()
    sys.stdout.flush()
    r = subprocess.run( cmd,
                        capture_output = True,
                        cwd = wd )
    print("MCTools TestLauncher done running command.")
    if r.stderr:
        #Todo support this, by merging them (probably needs subprocess.Popen not
        #subprocess.run):
        for line in r.stderr.splitlines():
            sys.stderr.buffer.write(b'stderr> '+line+'\n'.encode())
        sys.stderr.flush()
        raise SystemExit('Error: Process emitted output on'
                         ' stderr (not supported yet with ref logs)')
    output_raw = r.stdout
    newout = pathlib.Path('./output.log').absolute()
    newout.unlink(missing_ok=True)
    assert not newout.exists()
    newout.write_bytes(output_raw)

    if r.returncode != 0:
        raise SystemExit(f'Error: Command ended with exit code {r.returncode}')
    if reflogfile is None:
        return #Done!
    refoutput = reflogfile.read_bytes()
    if output_raw == refoutput:
        sys.stdout.flush()
        print("Reference log-files are exact byte-for-byte match!")
        sys.stdout.flush()
        return
    output = output_raw.decode().splitlines()
    refoutput = refoutput.decode().splitlines()
    if output == reflogfile:
        sys.stdout.flush()
        print("Reference log-files match!")
        return
    def qp( p ):
        return shlex.quote(str(p.absolute()))
    raise SystemExit(f"""
ERROR: Output does not match that of the reference log.
New output is at: {newout.absolute()}
Reference output is at: {reflogfile.absolute()}
Unix commands to diff and update:

    colordiff -y {qp(reflogfile)} {qp(newout)} | less -r
    diff {qp(reflogfile)} {qp(newout)}
    cp {qp(reflogfile)} {qp(newout)}

""")

def main( ):
    import sys
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
