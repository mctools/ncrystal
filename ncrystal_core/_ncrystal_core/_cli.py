

def ncrystal_cfg_wrapper():
    import subprocess
    import pathlib
    f = pathlib.Path(__file__).parent / 'ncrystal_pyinst_data' / 'bin' / 'ncrystal-config'
    import sys
    a = sys.argv[:]
    a[0] = f
    rv = subprocess.run( a )
    raise SystemExit(rv.returncode)
