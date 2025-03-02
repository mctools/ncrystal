
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

def script_ok_as_pure_python( name, fh ):
    for line in fh:
        line = line.split('#',1)[0]
        if 'NCTestUtils.loadlib' in line:
            return False
    return True

def filter_py_line( line ):
    if 'NCrystalDev' not in line:
        return line
    if 'import NCrystalDev.' in line:
        line = line.replace('import NCrystalDev.','import NCrystal.')
    elif 'from NCrystalDev.' in line:
        line = line.replace('from NCrystalDev.','from NCrystal.')
    elif 'import NCrystalDev' in line:
        line =  line.replace('import NCrystalDev','import NCrystal')
    assert 'NCrystalDev' not in line
    return line

def read_pyfile_filtered( f ):
    return ''.join( filter_py_line(e) for e in
                    f.read_text().splitlines(keepends=True) )

def handle_scripts( testinfo, tgtdir = None ):
    if tgtdir:
        import builtins
        print = builtins.print
        tgtdir.mkdir( parents = True, exist_ok = True )
    else:
        def print(*a,**kw):
            pass
    needs = set()
    for name, info in sorted(testinfo['scripts'].items()):
        if name=='ncdevtoolchecks':
            print(f' .. skipping script: {name}')
            continue
        with info['pyfile'].open('rt') as fh:
            if not script_ok_as_pure_python(name, fh):
                print(f' .. skipping non-pure script: {name}')
                continue
            print(f" .. using script: {name}")
            if not tgtdir:
                #Don't write anything, just gather needs
                for line in info['pyfile'].read_text().splitlines():
                    if line.startswith('# NEEDS: '):
                        needs.update(set(line[len('# NEEDS: '):].split()))
                        break
                continue
            tgtdir.joinpath(f'{name}.py').write_text(
                read_pyfile_filtered( info['pyfile'] )
            )
            if info['logfile']:
                tgtdir.joinpath(f'{name}.log').write_text(
                    info['logfile'].read_text()
                )
    return needs

def handle_pymods( testinfo, tgtdir ):
    for modname, pyfiles in sorted(testinfo['testmods'].items()):
        d = tgtdir.joinpath(modname)
        d.mkdir(parents=True)
        for f in sorted(pyfiles):
            dest = d.joinpath(f.name)
            assert f.name.endswith('.py')
            m = '%s.%s'%(modname,f.stem)
            if modname=='NCTestUtils' and f.name=='loadlib.py':
                print(f" .. ignoring module: {m}")
                continue#Ignore
            if modname=='NCTestUtils' and f.name=='enable_fpe.py':
                print(f" .. forcing empty module: {m}")
                dest.touch()#make this module noop
                continue
            print(f" .. copying module: {m}")
            dest.write_text( read_pyfile_filtered( f ) )

def handle_datadirs( testinfo, tgtdir ):
    for subdir, f in sorted( testinfo['datafiles'] ):
        subpath = '%s/%s'%(subdir,f.name) if subdir else f.name
        print(f" .. copying data file: {subpath}")
        d = tgtdir.joinpath(subdir) if subdir else tgtdir
        d.mkdir(parents=True,exist_ok=True)
        d.joinpath( f.name ).write_bytes( f.read_bytes() )

def safe_copy_file( f, targetdir ):
    tf = targetdir.joinpath(f.name)
    if tf.exists():
        raise SystemExit(f'ERROR: Multiple sources for {tf}')
    targetdir.mkdir(parents=True,exist_ok=True)
    tf.write_bytes(f.read_bytes())

def handle_extra_scripts( srcdir, targetdir ):
    for f in srcdir.glob('*.py'):
        safe_copy_file( f, targetdir )
    for f in srcdir.glob('*.log'):
        safe_copy_file( f, targetdir )

def handle_extra_data( srcdir, targetdir ):
    for f in srcdir.glob('*'):
        if '~' in f.name or '#' in f.name:
            continue
        safe_copy_file( f, targetdir )

def main():
    import sys
    import pathlib
    srcroot = pathlib.Path(__file__).parent
    reporoot = srcroot.parent
    sys.path.insert(0,str(reporoot.joinpath('devel','pypath')))
    from ncrystal_repo_tools.testinfo import load as testinfo_load
    testinfo = testinfo_load()

    if '--gather-scripts-deps' in sys.argv[1:]:
        #Special mode needed to verify dependency list in pyproject.toml
        needs = handle_scripts( testinfo, None )
        import json
        print( json.dumps( sorted(needs) ) )
        return

    tgt = srcroot.joinpath('skbld_autogen','_ncrystal_verify')
    assert not tgt.exists(), f"Assume cleaned target dir {tgt}"
    tgt.mkdir(parents=True)
    tgt.joinpath('__init__.py').touch()
    tgt.joinpath('_cli.py').write_text( srcroot.joinpath('cli.py').read_text() )

    handle_scripts( testinfo, tgt.joinpath('data','scripts') )
    handle_pymods( testinfo, tgt.joinpath('data','pypath') )
    handle_datadirs( testinfo, tgt.joinpath('data','data') )

    handle_extra_scripts( srcroot.joinpath('extra','scripts'),
                          tgt.joinpath('data','scripts') )
    handle_extra_data( srcroot.joinpath('extra','data'),
                       tgt.joinpath('data','data') )

if __name__ == '__main__':
    main()
