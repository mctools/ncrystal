#!/usr/bin/env python3

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

# NEEDS: cmake

def get_downstream_project_files():
    from NCTestUtils.dirs import test_data_dir as td
    return [
        ( 'CMakeLists.txt',
          td.joinpath('downstream_cmake_CMakeLists.txt').read_text() ),
        ( 'main.cc',
          td.joinpath('downstream_cmake_main.cc').read_text() ),
    ]

def get_ncrystal_shlibdir():
    import subprocess
    import pathlib
    rv = subprocess.run(['ncrystal-config','--show','shlibdir'],
                        check = True,
                        capture_output=True)
    assert rv.returncode == 0 and not rv.stderr
    ld = pathlib.Path(rv.stdout.decode().strip())
    assert ld.is_dir()
    return ld

def main():
    from NCTestUtils.common import work_in_tmpdir
    import shutil
    import subprocess
    import pathlib
    import shlex
    import platform
    import os
    is_osx = False
    is_win = False
    if platform.system() == 'Darwin':
        is_osx = True
    elif platform.system() == 'Windows':
        is_win = True

    print('starting testing of compiled downstream cmake-based projects')
    cmakecmd = shutil.which('cmake')
    if not cmakecmd:
        raise RuntimeError('cmake command not found')
    cmake_extra_args = shlex.split(os.environ.get('CMAKE_ARGS',''))

    print('Using cmake command: %s'%cmakecmd)
    ncrystal_libdir = get_ncrystal_shlibdir()

    def prepend_to_path( envdict, pathvar, entry, sep = ':' ):
        entry = str(entry)
        if not entry:
            return
        ll = [ entry ]
        orig = envdict.get(pathvar,'')
        if orig:
            ll.append( orig )
        envdict[pathvar] = sep.join(ll)
    runtime_env = os.environ.copy()
    if is_win:
        prepend_to_path( runtime_env, 'PATH', ncrystal_libdir, sep = ';' )
    elif is_osx:
        prepend_to_path( runtime_env, 'DYLD_LIBRARY_PATH', ncrystal_libdir )
    else:
        prepend_to_path( runtime_env, 'LD_LIBRARY_PATH', ncrystal_libdir )
    print('Testing with libdir added to relevant path variable')

    with work_in_tmpdir():
        td = pathlib.Path('.').absolute()
        blddir = td / 'bld'
        srcdir = td / 'src'
        instdir = td / 'install'
        srcdir.mkdir()
        blddir.mkdir()
        for fn, content in get_downstream_project_files():
            srcdir.joinpath(fn).write_text(content)
        cmd = ( [cmakecmd]
                + cmake_extra_args
                + ['-B', str(blddir),'-S',str(srcdir),
                   '-DCMAKE_INSTALL_PREFIX=%s'%instdir] )
        subprocess.run( cmd, check=True )
        subprocess.run( [cmakecmd,
                         '--build',str(blddir),
                         '--config','Release'], check=True )
        subprocess.run( [cmakecmd,
                         '--install',str(blddir),
                         '--config','Release'], check=True )

        app = instdir.joinpath('bin',
                               ('testapp.exe' if is_win else 'testapp') )
        if not app.exists():
            raise RuntimeError('Could not produce test app'
                               ' in downstream CMake project')
        print(f'Launching: {app}')
        if is_osx:
            #osx's system protection #$%#^ is so damn annoying:
            cmd=shlex.quote(str(app))
            v = runtime_env.get('DYLD_LIBRARY_PATH')
            if v is not None:
                v = shlex.quote(v)
                cmd = f'export DYLD_LIBRARY_PATH={v} && {cmd}'
            subprocess.run(cmd,shell=True,check=True,env=runtime_env)
        else:
            subprocess.run([str(app)],check=True,env=runtime_env)
        print('App ran OK')
    print('testing of compiled downstream cmake-based projects done')

if __name__ == '__main__':
    main()
