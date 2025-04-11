
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

from . import dirs
import subprocess
import os
import shlex
import shutil
import platform
from pathlib import Path as plPath

cmakerunner_modes = ('ctest','install','buildonly')

class CMakeRunner:

    def __init__( self, *,
                  blddir,
                  instdir,
                  mode,
                  cmake_flags = None,
                  build_types = None,
                  generator = None,
                  force = False,
                  cmake_cmd = None,
                  ctest_cmd = None,
                  nprocs = None,
                  nprocs_bld = None,
                  nprocs_ctest = None,
                 ):
        #either CTest or install+test:
        if generator:
            assert generator in ('multi','single')
        else:
            generator = 'multi' if platform.system()=='Windows' else 'single'
        self.generator = generator

        if nprocs is not None:
            assert nprocs_bld is None and nprocs_ctest is None
            nprocs_ctest = nprocs
            nprocs_bld = nprocs
        self.nprocs_bld = int(nprocs_bld) if nprocs_bld else None
        self.nprocs_ctest = int(nprocs_ctest) if nprocs_ctest else None

        build_types = build_types or ['rel']
        assert all( bt in ('rel','reldbg','dbg') for bt in build_types )

        self.build_types = sorted(set(build_types))
        if self.generator == 'single':
            assert len(self.build_types)==1
        assert mode in cmakerunner_modes
        self.mode = mode
        self.cmake_flags = ( shlex.split(os.environ.get('CMAKE_ARGS',''))
                             + ( cmake_flags or [] ) )
        self.blddir = plPath(blddir).absolute()
        self.instdir = plPath(instdir).absolute()
        from .util import path_is_relative_to
        assert not path_is_relative_to( self.instdir, self.blddir )
        assert not path_is_relative_to( self.blddir, self.instdir )
        self.force = force
        self.stage = 'none'
        self.cmake_cmd = cmake_cmd or shutil.which('cmake')
        if self.mode == 'ctest':
            self.ctest_cmd = ctest_cmd or shutil.which('ctest')

    def _invoke( self, cmdname, arg_list, cwd = None, env = None,
                 capture_output = False ):
        assert not isinstance(arg_list,str)
        assert plPath(self.cmake_cmd).is_file()
        cmd = [ str(cmdname) ] + [ str(e) for e in arg_list ]
        cmdstr = shlex.join( cmd )
        print()
        print(f'LAUNCHING: {cmdstr}')
        print()
        rv = subprocess.run( cmd,
                             cwd=cwd,
                             env=env,
                             capture_output=capture_output )
        if rv.returncode!=0:
            raise RuntimeError(f'Command failed: {cmdstr}')
        if capture_output:
            if rv.stderr:
                print(rv.stderr)
            return rv.stdout

    def _bt2cmakebt( self, buildtype ):
        return { 'rel' : 'Release',
                 'dbg' : 'Debug',
                 'reldbg' : 'RelWithDebInfo' }[buildtype]

    def _determine_nprocs( self, nprocs_requested, envvar ):
        if nprocs_requested is not None:
            return max(nprocs_requested, 1)
        ev = os.environ.get(envvar)
        if ev is not None:
            assert ev.isdigit()
            return max(int(ev), 1)
        from .utils import get_nprocs
        return get_nprocs()

    def _msg( self, *a, **kw ):
        print('>>>> CMakeRunner:',*a,**kw)

    def _prepare_dirs( self ):
        assert self.stage=='none'
        for dirname, thedir in ( ('Build',self.blddir),
                                 ('Install',self.instdir) ):
            if thedir.exists() and not dirs.is_empty_dir(thedir):
                if self.force:
                    if not dirs.is_empty_dir( thedir ):
                        self._msg(f'forcefully removing {thedir}')
                        shutil.rmtree(thedir)
                else:
                    raise RuntimeError(f'{dirname} dir already exists: {thedir}')

    def do_cfg( self ):
        assert self.stage=='none'
        self._prepare_dirs()
        args = [ '-S',dirs.reporoot,'-B', self.blddir]
        if self.generator=='multi' and platform.system()!='Windows':
            args += ['-G','Ninja Multi-Config']
        elif self.generator=='single' and platform.system()=='Windows':
            args += ['-G','Ninja']
        self._msg(f'using build dir {self.blddir}')
        if self.mode == 'install':
            args.append('-DCMAKE_INSTALL_PREFIX=%s'%self.instdir)
            self._msg(f'using install dir {self.blddir}')
        elif self.mode == 'buildonly':
            pass
        else:
            assert self.mode == 'ctest'
            args.append('-DNCRYSTAL_ENABLE_TESTING=ON')
            args.append('-DNCRYSTAL_SKIP_INSTALL=ON')
            args.append('-DMCTOOLS_REQUIRE_ALL_TEST_DEPS=ON')

        if self.generator=='single':
            assert len(self.build_types)==1
            args.append('-DCMAKE_BUILD_TYPE=%s'%self._bt2cmakebt(self.build_types[0]))

        args += self.cmake_flags
        self._invoke( self.cmake_cmd, args )
        self.stage='cfg'

    def do_build( self ):
        assert self.stage=='cfg'
        nprocs = self._determine_nprocs( self.nprocs_bld,
                                         'CMAKE_BUILD_PARALLEL_LEVEL' )
        args = [ '--build', self.blddir,
                 '--parallel', nprocs ]
        for bt in self.build_types:
            self._invoke( self.cmake_cmd,
                          args + ['--config',self._bt2cmakebt(bt)] )
        self.stage = 'bld'

    def do_ctest( self, extra_args = None ):
        assert self.mode == 'ctest'
        assert self.stage == 'bld'
        nprocs = self._determine_nprocs( self.nprocs_ctest,
                                         'CTEST_PARALLEL_LEVEL' )
        args = [ '--output-on-failure',
                 '--no-tests=error',
                 '--test-output-size-failed', '10000',
                 '--test-output-truncation', 'middle',
                 '--parallel',nprocs,
                ] + ( extra_args or [] )

        for bt in self.build_types:
            self._invoke( self.ctest_cmd,
                          args + ['--build-config',self._bt2cmakebt(bt)],
                          cwd = self.blddir )

    def do_install( self ):
        assert self.mode == 'install'
        assert self.stage == 'bld'
        args = [ '--install', self.blddir ]
        for bt in self.build_types:
            self._invoke( self.cmake_cmd,
                          args + ['--config',self._bt2cmakebt(bt)] )
        self.stage = 'install'

    def do_test_install( self ):
        #Minimal test.
        import pathlib
        assert self.stage == 'install'
        assert self.mode == 'install'
        self._invoke( self.instdir/'bin/ncrystal-config', [ '-s' ] )
        for p in ( 'libpath', 'shlibpath' ):
            o = self._invoke( self.instdir/'bin/ncrystal-config', [ '--show',p ],
                              capture_output = True)
            lp = pathlib.Path(o.decode().strip())
            if not lp.is_file():
                raise SystemExit(f'File not found: {lp}')
            else:
                print("Verified existence of file at: "
                      f"ncrystal-config  --show {p}")
