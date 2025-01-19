
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

def short_description():
    return 'Build and test code with CMake'

def main( parser ):
    from .cmake import cmakerunner_modes
    parser.init( 'Launch CMake to build the files in ncrystal_core.' )

    parser.add_argument(
        '-b','--build-dir', metavar='DIR',
        help="""CMake build directory to use. By default will work in temporary
        directory which is cleaned up afterwards."""
    )
    parser.add_argument(
        '-i','--install-dir', metavar='DIR',
        help="""CMake install directory to use. By default will work in temporary
        directory which is cleaned up afterwards."""
    )
    parser.add_argument(
        '-f','--force', action='store_true',
        help="""This flag will force any existing build and install dirs to be
        removed before running, instead of ending in an error."""
    )
    parser.add_argument(
        '--multi', action='store_true',
        help="""Will use multi generator mode even on unix (requires Ninja)."""
    )
    parser.add_argument(
        '--dbg', action='store_true',
        help="""Will use Debug mode rather than Release."""
    )
    assert 'ctest' in cmakerunner_modes
    parser.add_argument(
        '-m','--mode', type = str, choices = cmakerunner_modes,
        default = 'ctest',
        help="""Mode (default: ctest)."""
    )
    parser.add_argument(
        '-s','--strict', type = str,
        choices = ('OFF','ON','11','14','17','20','23'),
        default = '11',
        help="""BUILD_STRICT mode (default: '11')."""
    )
    parser.add_argument(
        '-j',type=int,default=0,
        help="""Use this many processes (default: auto detect)."""
    )
    args = parser.parse_args()

    from .util import get_nprocs
    nprocs = args.j or get_nprocs()

    needs_tmpdir = (not args.install_dir) or ( not args.build_dir )
    runner_args = dict( force = args.force,
                        mode = args.mode,
                        cmake_flags = None,
                        build_types = ['rel'],
                        nprocs_bld = nprocs,
                        nprocs_ctest = nprocs )
    if args.multi:
        runner_args['generator'] = 'multi'
    if args.dbg:
        runner_args['build_types'] = ['dbg']

    def do_work( blddir, instdir, mode ):
        from .cmake import CMakeRunner
        c = CMakeRunner( blddir = blddir,
                         instdir = instdir,
                         **runner_args )
        c.do_cfg()
        c.do_build()
        if mode == 'ctest':
            c.do_ctest()
        elif mode == 'install':
            c.do_install()
            c.do_test_install()

    if needs_tmpdir:
        from .util import work_in_tmpdir
        with work_in_tmpdir():
            do_work( args.build_dir or './build',
                     args.install_dir or '/install',
                     args.mode )
    else:
        do_work( args.build_dir, args.install_dir, args.mode )
