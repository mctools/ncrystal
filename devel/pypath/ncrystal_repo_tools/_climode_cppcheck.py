
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

def short_description():
    return 'Run cppcheck during a CMake build'

def main( parser ):
    parser.init( 'Run cppcheck during a CMake build.' )
    parser.add_argument(
        '-d','--debug', action='store_true',
        help="""Use CMAKE_BUILD_TYPE=Debug (default is Release)."""
    )
    parser.add_argument(
        '-j',type=int,default=0,
        help="""Use this many processes (default: auto detect)."""
    )
    args = parser.parse_args()

    from .util import work_in_tmpdir
    from .cmake import CMakeRunner
    from .util import get_nprocs
    import shutil
    if not shutil.which('cppcheck'):
        raise SystemExit('Error: cppcheck command missing')

    def cppcheck_arg(lang):
        assert lang in ('c','cxx')
        cmd = [ 'cppcheck',
                '--error-exitcode=3',
                '--inline-suppr',
                '--check-level=exhaustive',
                #'--quiet','--verbose',
                '--language=%s'%({'c':'c','cxx':'c++'}[lang])
               ]
        return ('-DCMAKE_%s_CPPCHECK='%(lang.upper())) + ';'.join(cmd)

    runner_args = dict( mode = 'ctest',
                        build_types = ['dbg' if args.debug else 'rel'],
                        generator = 'single',
                        nprocs_bld = args.j or get_nprocs(),
                        blddir = './build',
                        instdir = './install',
                        cmake_flags = [ cppcheck_arg('c'),
                                        cppcheck_arg('cxx'),
                                        '-DNCRYSTAL_ENABLE_TESTING=ON',
                                        '-DNCRYSTAL_ENABLE_EXAMPLES=ON',
                                        #Need no ncmat2cpp since cppcheck hangs
                                        #on the autogenerated .cc file:
                                        '-DNCRYSTAL_ENABLE_DATA=ON',
                                       ]
                       )

    with work_in_tmpdir():
        c = CMakeRunner( **runner_args )
        c.do_cfg()
        c.do_build()
