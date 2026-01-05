
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


def short_description():
    return 'Output number of processes available for multiprocessing'

def main( parser ):
    parser.init( short_description() )
    parser.add_argument(
        '--enable-github-parallel', action='store_true',
        help="""Appends CMAKE_BUILD_PARALLEL_LEVEL=<nprocs> to $GITHUB_ENV"""
    )

    args = parser.parse_args()
    from .util import get_nprocs
    n = get_nprocs( nice_factor = 1.0 )
    if args.enable_github_parallel:
        import os
        import pathlib
        p = pathlib.Path(os.environ['GITHUB_ENV'])
        kv = []
        kv.append( ('CMAKE_BUILD_PARALLEL_LEVEL',str(n) ) )
        #kv.append( ('CTEST_PARALLEL_LEVEL',str(n) ) )
        s = ''
        for k,v in kv:
            print(f"Appending {k}={v} to " "${GITHUB_ENV}")
            s += f'{k}={v}\n'
        with p.open('at') as fh:
            fh.write(s)
    else:
        print( n )
