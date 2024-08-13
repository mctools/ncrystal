"""

Module which is intended to trigger the built-in test from the command line by running python3 -m NCrystal.test

"""

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

from ._testimpl import test, test_cmdline, test_cmake, test_extra, test_all

if __name__ == '__main__':
    from . import _testimpl as _tests
    import sys
    args=set(sys.argv[1:])
    if '-h' in args or '--help' in args:
        print('Run tests. By default the standard tests if no keywords (std/cmake/cmdline/extra/all)\n'
              'are supplied (use "verbose"/"quiet" for printout control).')
        raise SystemExit
    do_quiet = ('quiet' in args)
    do_verbose = ( not do_quiet) and ('verbose' in args)
    do_cmake = ('cmake' in args)
    do_cmdline = ('cmdline' in args)
    do_extra = ('extra' in args)
    do_all = ('all' in args)
    do_std = ('std' in args)

    if not do_cmake and not do_cmdline and not do_extra and not do_all:
        do_std = True

    unknown = args - set(['cmake','cmdline','extra','all','std','verbose','quiet'])
    for u in unknown:
        raise SystemExit('Unknown keyword: %s'%u)

    test_kwargs = dict( verbose = ('quiet' if do_quiet else do_verbose) )
    if do_extra and do_cmake:
        do_all = True

    if ( do_cmake and do_std and do_cmdline ):
        do_all = True

    if do_all:
        _tests.test( **test_kwargs )
        _tests.test_cmdline( **test_kwargs )
        _tests.test_cmake( **test_kwargs )
        raise SystemExit

    if do_extra:
        _tests.test_extra( **test_kwargs )
        do_std = False
        do_cmdline = False

    if do_std:
        _tests.test( **test_kwargs )

    if do_cmdline:
        _tests.test_cmdline( **test_kwargs )

    if do_cmake:
        _tests.test_cmake( **test_kwargs )
