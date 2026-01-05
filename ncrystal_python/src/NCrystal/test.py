
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

"""Module which is intended to trigger the built-in test from the command line
by running python3 -m NCrystal.test"""

from ._testimpl import ( test, test_cmdline, test_all ) # noqa F401

if __name__ == '__main__':
    from . import _testimpl as _tests
    import sys
    args=set(sys.argv[1:])
    if '-h' in args or '--help' in args:
        print('Run tests. By default the standard tests if no keywords'
              ' (std/cmdline/all)\nare supplied (use "verbose"/"quiet" for'
              ' printout control).')
        raise SystemExit
    do_quiet = ('quiet' in args)
    do_verbose = ( not do_quiet) and ('verbose' in args)
    do_cmdline = ('cmdline' in args)
    do_all = ('all' in args)
    do_std = ('std' in args)

    if not do_cmdline and not do_all:
        do_std = True

    if do_all:
        do_std = True
        do_cmdline = True

    unknown = args - set(['cmdline','all','std','verbose','quiet'])
    for u in unknown:
        raise SystemExit('Unknown keyword: %s'%u)

    test_kwargs = dict( verbose = ('quiet' if do_quiet else do_verbose) )

    if do_std:
        _tests.test( **test_kwargs )

    if do_cmdline:
        _tests.test_cmdline( **test_kwargs )
