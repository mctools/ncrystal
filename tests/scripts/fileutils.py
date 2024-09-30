#!/usr/bin/env python3

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

from NCTestUtils.loadlib import Lib
lib = Lib('testfileutils')
lib.dump()
assert hasattr(lib,'nctest_file_exists')
assert hasattr(lib,'nctest_ncgetcwd')

import pathlib
d = pathlib.Path('.')
(d / 'foo.txt').write_text('foo')

assert lib.nctest_file_exists('./foo.txt')
assert lib.nctest_file_exists('foo.txt')
assert not lib.nctest_file_exists('bar.txt')
assert not lib.nctest_file_exists('./bar.txt')
assert not lib.nctest_file_exists('/some/file/that/does/not/exist')
print("Testing that two files are identical:")
_cwd = lib.nctest_ncgetcwd()
d_abs = d.resolve().absolute()
assert d.samefile(d_abs)
print("  1) %s"%d_abs)
print("  2) %s"%_cwd)
assert d_abs.samefile ( pathlib.Path(_cwd).resolve().absolute() )
assert d_abs.samefile ( _cwd )
assert d.samefile ( _cwd )

#FIXME : Step into dir with on-ascii name and check ncgetcwd ??

#print('cwd',lib.nctest_ncgetcwd())

#FIXME: Much more, including globbing, is_absolute_path("c:[\]bla.ncmat) and
#whatever is likely to cause issues on Windows.
