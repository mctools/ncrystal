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

import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.loadlib import Lib
from NCTestUtils.common import explicit_unicode_str
import pathlib

lib = Lib('testfileutils')
lib.dump()
assert hasattr(lib,'nctest_file_exists')
assert hasattr(lib,'nctest_ncgetcwd')
assert hasattr(lib,'nctest_readEntireFileToString')
assert hasattr(lib,'nctest_ncglob')

_raw_ncglob = lib.nctest_ncglob
def nctest_ncglob( pattern ):
    return _raw_ncglob(pattern).split('<<@>>')
lib.nctest_ncglob = nctest_ncglob

def test1():
    d = pathlib.Path('.')
    (d / 'foo.txt').write_text('foo')

    assert lib.nctest_file_exists('./foo.txt')
    assert lib.nctest_file_exists('foo.txt')
    assert not lib.nctest_file_exists('bar.txt')
    assert not lib.nctest_file_exists('./bar.txt')
    assert not lib.nctest_file_exists('/some/file/that/does/not/exist')
    print("Testing ncgetcwd...")
    _cwd = lib.nctest_ncgetcwd()
    d_abs = d.resolve().absolute()
    assert d.samefile(d_abs)
    assert d_abs.samefile ( pathlib.Path(_cwd).resolve().absolute() )
    assert d_abs.samefile ( _cwd )
    assert d.samefile ( _cwd )
    print("...ncgetcwd ok")

def test2():
    dirname = 'unicodedir_test\u4500abc'
    subdir = pathlib.Path('.') / dirname
    subdir.mkdir()
    for fn in [ 'a.txt', 'b.txt', 'b\u2030.ncmat' ]:
        (subdir / fn).write_text('dummy')
    testtext = 'some\nmulti\nline funny \u2030 text\n'
    print("Writing %i chars"%len(testtext))

    (subdir / 'b\u2030.ncmat').write_text(testtext,encoding='utf8')#NCrystal
                                                                   #only loads
                                                                   #utf8 encoded
                                                                   #text files!

    #NB: this fails, since cp1252 (windows Western) is incomplete:
    #'\u2030'.encode('cp1252') So on windows (where that encoding is common) we
    #can not create a filename with '\u2030' inside. Of course, anything we pick
    #might have trouble on some machine, if set to e.g. greek, japanese, or
    #whatever. For now, we focus on what works on github runners (cp1252).

    content = lib.nctest_readEntireFileToString(
        'unicodedir_test\u4500abc/b\u2030.ncmat'
    )
    print("Read %i chars"%len(content))
    for e in content.splitlines():
        print('READ>',repr(explicit_unicode_str(e)))
    assert testtext == content

    def testglob(pattern,nexpect):
        g = lib.nctest_ncglob(pattern)
        assert len(g) == nexpect
        g = list(explicit_unicode_str(e) for e in g)
        print(f'ncglob({repr(explicit_unicode_str(pattern))}) -> {repr(g)}')
    testglob('unicode*/*.ncmat',1)
    testglob('unicode*/*.txt',2)
    testglob('*abc/*\u2030*',1)

def main():


    test1()
    test2()

if __name__=='__main__':
    main()


#FIXME: Much more, including globbing, is_absolute_path("c:[\]bla.ncmat) and
#whatever is likely to cause issues on Windows.
