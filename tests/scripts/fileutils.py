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

import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.loadlib import Lib
from NCTestUtils.common import ( ensure_error,
                                 explicit_unicode_str,
                                 is_windows )
import pathlib

lib = Lib('testfileutils')
lib.dump()
assert hasattr(lib,'nctest_file_exists')
assert hasattr(lib,'nctest_ncgetcwd')
assert hasattr(lib,'nctest_readEntireFileToString')
assert hasattr(lib,'nctest_ncglob')

_raw_ncglob = lib.nctest_ncglob
def nctest_ncglob( pattern ):
    res = _raw_ncglob(pattern)
    return res.split('<<@>>') if res else []
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
    testpaths = [ r'C:\somewhere\bla.txt',
                  '/somewhere/bla.txt',
                  r'C:somewhere\bla.txt',
                  r'somewhere\bla.txt',
                  r'\somewhere\bla.txt',
                  'somewhere/bla.txt',
                  '.',
                  '',
                  '/',
                  './'#os.path.basename of this is empty
                  '\\',
                  '/some/where//\\bla.txt',
                  'some\\where/bla.txt',
                  '/some/where/../bla.txt',
                  'some/where/../bla.txt',
                  '../bla.txt',
                  '/test/bla/'#pathlib.Path(..).name is 'bla' but
                              #os.path.basename is ''
                 ]

    fctnames = ['dirname','basename','getfileext',
                'path_is_absolute','normalise']
    for i,p in enumerate(testpaths):
        print(f'{i} Testing path "{p}":')
        for fct in fctnames:
            res = getattr(lib,f'nctest_{fct}')(p)
            if is_windows() and fct in ('dirname','normalise'):
                #ensure test reproducibility:
                assert '/' not in res, "normalised paths should not contain / on windows"
                is_windows_path = (len(p)>1 and p[1]==':')
                if not is_windows_path:
                    res = res.replace('\\','/')
            print(f'   NCrystal::{fct} = "{res}"')
        #test versus refs, but pathlib and os.path does not like back slashes on
        #unix:

        nc_basename = lib.nctest_basename(p)
        def decode_refbn(p):
            is_windows_path = '\\' in p or (len(p)>1 and p[1]==':')
            if is_windows_path:
                p = p.replace('/','\\')
                return dict(pathlib=pathlib.PureWindowsPath(p).name)
            else:
                p = p.replace('\\','/')
                return dict(pathlib=pathlib.PurePosixPath(p).name)
        for refsrc, refbasename in decode_refbn(p).items():
            if refbasename != nc_basename:
                raise SystemExit(f"basename({repr(p)}) mismatch:"
                                 f" ncrystal={repr(nc_basename)}"
                                 f" {refsrc}={repr(refbasename)}")

        assert ( lib.nctest_path_is_absolute( lib.nctest_normalise( p ) )
                 == lib.nctest_path_is_absolute( p ) ), "normalisation alters is_absolute"
        print()


def test3():
    dirname_simple = 'some_sub_dir'
    subdir_simple = pathlib.Path('.') / dirname_simple
    subdir_simple.mkdir()
    for fn in [ 'a.txt', 'b.txt', 'b_bla.ncmat' ]:
        (subdir_simple / fn).write_text('dummy')

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
    content_newlines_norm = content.replace('\r\n','\n')
    print("Read %i chars (after newline normalisation)"%len(content_newlines_norm))
    for e in content.splitlines():
        print('READ>',repr(explicit_unicode_str(e)),flush=True)
    assert testtext == content_newlines_norm
    if not is_windows():
        #this breaks on windows since readEntireFileToString does not discard
        #extra newline chars (they are anyway supported in .ncmat data, so not
        #really a deal breaker).
        assert testtext == content

    def testglob(pattern,nexpect,is_absolute=False):
        def fmtpat( pattern ):
            if is_absolute:
                assert pathlib.Path(pattern).is_absolute()
                s = '/SOMEWHERE/'+pathlib.Path(pattern).name
            else:
                s = str(pattern)
            if is_windows():
                s=s.replace('\\','/')#for test reproducibility
            return repr(explicit_unicode_str(s))
        print(f'Testing ncglob({fmtpat(pattern)}):')
        import glob
        g = sorted(glob.glob(pattern))
        print('  --> Python glob got %i results:'%len(g))
        for e in g:
            print(f'   *: {fmtpat(e)}')
        g = lib.nctest_ncglob(pattern)
        print('  --> got %i results:'%len(g))
        all_ok = True
        for e in g:
            badstr = ''
            if not pathlib.Path(e).exists():
                badstr = ' (ERROR NOT ACTUALLY FOUND)'
                all_ok = False
            print(f'   *: {fmtpat(e)}{badstr}')
        if len(g) != nexpect:
            raise SystemExit('Error: did not get expected %i entries'%nexpect)
        if not all_ok:
            raise SystemExit('Error: some hits were invalid')

    #windows: testglob(r'some_sub_dir\*.txt',2)
    testglob('some_sub_dir/*.txt',2)
    testglob('some_sub_dir/*',3)
    testglob('some_sub_dir/***',3)
    testglob('some_sub_dir/*txt*',2)
    testglob('some_sub_dir/?.txt*',2)
    testglob('some_sub_dir/*.ncmat',1)
    testglob('some_sub_dir/*b*',2)
    testglob('some_sub_dir/*a*',2)
    testglob('some_sub_dir/*.ncm',0)#could catch the .ncmat file if using wrong
                                    #windows API, where only 3 characters in
                                    #extensions are compared.
    testglob('unicodedir_test\u4500abc/*.ncmat',1)
    testglob('unicodedir_test\u4500abc/*.txt',2)
    testglob('unicodedir_test\u4500abc/*\u2030*',1)

    testglob('unicodedir_test\u4500ab*',1)
    testglob('unicodedir_test\u4500abc*',1)
    testglob('unicodedir_*',1)

    import os
    testglob(os.path.join(
        str(pathlib.Path('unicodedir_test\u4500abc').absolute()),
        '*') ,3,is_absolute = True)

    print("Triggering expected error (hopefully):",flush=True)
    with ensure_error(RuntimeError,
                      'ncglob only supports wildcards in the'
                      ' last file or directory name'):
        lib.nctest_ncglob('some*sub_dir/a.txt')

    #Finally, test with paths longer than 260 chars:
    dummy_path = pathlib.Path('a'*150).joinpath('b'*150).joinpath('bla.txt')
    dummy_path.parent.parent.mkdir()
    dummy_path.parent.mkdir()
    dummy_path.write_text('hello')
    testglob( str(dummy_path.parent)+'/bl*.txt',1)

def test4():
    def testbn(path,expected):
        res=lib.nctest_basename(path)
        if res != expected:
            raise SystemExit(f'NCrystal::basename({repr(path)}) gave'
                             f' {repr(res)} and not the'
                             f' expected {repr(expected)}.')
    testbn("","")
    testbn("hej","hej")
    testbn("hej.txt","hej.txt")
    testbn("hej.txt.0","hej.txt.0")
    testbn("./","")
    testbn("./hej","hej")
    testbn("./hej.txt","hej.txt")
    testbn("./hej.txt.0","hej.txt.0")
    testbn("/lala/bla/","bla")#NB: pre NCrystal 4.0 this gave ""
    testbn("/lala/bla/hej","hej")
    testbn("/lala/bla/hej.txt","hej.txt")
    testbn("/lala/bla/hej.txt.0","hej.txt.0")
    testbn("~lala/bla/","bla")#NB: pre NCrystal 4.0 this gave ""
    testbn("~lala/bla/hej","hej")
    testbn("~lala/bla/hej.txt","hej.txt")
    testbn("~lala/bla/hej.txt.0","hej.txt.0")
    testbn("../bla/","bla")#NB: pre NCrystal 4.0 this gave ""
    testbn("../bla/hej","hej")
    testbn("../bla/hej.txt","hej.txt")
    testbn("../bla/hej.txt.0","hej.txt.0")

    def testext(path,expected):
        res=lib.nctest_getfileext(path)
        if res != expected:
            raise SystemExit(f'NCrystal::getfileext({repr(path)}) gave'
                             f' {repr(res)} and not the'
                             f' expected {repr(expected)}.')

    testext("../bla/hej.txt","txt")
    testext("../bla/hej.lala.txt","txt")
    testext("../bla/","")
    testext("../bla/lala","")
    testext("../bla/.txt","txt")

def main():
    test1()
    test2()
    test3()
    test4()

if __name__=='__main__':
    main()
