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
from NCTestUtils.common import ( work_in_tmpdir,
                                 explicit_unicode_str,
                                 is_windows )
import pathlib
import os

#import builtins
#def print(*a,**kw):
#    sys.stdout.flush()
#    builtins.print(*a,**kw)
#    sys.stdout.flush()

native_sep = '\\' if is_windows() else '/'
nonnative_sep = '/' if is_windows() else '\\'

pPath = pathlib.Path
lib = Lib('testcfileutils')
lib.dump()
assert hasattr(lib,'nctest_get_current_working_dir')
assert hasattr(lib,'nctest_real_path')
assert hasattr(lib,'nctest_absolute_path')
assert hasattr(lib,'nctest_fopen_and_read_text')
assert hasattr(lib,'nctest_pathseps_platform')
assert hasattr(lib,'nctest_pathseps_generic')
assert hasattr(lib,'nctest_path_is_absolute')
assert hasattr(lib,'nctest_path_is_relative')
assert hasattr(lib,'nctest_drive_letter')
assert hasattr(lib,'nctest_basename')
assert hasattr(lib,'nctest_basename_view')
assert hasattr(lib,'nctest_dirname')
assert hasattr(lib,'nctest_fileextension')
assert hasattr(lib,'nctest_fileextension_view')
assert hasattr(lib,'nctest_path_join')
assert hasattr(lib,'nctest_is_same_file')
assert hasattr(lib,'nctest_is_file')
assert hasattr(lib,'nctest_is_dir')
assert hasattr(lib,'nctest_exists')
assert hasattr(lib,'nctest_expand_path')

def test1():
    d = pPath('.')
    (d / 'foo.txt').write_text('foo')

    assert lib.nctest_is_file('./foo.txt')
    assert lib.nctest_is_file(r'.\foo.txt')
    assert lib.nctest_is_file('foo.txt')
    assert not lib.nctest_is_file('bar.txt')
    assert not lib.nctest_is_file('./bar.txt')
    assert not lib.nctest_is_file('/some/file/that/does/not/exist')
    print("Testing get_current_working_dir...")
    _cwd = lib.nctest_get_current_working_dir()
    d_abs = d.resolve().absolute()
    assert d.samefile(d_abs)
    assert d_abs.samefile ( pPath(_cwd).resolve().absolute() )
    assert d_abs.samefile ( _cwd )
    assert d.samefile ( _cwd )
    print("...get_current_working_dir ok")

def test2():
    testpaths = [ 'hello.txt',
                  r'C:\somewhere\bla.txt',
                  '/somewhere/bla.txt',
                  r'c:somewhere\bla.txt',
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
                 ]

    fctnames = ['dirname','basename','basename_view',
                'fileextension','fileextension_view',
                'is_file','is_dir','exists',
                'path_is_absolute','path_is_relative',
                'pathseps_platform','pathseps_generic',
                'drive_letter','expand_path',
                ]

    fctnames_same_is_absolute = ['dirname','pathseps_platform',
                                 'pathseps_generic','expand_path']

    for i,p in enumerate(testpaths):
        print(f'{i} Testing path "{p}":')
        for fct in fctnames:
            res = getattr(lib,f'nctest_{fct}')(p)
            res_print = res
            if isinstance(res,str):
                res_print = res_print.replace(native_sep,'/')
                if len(res)>=2 and res[1]==':':
                    assert res[0] == res[0].upper()
                forbid_sep, allow_sep = nonnative_sep, native_sep
                if fct=='pathseps_generic':
                    forbid_sep, allow_sep = '\\', '/'
                    assert forbid_sep not in res
                    res = res.replace(allow_sep,'/')
            print(f'   mctools_{fct} = "{res_print}"')

        #test versus refs, but pathlib and os.path does not like back slashes on
        #unix:
        nc_basename = lib.nctest_basename(p)
        def decode_refbn(p):
            is_windows_path = '\\' in p or (len(p)>1 and p[1]==':')
            if is_windows_path:
                p = p.replace('/','\\')
                pobj = pathlib.PureWindowsPath(p)
            else:
                p = p.replace('\\','/')
                pobj=pathlib.PurePosixPath(p)
            return dict(pathlib=pobj.name)
        for refsrc, refbasename in decode_refbn(p).items():
            if refbasename != nc_basename:
                raise SystemExit(f"basename({repr(p)}) mismatch:"
                                 f" ncrystal={repr(nc_basename)}"
                                 f" {refsrc}={repr(refbasename)}")



        _isabs = lib.nctest_path_is_absolute( p )
        for fct in fctnames_same_is_absolute:
            pmod = getattr(lib,f'nctest_{fct}')(p)
            newisabs = lib.nctest_path_is_absolute( pmod )
            if _isabs != newisabs:
                raise SystemExit(f"ERROR: {fct}({p})={pmod} "
                                 "has changed is_absolute_path"
                                 f" from {_isabs} to {newisabs}")
        print()


def test3():
    dirname_simple = 'some_sub_dir'
    subdir_simple = pPath('.') / dirname_simple
    subdir_simple.mkdir()
    for fn in [ 'a.txt', 'b.txt', 'b_bla.ncmat' ]:
        (subdir_simple / fn).write_text('dummy')

    dirname = 'unicodedir_test\u4500abc'
    subdir = pPath('.') / dirname
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

    content = lib.nctest_fopen_and_read_text(
        'unicodedir_test\u4500abc/b\u2030.ncmat'
    )
    content_newlines_norm = content.replace('\r\n','\n')
    print("Read %i chars (after newline normalisation)"%len(content_newlines_norm))
    for e in content.splitlines():
        print('READ>',repr(explicit_unicode_str(e)),flush=True)
    assert testtext == content_newlines_norm
    if not is_windows():
        #this breaks on windows since file reading with fread does not discard
        #extra newline chars (they are anyway supported in .ncmat data, so not
        #really a deal breaker).
        assert testtext == content

def test4():

    def testdirname(path,expected):
        res=lib.nctest_dirname(path)
        assert nonnative_sep not in res
        expected = expected.replace('/',native_sep)
        if res != expected:
            raise SystemExit(f'NCrystal::dirname({repr(path)}) gave'
                             f' {repr(res)} and not the'
                             f' expected {repr(expected)}.')

    def testbn(path,expected):
        res=lib.nctest_basename(path)
        res_view=lib.nctest_basename_view(path)
        if res != res_view:
            raise SystemExit(f'NCrystal::basename({repr(path)})={repr(res)}'
                             ' and '
                             f' NCrystal::basename_view({repr(path)})'
                             f'={repr(res_view)} differ!!')
        if res != expected:
            raise SystemExit(f'NCrystal::basename({repr(path)}) gave'
                             f' {repr(res)} and not the'
                             f' expected {repr(expected)}.')

    testdirname(".",".")
    testbn(".","")

    testdirname("hej","")
    testbn("hej","hej")

    testdirname("hej/","hej")
    testbn("hej/","")

    testdirname("/hej","/")
    testbn("/hej","hej")

    testdirname("/","/")
    testbn("/","")

    testdirname("D:.","D:")
    testbn("D:./","")
    testdirname("D:./","D:")
    testbn("D:.","")
    testbn("D:","")
    testdirname("D:","D:")

    testbn("hej","hej")
    testbn("hej.txt","hej.txt")
    testbn("hej.txt.0","hej.txt.0")
    testbn("./","")
    testbn("./hej","hej")
    testbn("./hej.txt","hej.txt")
    testbn("./hej.txt.0","hej.txt.0")
    testbn("/lala/bla/","")
    testbn("/lala/bla/hej","hej")
    testbn("/lala/bla/hej.txt","hej.txt")
    testbn("/lala/bla/hej.txt.0","hej.txt.0")
    testbn("~lala/bla/","")
    testbn("~lala/bla/hej","hej")
    testbn("~lala/bla/hej.txt","hej.txt")
    testbn("~lala/bla/hej.txt.0","hej.txt.0")
    testbn("../bla/","")
    testbn("../bla/hej","hej")
    testbn("../bla/hej.txt","hej.txt")
    testbn("../bla/hej.txt.0","hej.txt.0")
    testbn("","")
    testbn("/","")
    testbn("//","")
    testbn("\\","")
    testbn("\\\\","")
    testbn('/test/bla/','')
    testbn('/test/bla','bla')

    def testext(path,expected):
        res=lib.nctest_fileextension(path)
        res_view=lib.nctest_fileextension_view(path)
        assert res == res_view
        if res != expected:
            raise SystemExit(f'NCrystal::fileextension({repr(path)}) gave'
                             f' {repr(res)} and not the'
                             f' expected {repr(expected)}.')

    testext("../bla/hej.txt","txt")
    testext("../bla/hej.lala.txt","txt")
    testext("../bla/","")
    testext("../bla/lala","")
    testext("../bla/.txt","txt")
    testext("~lala/bla/hej.txt.0","0")

def test5( workdir ):
    #exercise a bunch of things, including is_same_file, is_dir, real_path,
    #path_join.

    #Damn... we have to stop doing any symlink testing on windows...
    td = workdir
    files = [
        #We only test symlink to files, not directories here!
        dict(name='asimplefile.txt',content='hello'),
        dict(name='sd 1',is_dir=True),
        dict(name='sd 1/lala.txt',content='bla'),
        dict(name='sd 1/\u4500abc',is_dir=True),
        dict(name='sd 1/\u4500abc/.\u4500',content='foobar'),
        dict(name='sd 1/\u4500abc/.foo~',content='foobar'),
        #dict(name='asimplefilelinked.txt',symlink='asimplefile.txt'),
        #dict(name='sd 1/asimplefilelinked.txt',symlink='../asimplefile.txt'),
        dict(name='sd2',is_dir=True),
        #dict(name='sd2/sl2',symlink=str(td.joinpath('sd 1').absolute()),is_dir=True),
        dict(name='sd  __   3',is_dir=True),
        dict(name='sd  __   3/bla.ncmat',content='NCMAT'),
    ]
    #Create layout:
    for finfo in files:
        f = finfo['name']
        fabs = td.joinpath(f)
        assert not fabs.exists()
        symlink = finfo.get('symlink','')
        is_dir = bool(finfo.get('is_dir'))
        print(f'-->  Creating test file {repr(f)}'
              f' (is_dir={is_dir}, is_symlink={bool(symlink)})')
        if symlink:
            #might have is_dir True
            assert not is_dir, "this will almost certainly not work in windows"
            fabs.symlink_to(symlink,target_is_directory=is_dir)
            if not fabs.exists():
                raise SystemExit("Failed to create symlink!")
            assert fabs.exists()
            assert fabs.is_dir() == is_dir
            assert fabs.is_symlink()
            if pPath(symlink).is_absolute():
                assert fabs.samefile( pPath(symlink) )
            else:
                expected_linkdest = pPath(fabs.parent).joinpath(symlink)
                assert fabs.samefile( expected_linkdest )
        elif is_dir:
            fabs.mkdir()
            continue
        else:
            content = finfo.get('content')
            if content is not None:
                fabs.write_text(content)
            else:
                fabs.touch()
    #Analyse it:
    for ifile,finfo in enumerate(files):
        frawname = finfo['name']
        print( f'Testing file {ifile}: "{frawname}"' )
        f = pPath(frawname)
        fabs = f.absolute()
        fabs_resolved = f.resolve()

        abs_path0 = lib.nctest_absolute_path( fabs_resolved )#fixme test more?
        assert abs_path0
        assert lib.nctest_path_is_absolute(abs_path0)

        real_path0 = lib.nctest_real_path( fabs_resolved )

        if not real_path0:
            print( 'ERROR: real_path returns empty (fails) for:')
            print( f"  arg: {repr(str(fabs_resolved))} (fabs_resolved)")
            print( f"  fabs= {repr(str(fabs))}")
            print( f"  cwd= {repr(str(os.getcwd()))}")
            raise SystemExit(1)

        transforms = [ str,
                       lambda fin : str(fin).replace('/','\\'),
                       lambda fin : str(fin).replace('\\','/') ]
        is_dir = bool( fabs_resolved.is_dir() )
        assert is_dir == bool(finfo.get('is_dir'))
        fabs_prev = td.parent
        for is_abs, fin in [ (False,f), (True,fabs) ]:
            assert not lib.nctest_is_same_file( str(fabs_prev), str(fin) )
            for ftest in sorted(set([t(fin) for t in transforms])):
                ftest_shown = ftest
                if is_abs:
                    ftest_shown = '<hiddenabspart>'+ftest[-len(str(f)):]
                print( f'  Testing with form: "{ftest_shown}" (is_abs={is_abs})' )
                assert bool(lib.nctest_is_dir(ftest)) == is_dir
                assert bool(lib.nctest_path_is_absolute(ftest)) == bool(is_abs)
                assert bool(lib.nctest_path_is_relative(ftest)) == bool(not is_abs)

                assert lib.nctest_exists(ftest)
                assert lib.nctest_is_dir(ftest) == is_dir
                assert lib.nctest_is_file(ftest) == bool(not is_dir)

                if lib.nctest_is_same_file(ftest,str(fabs_resolved)) != bool(not is_dir):
                    raise SystemExit("is_same_file does not return True"
                                     f" for {repr(ftest)} vs "
                                     f"{repr(str(fabs_resolved))}")
                assert lib.nctest_is_same_file(ftest,str(fabs_resolved)) == bool(not is_dir)
                assert lib.nctest_is_same_file(str(fabs_resolved),ftest) == bool(not is_dir)

                rp_ftest = lib.nctest_real_path( ftest )
                if not real_path0 == rp_ftest:
                    print( 'ERROR: real_path(ftest) results unexpected for:')
                    print( f"  ftest: {repr(str(ftest))}")
                    print( f"  cwd= {repr(str(os.getcwd()))}")
                    print( f"  expected= {repr(str(real_path0))}")
                    print( f"  got     = {repr(str(rp_ftest))}")
                    raise SystemExit(1)
                assert real_path0 == rp_ftest

                abs_ftest = lib.nctest_absolute_path( ftest )
                assert pPath(abs_ftest).samefile(pPath(fabs))

                if not is_abs:
                    pj = lib.nctest_path_join(str(td),ftest)
                    pj2 = lib.nctest_path_join(
                        lib.nctest_dirname(str(fabs)),
                        lib.nctest_basename(str(fabs))
                    )
                    pj_native = lib.nctest_pathseps_platform(pj)
                    assert pPath(fabs_resolved).samefile(pPath(pj_native))
                    assert pPath(fabs_resolved).exists()
                    assert pPath(fabs_resolved).samefile( f )
                    assert pPath(fabs_resolved).samefile( pPath(pj2) )
                    assert pPath(fabs_resolved).samefile( fabs )
                    if not pPath(fabs_resolved).samefile( pPath(real_path0) ):
                        print( 'ERROR: pathlib.Path.samefile_file False for:')
                        print( f"  arg1: {repr(str(fabs_resolved))} (fabs_resolved)")
                        print( f"  arg2: {repr(real_path0)} (real_path0)")
                        print( f"  fabs= {repr(str(fabs))}")
                        print( f"  fin= {repr(str(fin))}")
                        print( f"  ftest= {repr(str(ftest))}")
                        print( f"  cwd= {repr(str(os.getcwd()))}")
                        raise SystemExit(1)
                    assert pPath(fabs_resolved).samefile( real_path0 )
        fabs_prev = fabs

def test_dirsymlinks( workdir ):
    #FIXME: Test symlink to dir here (with no printouts)
    f1 = workdir / 'bla.txt'
    f1.write_text('bla')
    fd = workdir / 'somedir'
    fd.mkdir()
    (workdir / 'yihadir').mkdir()
    f2 = fd / 'foo.bar'
    f2.symlink_to('../bla.txt')
    assert f2.name == 'foo.bar'
    assert f2.resolve().name == 'bla.txt'
    assert lib.nctest_basename(f2) == 'foo.bar'
    assert lib.nctest_basename(lib.nctest_real_path(f2)) == 'bla.txt'

    def test_symlink( f ):
        if f.exists():
            return True
        if is_windows():
            #Symlinks are not always available Windows, so simply skip this
            #test.
            return False
        else:
            raise SystemExit(f'failed to create symlink: {repr(f)}')

    def ensure_is_file(f):
        assert not lib.nctest_is_dir(f)
        assert lib.nctest_is_file(f)
        assert lib.nctest_exists(f)

    def ensure_is_dir(f):
        assert lib.nctest_is_dir(f)
        assert not lib.nctest_is_file(f)
        assert lib.nctest_exists(f)

    #a relative symlink to a file:
    f3 = fd / 'muahaha.bar'
    f3.symlink_to('../bla.txt')
    if test_symlink( f3 ):
        ensure_is_file(f3)
        assert f3.name == 'muahaha.bar'
        assert f3.resolve().name == 'bla.txt'
        assert lib.nctest_basename(f3) == 'muahaha.bar'
        rp=lib.nctest_real_path(f3)
        ensure_is_file(rp)
        assert lib.nctest_basename(rp) == 'bla.txt'
        assert lib.nctest_is_same_file(rp,f3)

    #a relative symlink to a dir:
    fr1 = (fd / 'a real dir1')
    fr2 = (fd / 'a real dir2')
    fr1_file = (fr1/'foo.txt')
    fr1.mkdir()
    fr2.mkdir()
    fr1_file.write_text('yo')

    fl1 = fr2 / 'a relsymlink to dir'
    fl1.symlink_to('../a real dir1')
    if test_symlink( fl1 ):
        ensure_is_dir(fl1)
        assert fl1.name == 'a relsymlink to dir'
        assert fl1.resolve().name == 'a real dir1'
        assert lib.nctest_basename(fl1) == 'a relsymlink to dir'
        rp=lib.nctest_real_path(fl1)
        ensure_is_dir(rp)
        assert lib.nctest_basename(rp) == 'a real dir1'
        assert pPath(rp).samefile(fl1)#pathlib's samefile works for dirs
        assert not lib.nctest_is_same_file(rp,fl1)#our ALWAYS return false for
                                                  #dirs
    #an absolute symlink to a dir:
    fl2 = fr2 / 'an abssymlink to dir'
    fl2.symlink_to(fr1.absolute())
    if test_symlink( fl2 ):
        ensure_is_dir(fl2)
        assert fl2.name == 'an abssymlink to dir'
        assert fl2.resolve().name == 'a real dir1'
        assert lib.nctest_basename(fl2) == 'an abssymlink to dir'
        rp=lib.nctest_real_path(fl2)
        ensure_is_dir(rp)
        assert lib.nctest_basename(rp) == 'a real dir1'
        assert pPath(rp).samefile(fl2)#pathlib's samefile works for dirs
        assert not lib.nctest_is_same_file(rp,fl2)#our ALWAYS return false for
                                                  #dirs

    #an absolute symlink to a file:
    fl3 = fr2 / 'an abssymlink to file'
    fl3.symlink_to(fr1_file.absolute())
    if test_symlink( fl3 ):
        ensure_is_file(fl3)
        assert fl3.name == 'an abssymlink to file'
        assert fl3.resolve().name == 'foo.txt'
        assert lib.nctest_basename(fl3) == 'an abssymlink to file'
        rp=lib.nctest_real_path(fl3)
        ensure_is_file(rp)
        assert lib.nctest_basename(rp) == 'foo.txt'
        assert lib.nctest_is_same_file(rp,fl3)

    #a relative symlink to a symlink to a file:
    fss1 = fd / 'ss1'
    fss1.symlink_to('../bla.txt')
    if test_symlink( fss1 ):
        ensure_is_file(fss1)
        assert fss1.name == 'ss1'
        assert fss1.resolve().name == 'bla.txt'
        assert lib.nctest_basename(fss1) == 'ss1'
        rp=lib.nctest_real_path(fss1)
        ensure_is_file(rp)
        assert lib.nctest_basename(rp) == 'bla.txt'
        assert lib.nctest_is_same_file(rp,fss1)

        fss2 = fd / 'ss2'
        fss2.symlink_to('ss1')
        if test_symlink( fss2 ):
            ensure_is_file(fss2)
            assert fss2.name == 'ss2'
            assert fss2.resolve().name == 'bla.txt'
            assert lib.nctest_basename(fss2) == 'ss2'
            rp=lib.nctest_real_path(fss2)
            ensure_is_file(rp)
            assert lib.nctest_basename(rp) == 'bla.txt'
            assert lib.nctest_is_same_file(rp,fss2)



def test6():
    pj = lib.nctest_path_join('/some/where','bla.txt')
    assert nonnative_sep not in pj
    expected = '/some/where/bla.txt'.replace('/',native_sep)
    assert pj == expected


def main():
    tempdir1 = pPath().joinpath('td1').absolute()
    tempdir2 = pPath().joinpath('td2').absolute()
    test1()
    test2()
    test3()
    test4()
    tempdir1.mkdir()
    os.chdir(tempdir1)
    test5(tempdir1)
    tempdir2.mkdir()
    os.chdir(tempdir2)
    test_dirsymlinks(tempdir2)
    test6()

if __name__=='__main__':
    import sys
    if '--no-tmpdir' in sys.argv[1:]:
        main()
    else:
        with work_in_tmpdir():
            main()

#fixme check with hardlinks + files without read permission
