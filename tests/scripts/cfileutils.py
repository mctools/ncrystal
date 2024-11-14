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
                                 #ensure_error,
                                 explicit_unicode_str,
                                 is_windows )
import pathlib
pPath = pathlib.Path
import os

lib = Lib('testcfileutils')
lib.dump()
assert hasattr(lib,'nctest_file_exists_and_readable')
assert hasattr(lib,'nctest_get_current_working_dir')
assert hasattr(lib,'nctest_real_path')
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
assert hasattr(lib,'nctest_is_dir')

def test1():
    d = pPath('.')
    (d / 'foo.txt').write_text('foo')

    assert lib.nctest_file_exists_and_readable('./foo.txt')
    assert lib.nctest_file_exists_and_readable(r'.\foo.txt')
    assert lib.nctest_file_exists_and_readable('foo.txt')
    assert not lib.nctest_file_exists_and_readable('bar.txt')
    assert not lib.nctest_file_exists_and_readable('./bar.txt')
    assert not lib.nctest_file_exists_and_readable('/some/file/that/does/not/exist')
    print("Testing get_current_working_dir...")
    _cwd = lib.nctest_get_current_working_dir()
    d_abs = d.resolve().absolute()
    assert d.samefile(d_abs)
    assert d_abs.samefile ( pPath(_cwd).resolve().absolute() )
    assert d_abs.samefile ( _cwd )
    assert d.samefile ( _cwd )
    print("...get_current_working_dir ok")

def test2():
    #FIXME: Replicate the stuff below!

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
                  #FIXME: What about '/some/where/../bla.txt'? Should be
                  #'/some/bla.txt'. But '../bla.txt' should keep the relative
                  #path.
                  #'/test/bla/'#pPath(..).name and unix basename gives
                  #            #'bla' but os.path.basename is ''
                 ]

    fctnames = ['dirname','basename','basename_view',
                'fileextension','fileextension_view',
                'file_exists_and_readable',
                'path_is_absolute','path_is_relative',
                'pathseps_platform','pathseps_generic',
                'drive_letter'
                ]


    fctnames_same_is_absolute = ['dirname','pathseps_platform',
                                 'pathseps_generic']

    native_sep = '\\' if is_windows() else '/'

    for i,p in enumerate(testpaths):
        print(f'{i} Testing path "{p}":')
        for fct in fctnames:
            res = getattr(lib,f'nctest_{fct}')(p)
            if fct == 'pathseps_platform':
                #ensure test reproducibility:
                res = res.replace(native_sep,'@NATIVESEP@')
            print(f'   mctools_{fct} = "{res}"')

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
            #return dict(pathlib=os.path.basename(pobj))
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

    #def testglob(pattern,nexpect,is_absolute=False):
    #    def fmtpat( pattern ):
    #        if is_absolute:
    #            assert pPath(pattern).is_absolute()
    #            s = '/SOMEWHERE/'+pPath(pattern).name
    #        else:
    #            s = str(pattern)
    #        if is_windows():
    #            s=s.replace('\\','/')#for test reproducibility
    #        return repr(explicit_unicode_str(s))
    #    print(f'Testing ncglob({fmtpat(pattern)}):')
    #    import glob
    #    g = sorted(glob.glob(pattern))
    #    print('  --> Python glob got %i results:'%len(g))
    #    for e in g:
    #        print(f'   *: {fmtpat(e)}')
    #    g = lib.nctest_ncglob(pattern)
    #    print('  --> got %i results:'%len(g))
    #    all_ok = True
    #    for e in g:
    #        badstr = ''
    #        if not pPath(e).exists():
    #            badstr = ' (ERROR NOT ACTUALLY FOUND)'
    #            all_ok = False
    #        print(f'   *: {fmtpat(e)}{badstr}')
    #    if len(g) != nexpect:
    #        raise SystemExit('Error: did not get expected %i entries'%nexpect)
    #    if not all_ok:
    #        raise SystemExit('Error: some hits were invalid')
    #
    ##windows: testglob(r'some_sub_dir\*.txt',2)
    #testglob('some_sub_dir/*.txt',2)
    #testglob('some_sub_dir/*',3)
    #testglob('some_sub_dir/***',3)
    #testglob('some_sub_dir/*txt*',2)
    #testglob('some_sub_dir/?.txt*',2)
    #testglob('some_sub_dir/*.ncmat',1)
    #testglob('some_sub_dir/*b*',2)
    #testglob('some_sub_dir/*a*',2)
    #testglob('some_sub_dir/*.ncm',0)#could catch the .ncmat file if using wrong
    #                                #windows API, where only 3 characters in
    #                                #extensions are compared.
    #testglob('unicodedir_test\u4500abc/*.ncmat',1)
    #testglob('unicodedir_test\u4500abc/*.txt',2)
    #testglob('unicodedir_test\u4500abc/*\u2030*',1)
    #
    #testglob('unicodedir_test\u4500ab*',1)
    #testglob('unicodedir_test\u4500abc*',1)
    #testglob('unicodedir_*',1)
    #
    #import os
    #testglob(os.path.join(
    #    str(pPath('unicodedir_test\u4500abc').absolute()),
    #    '*') ,3,is_absolute = True)
    #
    ##Fixme: we don't propagate C++ exceptions nicely:
    ##with ensure_error(RuntimeError,
    ##                  'ncglob only supports wildcards in the'
    ##                  ' last file or directory name'):
    ##    testglob('*ome_sub*/*.txt',2)
    ##so we do instead:
    #print("Triggering expected error (hopefully):",flush=True)
    #lib.nctest_ncglob('some*sub_dir/a.txt')
    #
    ##Finally, test with paths longer than 260 chars:
    #dummy_path = pPath('a'*150).joinpath('b'*150).joinpath('bla.txt')
    #dummy_path.parent.parent.mkdir()
    #dummy_path.parent.mkdir()
    #dummy_path.write_text('hello')
    #testglob( str(dummy_path.parent)+'/bl*.txt',1)

def test4():

    def testdirname(path,expected):
        res=lib.nctest_dirname(path)
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

def test5():
    #exercise a bunch of things, including is_same_file, is_dir, real_path,
    #path_join.
    td = pPath().joinpath('rp').absolute()
    td.mkdir()
    os.chdir(td)
    files = [
        dict(name='asimplefile.txt',content='hello'),
        dict(name='sd 1',is_dir=True),
        dict(name='sd 1/lala.txt',content='bla'),
        dict(name='sd 1/\u4500abc/.\u4500',content='foobar'),
        dict(name='sd 1/\u4500abc/.foo~',content='foobar'),
        dict(name='sd2',is_dir=True),
        dict(name='sd2/sl1',symlink='../sd 1'),
        dict(name='sd2/sl2',symlink=str(td.joinpath('sd 1').absolute())),
        dict(name='sd    \t 3/bla.ncmat',content='NCMAT'),
    ]
    #Create layout:
    for finfo in files:
        f = finfo['name']
        fabs = td.joinpath(f)
        symlink=finfo.get('symlink')
        if finfo.get('is_dir'):
            fabs.mkdir()
            continue
        elif symlink:
            fabs.symlink_to(symlink)
        else:
            fabs.parent.mkdir(exist_ok=True)
            content = finfo.get('content')
            if content is not None:
                fabs.write_text(content)
            else:
                fabs.touch()
    #Analyse it:

    for ifile,f in enumerate(sorted(pPath('.').rglob('**/*'))):
        fabs = f.absolute()
        fabs_resolved = f.resolve()
        real_path0 = lib.nctest_real_path( fabs_resolved )
        transforms = [ str,
                       lambda fin : str(fin).replace('/','\\'),
                       lambda fin : str(fin).replace('\\','/') ]
        is_dir = fabs.is_dir()
        fabs_prev = td.parent
        for is_abs, fin in [ (True,fabs), (False,f) ]:
            print( f'Testing file {ifile}: "{f}"' )
            assert not lib.nctest_is_same_file( str(fabs_prev), str(fin) )
            for ftest in sorted(set([t(fin) for t in transforms])):
                assert bool(lib.nctest_is_dir(ftest)) == is_dir
                assert bool(lib.nctest_path_is_absolute(ftest)) == bool(is_abs)
                assert bool(lib.nctest_path_is_relative(ftest)) == bool(not is_abs)
                assert lib.nctest_file_exists_and_readable(ftest)
                assert lib.nctest_is_same_file(ftest,str(fabs_resolved))
                assert lib.nctest_is_same_file(str(fabs_resolved),ftest)
                assert real_path0 == lib.nctest_real_path( ftest )
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
                    assert pPath(fabs_resolved).samefile( real_path0 )
        fabs_prev = fabs

def main():
    test1()
    test2()
    test3()
    test4()
    test5()#nb changes dir

if __name__=='__main__':
    import sys
    if '--no-tmpdir' in sys.argv[1:]:
        main()
    else:
        with work_in_tmpdir():
            main()
