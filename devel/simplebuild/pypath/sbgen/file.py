
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

from . import dirs
import os
import stat

def chmod_x( path ):
    st = os.stat(path)
    if not st.st_mode & stat.S_IEXEC:
        os.chmod(path, st.st_mode | stat.S_IEXEC)

class LinkFile:
    #A file which is actually a symlink to another file
    def __init__( self, real_path, name = None ):
        self.__name = name or real_path.name
        self.__path = real_path

    @property
    def name(self):
        return self.__name

    def create_in_dir( self, dir_path ):
        if not self.__path.is_file():
            raise RuntimeError(f'File not found: {self.__path}')
        f = dir_path.joinpath(self.__name)
        if f.is_file():
            #Check if already valid:
            if f.is_symlink() and f.samefile( self.__path ):
                return#already fine
            print( "Removing outdated link:",f)
            f.unlink()
        elif f.is_symlink():
            #must be broken link
            print( "Removing broken link:",f)
            f.unlink()
        dir_path.mkdir(parents=True, exist_ok=True)
        f.symlink_to( self.__path )
        print( "Created link:",f)

class ContentFile:

    def __init__( self, name, content ):
        self.__name = name
        self.__content = content

    @property
    def name(self):
        return self.__name

    def create_in_dir( self, dir_path ):
        f = dir_path.joinpath(self.__name)
        if f.is_file():
            #Check if already valid:
            if f.is_symlink():
                print( "Removing outdated link:",f)
                f.unlink()
            elif f.read_text() != self.__content:
                print( "Updating contents of:",f)
                f.write_text(self.__content)
                return
            else:
                return#already fine
        elif f.is_symlink():
            #Must be broken link
            print( "Removing broken link:",f)
            f.unlink()
        dir_path.mkdir(parents=True, exist_ok=True)
        f.write_text(self.__content)
        print( "Created file:",f)

class File:

    def __init__( self,
                  path_in_genroot,
                  content = None,
                  link_target = None,
                  make_executable = False ):
        f = dirs.genroot.joinpath(path_in_genroot)
        assert '#' not in f.name
        assert '!' not in f.name
        self._dir = f.parent
        self._exe = make_executable
        if content is not None:
            assert link_target is None
            self._file = ContentFile( name = f.name,
                                      content = content )
        else:
            assert link_target is not None
            self._file = LinkFile( real_path = link_target,
                                   name = f.name )

    def full_path( self ):
        return self._dir.joinpath( self._file.name )

    def create(self):
        self._file.create_in_dir( self._dir )
        if self._exe:
            chmod_x( self.full_path() )


all_files = []
all_files_outpaths = set()
def add_file( *a, **kw ):
    f = File( *a, **kw )
    all_files.append( f )
    fp = f.full_path()
    if fp in all_files_outpaths:
        raise RuntimeError(f'Created {fp} more than once')
    all_files_outpaths.add(fp)
    return f

def path_is_relative_to( p, pother ):
    #Path.is_relative_to(..) was introduced in Python 3.9, this function lets us
    #support python 3.8.
    import pathlib
    assert isinstance( p, pathlib.Path )
    if hasattr(p,'is_relative_to'):
        #Python 3.9+:
        return p.is_relative_to(pother)
    else:
        #Python 3.8:
        try:
            p.relative_to(pother)
            return True
        except ValueError:
            return False

def create_files():
    #First remove anything that disappared:
    import shutil
    flist = set()
    for p in all_files_outpaths:
        is_dir = False
        while True:
            flist.add( ( str(p), is_dir ) )
            is_dir = True
            p = p.parent
            if not path_is_relative_to(p,dirs.genroot):
                break
    for f in dirs.genroot.rglob('**/*'):
        key=(str(f),f.is_dir())
        if key not in flist:
            if f.is_symlink():
                print( "Removing link:",f)
                f.unlink()
            elif f.is_file():
                print( "Removing file:",f)
                f.unlink()
            else:
                print( "Removing directory:",f)
                shutil.rmtree(f)
    #Now create all the necessary files, dirs, and symlinks:
    for f in all_files:
        f.create()

