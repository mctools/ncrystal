
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

"""

Utilities for customising NCrystal's data sources, for instance by adding
directories with data files to the search path, or adding in-memory files.

"""

from ._chooks import _get_raw_cfcts,_str2cstr
import pathlib as _pathlib
_rawfct = _get_raw_cfcts()

def registerInMemoryFileData(virtual_filename,data):
    """Register in-memory file data. This needs a "filename" and the content of this
       virtual file. After registering such in-memory "files", they can be used
       as file names in cfg strings or MatCfg objects. Registering the same
       filename more than once, will simply override the content.

       As a special case data can specified as "ondisk://<path>",
       which will instead create a virtual alias for an on-disk file.

       Technically, the registered file data will be delivered to users from a
       factory named "virtual", so prefixing the filename with "virtual::" in
       future requests, will guarantee that only data registered in this manner
       is returned.
    """
    if ( isinstance(data,str) and data.startswith('ondisk://')):
        data = 'ondisk://'+str(_pathlib.Path(data[9:]).resolve())
    _rawfct['ncrystal_register_in_mem_file_data'](virtual_filename,data)


def addCustomSearchDirectory(dirpath):
    """Register custom directories to be monitored for data files."""
    _rawfct['ncrystal_add_custom_search_dir'](_str2cstr(str(_pathlib.Path(dirpath).resolve())))

def removeCustomSearchDirectories():
    """Remove all search directories added with addCustomSearchDirectory."""
    _rawfct['ncrystal_remove_custom_search_dirs']()

def removeAllDataSources():
    """Disable all standard data sources, remove all TextData factories as well,
       clear all registered virtual files and custom search directories. Finish
       by calling global clearCaches function ("Ripley: I say we take off and
       nuke the entire site from orbit. It's the only way to be sure.").
    """
    _rawfct['ncrystal_remove_all_data_sources']()

def enableAbsolutePaths( enable = True ):
    """Whether or not absolute file paths are allowed."""
    _rawfct['ncrystal_enable_abspaths'](1 if enable else 0)

def enableRelativePaths( enable = True ):
    """Whether or not paths relative to current working directory are allowed."""
    _rawfct['ncrystal_enable_relpaths'](1 if enable else 0)

def enableStandardSearchPath( enable = True ):
    """Whether or not the standard search path should be searched. This standard
      search path is is by default searched *after* the standard data library,
      and is built by concatenating entries in the NCRYSTAL_DATA_PATH
      environment variables with entries in the compile time definition of the
      same name (in that order). Note that by default the standard search path
      is searched *after* the standard data library.
    """
    _rawfct['ncrystal_enable_stdsearchpath'](1 if enable else 0)

def enableStandardDataLibrary( enable = True, dirpath_override = None ):
    """Whether or not the standard data library shipped with NCrystal should be
       searched.

       Unless NCrystal is configured to have the standard data library embedded
       into the binary at compilation time, the location (directory path) of the
       standard data library is taken from the NCRYSTAL_DATADIR environment
       variable. If the environment variable is not set, the location is taken
       from the compile time definition of the same name. If neither is set, and
       data was not embedded at compilation time, the standard data library will
       be disabled by default and the location must be provided before it can be
       enabled. In all cases, the location can be overridden if explicitly
       provided by the user as the second parameter to this function.
    """
    import ctypes
    d = _str2cstr(str(_pathlib.Path(dirpath_override).resolve())) if dirpath_override else ctypes.cast(None, ctypes.c_char_p)
    _rawfct['ncrystal_enable_stddatalib'](1 if enable else 0, d)

class FileListEntry:
    """Entry in list returned by browseFiles."""
    def __init__(self,*,name,source,factName,priority):
        self.__n = name or None
        self.__f = factName or None
        self.__p = int(priority) if priority.isdigit() else priority
        self.__s = source or None

    @property
    def name(self):
        """The (possibly virtual) filename needed to select this entry"""
        return self.__n

    @property
    def source(self):
        """Description (such as the parent directory in case of on-disk files)"""
        return self.__s

    @property
    def factName(self):
        """Name of the factory delivering entry."""
        return self.__f

    @property
    def priority(self):
        """The priority value of the entry (important in case multiple factories
        delivers content with the the same name). Can be 'Unable',
        'OnlyOnExplicitRequest' or an integer priority value (entries with
        higher values will be preferred).
        """
        return self.__p

    @property
    def fullKey(self):
        """The string '%s::%s'%(self.factName,self.name), which can be used to
           explicitly request this entry without interference from similarly
           named entries in other factories.
        """
        return '%s::%s'%(self.__f,self.__n)

    def __str__(self):
        ll=[]
        if self.__n:
            ll+=['name=%s'%self.__n]
        if self.__s:
            ll+=['source=%s'%self.__s]
        if self.__f:
            ll+=['factory=%s'%self.__f]
        ll+=['priority=%s'%self.__p]
        return 'FileListEntry(%s)'%(', '.join(ll))

    def __lt__(self, other):
        if not isinstance(other, FileListEntry):
            return False
        return ( (self.__f,self.__n,self.__s,self.__p)
                 < (other.__f,other.__n,other.__s,other.__p) )

def browseFiles(dump=False,factory=None):
    """Browse list of available input files (virtual or on-disk). The list is not
       guaranteed to be exhaustive, but will usually include all files in
       supported files in the most obvious locations (the NCrystal data
       directory and other directories of the standard search path, the current
       working directory, virtual files embedded in the NCrystal library or
       registered dynamically.

       Returns a list of FileListEntry objects. If the dump flag is set to True,
       the list will also be printed to stdout in a human readable form.

       Setting factory parameter will only return / print entries from the
       factory of that name.

    """
    res=[]
    def sortkey(e):
        praw = e.priority
        if praw=='Unable':
            p=-2
        elif isinstance(praw,int):
            p=praw
        else:
            assert praw=='OnlyOnExplicitRequest'
            p=-1
        return (-p, e.factName,e.source,e.name)
    for n,s,f,p in _rawfct['ncrystal_get_filelist']():
        res.append( FileListEntry(name=n,source=s,factName=f,priority=p) )
    res.sort(key=sortkey)
    if dump:
        seen_names=set()
        def groupfct( e ):
            return (e.factName,e.source,e.priority)
        lastgroup = None
        pending=[]
        def print_pending():
            if not pending:
                return
            if factory is not None and lastgroup[0]!=factory:
                pending.clear()
                return
            n=len(pending) - 1
            pending[0] = pending[0]%('%s files'%n if n!=1 else '%s file'%n )
            for line in pending:
                print (line)
            pending.clear()
        for e in res:
            group = groupfct(e)
            if lastgroup != group:
                print_pending()
                lastgroup = group
                pending.append('==> %%s from "%s" (%s, priority=%s):'%group)
            hidden = e.name in seen_names
            seen_names.add(e.name)
            extra=''
            prname=e.name
            if e.priority=='OnlyOnExplicitRequest':
                prname='%s::%s'%(e.factName,e.name)
            elif hidden:
                extra=' <--- Hidden by higher priority entries (select as "%s::%s")'%(e.factName,e.name)
            pending.append(    '    %s%s'%(prname,extra))
        print_pending()
        return #return None in this case, to avoid spurious printouts in an interactive session
    if factory is None:
        return res
    return [e for e in res if e.factName==factory]
