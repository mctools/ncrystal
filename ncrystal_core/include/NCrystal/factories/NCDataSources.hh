#ifndef NCrystal_DataSources_hh
#define NCrystal_DataSources_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/factories/NCFactImpl.hh"

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Utilities for customising NCrystal's data sources, for instance by adding  //
// directories with data files to the search path, or adding in-memory files. //
//                                                                            //
// The functions here are merely convenience functions dealing with the most  //
// typical use cases. Behind the scenes, the functions simply modify the list //
// of registered TextData factories (cf. NCFactImpl.hh). For unusual          //
// use-cases not served by the functions below, more flexiblity might thus    //
// be achieved by implementing and registering custom TextData factories.     //
//                                                                            //
// Most users are likely only interested in the first few functions near the  //
// top of the file, allowing the registry of new file content directly in     //
// memory.                                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace DataSources {

    //For reference the default priorities used by the various NCrystal TextData
    //factories are:

    constexpr unsigned default_priority_abspath = 150;
    constexpr unsigned default_priority_relpath = 140;
    constexpr unsigned default_priority_virtual = 130;
    constexpr unsigned default_priority_stdlib  = 120;
    constexpr unsigned default_priority_stdpath = 110;

    ////////////////////////////////////////////////////////////////////////////
    //Register custom virtual file data, consisting of in-memory files or a
    //virtual alias for an on-disk file. This needs a "filename" and the
    //associated content or location. After registering such virtual "files",
    //they can be used as file names in cfg strings or MatCfg objects (for input
    //types which support it, which certainly includes NCMAT file data, for
    //which even the virtual "filename" should end with ".ncmat"). Registering
    //the same virtual filename more than once, will simply override the
    //previous registration.
    //
    //Technically, the registered file data will be delivered to users from a
    //factory named "virtual", so prefixing the filename with "virtual::" will
    //guarantee that only the just registered file is returned.

    //Version which takes ownership of provided data:
    NCRYSTAL_API void registerInMemoryFileData( std::string virtualFileName,
                                                std::string&& data,
                                                Priority = Priority{default_priority_virtual} );

    //Version which just registers the address of data and does not copy
    //it. Naturally the lifetime of the static_data should be longer than any
    //direct or indirect calls to FactImpl::createTextData (it is intended for
    //efficiently hard-coding a large database directly in C/C++ code):
    NCRYSTAL_API void registerInMemoryStaticFileData( std::string virtualFileName,
                                                      const char* static_data,
                                                      Priority = Priority{default_priority_virtual} );

    //A somewhat more unusual use-case perhaps, setting up a virtual alias for
    //on-disk file (will also be served the "virtual" factory):
    NCRYSTAL_API void registerVirtualFileAlias( std::string virtualFileName,
                                                std::string realAbsPathFileName,
                                                Priority = Priority{default_priority_virtual} );

    //Multiple virtual files can be registered in dedicated factory (e.g. to
    //register all files from a given analysis or NCrystal plugin in a nice
    //group):
    NCRYSTAL_API void registerNamedVirtualDataSource( const std::string& factoryName,
                                                      std::map<std::string,TextDataSource>&& virtualfiles,
                                                      Priority );

    ////////////////////////////////////////////////////////////////////////////
    // Register custom directories to be monitored for data files.
    //
    // Directories do not need to exist at the time of this call, and
    // registering the same directory more than once, will simply override the
    // previous registration.

    NCRYSTAL_API void addCustomSearchDirectory( std::string dirname,
                                                Priority = Priority{default_priority_stdpath+1} );

    //Remove all search directories added with addCustomSearchDirectory:
    NCRYSTAL_API void removeCustomSearchDirectories();

    ////////////////////////////////////////////////////////////////////////////
    // Manipulate standard factory behaviour. Most applications will not do
    // this, but it might be useful for specialised applications.
    //

    //Whether or not absolute file paths are allowed (affects the "abspath"
    //factory):
    NCRYSTAL_API void enableAbsolutePaths(bool = true);

    //Whether or not paths relative to current working directory are allowed
    //(affects the "relpath" factory):
    NCRYSTAL_API void enableRelativePaths(bool = true);

    //Whether or not the standard data library shipped with NCrystal should be
    //searched (affects the "stdlib" factory).
    //
    //Unless NCrystal is configured to have the standard data library embedded
    //into the binary at compilation time, the location (directory path) of the
    //standard data library is taken from the NCRYSTAL_DATADIR environment
    //variable. If the environment variable is not set, the location is taken
    //from the compile time definition of the same name. If neither is set, and
    //data was not embedded at compilation time, the standard data library will
    //be disabled by default and the location must be provided before it can be
    //enabled. In all cases, the location can be overridden if explicitly
    //provided by the user.
    NCRYSTAL_API void enableStandardDataLibrary( bool = true,
                                                 Optional<std::string> path = NullOpt );

    //Whether or not the standard search path should be searched (affects the
    //"stdpath" factory). This standard search path is built by concatenating
    //entries in the NCRYSTAL_DATA_PATH environment variables with entries in
    //the compile time definition of the same name (in that order):
    NCRYSTAL_API void enableStandardSearchPath(bool = true);

    //Whether or not to search data directories provided by plugins. These
    //directories will in any case only be used when explicitly requested
    //(e.g. "plugins::SomePlugin/somefile.ncmat").
    NCRYSTAL_API void enablePluginSearchPaths(bool = true);

    // Disable all standard data sources, remove all TextData factories as well,
    // clear all registered virtual files and custom search directories. Finish
    // by calling global clearCaches function ("Ripley: I say we take off and
    // nuke the entire site from orbit. It's the only way to be sure."):
    NCRYSTAL_API void removeAllDataSources();

    ////////////////////////////////////////////////////////////////////////////
    //For "browsing" available files. This information is ultimately provided by
    //factories, and might not be exhaustive. For instance, absolute paths are
    //never reported (as it would mean NCrystal would have to search the entire
    //filesystem!). For monitored directories, only files with a limited number
    //of recognised file extensions will be reported.

    struct NCRYSTAL_API FileListEntry {
      std::string name;//the filename (virtual or real)
      std::string source;//directory path or description
      std::string factName;//the factory name
      Priority priority;//the priority value it would be delivered with (used to
                        //disentagle conflicting file names in different
                        //factories).
    };
    NCRYSTAL_API std::vector<FileListEntry> listAvailableFiles();

    //Augment the list of recognised file extensions (it is recommended that
    //plugins enabling support for a given file type also use this method to
    //register their associated extension here). This list only affects the
    //ability to browse available files, not which files can actually be
    //provided by e.g. createTextData:
    NCRYSTAL_API void addRecognisedFileExtensions(std::string);
    NCRYSTAL_API VectS recognisedFileExtensions();
  }

}

#endif
