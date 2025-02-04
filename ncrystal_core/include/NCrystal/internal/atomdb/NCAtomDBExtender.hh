#ifndef NCrystal_AtomDBExtender_hh
#define NCrystal_AtomDBExtender_hh

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

#include "NCrystal/interfaces/NCAtomData.hh"

namespace NCRYSTAL_NAMESPACE {

  //Helper class which can be used to implement an extension to the hard-wired
  //atom database (from NCAtomDB.hh), with the possibility to add new data,
  //override existing data, and add custom mixtures of isotopes and
  //elements. The interface is string-based, as the primary purpose is to
  //support @ATOMDB sections in NCMAT files, or the atomdb configuration
  //keyword.
  //
  //In order to ensure the same AtomData instances are created when the data is
  //identical (e.g. when two different NCMAT files define the same entry), a
  //global database is kept of any AtomData instance provided by any extended
  //database. This is important since downstream caches can then rely on the
  //unique id of the AtomData instances for caching purposes. The static
  //clearGlobalCache method can be used to clear it specifically (it will also
  //be cleared by a call to the global clearCaches() function from NCMem.hh);


  class AtomDBExtender final : private MoveOnly {
  public:
    //Constructor. If allowInbuiltDB=false, queries to the DB will not fall back
    //to the in-built database of natural elements and isotopes (from
    //NCAtomDB.hh):
    AtomDBExtender( bool allowInbuiltDB = true )
      : m_allowInbuiltDB(allowInbuiltDB) {}
    ~AtomDBExtender() = default;

    //Add data. A format_version must be provided. This should be the NCMAT
    //format version when data comes from an NCMAT file, it should be the most
    //recent NCMAT version (known to this installation of NCrystal) when data
    //comes from a configuration string. Throws BadInput exceptions in case of
    //problems. Each call to addData takes a single "line", consisting of
    //multiple words. In general only simple ASCII characters are allowed, and
    //whitespace is normalised and trimmed.
    constexpr static int latest_version = static_cast<int>(supported_ncmat_format_version_max);
    void addData( const VectS& words, unsigned format_version = latest_version );
    void addData( const std::string& line, unsigned format_version = latest_version );

    //Query DB (returns nullptr if not available):
    OptionalAtomDataSP lookupAtomDataAllowMissing(const std::string&);

    //Query DB (throws error if not available):
    AtomDataSP lookupAtomData(const std::string&);

    //Clear global cache as explained above:
    static void clearGlobalCache();

    //Plumbing:
    AtomDBExtender( const AtomDBExtender& ) = delete;
    AtomDBExtender& operator=( const AtomDBExtender& ) = delete;
    AtomDBExtender( AtomDBExtender&& ) = default;
    AtomDBExtender& operator=( AtomDBExtender&& ) = default;

  private:
    bool m_allowInbuiltDB;
    void populateDB(const std::string&,AtomDataSP);
    std::map<std::string,AtomDataSP> m_db;
  };

}

#endif
