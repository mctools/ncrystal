#ifndef NCrystal_TestFindData_hh
#define NCrystal_TestFindData_hh

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

#include "NCrystal/internal/utils/NCFileUtils.hh"
#include <cstdlib>

// Utility functions for locating the test data in <reporoot>/tests/data

namespace nctest {

  inline std::string find_test_data( const char * subdir, const char * filename )
  {
    std::string datadir;
    {
#ifdef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
      const char * sbld_datadir = std::getenv("SBLD_DATA_DIR" );
      nc_assert_always(sbld_datadir);
      std::string pkgname;
      if ( subdir ) {
        pkgname = "NCTestData_";
        pkgname += subdir;
      } else {
        pkgname = "NCTestUtils";
      }
      datadir = NCrystal::path_join(sbld_datadir,pkgname);
#else
      const char * datadir_base = std::getenv("NCRYSTALTEST_DATADIR");
      nc_assert_always(datadir_base);
      datadir = ( subdir ?
                  NCrystal::path_join(datadir_base, subdir)
                  : std::string(datadir_base) );
#endif
    }
    std::string f = NCrystal::path_join(datadir, filename);
    nc_assert_always(NCrystal::file_exists(f));
    return f;
  }

  inline std::string find_test_data( const char * filename ){
    return find_test_data( nullptr, filename );
  }
}

#endif
