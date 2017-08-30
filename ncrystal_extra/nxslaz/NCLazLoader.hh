#ifndef NCrystal_LazLoader_hh
#define NCrystal_LazLoader_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include <string>
#include <vector>

#include "NCNXSLib.hh"

namespace NCrystal {

  class Info;

  class LazLoader{
  public:
    LazLoader(std::string laz_file, double dcutlow, double dcutup, double temp);

    ~LazLoader();
    void read();
    Info* getCrystalInfo();
  protected:
    std::string m_full_path;

    Info *m_cinfo;
    double m_dcutlow;
    double m_dcutup;
    double m_temp;
    typedef std::vector<std::string>::const_iterator StrVecItr;
    typedef std::vector<std::vector<std::string> >::const_iterator RawItr;
    std::vector<std::vector<std::string> > m_raw_header;
    std::vector<std::vector<std::string> > m_raw_data;
    bool search_parameter(std::string attr, double &result);
    bool search_index(std::string attr, unsigned &result);
    bool search_spacegroup(unsigned &result);
    bool search_multiplicity(unsigned &result );
    unsigned countAtom(std::string formula);
    bool setupSgInfo(unsigned spaceGroup, nxs::T_SgInfo& sgInfo);
  };
}

#endif
