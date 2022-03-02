#ifndef NCrystal_LazLoader_hh
#define NCrystal_LazLoader_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2022 NCrystal developers                                   //
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

#include "NCrystal/NCTextData.hh"
#include "NCrystal/NCInfoBuilder.hh"
#include "NCNXSLib.hh"

namespace NCrystal {

  class LazLoader{
  public:
    LazLoader(const TextData& laz_file, double dcutlow, double dcutup, Temperature);
    ~LazLoader() = default;
    InfoBuilder::SinglePhaseBuilder read();
  protected:

    void preParse(const TextData&);
    DataSourceName m_inputDescription;
    std::vector<VectS> m_raw_header;
    std::vector<VectS> m_raw_data;

    double m_dcutlow;
    double m_dcutup;
    Temperature m_temp;
    typedef VectS::const_iterator StrVecItr;
    typedef std::vector<VectS >::const_iterator RawItr;
    bool search_parameter(const std::string& attr, double &result) const;
    bool search_index(const std::string& attr, unsigned &result) const;
    bool search_spacegroup(unsigned &result) const;
    bool search_multiplicity(unsigned &result ) const;
    unsigned countAtom(const std::string& formula) const;
    bool setupSgInfo(unsigned spaceGroup, nxs::T_SgInfo& sgInfo) const;
  };
}

#endif
