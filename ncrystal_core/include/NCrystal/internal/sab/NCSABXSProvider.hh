#ifndef NCrystal_SABXSProvider_hh
#define NCrystal_SABXSProvider_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/sab/NCSABExtender.hh"

namespace NCRYSTAL_NAMESPACE {

  class SABXSProvider final : private MoveOnly {
  public:
    SABXSProvider( VectD&& egrid, VectD&& xsvals,
                   std::shared_ptr<const SAB::SABExtender> );
    void setData( VectD&& egrid, VectD&& xsvals,
                  std::shared_ptr<const SAB::SABExtender> );
    SABXSProvider() = default;//default constructs invalid instance
    ~SABXSProvider();
    CrossSect crossSection(NeutronEnergy) const;

    //Move ok:
    SABXSProvider( SABXSProvider&& ) = default;
    SABXSProvider& operator=( SABXSProvider&& ) = default;

    //For reference:
    const VectD & internalEGrid() const { return m_egrid; }
    const VectD & internalXSGrid() const { return m_xs; }
  private:
    VectD m_egrid, m_xs;
    std::shared_ptr<const SAB::SABExtender> m_extender;
    double m_kExtension;
  };

}

#endif
