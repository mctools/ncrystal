#ifndef NCrystal_AbsOOV_hh
#define NCrystal_AbsOOV_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/NCProcImpl.hh"
#include "NCrystal/NCMatInfo.hh"

namespace NCrystal {

  class AbsOOV : public ProcImpl::AbsorptionIsotropicMat {
  public:

    // Provide absorption cross section based on simple 1/velocity
    // (=OneOverVelocity=OOV) scaling. This is non-oriented.


    AbsOOV( SigmaAbsorption );
    AbsOOV( const MatInfo& );

    const char * name() const noexcept final { return "AbsOOV"; }
    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const final;

    std::shared_ptr<Process> createMerged( const Process& ) const override;

  private:
    double m_c;
  };
}

#endif
