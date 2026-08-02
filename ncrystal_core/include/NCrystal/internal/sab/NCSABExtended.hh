#ifndef NCrystal_SABExtended_hh
#define NCrystal_SABExtended_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/sab/NCSABExtender.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    class SABExtended final : private MoveOnly {

    public:

      //Class combining SABProcessor and SABExtender instances, in order to
      //extend cross-sections and sampling capabilities of an S(alpha,beta)
      //table to large energies.

      //Constructor initialises based on an extender and processor. IMPORTANT
      //NOTE: It is very important that the SABExtender was created with
      //SigmaBound=1 (fixme: redesign so this is always guaranteed?).
      SABExtended( shared_obj<const SABProcessor>,
                   shared_obj<const SAB::SABExtender> );

      //Convenience function for initialising directly from SABData and using a
      //free-gas extension:
      static shared_obj<const SABExtended>
      createWithFGExtender( const SABCfg::Cfg&,
                            shared_obj<const SABData>,
                            std::shared_ptr<const VectD>
                            energyGrid = nullptr );

      //Access cross sections, always assuming SigmaBound=1barn.
      CrossSect crossSectionUnitSigmaBound( NeutronEnergy ekin ) const;

      //Sample scatterings:
      ScatterOutcomeIsotropic sampleScatter( RNG&, NeutronEnergy ) const;
      PairDD sampleScatterAlphaBeta( RNG& rng, NeutronEnergy ekin ) const;

    private:
      shared_obj<const SABProcessor> m_p;
      shared_obj<const SAB::SABExtender> m_e;
      NeutronEnergy m_emax;
      double m_kT;
      double m_invkT;
      double m_extSAtEmax;
      double m_tableSAtEmax;
      double m_kExtension;
      CrossSect xsHighE( NeutronEnergy ) const;
      PairDD scatABHighE( RNG&, NeutronEnergy ) const;
    };
  }
}

////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::CrossSect NCrystal::SABUtils::
SABExtended::crossSectionUnitSigmaBound( NeutronEnergy ekin ) const
{
  if ( ekin <= m_emax )
    return m_p->crossSectionUnitSigmaBound(ekin);
  return xsHighE(ekin);
}

inline NCrystal::PairDD NCrystal::SABUtils::
SABExtended::sampleScatterAlphaBeta( RNG& rng, NeutronEnergy ekin ) const
{
  if ( ekin <= m_emax ) {
    auto ab = m_p->sampleScatterAlphaBeta( rng, ekin );
    return { ab.alpha, ab.beta };
  }
  return scatABHighE( rng, ekin );
}

#endif
