#ifndef NCrystal_MMC_Tally_hh
#define NCrystal_MMC_Tally_hh

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

#include "NCrystal/internal/minimc/NCMMC_Defs.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/minimc/NCMMC_UBView.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //All tallies must implement the following interface.

    class TallyBase : NoCopyMove {//fixme rename as Tally.
    public:
      virtual ~TallyBase() = default;

      //To support multi-threading, tallies can be cloned from a template object
      //(without any data from simulations), and later merged:
      virtual shared_obj<TallyBase> cloneSetup() const = 0;
      virtual void merge(TallyBase&&) = 0;

      //Ultimately, a tally provides is data as JSON (this is usually done
      //post-merge, so only happens once per simulation):
      virtual void tallyItemToJSON( std::ostream&, StrView itemName ) const = 0;
      virtual VectS tallyItemNames() const = 0;

      //Results from simulations are delivered to the tally in the form of
      //baskets of neutrons. These are usually neutrons exiting (or having
      //missed) the geometry.
      virtual void registerResultsUB( const UniversalBasket& ) = 0;//fixme: rename without UB
    };
    using TallyPtr = shared_obj<TallyBase>;

    class TallyMgr final : NoCopyMove {
    public:
      //Manager responsible for dishing out tally objects for all threads, and
      //for later merging them upon work completion.
      TallyMgr( TallyPtr tally_template );
      TallyPtr getIndependentTallyPtr() const;
      void addResult( TallyPtr );
      TallyPtr getFinalResult();
    private:
      TallyPtr m_template;
      std::shared_ptr<TallyBase> m_final;
      std::mutex m_final_mutex;
    };

  }
}

#endif
