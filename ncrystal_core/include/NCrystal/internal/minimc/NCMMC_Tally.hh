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

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //All tallies must implement the following interface, split into two classes
    //to allow templating on the basket type, but keeping most of the interface
    //non-templated.

    class TallyBase : NoCopyMove {
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
    };
    using TallyPtr = shared_obj<TallyBase>;

    template<class TBasket>
    class Tally : public TallyBase {
    public:
      //Results from simulations are delivered to the tally in the form of
      //baskets of neutrons. These are usually neutrons exiting (or having
      //missed) the geometry, but for e.g. 3D visualisation purposes we could
      //envision other schemes. Naturally, the simulation engine and tally
      //implementation must be in agreement on this:
      virtual void registerResults( const TBasket& ) = 0;
    };

    //The tally manager is a helper object responsible for dishing out tally
    //objects for all threads, and for later merging them upon work completion:

    class TallyMgr final : NoCopyMove {
    public:

      TallyMgr( TallyPtr tally_template )
        : m_template(std::move(tally_template)) {}

      TallyPtr getIndependentTallyPtr() const
      {
        return m_template->cloneSetup();
      }

      void addResult( TallyPtr res_SO )
      {
        NCRYSTAL_DEBUGMMCMSG("TallyMgr::addResult");
        std::shared_ptr<TallyBase> res = std::move(res_SO);
        std::shared_ptr<TallyBase> to_merge;
        {
          NCRYSTAL_LOCK_GUARD(m_final_mutex);//hold it only briefly!
          if ( m_final == nullptr ) {
            //first time, or someone else is currently merging.
            m_final = std::move(res);
            return;
          } else {
            to_merge.swap(m_final);
          }
        }
        if ( to_merge != nullptr ) {
          //perform the merging in the current thread without holding a lock,
          //and then put the merged result back:
          to_merge->merge( std::move( *res.get() ) );
          //Put result back:
          res = nullptr;//to be safe (not strictly needed)
          this->addResult( std::move(to_merge) );
        }
      }
      TallyPtr getFinalResult()
      {
        NCRYSTAL_LOCK_GUARD(m_final_mutex);//should not really be needed if used correctly
        nc_assert_always(m_final!=nullptr);
        return std::move(m_final);
      }
    private:
      TallyPtr m_template;
      std::shared_ptr<TallyBase> m_final;
      std::mutex m_final_mutex;
    };

  }
}

#endif
