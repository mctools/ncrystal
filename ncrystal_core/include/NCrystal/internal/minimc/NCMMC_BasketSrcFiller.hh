#ifndef NCrystal_MMC_BasketSrcFiller_hh
#define NCrystal_MMC_BasketSrcFiller_hh

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
#include "NCrystal/internal/minimc/NCMMC_BasketMgr.hh"
#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/minimc/NCMMC_Source.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //A helper class, extending the BasketMgr with the ability to also "top up"
    //pending baskets with new baskets from a source, and to handle the
    //propagation of src particles to the modelled geometry. It can also be used
    //to signal the end of simulations (i.e. "please do not provide any more
    //fresh source particles"), even with a handy "button". It will be MT-safe
    //to use, unless told not to be.

    enum class ThreadedUsage { Single, Multi };

    namespace detail {
      void propagateDistance( NeutronBasket& nb, Span<const double> distances,
                              std::size_t offset = 0 ) ncnoexceptndebug;
    }

    template<class TBasket>
    class BasketSrcFiller final : NoCopyMove {
    public:
      using basket_t = TBasket;
      using basket_holder_t = BasketHolder<TBasket>;
      using basketmgr_t = BasketMgr<TBasket>;

      BasketSrcFiller( GeometryPtr geom,
                       SourcePtr src,
                       shared_obj<basketmgr_t> bm,
                       ThreadedUsage tu = ThreadedUsage::Multi )
        : m_geom(std::move(geom)),
          m_src(std::move(src)),
          m_basketmgr(std::move(bm)),
          m_srcParticlesMightBeOutside(m_src->particlesMightBeOutside(m_geom))
      {
        m_srcHalted.store(false);
        if ( tu == ThreadedUsage::Multi && !m_src->metaData().concurrent )
          m_srcmutex.emplace();
      }

      const Geometry& geometry() const { return m_geom; }
      const basketmgr_t& basketMgr() const { return m_basketmgr; }
      basketmgr_t& basketMgr() { return m_basketmgr; }

      void haltSource()
      {
        m_srcHalted.store( true );
      }

      basket_holder_t getPendingBasket( ThreadCount nthreads,
                                        RNG& rng,
                                        const std::function<void(const basket_t&)>& resultFct )
      {
        return this->getPendingBasketImpl(nthreads,rng,10,resultFct);
      }

    private:
      void propagateToVolume( basket_t& b, std::size_t offset,
                              const std::function<void(const basket_t&)>& resultFct ) {
        //We must edit the basket, propagating all the neutrons to the volume if
        //they can. Those that can NOT do that, should be marked as having
        //missed the volume and served up as a result.

        //First do the geometry distance calculations:
        nc_assert(m_srcParticlesMightBeOutside);
        double dist_results[NeutronBasket::N];
        m_geom->distToVolumeEntry( b.neutrons, dist_results );

        //Next, reserve a basket for results (those that missed):
        auto bh_results = m_basketmgr->allocateBasket();
        auto& bmiss = bh_results.basket();

        std::size_t i_first_hole = NeutronBasket::N;
        const std::size_t b_size = b.size();
        for ( std::size_t i = offset; i < b_size; ++i ) {
          double dist = dist_results[i];
          if ( dist < 0.0 ) {
            //Missed the volume, register as result.
            bmiss.appendEntryFromOther( b, i );
            if ( i_first_hole == NeutronBasket::N )
              i_first_hole = i;
          } else {
            //Keep this. In case of holes, we also have to move it.
            std::size_t inew = i;
            if ( i_first_hole < i ) {
              b.copyEntry( (inew=i_first_hole++), i );
              dist_results[inew] = dist;
            }
          }
        }

        //shrink-to-fit:
        if ( i_first_hole != NeutronBasket::N )
          b.neutrons.nused = i_first_hole;

        //Now we should propagate all the neutrons that were not already inside
        //the volume:
        detail::propagateDistance( b.neutrons, dist_results, offset );

        //And finally, return any neutrons that missed as a result, after
        //marking them as having missed the target.
        if ( bh_results.basket().empty() ) {
          m_basketmgr->deallocateBasket( std::move(bh_results) );
        } else {
          for ( auto i : ncrange( bmiss.size() ) )
            bmiss.cache.markAsMissedTarget(i);

          resultFct( bh_results.basket() );
          m_basketmgr->deallocateBasket( std::move(bh_results) );
        }
      }

      basket_holder_t getPendingBasketImpl( ThreadCount nthreads,
                                            RNG& rng,
                                            unsigned nretry,
                                            const std::function<void(const basket_t&)>& resultFct ) {

        //Get via basket mgr:
        nc_assert_always(nthreads.get()>=1);
        //Pass nthreads to basketmgr call in the next line, to avoid merging
        //baskets if there are anyway less pending baskets than the number of
        //threads:
        auto bh = m_basketmgr->getPendingBasketOrAllocateEmpty( nthreads );

        nc_assert(bh.valid());
        const std::size_t size_orig = bh.basket().size();
        if ( size_orig >= basket_N_almost_Full )
          return bh;

        //Ok, we need to try to fill the basket from the source (while locking
        //m_srcmutex if source does not support concurrent access).
        bool src_has_more = !m_srcHalted.load();
        if ( src_has_more ) {
          if ( !m_srcmutex.has_value() ) {
            m_src->fillBasket( rng, bh.basket().neutrons );
          } else {
            NCRYSTAL_LOCK_GUARD(m_srcmutex.value());
            m_src->fillBasket( rng, bh.basket().neutrons );
          }
          //Source provided neutron parameters, but we still need to initialise
          //any extra cache parameters as well.
          for( std::size_t i = size_orig; i < bh.basket().size(); ++i)
            bh.basket().cache.init( i );
          const bool src_ran_out = !bh.basket().full();
          if ( src_ran_out ) {
            m_srcHalted.store( true );
            src_has_more = false;
          }
        }

        //If source particles might be outside volume, we have to propagate them
        //to the volume (and record the rest as results already):
        if ( m_srcParticlesMightBeOutside ) {
          propagateToVolume( bh.basket(), size_orig, resultFct );
          if ( bh.basket().empty() ) {
            //Original basket was empty, and they *all* missed!
            //Try again, unless src ran out:
            m_basketmgr->deallocateBasket(std::move(bh));
            if ( !src_has_more )
              return basket_holder_t{ no_init };
            if ( nretry==0 )
              throw std::runtime_error("Source particles consistently seems to miss the geometry.");
            return this->getPendingBasketImpl( nthreads, rng, nretry-1, resultFct );
          }
        }

        //NB: After the propagation step, we still might have most particles
        //missing the volume! We could try to "top off" the basket again in that
        //case? (TODO)

        if ( bh.basket().empty() ) {
          //Apparently we are done.
          m_basketmgr->deallocateBasket(std::move(bh));
          return basket_holder_t{ no_init };
        } else {
          return bh;
        }
      }

      GeometryPtr m_geom;
      SourcePtr m_src;
      shared_obj<basketmgr_t> m_basketmgr;
      Optional<std::mutex> m_srcmutex;
      std::atomic<bool> m_srcHalted;
      bool m_srcParticlesMightBeOutside;


    };
  }
}

#endif
