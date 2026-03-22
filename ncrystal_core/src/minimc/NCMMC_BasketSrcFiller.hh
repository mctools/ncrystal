#ifndef NCrystal_MMC_BasketSrcFiller_hh
#define NCrystal_MMC_BasketSrcFiller_hh

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

#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/minimc/NCMMC_Source.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //A helper class for "topping up" pending baskets with new baskets from a
    //source, and to handle the propagation of src particles to the modelled
    //geometry.

    namespace detail {
      //propagator fct which supports offset:
      void propagateDistance( NeutronBasket& nb,
                              std::size_t basket_offset,
                              const BasketValBufDbl& ) ncnoexceptndebug;
    }

    template<class TBasket>
    class BasketSrcFiller final : NoCopyMove {
    public:
      using basket_t = TBasket;

      BasketSrcFiller( GeometryPtr geom,
                       SourcePtr src )
        : m_geom(std::move(geom)),
          m_src(std::move(src)),
          m_srcParticlesMightBeOutside(m_src->particlesMightBeOutside(m_geom))
      {
      }

      //Main function. returns false if source ran out (or has run out earier):
      bool fillFromSource( basket_t& b,
                           RNG& rng,
                           const std::function<void(const basket_t&)>& resultFct,
                           ParticleCountSum& missCount,
                           basket_t& extra_basket_buffer,
                           unsigned nretry = 10 ) {

        b.validateIfDbg();
        const std::size_t size_orig = b.size();

        //Ok, we need to try to fill the basket from the source:
        bool src_has_more = true;
        m_src->fillBasket( rng, b.get_neutrons() );

        //Source provided neutron parameters, but we still need to initialise
        //any other parameters:
        for( std::size_t i = size_orig; i < b.size(); ++i)
          b.init_extra( i );

        b.validateIfDbg();
        if ( !b.full() )
          src_has_more = false;

        //If source particles might be outside volume, we have to propagate them
        //to the volume (and record the rest as results already):
        b.validateIfDbg();
        if ( m_srcParticlesMightBeOutside ) {
          propagateToVolume( b, size_orig, resultFct, extra_basket_buffer,
                             missCount );
          b.validateIfDbg();
        }

        if ( src_has_more && b.size() < basket_N_almost_Full ) {
          //Not quite full, lots of misses. We can try again a few times:
          if ( nretry > 0 ) {
            return fillFromSource( b, rng, resultFct, missCount,
                                   extra_basket_buffer, nretry - 1 );
          } else {
            NCRYSTAL_THROW(CalcError,"Source particles consistently "
                           "seem to miss the geometry.");
          }
        }
        return src_has_more;
      }

    private:
      void propagateToVolume( basket_t& b, std::size_t offset,
                              const std::function<void(const basket_t&)>& resultFct,
                              basket_t& extra_basket_buffer,
                              ParticleCountSum& missCount ) {
        //We must edit the basket, propagating all the neutrons to the volume if
        //they can. Those that can NOT do that, should be marked as having
        //missed the volume and served up as a result (unless ignoring misses
        //completely, indicated by resultFct == nullptr).

        if ( offset >= b.size() )
          return;

        b.validateIfDbg();

        const bool ignore_miss = (resultFct == nullptr);

        //First do the geometry distance calculations:
        nc_assert(m_srcParticlesMightBeOutside);
        BasketValBufDbl dist_results;
        m_geom->distToVolumeEntry( b.get_neutrons(), dist_results, offset );

        //Next, prepare a basket for results (those that missed):
        basket_t * bmiss = nullptr;
        if ( !ignore_miss ) {
          bmiss = &extra_basket_buffer;
          bmiss->get_neutrons().nused = 0;
        }

        std::size_t i_first_hole = basket_N;
        const std::size_t b_size = b.size();
        ParticleCountSum missStat;

        for ( std::size_t i = offset; i < b_size; ++i ) {
          double dist = dist_results[i];
          if ( dist < 0.0 ) {
            //Missed the volume, register as result (or just ignore). In any
            //case, record the miss for statistics.
            missStat.weight += b.get_neutrons().fields.w[i];
            ++(missStat.count);
            if (bmiss)
              bmiss->append1( b, i );
            if ( i_first_hole == basket_N )
              i_first_hole = i;
          } else {
            //Keep this. In case of holes, we also have to move it.
            std::size_t inew = i;
            if ( i_first_hole < i ) {
              b.copy1( b, i, (inew=i_first_hole++) );
              dist_results.data[inew] = dist;
            }
          }
        }
        missCount.weight += missStat.weight;
        missCount.count += missStat.count;

        //shrink-to-fit:
        if ( i_first_hole != basket_N )
          b.get_neutrons().nused = i_first_hole;

        nc_assert( offset <= b.get_neutrons().nused );

        //Now we should propagate all the neutrons that were not already inside
        //the volume (those already inside have dist_results[i]=0, so the
        //propagation won't do anything)::
        detail::propagateDistance( b.get_neutrons(), offset, dist_results );

        //And finally, return any neutrons that missed as a result, after
        //marking them as having missed the target.
        if ( bmiss && !bmiss->empty() ) {
          for ( auto i : ncrange( bmiss->size() ) )
            bmiss->markAsMissedTarget(i);
          resultFct( *bmiss );
        }
      }

      GeometryPtr m_geom;
      SourcePtr m_src;
      bool m_srcParticlesMightBeOutside;
    };
  }
}

#endif
