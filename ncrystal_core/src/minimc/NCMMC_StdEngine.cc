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

#include "NCrystal/internal/minimc/NCMMC_StdEngine.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/extd_utils/NCABIUtils.hh"
#include "NCrystal/internal/minimc/NCMMC_Utils.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
      class StdSimEngine final : public SimEngine {
        //The input:
        EngineOpts m_opt;
        GeometryPtr m_geom;
        MatDef m_mat;

        //Derived values and caches:
        double m_opt_roulette_survivor_boost;
        CachePtr m_sct_cacheptr;
        CachePtr m_abs_cacheptr;

        //We need a few buffers (careful with memory usage, but better here on
        //the heap than on the stack):
        BasketValBufDbl m_buf_disttoexit;
        BasketValBufDbl m_buf_xs_abs;
        BasketValBufDbl m_buf_ptransm;
        BasketValBufDbl m_buf_disttoscat;

        static void checkBasketFields( const UniversalBasket& b )
        {
          if ( ! ( b.buf1 && b.nscat && b.nscat_inelas ) )
            NCRYSTAL_THROW(LogicError,"StdSimEngine requires a basket with"
                           " nscat, nscat_inelas, and buf1 available.");
        }

      public:

        shared_obj<SimEngine> clone() const override
        {
          return makeSO<StdSimEngine>( m_opt, m_geom, m_mat );
        }

        StdSimEngine( const EngineOpts& eopts,
                      GeometryPtr geom,
                      MatDef mat )
          : m_opt( eopts ),
            m_geom( std::move(geom) ),
            m_mat( std::move(mat) )
        {
          if ( ! ( m_opt.roulette.survival_probability > 1e-20 ) )
            NCRYSTAL_THROW(BadInput,"roulette_survival_probability must be >1e-20");
          if ( ! ( m_opt.roulette.survival_probability < 1.0 ) )
            NCRYSTAL_THROW(BadInput,"roulette_survival_probability must be <1.0");
          if ( ! ( m_opt.roulette.weight_threshold > 0.0 ) )
            NCRYSTAL_THROW(BadInput,"roulette_weight_threshold must be >0.0");
          m_opt_roulette_survivor_boost
            = 1.0 / m_opt.roulette.survival_probability;
          nc_assert_always( m_opt.nScatLimit.value_or(0)
                            <= (std::numeric_limits<int>::max()-1) );
        }

        void step( UniversalBasket inbasket,
                   RNG& rng,
                   UniversalBasketMgr& mgr,
                   const TallyFct& tallyfct ) override
        {
          //In this engine, each step produce 1 basket for tallying and 1
          //pending basket for further processing in a new step. In case the
          //pending basket does not need filling, we immediately process it
          //again, thus reducing the need to communicate with the global
          //thread-safe queue in the mgr.
          //
          //Additionally, and for the same reason, we provide a single buffer
          //basket to be reused between these calls.
          UniversalBasket buf;
          do {
            inbasket = processBasket( std::move(inbasket), buf,
                                      rng, mgr, tallyfct );
          } while ( inbasket.valid()
                    && inbasket.size() >= basket_N_almost_Full );
          if ( inbasket.valid() ) {
              //transfer to global queue:
              mgr.addPendingBasket( std::move(inbasket) );
          }
          if ( buf.valid() )
            mgr.deallocateBasket( std::move(buf) );
        }

      private:
        UniversalBasket processBasket( UniversalBasket inbasket,
                                       UniversalBasket& basket_buffer,
                                       RNG& rng,
                                       UniversalBasketMgr& mgr,
                                       const TallyFct& tallyfct )
        {
          NCRYSTAL_DEBUGMMCMSG("StdSimEngine::processBasket inbasket.size()="
                               <<inbasket.size());

          UniversalBasket result_basket;
          nc_assert(!result_basket.valid());

          checkBasketFields( inbasket );

          const int nscatlimit = ( m_opt.nScatLimit.has_value()
                                   ? static_cast<int>(m_opt.nScatLimit.value() )
                                   : -1 );
          nc_assert( tallyfct != nullptr );
          nc_assert( inbasket.valid() );
          nc_assert( !inbasket.empty() );
          const bool has_scat = !( m_mat.scatter == nullptr
                                   || m_mat.scatter->isNull() );
          const bool has_abs = !( m_mat.absorption == nullptr
                                  || m_mat.absorption->isNull() );
          const bool scatter_is_isotropic = !(has_scat&&m_mat.scatter->isOriented());
          const bool absorption_is_isotropic = !(has_abs&&m_mat.absorption->isOriented());
          const bool geom_is_unbounded = m_geom->hasUnboundedDistToVolExit();

          //Get distances out for all the particles (note, if geom_is_unbounded,
          //this might include infinities):
          m_geom->distToVolumeExit( *inbasket.neutrons, m_buf_disttoexit.data );//fixme: migrate fully to buffers

          //Get absorption cross sections:
          const double * values_abs_xs_or_nullptr = nullptr;
          if ( has_abs ) {
            if ( absorption_is_isotropic ) {
              ProcImpl::NewABI::evalManyXSIsotropic( *m_mat.absorption,
                                                     m_abs_cacheptr,
                                                     inbasket.neutrons->fields.ekin.data,
                                                     inbasket.size(),
                                                     m_buf_xs_abs.data );
            } else {
              for ( auto i : ncrange( inbasket.size() ) ) {
                m_buf_xs_abs[i]
                  = m_mat.absorption->crossSection( m_abs_cacheptr,
                                                    inbasket.neutrons->ekin_obj(i),
                                                    inbasket.neutrons->dir_obj(i) ).dbl();
              }
            }
            values_abs_xs_or_nullptr = &m_buf_xs_abs[0];
          }

          if (!has_scat) {
            //Special case of no scattering, just transmit (in-place) and return
            //an empty basket (i.e. nothing needs further processing):
            MiniMC::Utils::propagateAndAttenuate( *inbasket.neutrons,
                                                  m_mat.numDens,
                                                  geom_is_unbounded,
                                                  m_buf_disttoexit.data,
                                                  values_abs_xs_or_nullptr );

            tallyfct( inbasket );
            nc_assert_always(!basket_buffer.valid());
            basket_buffer = std::move(inbasket);
            return result_basket;
          }

          //Get scattering cross sections (potentially with cached
          //values). Also, if nscatlimit is enabled, anything already scattered
          //nscatlimit times will get a phony scattering cross section value of
          //0:

          auto buf_scatxsval = [](const UniversalBasket&b) -> BasketValBufDbl&
          {
            nc_assert(b.buf1!=nullptr);
            return *b.buf1;
          };

          if ( nscatlimit == 0 ) {
            //special case, completely disable scattering.
            for ( auto i : ncrange( inbasket.size() ) )
              buf_scatxsval(inbasket)[i] = 0.0;
          } else {
            if ( scatter_is_isotropic ) {
              if ( nscatlimit>=0 ) {
                for ( auto i : ncrange( inbasket.size() ) )
                  if ( inbasket.nscat->data[i] >= nscatlimit )
                    buf_scatxsval(inbasket)[i] = 0.0;
              }
              for ( auto i : ncrange(nscatlimit==0?0:inbasket.size()) ) {
                if ( inbasket.nscat->data[i] <= 0
                     || buf_scatxsval(inbasket)[i] < 0.0 )
                  buf_scatxsval(inbasket)[i]
                    = m_mat.scatter->crossSectionIsotropic( m_sct_cacheptr,
                                                            inbasket.neutrons->ekin_obj(i) ).dbl();
              }
            } else {
              //not isotropic, always recalculate all xs values (except if
              //nscatlimit is exceeded of course):
              for ( auto i : ncrange( inbasket.size() ) ) {
                if ( nscatlimit >= 0 && inbasket.nscat->data[i] >= nscatlimit )
                  buf_scatxsval(inbasket)[i] = 0.0;
                else
                  buf_scatxsval(inbasket)[i]
                    = m_mat.scatter->crossSection( m_sct_cacheptr,
                                                   inbasket.neutrons->ekin_obj(i),
                                                   inbasket.neutrons->dir_obj(i) ).dbl();
              }
            }
          }

          //Transmission probability (note, as a special case this will be set
          //to 0 if disttoexit=inf):
          MiniMC::Utils::calcProbTransm( m_mat.numDens,
                                         inbasket.size(),
                                         geom_is_unbounded,
                                         buf_scatxsval(inbasket).data,
                                         m_buf_disttoexit.data,
                                         m_buf_ptransm.data );

          //FIXME: Make sure we handle xs_scat=0 correctly!!

          //fixme: stop recalculating macro_xs, but immediately convert xs
          //values to macroscopic values right after query!

          //Pick scattering points (note returns kInfinity if macro scat xs=0):

          MiniMC::Utils::sampleRandDists(rng,
                                         m_mat.numDens,
                                         m_buf_disttoexit.data,
                                         buf_scatxsval(inbasket).data,
                                         inbasket.size(),
                                         m_buf_disttoscat.data );

          {
            //Add a pending basket with scattered particles for further
            //simulations. We also implement a russian roulette particle-killing
            //scheme (otherwise the simulations would never terminate with our
            //forced-collision scheme).
            //
            //We also skip any particle with weight=0 or vanishing scattering
            //cross section.

            UniversalBasket outb = std::move(basket_buffer);
            if ( outb.valid() )
              outb.neutrons->nused = 0;
            else
              outb = mgr.allocateBasket();
            nc_assert( outb.valid() && outb.empty() );
            checkBasketFields( outb );
            auto& inb = inbasket;
            nc_assert( m_opt.roulette.survival_probability < 1.0 );

            for ( auto i : ncrange( inb.size() ) ) {
              if ( inb.neutrons->fields.w[i]==0.0 )
                continue;//no actual flux here, killed!

              if ( std::isinf(m_buf_disttoscat[i] ) )
                continue;//this comes from no scattering cross section (fixme: assert
              //that this is so).

              const double macro_scat_xs
                = macroXS( m_mat.numDens,
                           CrossSect{ buf_scatxsval(inbasket)[i] } );

              if (!(macro_scat_xs > 0.0))
                continue;//can't scatter (all relevant weight should have been
                         //delivered already).

              //to the tally already via ptransm).
              //Implement russian-roulette:
              double roulette_weight_factor = 1.0;
              double roulette = ( ( inb.nscat->data[i] >= m_opt.roulette.nscat_threshold
                                    && inb.neutrons->fields.w[i] < m_opt.roulette.weight_threshold )
                                  ? m_opt.roulette.survival_probability
                                  : 1.0 );
              if ( roulette < 1.0 ) {
                if ( rng.generate() > roulette ) {
                  continue;//killed!
                } else {
                  nc_assert( m_opt.roulette.survival_probability > 0.0 );
                  nc_assert( m_opt.roulette.survival_probability <= 1.0 );
                  roulette_weight_factor = m_opt_roulette_survivor_boost;
                }
              }

              //Get macroscopic scattering cross section:
              //fixme cache factor macroXS( m_mat.numDens, CrossSect{1.0}).
              const double macro_abs_xs
                = ( values_abs_xs_or_nullptr
                    ? macroXS( m_mat.numDens,
                               CrossSect{ *std::next(values_abs_xs_or_nullptr,i) } )
                    : 0.0 );

              nc_assert( !std::isinf(m_buf_disttoscat[i]));
              const double weight_reduction_factor
                = std::exp( -macro_abs_xs * m_buf_disttoscat[i] );

              if (!(weight_reduction_factor>0.0))
                continue;//did not survive to scatter point, kill!

              //We are now in the "standard" situation:
              //  macro_scat_xs > 0.0
              //  weight_reduction_factor > 0.0
              //  disttoscat < infinity
              //
              //So we process this particle further, i.e. copy it to the pending
              //basket and update the weight (and then scatter it below):
              nc_assert( outb.size() < basket_N );
              std::size_t j = outb.append1( inbasket, i );
              auto& outb_fields = outb.neutrons->fields;
              outb_fields.w[j] *= roulette_weight_factor;

              //Move to scattering point and attenuate: (fixme: can disttoscat
              //be inf if xs=0??)
              outb_fields.x[j] += m_buf_disttoscat[i] * outb_fields.ux[j];
              outb_fields.y[j] += m_buf_disttoscat[i] * outb_fields.uy[j];
              outb_fields.z[j] += m_buf_disttoscat[i] * outb_fields.uz[j];
              outb_fields.w[j] *= weight_reduction_factor;

              //Scatter:
              nc_assert( has_scat );
              auto outcome = m_mat.scatter->sampleScatter( m_sct_cacheptr,
                                                           rng,
                                                           outb.neutrons->ekin_obj(j),
                                                           outb.neutrons->dir_obj(j));
              nc_assert( ncabs(outcome.direction.as<Vector>().mag()-1) < 1e-9 );
              outb_fields.ux[j] = outcome.direction[0];
              outb_fields.uy[j] = outcome.direction[1];
              outb_fields.uz[j] = outcome.direction[2];
              const bool was_elastic = (outb_fields.ekin[j] == outcome.ekin.dbl());
              outb_fields.ekin[j] = outcome.ekin.dbl();
              ++( outb.nscat->data[j] );
              if ( !was_elastic ) {
                ++( outb.nscat_inelas->data[j] );
                buf_scatxsval(outb)[j] = -1.0;//xs might have changed
              }

              //Fix weights for the forced collision. Note that when
              //disttoexit=inf, we have ptransm=0 (to avoid anything transmitted
              //cropping up in the tallies), so it is right that the final
              //factor is (1-ptransm)=1. Nothing made it to the tally, we kept
              //all for the scattering. In case scattering XS was also 0, that
              //particle was simply lost
              outb_fields.w[j] *= ( 1.0 - m_buf_ptransm[i] );
            }

            result_basket = std::move(outb);
          }

          {
            //Add a result basket with directly transmitted particles. For
            //efficiency, we simply reuse the input basket (in case of infinite
            //disttoexit we end up with w=0, but these will be ignored later):
            auto& outb = inbasket;
            MiniMC::Utils::propagateAndAttenuate( *outb.neutrons,
                                                  m_mat.numDens,
                                                  geom_is_unbounded,
                                                  m_buf_disttoexit.data,
                                                  values_abs_xs_or_nullptr );
            //We also reduce with the transmission probability
            //(i.e. scatter-process attenuation, the above propagateAndAttenuate
            //only took care of the absorption-process attenuation:
            for ( auto i : ncrange(outb.size()) )
              outb.neutrons->fields.w[i] *= m_buf_ptransm[i];

            tallyfct( outb );
            //store the buffer for future use:
            nc_assert( !basket_buffer.valid() );
            basket_buffer = std::move(inbasket);
          }

          return result_basket;
        }
      };
    }
  }
}

NC::shared_obj<NCMMC::SimEngine>
NCMMC::createStdSimEngine( GeometryPtr geom,
                           MatDef mat,
                           const EngineOpts& opts )
{
  return makeSO<StdSimEngine>( opts, std::move(geom), std::move(mat) );
}
