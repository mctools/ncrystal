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

NCMMC::StdEngine::StdEngine( matdef_t md, StdEngineOptions opts )
  : m_opt( std::move(opts) ),
    m_mat( std::move(md) )
{
  if ( ! ( m_opt.roulette_survival_probability > 1e-20 ) )
    NCRYSTAL_THROW(BadInput,"roulette_survival_probability must be >1e-20");

  if ( ! ( m_opt.roulette_survival_probability < 1.0 ) )
    NCRYSTAL_THROW(BadInput,"roulette_survival_probability must be <1.0");

  if ( ! ( m_opt.roulette_weight_threshold > 0.0 ) )
    NCRYSTAL_THROW(BadInput,"roulette_weight_threshold must be >0.0");

  //derived values:
  m_opt_roulette_survivor_boost = 1.0 / m_opt.roulette_survival_probability;
}

void NCMMC::StdEngine::advanceSimulation( RNG& rng,
                                         const Geometry& geom,
                                         basket_holder_t&& inbasket_holder,
                                         basketmgr_t& mgr,
                                         const resultfct_t& resultFct )
{
  NCRYSTAL_DEBUGMMCMSG("StdEngine::advanceSimulation inbasket.size()="
                       <<inbasket_holder.basket().size());

  const auto& eopts = m_opt.general_options;
  nc_assert_always( eopts.nScatLimit.value_or(0)
                    <= (std::numeric_limits<int>::max()-1) );
  const int nscatlimit = ( eopts.nScatLimit.has_value()
                         ? static_cast<int>(eopts.nScatLimit.value() )
                         : -1 );
  nc_assert( resultFct != nullptr );
  nc_assert( inbasket_holder.valid() );
  nc_assert( !inbasket_holder.basket().empty() );
  auto& inbasket = inbasket_holder.basket();//Not const, since we will update
                                            //xsects in-place.
  const bool has_scat = !( m_mat.scatter == nullptr
                           || m_mat.scatter->isNull() );
  const bool has_abs = !( m_mat.absorption == nullptr
                          || m_mat.absorption->isNull() );
  const bool scatter_is_isotropic = !(has_scat&&m_mat.scatter->isOriented());
  const bool absorption_is_isotropic = !(has_abs&&m_mat.absorption->isOriented());
  const bool geom_is_unbounded = geom.hasUnboundedDistToVolExit();

  //Get distances out for all the particles (note, if geom_is_unbounded, this
  //might include infinities):
  geom.distToVolumeExit( inbasket.neutrons, m_buf_disttoexit );

  //Get absorption cross sections:
  const double * values_abs_xs_or_nullptr = nullptr;
  if ( has_abs ) {
    if ( absorption_is_isotropic ) {
      ProcImpl::NewABI::evalManyXSIsotropic( *m_mat.absorption,
                                             m_abs_cacheptr,
                                             inbasket.neutrons.ekin,
                                             inbasket.size(),
                                             m_buf_xs_abs );
    } else {
      for ( auto i : ncrange( inbasket.size() ) ) {
        m_buf_xs_abs[i]
          = m_mat.absorption->crossSection( m_abs_cacheptr,
                                            inbasket.neutrons.ekin_obj(i),
                                            inbasket.neutrons.dir_obj(i) ).dbl();
      }
    }
    values_abs_xs_or_nullptr = &m_buf_xs_abs[0];
  }

  if (!has_scat) {
    //Special case of no scattering, just transmit (in-place) and return:
    MiniMC::Utils::propagateAndAttenuate( inbasket_holder.basket().neutrons,
                                          m_mat.numDens,
                                          geom_is_unbounded,
                                          m_buf_disttoexit,
                                          values_abs_xs_or_nullptr );

    resultFct( inbasket_holder.basket() );
    deallocateBasket( mgr, std::move(inbasket_holder) );
    return;
  }

  //Get scattering cross sections (potentially with cached values). Also, if
  //nscatlimit is enabled, anything already scattered nscatlimit times will get
  //a phony scattering cross section value of 0:

  if ( scatter_is_isotropic ) {
    if ( nscatlimit>=0 ) {
      for ( auto i : ncrange( inbasket.size() ) )
        if ( inbasket.cache.nscat[i] >= nscatlimit )
          inbasket.cache.scatxsval[i] = 0.0;
    }
    for ( auto i : ncrange( inbasket.size() ) ) {
      if ( inbasket.cache.scatxsval[i] < 0.0 )
        inbasket.cache.scatxsval[i]
          = m_mat.scatter->crossSectionIsotropic( m_sct_cacheptr,
                                                  inbasket.neutrons.ekin_obj(i) ).dbl();
    }
  } else {
    //not isotropic, always recalculate all xs values (except if nscatlimit is
    //exceeded of course):
    for ( auto i : ncrange( inbasket.size() ) ) {
      if ( nscatlimit >= 0 && inbasket.cache.nscat[i] >= nscatlimit )
        inbasket.cache.scatxsval[i] = 0.0;
      else
        inbasket.cache.scatxsval[i]
          = m_mat.scatter->crossSection( m_sct_cacheptr,
                                         inbasket.neutrons.ekin_obj(i),
                                         inbasket.neutrons.dir_obj(i) ).dbl();
    }
  }

  //Transmission probability (note, as a special case this will be set to 0 if
  //disttoexit=inf):
  MiniMC::Utils::calcProbTransm( m_mat.numDens,
                                 inbasket.size(),
                                 geom_is_unbounded,
                                 inbasket.cache.scatxsval,
                                 m_buf_disttoexit,
                                 m_buf_ptransm );

  //FIXME: Make sure we handle xs_scat=0 correctly!!

  //fixme: stop recalculating macro_xs, but immediately convert xs values to
  //macroscopic values right after query!

  //Pick scattering points (note returns kInfinity if macro scat xs=0):

  MiniMC::Utils::sampleRandDists(rng,
                                 m_mat.numDens,
                                 m_buf_disttoexit,
                                 inbasket.cache.scatxsval,
                                 inbasket.size(),
                                 m_buf_disttoscat );

  {
    //Add a pending basket with scattered particles for further simulations. We
    //also implement a russian roulette particle-killing scheme (otherwise the
    //simulations would never terminate with our forced-collision scheme).
    //
    //Finally, we skip any particle with weight=0 or vanishing scattering cross
    //section.

    basket_holder_t pending = allocateBasket(mgr);
    nc_assert( pending.valid() && pending.basket().empty() );
    auto& outb = pending.basket();

    auto& inb = inbasket;
    nc_assert( m_opt.roulette_survival_probability < 1.0 );

    for ( auto i : ncrange( inb.size() ) ) {
      if ( inb.neutrons.w[i]==0.0 )
        continue;//no actual flux here, killed!

      if ( std::isinf(m_buf_disttoscat[i] ) )
        continue;//this comes from no scattering cross section (fixme: assert
                 //that this is so).

      const double macro_scat_xs
        = macroXS( m_mat.numDens,
                   CrossSect{ inbasket.cache.scatxsval[i] } );
      if (!(macro_scat_xs > 0.0))
        continue;//can't scatter (all relevant weight should have been delivered
                 //to the tally already via ptransm).
      //Implement russian-roulette:
      double roulette_weight_factor = 1.0;
      double roulette = ( ( inb.cache.nscat[i] >= m_opt.roulette_nscat_threshold
                            && inb.neutrons.w[i] < m_opt.roulette_weight_threshold )
                          ? m_opt.roulette_survival_probability
                          : 1.0 );
      if ( roulette < 1.0 ) {
        if ( rng.generate() > roulette ) {
          continue;//killed!
        } else {
          nc_assert( m_opt.roulette_survival_probability > 0.0 );
          nc_assert( m_opt.roulette_survival_probability <= 1.0 );
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
      //So we process this particle further, i.e. copy it to the pending basket
      //and update the weight (and then scatter it below):
      nc_assert( outb.size() < basket_N );
      std::size_t j = outb.appendEntryFromOther( inbasket, i );
      outb.neutrons.w[j] *= roulette_weight_factor;

      //Move to scattering point and attenuate: (fixme: can disttoscat be inf if xs=0??)
      outb.neutrons.x[j] += m_buf_disttoscat[i] * outb.neutrons.ux[j];
      outb.neutrons.y[j] += m_buf_disttoscat[i] * outb.neutrons.uy[j];
      outb.neutrons.z[j] += m_buf_disttoscat[i] * outb.neutrons.uz[j];
      outb.neutrons.w[j] *= weight_reduction_factor;

      //Scatter:
      nc_assert( has_scat );
      auto outcome = m_mat.scatter->sampleScatter( m_sct_cacheptr,
                                                   rng,
                                                   outb.neutrons.ekin_obj(j),
                                                   outb.neutrons.dir_obj(j));
      nc_assert( ncabs(outcome.direction.as<Vector>().mag()-1) < 1e-9 );
      outb.neutrons.ux[j] = outcome.direction[0];
      outb.neutrons.uy[j] = outcome.direction[1];
      outb.neutrons.uz[j] = outcome.direction[2];
      const bool was_elastic = (outb.neutrons.ekin[j] == outcome.ekin.dbl());
      outb.neutrons.ekin[j] = outcome.ekin.dbl();
      if ( was_elastic ) {
        outb.cache.markScatteredElastic(j);
      } else {
        outb.cache.markScatteredInelastic(j);
        outb.cache.scatxsval[j] = -1.0;//xs might have changed
      }

      //Fix weights for the forced collision. Note that when disttoexit=inf, we
      //have ptransm=0 (to avoid anything transmitted cropping up in the
      //tallies), so it is right that the final factor in (1-ptransm)=1. Nothing
      //made it to the tally, we kept all for the scattering. In case scattering
      //XS was also 0, that particle was simply lost (FIXME: This will have been
      //counted as "absorbed" in our plots. We should calculate the "lost"
      //separately? Or at least think carefully about what the correct behaviour
      //is and then test that we get it)
      outb.neutrons.w[j] *= ( 1.0 - m_buf_ptransm[i] );
    }
    if ( !outb.empty() ) {
      mgr.addPendingBasket( std::move(pending) );
    } else {
      deallocateBasket( mgr, std::move(pending) );
    }
  }

  {
    //Add a result basket with directly transmitted particles. For efficiency,
    //we simply reuse the input basket (in case of infinite disttoexit we end up
    //with w=0, but these will be ignored later):
    auto& outb = inbasket_holder.basket();
    MiniMC::Utils::propagateAndAttenuate( outb.neutrons,
                                          m_mat.numDens,
                                          geom_is_unbounded,
                                          m_buf_disttoexit,
                                          values_abs_xs_or_nullptr );
    //We also reduce with the transmission probability (i.e. scatter-process
    //attenuation, the above propagateAndAttenuate only took care of the
    //absorption-process attenuation:
    for ( auto i : ncrange(outb.size()) )
      outb.neutrons.w[i] *= m_buf_ptransm[i];

    resultFct( inbasket_holder.basket() );
    deallocateBasket( mgr, std::move(inbasket_holder) );
    return;
  }
}
