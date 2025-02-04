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
  NCRYSTAL_DEBUGMMCMSG("StdEngine::advanceSimulation inbasket.size()="<<inbasket_holder.basket().size());

  nc_assert( inbasket_holder.valid() );
  nc_assert( !inbasket_holder.basket().empty() );
  auto& inbasket = inbasket_holder.basket();//Not const, since we will update
                                            //xsects in-place.
  const bool has_scat = !m_mat.scatter->isNull();
  const bool has_abs = !m_mat.absorption->isNull();
  const bool scatter_is_isotropic = !m_mat.scatter->isOriented();
  const bool absorption_is_isotropic = !m_mat.absorption->isOriented();


  //Get distances out for all the particles:
  geom.distToVolumeExit( inbasket.neutrons, m_buf_disttoexit );

  //Get absorption cross sections:
  const double * values_abs_xs_or_nullptr = nullptr;
  if ( has_abs ) {
    if ( absorption_is_isotropic ) {
      ProcImpl::NewABI::evalManyXSIsotropic( m_mat.absorption,
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
                                          m_buf_disttoexit,
                                          values_abs_xs_or_nullptr );

    resultFct( inbasket_holder.basket() );
    deallocateBasket( mgr, std::move(inbasket_holder) );
    return;
  }

  //Get scattering cross sections (potentially with cached values):
  if ( scatter_is_isotropic ) {
    for ( auto i : ncrange( inbasket.size() ) ) {
      if ( inbasket.cache.scatxsval[i] < 0.0 )
        inbasket.cache.scatxsval[i]
          = m_mat.scatter->crossSectionIsotropic( m_sct_cacheptr,
                                                  inbasket.neutrons.ekin_obj(i) ).dbl();
    }
  } else {
    //not isotropic, always recalculate all xs values:
    for ( auto i : ncrange( inbasket.size() ) ) {
      inbasket.cache.scatxsval[i]
        = m_mat.scatter->crossSection( m_sct_cacheptr,
                                       inbasket.neutrons.ekin_obj(i),
                                       inbasket.neutrons.dir_obj(i) ).dbl();
    }
  }

  //Transmission probability:
  MiniMC::Utils::calcProbTransm( m_mat.numDens,
                                 inbasket.size(),
                                 inbasket.cache.scatxsval,
                                 m_buf_disttoexit,
                                 m_buf_ptransm );

  //Pick scattering points:
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
    basket_holder_t pending = allocateBasket(mgr);
    nc_assert( pending.valid() && pending.basket().empty() );
    auto& outb = pending.basket();

    auto& inb = inbasket;
    nc_assert( m_opt.roulette_survival_probability < 1.0 );

    for ( auto i : ncrange( inb.size() ) ) {
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

      //Process this particle further, i.e. copy it to the pending
      //basket and update the weight (and then scatter it below):
      nc_assert( outb.size() < basket_N );
      std::size_t j = outb.appendEntryFromOther( inbasket, i );
      outb.neutrons.w[j] *= roulette_weight_factor;

      //Move to scattering point and attenuate:
      {
        const double disttoscat = m_buf_disttoscat[i];
        outb.neutrons.x[j] += disttoscat * outb.neutrons.ux[j];
        outb.neutrons.y[j] += disttoscat * outb.neutrons.uy[j];
        outb.neutrons.z[j] += disttoscat * outb.neutrons.uz[j];
        if ( values_abs_xs_or_nullptr ) {
          const double xsval_abs = values_abs_xs_or_nullptr[i];
          outb.neutrons.w[j] *= std::exp( -macroXS( m_mat.numDens, CrossSect{ xsval_abs } ) * disttoscat );
        }
      }

      //Scatter:
      nc_assert( has_scat );
      auto outcome = m_mat.scatter->sampleScatter( m_sct_cacheptr,
                                                   rng,
                                                   outb.neutrons.ekin_obj(j),
                                                   outb.neutrons.dir_obj(j));
      outb.neutrons.ux[j] = outcome.direction[0];
      outb.neutrons.uy[j] = outcome.direction[1];
      outb.neutrons.uz[j] = outcome.direction[2];
      bool was_elastic = (outb.neutrons.ekin[j] == outcome.ekin.dbl());
      outb.neutrons.ekin[j] = outcome.ekin.dbl();
      if ( was_elastic ) {
        outb.cache.markScatteredElastic(j);
      } else {
        outb.cache.markScatteredInelastic(j);
        outb.cache.scatxsval[j] = -1.0;//xs might have changed
      }
      //Needless since we never use the cached xs in case of oriented materials:
      //if (!scatter_is_isotropic)
      //  outb.cache.scatxsval[j] = -1.0;

      //Fix weights for the forced collision:
      outb.neutrons.w[j] *= ( 1.0 - m_buf_ptransm[i] );

    }
    if ( !outb.empty() ) {
      mgr.addPendingBasket( std::move(pending) );
    } else {
      deallocateBasket( mgr, std::move(pending) );
    }
  }

  {
    //Add a result basket with directly transmitted particles. For
    //efficiency, we simply reuse the input basket:
    auto& outb = inbasket_holder.basket();
    MiniMC::Utils::propagateAndAttenuate( outb.neutrons,
                                          m_mat.numDens,
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
