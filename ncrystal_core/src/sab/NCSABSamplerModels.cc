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

#include "NCrystal/internal/sab/NCSABSamplerModels.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
namespace NC = NCrystal;

NC::SAB::SABSamplerAtE_Alg1::SABSamplerAtE_Alg1( std::shared_ptr<const CommonCache> common,
                                                 VectD&& betaVals,
                                                 VectD&& betaWeights,
                                                 std::vector<AlphaSampleInfo>&& alphaSamplerInfos,
                                                 std::size_t ibetaOffset,
                                                 double firstBinKinematicEndpointValue )
  : m_common( std::move(common) ),
    m_betaSampler(VectD(betaVals.begin(),betaVals.end()),//todo: in principle no need to copy here.
                  VectD(betaWeights.begin(),betaWeights.end()) ),
    m_alphaSamplerInfos( std::move(alphaSamplerInfos) ),
    m_ibetaOffset( ibetaOffset ),
    m_firstBinKinematicEndpointValue(firstBinKinematicEndpointValue)
{
  nc_assert( !!m_common );
  nc_assert( betaWeights.size() == betaVals.size() );

  //+1 in the next two asserts since vals,weights starts with (beta_lower,0.0):
  nc_assert( m_alphaSamplerInfos.size()+1 == betaVals.size() );
  nc_assert( ibetaOffset+betaVals.size() == m_common->data->betaGrid().size()+1 );
}

NC::PairDD NC::SAB::SABSamplerAtE_Alg1::sampleAlphaBeta(double ekin_div_kT, RNG&rng) const
{
  nc_assert(!!m_common);
  const auto& betaGrid = m_common->data->betaGrid();
  //  const auto& alphaGrid = m_common->data->alphaGrid();
  nc_assert(m_ibetaOffset<betaGrid.size());

  //Allow only loopmax sample attempts, to make sure we detect if code gets too
  //inefficient. However, make sure users can override this if needed.
  static const unsigned s_loopmax = ncgetenv_int("SABSAMPLE_LOOPMAX", 100 );

  unsigned iloopmax(s_loopmax+1);
  while (--iloopmax) {
    double beta;
    unsigned ibetaSampled;
    std::tie(beta,ibetaSampled) = m_betaSampler.percentileWithIndex( rng() );

    nc_assert( !ncisnan(beta) );
    nc_assert( beta <= betaGrid.back() );
    nc_assert( ibetaSampled < m_betaSampler.getXVals().size() );

    if ( ibetaSampled == 0 && m_firstBinKinematicEndpointValue <= 0.0 ) {
      //Resample the first starting at m_firstBinKinematicEndpointValue rather
      //than the fake initial value in the beta sampler (see also comment in
      //NCSABIntegrator.cc). Assuming that the first actual alpha sampler
      //contains the best bet for the alpha-dependency over the entire
      //rectangle, we first sample alpha according to the first alpha
      //sampler. Beta is then sampled uniformly, and we discard values outside
      //the kinematic boundary.
      const double b0 = m_firstBinKinematicEndpointValue;
      const double b1 = vectAt(m_betaSampler.getXVals(),1);
      if ( b1 < -ekin_div_kT )
        continue;//reject no matter what
      nc_assert( b0 > m_betaSampler.getXVals().at(0) );//because we moved m_betaSampler.getXVals()[0] down by 4/3*(b1-b0)
      const double delta_beta = b1 - b0;
      nc_assert(delta_beta>0.0);
      double alphaval;
      constexpr auto nsampletries = 30;
      for ( auto iii : ncrange(nsampletries) ) {
        (void)iii;
        beta = ncmax(m_firstBinKinematicEndpointValue, b0 + delta_beta*rng.generate());
        if ( beta < -ekin_div_kT )
          break;//reject
        alphaval = sampleAlpha(m_ibetaOffset, rng.generate());
        auto alims = getAlphaLimits(-m_firstBinKinematicEndpointValue,beta);
        if ( valueInInterval(alims.first,alims.second,alphaval) )
          break;
        constexpr auto nsampletriesm1 = nsampletries - 1;
        if ( iii == nsampletriesm1 ) {
          //Too inefficient. Fall back to isotropic alpha at the given beta.
          static std::atomic<unsigned> s_nfail{0};
          static constexpr unsigned nfail_max_warn = 20;
          auto nfail = ++s_nfail;
          if ( nfail <= nfail_max_warn ) {
            NCRYSTAL_WARN("SABSampler reverts to isotropic model"
                          " after "<<nsampletries<<" rejected attempts"
                          << ( nfail == nfail_max_warn ? " (suppressing"
                               " further warnings of this type)" : "" ));
          }
          alphaval = 0.5*(alims.first+alims.second);
          break;
        }
      }
      if ( beta < -ekin_div_kT )
        continue;//reject
      auto alimits = getAlphaLimits( ekin_div_kT, beta );
      if ( valueInInterval( alimits.first, alimits.second, alphaval ) )
        return { alphaval, beta };//accept
      continue;//reject
    }

    if (beta <= ncmax(-ekin_div_kT,betaGrid.front()))
      continue;//reject

    double alphal(-1.0), bl;
    std::size_t ibeta;

    double rand_percentile = rng.generate();
    nc_assert ( beta >= m_betaSampler.getXVals().at(1) );
    ibeta = m_ibetaOffset + ibetaSampled;
    nc_assert( ibeta>0 );
    bl = betaGrid.at(ibeta-1);
    alphal = sampleAlpha(ibeta-1, rand_percentile);
    nc_assert( valueInInterval(betaGrid.at(ibeta-1),betaGrid.at(ibeta),beta) );
    nc_assert( alphal>=0.0 );
    nc_assert( ibeta < betaGrid.size() );

    //Sample alphas (with same "random percentile" at two neighbouring beta grid
    //points and combine with linear interpolation, as in line 9 in Algorithm 1
    //of the paper:
    const double bh = betaGrid.at(ibeta);
    const double alphah = sampleAlpha(ibeta, rand_percentile);
    nc_assert(bh-bl>0.0);
    nc_assert(alphal>=0.0);
    nc_assert(alphah>=0.0);
    const double alpha = alphal + (alphah-alphal) * (beta-bl)/(bh-bl);

    //Check if we can accept this:
    auto alimits = getAlphaLimits( ekin_div_kT, beta );
    if ( valueInInterval( alimits.first, alimits.second, alpha ) )
      return { alpha, beta };//accept
  }
  NCRYSTAL_THROW2(CalcError,"Rejection method failed to sample kinematically valid (alpha,beta) point after "
                  <<s_loopmax<<" attempts. Perhaps energy grid is too sparse?"
                  " As a workaround it is possible to increase the allowed number of sampling attempts"
                  " by setting the NCRYSTAL_SABSAMPLE_LOOPMAX variable to a higher number"
                  " (but please consider reporting the issue to the NCrystal developers nonetheless).");
}

double NC::SAB::SABSamplerAtE_Alg1::sampleAlpha(std::size_t ibeta, double rand_percentile) const
{
  nc_assert( ibeta >= m_ibetaOffset );
  const auto& info = vectAt(m_alphaSamplerInfos,ibeta-m_ibetaOffset);

  const auto& cd = m_common->data;
  auto nalpha = cd->alphaGrid().size();
  auto cumul = SABUtils::sliceSABAtBetaIdx_const(m_common->alphaintegrals_cumul,nalpha,ibeta);
  auto sab = SABUtils::sliceSABAtBetaIdx_const(cd->sab(),nalpha,ibeta);
  auto logsab = SABUtils::sliceSABAtBetaIdx_const(m_common->logsab,nalpha,ibeta);
  auto clampRandNum = [](double r) { return ncclamp( r, std::numeric_limits<double>::min(), 1.0 ); };//ensure r is in (0,1]
  auto clampUnitInterval = [](double r) { return ncclamp( r, 0.0, 1.0 ); };//ensure r is in [0,1]

  if ( rand_percentile <= info.prob_front) {
    if ( info.prob_front == 2.0) {
      //special value indicating 0 cross-section at value, sample linearly in
      //[pt_front.alpha,pt_back.alpha] for lack of better options..
      double da = info.pt_back.alpha-info.pt_front.alpha;
      nc_assert(da>=0.0);
      return info.pt_front.alpha + rand_percentile*da;
    } else if ( info.prob_front == 1.0) {
      //Valid alpha range is narrow and contained within a single bin.
      return SABUtils::sampleLogLinDist_fast( info.pt_front.alpha, info.pt_front.sval,
                                              info.pt_back.alpha, info.pt_back.sval,
                                              rand_percentile,
                                              info.pt_front.logsval, info.pt_back.logsval );
    } else {
      //Sample front tail
      double percentile2 = clampRandNum( rand_percentile / info.prob_front );
      return SABUtils::sampleLogLinDist_fast( info.pt_front.alpha, info.pt_front.sval,
                                              vectAt(cd->alphaGrid(),info.pt_front.alpha_idx), sab[info.pt_front.alpha_idx],
                                              percentile2,
                                              info.pt_front.logsval, logsab[info.pt_front.alpha_idx] );
    }
  } else if ( rand_percentile <= info.prob_notback ) {
    //Middle section - sample over entire alpha bins.
    nc_assert( info.prob_notback - info.prob_front > 0.0 );
    double percentile2 = clampUnitInterval( ( rand_percentile - info.prob_front ) / ( info.prob_notback - info.prob_front ) );
    unsigned alphaidx_low(info.pt_front.alpha_idx), alphaidx_upp(info.pt_back.alpha_idx);
    auto itCumul_low = std::next(cumul.begin(),alphaidx_low);
    auto itCumul_upp = std::next(cumul.begin(),alphaidx_upp);
    nc_assert( itCumul_upp > itCumul_low && itCumul_upp<cumul.end() );
    double selectedArea = *itCumul_low + percentile2 * ( *itCumul_upp - *itCumul_low );
    auto itCumul_selected_edgeupp = std::upper_bound(itCumul_low, std::next(itCumul_upp), selectedArea);
    if ( itCumul_selected_edgeupp > itCumul_upp )
      return vectAt( cd->alphaGrid(), alphaidx_upp );
    if ( itCumul_selected_edgeupp <= itCumul_low )
      return vectAt( cd->alphaGrid(), alphaidx_low );

    nc_assert(itCumul_selected_edgeupp>itCumul_low);
    auto itCumul_selected_edgelow = std::prev(itCumul_selected_edgeupp);
    nc_assert( *itCumul_selected_edgelow <= selectedArea );
    nc_assert( *itCumul_selected_edgeupp >= selectedArea );
    double binArea = *itCumul_selected_edgeupp - *itCumul_selected_edgelow;
    nc_assert( binArea > 0.0);
    //rescale leftover parts of rand_percentile back to the unit interval:
    double rand_rescaled = clampRandNum((selectedArea-*itCumul_selected_edgelow)/binArea);
    nc_assert( rand_rescaled >= 0.0 );
    nc_assert( rand_rescaled <= 1.0 );
    auto a0 = itCumul_selected_edgelow - cumul.begin();
    auto a1 = a0 + 1;
    //Interpolate in selected bin for alpha value:
    return SABUtils::sampleLogLinDist_fast( vectAt(cd->alphaGrid(),a0), sab[a0],
                                            vectAt(cd->alphaGrid(),a1), sab[a1],
                                            rand_rescaled,
                                            logsab[a0], logsab[a1] );
  } else {
    //Sample back tail
    nc_assert( 1.0 - info.prob_notback > 0.0 );
    double percentile2 = clampRandNum ( ( rand_percentile - info.prob_notback ) / ( 1.0 - info.prob_notback ) );
    return SABUtils::sampleLogLinDist_fast( vectAt(cd->alphaGrid(),info.pt_back.alpha_idx), sab[info.pt_back.alpha_idx],
                                            info.pt_back.alpha, info.pt_back.sval,
                                            percentile2,
                                            logsab[info.pt_back.alpha_idx], info.pt_back.logsval );

  }

}
