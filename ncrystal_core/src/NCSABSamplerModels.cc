////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/internal/NCSABSamplerModels.hh"
#include "NCrystal/internal/NCSABUtils.hh"
#include "NCrystal/internal/NCString.hh"
namespace NC = NCrystal;

NC::SAB::SABSamplerAtE_Alg1::SABSamplerAtE_Alg1( std::shared_ptr<const CommonCache> common,
                                                 VectD&& betaVals,
                                                 VectD&& betaWeights,
                                                 std::vector<AlphaSampleInfo>&& alphaSamplerInfos,
                                                 std::size_t ibetaOffset )
  : m_common( std::move(common) ),
    m_betaSampler(VectD(betaVals.begin(),betaVals.end()),
                  VectD(betaWeights.begin(),betaWeights.end()) ),
    m_alphaSamplerInfos( std::move(alphaSamplerInfos) ),
    m_ibetaOffset( ibetaOffset )
{
  nc_assert( !!m_common );
  nc_assert( betaWeights.size() == betaVals.size() );

  //+1 in the next two asserts since vals,weights starts with (beta_lower,0.0):
  nc_assert( m_alphaSamplerInfos.size()+1 == betaVals.size() );
  nc_assert( ibetaOffset+betaVals.size() == m_common->data->betaGrid().size()+1 );
}

NC::PairDD NC::SAB::SABSamplerAtE_Alg1::sampleAlphaBeta(double ekin_div_kT, RandomBase&rng) const
{
  nc_assert(!!m_common);
  const auto& betaGrid = m_common->data->betaGrid();
  const auto& alphaGrid = m_common->data->alphaGrid();
  nc_assert(m_ibetaOffset<betaGrid.size());

  //Allow only loopmax sample attempts, to make sure we detect if code gets too
  //inefficient. However, make sure users can override this if needed.
  static const unsigned s_loopmax = []() -> unsigned
                                    {
                                      auto envstr = getenv("NCRYSTAL_SABSAMPLE_LOOPMAX");
                                      return envstr ? str2int(envstr) : 100;
                                    }();
  unsigned iloopmax(s_loopmax+1);
  while (--iloopmax) {
    double beta;
    unsigned ibetaSampled;
    std::tie(beta,ibetaSampled) = m_betaSampler.sampleWithIndex( rng );

    nc_assert( !ncisnan(beta) );
    nc_assert( beta <= betaGrid.back() );
    nc_assert( ibetaSampled < m_betaSampler.getXVals().size() );

    if (beta <= ncmax(-ekin_div_kT,betaGrid.front()))
      continue;//reject
    double alphal(-1.0), bl;
    std::size_t ibeta;

    double rand_percentile = rng.generate();

    if ( beta <= m_betaSampler.getXVals().at(1) ) {
      //Special case, beta is before the first actual betaGrid point used to
      //construct this instance of SABSamplerAtE_Alg1. Simply pick uniformly in
      //allowed region (nb: we could of course instead query the actual S-values
      //and interpolate among them).
      bl = m_betaSampler.getXVals().front();
      auto alimits_bl = getAlphaLimits( ekin_div_kT, bl );
      double alphal_low = ncclamp(alimits_bl.first,alphaGrid.front(),alphaGrid.back());
      double alphal_up = ncclamp(alimits_bl.second,alphaGrid.front(),alphaGrid.back());
      alphal = alphal_low + rng.generate()*(alphal_up-alphal_low);
      ibeta = m_ibetaOffset;
      //Double-check that chosen ibeta value gives correct bh=betaGrid.at(ibeta)
      //value, i.e. m_betaSampler.getXVals().at(1):
      nc_assert(floateq(m_betaSampler.getXVals().at(1),betaGrid.at(ibeta)));
    } else {
      ibeta = m_ibetaOffset + ibetaSampled;
      nc_assert( ibeta>0 );
      bl = betaGrid.at(ibeta-1);
      alphal = sampleAlpha(ibeta-1, rand_percentile);
    }

    nc_assert( valueInInterval(betaGrid.at(ibeta-1),betaGrid.at(ibeta),beta) );
    nc_assert( alphal>=0.0 );
    nc_assert( ibeta < betaGrid.size() );

    //Sample alphas (with same "random percentile" at two neighbouring beta grid
    //points and combine with linear interpolation, as in line 9 in Algorithm 1
    //of the paper:
    double bh = betaGrid.at(ibeta);
    double alphah = sampleAlpha(ibeta, rand_percentile);

    nc_assert(bh-bl>0.0);
    nc_assert(alphal>=0.0);
    nc_assert(alphah>=0.0);
    double alpha = alphal + (alphah-alphal) * (beta-bl)/(bh-bl);

    //Check if we can accept this:
    auto alimits = getAlphaLimits( ekin_div_kT, beta );
    if ( valueInInterval( alimits.first, alimits.second, alpha ) )
      return { alpha, beta };
  }
  NCRYSTAL_THROW2(CalcError,"Rejection method failed to sample kinematically valid (alpha,beta) point after "
                  <<s_loopmax<<" attempts. Perhaps energy grid is too sparse?"
                  " As a workaround it is possible to increase the allowed number of sampling attempts"
                  " by setting the NCRYSTAL_SABSAMPLE_LOOPMAX variable to a higher number"
                  " (but please consider reporting the issue to the NCrystal developers nonetheless).");
}

double NC::SAB::SABSamplerAtE_Alg1::sampleBeta(RandomBase& rng) const
{
  return m_betaSampler.sample(rng);
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
                                              vectAt(cd->alphaGrid(),info.pt_front.alpha_idx), span_at(sab,info.pt_front.alpha_idx),
                                              percentile2,
                                              info.pt_front.logsval, span_at(logsab,info.pt_front.alpha_idx) );
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
    return SABUtils::sampleLogLinDist_fast( vectAt(cd->alphaGrid(),a0), span_at(sab,a0),
                                            vectAt(cd->alphaGrid(),a1), span_at(sab,a1),
                                            rand_rescaled,
                                            span_at(logsab,a0), span_at(logsab,a1) );
  } else {
    //Sample back tail
    nc_assert( 1.0 - info.prob_notback > 0.0 );
    double percentile2 = clampRandNum ( ( rand_percentile - info.prob_notback ) / ( 1.0 - info.prob_notback ) );
    return SABUtils::sampleLogLinDist_fast( vectAt(cd->alphaGrid(),info.pt_back.alpha_idx), span_at(sab,info.pt_back.alpha_idx),
                                            info.pt_back.alpha, info.pt_back.sval,
                                            percentile2,
                                            span_at(logsab,info.pt_back.alpha_idx), info.pt_back.logsval );

  }

}
