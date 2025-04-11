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

#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {
      std::size_t detail_trimZeroEdgesFromKernel(ScatKnlData& input)
      {
        //trims edges but does NOT validate in this internal function
        const auto nalpha = input.alphaGrid.size();
        const auto nbeta = input.betaGrid.size();

        auto betaRowIsZero = [&input,nalpha] ( std::size_t ibeta )
        {
          auto idx0 = ibeta * nalpha;
          auto idxE = idx0 + nalpha;
          for ( auto idx = idx0; idx != idxE; ++idx )
            if ( vectAt( input.sab, idx ) )
              return false;
          return true;
        };
        auto alphaColIsZero = [&input,nalpha] ( std::size_t ialpha )
        {
          auto idx0 = ialpha;
          auto idxE = input.sab.size();
          for ( auto idx = idx0; idx < idxE; idx += nalpha )
            if ( vectAt( input.sab, idx ) )
              return false;
          return true;
        };

        std::size_t nTrimBetaUpper(0);
        std::size_t nTrimBetaLower(0);
        std::size_t nTrimAlphaUpper(0);

        for ( auto i : ncrange( nalpha ) ) {
          const auto ialpha = nalpha - i - 1;
          if ( vectAt(input.alphaGrid,ialpha) > 0.0 && alphaColIsZero( ialpha ) ) {
            ++nTrimAlphaUpper;
          } else {
            break;
          }
        }

        for ( auto i : ncrange( nbeta ) ) {
          const auto ibeta = nbeta - i - 1;
          if ( vectAt(input.betaGrid,ibeta) > 0.0 && betaRowIsZero( ibeta ) ) {
            ++nTrimBetaUpper;
          } else {
            break;
          }
        }

        for ( auto ibeta : ncrange( nbeta ) ) {
          if ( vectAt(input.betaGrid,ibeta) < 0.0 && betaRowIsZero( ibeta ) ) {
            ++nTrimBetaLower;
          } else {
            break;
          }
        }

        //Abort if too agressive (e.g. S=0 everywhere)
        if (  nTrimAlphaUpper >= nalpha ) {
          nTrimAlphaUpper = 0;
          nTrimBetaUpper = 0;
          nTrimBetaLower = 0;
        }

        const auto ntrimtot = nTrimBetaUpper + nTrimBetaLower + nTrimAlphaUpper;

        if ( ntrimtot == 0 )
          return 0;//nothing to trim

        //Trim:
        VectD new_sab;
        new_sab.reserve( ( nalpha - nTrimAlphaUpper ) * ( nbeta - (nTrimBetaUpper + nTrimBetaLower) ) );
        const auto nbetalim = nbeta - nTrimBetaUpper;
        for ( auto ibeta : ncrange( nbeta ) ) {
          if ( ! ( ibeta >= nTrimBetaLower && ibeta < nbetalim ) )
            continue;
          const auto offset = nalpha * ibeta;
          for ( auto ialpha : ncrange( nalpha-nTrimAlphaUpper ) )
            new_sab.push_back( vectAt( input.sab, ialpha + offset ) );
        }
        nc_assert_always( new_sab.size() == ( nalpha - nTrimAlphaUpper ) * ( nbeta - (nTrimBetaUpper + nTrimBetaLower) ) );
        std::swap( input.sab, new_sab );
        if ( nTrimAlphaUpper ) {
          VectD new_alpha( input.alphaGrid.begin(), std::next(input.alphaGrid.begin(),nalpha-nTrimAlphaUpper) );
          std::swap(input.alphaGrid,new_alpha);
        }
        if ( nTrimBetaUpper || nTrimBetaLower ) {
          VectD new_beta( std::next(input.betaGrid.begin(),nTrimBetaLower), std::next(input.betaGrid.begin(),nbeta-nTrimBetaUpper) );
          std::swap(input.betaGrid,new_beta);
        }
        nc_assert_always ( input.sab.size() == input.alphaGrid.size() * input.betaGrid.size() );
        return ntrimtot;
      }
    }
  }
}


NC::ScatKnlData NC::SABUtils::trimZeroEdgesFromKernel(ScatKnlData&& input)
{
  validateScatKnlData(input);
  NC::ScatKnlData res{ std::move(input) };
  detail_trimZeroEdgesFromKernel(res);
  return res;
}

NC::SABData NC::SABUtils::transformKernelToStdFormat( NC::ScatKnlData&& input_orig )
{
  validateScatKnlData(input_orig);

  //////////////////////
  // Trim null edges: //
  //////////////////////

  NC::ScatKnlData input{ std::move(input_orig) };
  auto ntrimmed = detail_trimZeroEdgesFromKernel(input);
  if ( ntrimmed )
    NCRYSTAL_WARN("Discarding "<<ntrimmed<<" edges of provided kernel"
                  " data due to missing S values.");

  ////////////////////////////////////////////////////
  // First convert ScatKnlData object to type "SAB" //
  ////////////////////////////////////////////////////

  if ( input.knltype == ScatKnlData::KnlType::SCALED_SYM_SAB ) {
    //Remove symmetry and actually specify full table :
    VectD complete_betagrid, complete_sab;
    expandBetaAndSABToAllBetas( input.betaGrid, input.alphaGrid, input.sab, complete_betagrid, complete_sab );
    std::swap( input.betaGrid, complete_betagrid );
    std::swap( input.sab, complete_sab );
    input.knltype = ScatKnlData::KnlType::SCALED_SAB;
  }

  if ( input.knltype == ScatKnlData::KnlType::SCALED_SAB ) {

    //Unscale, i.e. calculate S=S_scaled*exp(-beta/2):

    for ( auto&& beta : enumerate(input.betaGrid) ) {
      auto slice = sliceSABAtBetaIdx( input.sab, input.alphaGrid.size(), beta.idx );
      const double exparg = -0.5*beta.val;
      if ( exparg < 700.0 ) {
        //Straightforward, when beta is not below -1400, we can safely and
        //accurately calculate exp(-beta/2):
        const double expfact = std::exp( exparg );
        for ( double& sabEntry : slice )
          sabEntry *= expfact;
      } else {
        //The factor of exp(-beta/2) overflows at double precision. Fortunately,
        //this normally happens when S_scaled is itself a very small number,
        //thus balancing out the extreme factor (otherwise the physics would
        //probably not make much sense anyway). The numerically safe way to
        //calculate S is in this case to carry out the cancellation in
        //log-space, i.e. by using the more expensive rewrite:
        //
        // S = S_scaled*exp(-beta/2) = exp(-beta/2+log(S_scaled))
        //
        //Of course, when S_scaled=0, S is trivially also 0.
        //
        //Since log(eps), where eps is the smallest positive floating point
        //number in double precision, is roughly -744, it follows that any
        //scaled sab grid specified must have strictly S_scaled=0 for all points
        //at beta less than approximately -2900. Accurate specification of
        //S(alpha,beta) at large negative beta-values thus requires unscaled
        //S(alpha,beta) values to be specified directly (assuming of course that
        //whoever provided the table were able to circumvent numerical issues in
        //their calculations).
        for ( double& sabEntry : slice ) {
          if (sabEntry==0.0)
            continue;//0 is still 0 when multiplied with exp(-beta/2)
          double expargcombined = exparg + std::log(sabEntry);
          if ( expargcombined < 700.0 ) {
            sabEntry = std::exp(expargcombined);
          } else {
            NCRYSTAL_THROW2(BadInput,"Problems unscaling of S(alpha,beta), at point where S_scaled="
                            <<sabEntry<<" and beta="<<beta.val<<", since it requires evaluation"
                            " of exp("<<expargcombined<<") which is infinity at double precision."
                            " Most likely this indicates a problem with the input data.");
          }
        }
      }
    }//beta idx
    input.knltype = ScatKnlData::KnlType::SAB;
  }

  if ( input.knltype == ScatKnlData::KnlType::SQW ) {
    //Convert  S(q,w) -> S(alpha,beta).
    //Postponed for later release.
    NCRYSTAL_THROW2(LogicError,"Support for kernels in S(q,w) format is planned, but not yet implemented.");
#if 0
    //Leftover notes from previous attempt:

    //Still not working/verified (and we also need an SQW_SCALED version...):

    //First thing to note is that w=omega in the input is actually W=hbar*omega
    //(the unit is eV). Because of that, the conversion is (the beta and S
    //conversion would have had a factor of hbar if omega was actually omega
    //rather than hbar*omega):
    //
    // alpha = hbar^2/2*M*kT * qval^2
    // beta = -1/kT * W
    // The S-values themselves needs to be scaled with 1/kT.

    //other miscellaneous leftover notes:
    //        #beta = -B*w (if "w" is actually hbar*w, then B=1/kY)
    //        #alpha = A*q^2
    //        #dw = dbeta/B
    //        #int(S(q,w),dw) = (1/B) * int(S(alpha,beta)dbeta)
    //        #SO: If I get 1 for int(S(alpha,beta)dbeta), it must mean that S_ab = S_qw * B

    const double kT = constant_boltzmann*input.temperature;//unit is eV
    const double invkT = 1.0/kT;
    const double M = input.elementMassAMU*constant_dalton2eVc2;//Unit is eV/(Aa/s)^2 = eV*s^2/Aa^2
    const double hbar = constant_hbar;//eV*s
    const double c_qsq2a = 1e40*hbar*hbar/ (2.0 * M * kT);//Unit is Aa^2 [NOTICE: !!!!!!! 1e40 here is just an unvalidated try!]
    for (auto& e : input.alphaGrid)
      e = c_qsq2a * (e * e);
    for (auto& e : input.betaGrid)
      e *= -invkT;
    //reverse order, since conversion factor was negative.:
    std::reverse(std::begin(input.betaGrid), std::end(input.betaGrid));
    for (auto& e : input.sab)
      e *= invkT;
#endif
    input.knltype = ScatKnlData::KnlType::SAB;
  }

  nc_assert_always( input.knltype == ScatKnlData::KnlType::SAB );

  ///////////////////////////////////////////////////////////////////////////////
  // Thicken beta grid to reduce numerical artifacts in subsequent processing. //
  // This is obviously just a workaround and not the correct way to solve this //
  // (TODO: Remove once processing handles sparse grids better)!               //
  ///////////////////////////////////////////////////////////////////////////////

  int nminbeta_int = ncgetenv_int("SAB_BETATHICKENING_MINNBETA",500);
  nc_assert_always( nminbeta_int >=0 && nminbeta_int < 20000 );
  const unsigned nminbeta = static_cast<unsigned>(nminbeta_int);

  if ( !input.betaGridOptimised && input.betaGrid.size() < nminbeta ) {
    const unsigned nextra = nminbeta / input.betaGrid.size();//Number of extra beta points
                                                             //to insert between each point
    nc_assert_always( nextra >= 1 );
    auto newNBeta = 1 + ( input.betaGrid.size() -1 ) * ( 1 + nextra );
    VectD newb;
    newb.reserve( newNBeta );
    VectD newS;
    const std::size_t nalpha = input.alphaGrid.size();
    newS.reserve( newNBeta * nalpha );
    auto itB = input.betaGrid.begin();
    auto itBLast = std::prev(input.betaGrid.end());
    std::size_t ibeta(0);
    for ( ; itB != itBLast; ++itB, ++ibeta ) {
      auto itBNext = std::next(itB);
      auto alphaSlice = sliceSABAtBetaIdx_const( input.sab, nalpha, ibeta );
      auto alphaSliceNext = sliceSABAtBetaIdx_const( input.sab, nalpha, ibeta+1 );
      //First copy over existing row:
      newb.push_back(*itB);
      std::copy(alphaSlice.begin(), alphaSlice.end(), std::back_inserter(newS));
      //Then insert the new rows, one by one:
      const double dBeta = (*itBNext - *itB)/(nextra+1);
      for ( auto iextra : ncrange(nextra) ) {
        double beta = *itB + (iextra+1)*dBeta;
        newb.push_back(beta);
        for ( auto ialpha : ncrange(nalpha) ) {
          double S = alphaSlice.at(ialpha);
          double SNext = alphaSliceNext.at(ialpha);
          double kkk = ( beta - *itB) / (*itBNext - *itB);
          //newS.push_back( S + (SNext - S) * kkk  );
          newS.push_back( S * ( 1.0 - kkk) + kkk * SNext );
        }
      }
    }
    //Add last slide as well:
    auto alphaSliceLast = sliceSABAtBetaIdx_const( input.sab, nalpha, ibeta );
    std::copy(alphaSliceLast.begin(), alphaSliceLast.end(), std::back_inserter(newS));
    newb.push_back(*itBLast);
    nc_assert_always( newS.size() == newNBeta * nalpha );
    nc_assert_always( newb.size() == newNBeta );

    //Apply:
    std::swap(newS,input.sab);
    std::swap(newb,input.betaGrid);
  }

  ///////////////////////////////////////////
  // Transfer to SABData object and return //
  ///////////////////////////////////////////

  SABData out{ std::move(input.alphaGrid),
               std::move(input.betaGrid),
               std::move(input.sab),
               input.temperature,
               input.boundXS,
               input.elementMassAMU,
               input.suggestedEmax };

#ifndef NDEBUG
  validateScatKnlData(out);
#endif

  return out;
}

void NC::SABUtils::expandBetaAndSABToAllBetas( NC::Span<const double> halfbetagrid,
                                               NC::Span<const double> alphagrid,
                                               NC::Span<const double> sab_for_halfbetagrid,
                                               VectD& complete_betagrid,
                                               VectD& complete_sab )
{
  //Prepare and check:
  complete_betagrid.clear(  );
  complete_sab.clear();
  const std::size_t nalpha = alphagrid.size();
  const std::size_t nbeta_old = halfbetagrid.size();
  const std::size_t nbeta_positive = nbeta_old - 1;
  const std::size_t nbeta_new = nbeta_positive*2 + 1;
  nc_assert_always(!halfbetagrid.empty());
  nc_assert_always(halfbetagrid.front()==0.0);
  nc_assert_always( nbeta_old * nalpha == static_cast<std::size_t>(sab_for_halfbetagrid.size()) );

  //Step 1. Create complete beta-grid:
  // -> the negative values and zero:
  complete_betagrid.reserve( nbeta_new );
  for (auto it = halfbetagrid.rbegin(); it != halfbetagrid.rend(); ++it)
    complete_betagrid.emplace_back( - (*it) );

  // -> avoid signed negative zero, for aesthetic reasons :
  nc_assert(complete_betagrid.back()==0.0);
  complete_betagrid.back() = 0.0;
  // -> add the positive values (the +1 below avoids adding the zero again):
  for (auto e: Span<const double>(halfbetagrid.begin()+1,halfbetagrid.end()) )
    complete_betagrid.emplace_back(e);
  nc_assert_always( complete_betagrid.size() == nbeta_new );

  //Step 2. Expand the sab kernel to complete beta range

  //Step 2.1: fill out the parts for beta<0 using S(alpha,-beta) := S(alpha,beta):
  complete_sab.reserve( nbeta_new * alphagrid.size() );
  complete_sab.resize( alphagrid.size() * nbeta_positive, 0. );
  const auto srcbegin = sab_for_halfbetagrid.begin();
  auto targetrowbegin = complete_sab.begin();

  for ( std::size_t i = 0; i < nbeta_positive; ++i )  {
    auto srcrowbegin = srcbegin + (nbeta_positive-i)*nalpha;
    nc_assert( targetrowbegin + (nalpha-1) < complete_sab.end() );
    std::copy( srcrowbegin, srcrowbegin + nalpha, targetrowbegin );
    targetrowbegin += nalpha;
  }

  //Step 2.2: copy over the original values for beta>=0:
  complete_sab.insert( complete_sab.end(),
                       sab_for_halfbetagrid.begin(),
                       sab_for_halfbetagrid.end() );

  //Validate:
  nc_assert_always( complete_betagrid.size() == nbeta_new );
  nc_assert_always( complete_sab.size() == nalpha * nbeta_new );

}

void NC::SABUtils::activeGridCells( const NC::SABData& data,
                                    double ekin_div_kT,
                                    std::vector<std::pair<std::uint16_t,std::uint16_t>>& out_alpharanges,
                                    std::size_t& ibeta_low  )
{
  out_alpharanges.clear();
  ibeta_low = data.betaGrid().size();

  //We want to find the active 2D cells. First we find the 1D active ranges:
  std::vector<std::pair<std::uint16_t,std::uint16_t>> alpha1dranges;
  std::size_t ibeta_low1d;
  activeGridRanges( data, ekin_div_kT, alpha1dranges, ibeta_low1d );

  //Special case, nothing active:
  if ( alpha1dranges.empty() )
    return;

  nc_assert( ibeta_low1d < data.betaGrid().size() );

  ibeta_low = ibeta_low1d;
  std::size_t nexpected_cell_ranges = alpha1dranges.size();
  nc_assert( nexpected_cell_ranges > 0 );
  if ( ibeta_low1d > 0 ) {
    //Cells start 1 bin below ranges in this case, the first cell has only
    //ranges on the right.
    --ibeta_low;
    out_alpharanges.reserve(nexpected_cell_ranges);
    out_alpharanges.push_back(alpha1dranges.front());
    //Special case: Cells including beta=0.0 must always extend to lowest alpha
    //cell (see also below where we treat the same special case).
    nc_assert(ibeta_low+1<data.betaGrid().size());
    if ( valueInInterval( vectAt(data.betaGrid(),ibeta_low),
                          vectAt(data.betaGrid(),ibeta_low+1),
                          0.0 ) )
      out_alpharanges.back().first = 0;
  } else {
    --nexpected_cell_ranges;
    out_alpharanges.reserve(nexpected_cell_ranges);
  }

  nc_assert(alpha1dranges.size()>0);
  const std::size_t nnn = alpha1dranges.size()-1;

  auto itBeta = std::next(data.betaGrid().begin(),ibeta_low1d);

  const uint16_t nalpha_uint16 = data.alphaGrid().size();

  for ( std::size_t i = 0; i < nnn; ++i, ++itBeta ) {
    const auto& r0 = vectAt(alpha1dranges,i);
    const auto& r1 = vectAt(alpha1dranges,i+1);
    const bool r0Empty = r0.first >= nalpha_uint16;
    const bool r1Empty = r1.first >= nalpha_uint16;
    if ( r0Empty ) {
      out_alpharanges.emplace_back(r1);//even if r1 is also empty, this is correct (both empty => empty cells)
    } else if ( r1Empty ) {
      out_alpharanges.emplace_back(r0);
    } else {
      out_alpharanges.emplace_back(std::min<uint16_t>(r0.first,r1.first),
                                   std::max<uint16_t>(r0.second,r1.second));
    }
    //Special case: Cell including beta=0.0 must always extend to lowest alpha
    //cell (the code above might have missed exactly this, because the slope of
    //alpha-(beta) changes sign at beta=0.
    nc_assert(std::next(itBeta)<data.betaGrid().end());
    if ( valueInInterval( *itBeta, *std::next(itBeta), 0.0 ) )
      out_alpharanges.back().first = 0;
  }

  nc_assert( out_alpharanges.size() == nexpected_cell_ranges );
  return;
}

void NC::SABUtils::activeGridRanges( const NC::SABData& data,
                                     double ekin_div_kT,
                                     std::vector<std::pair<std::uint16_t,std::uint16_t>>& out_alpharanges,
                                     std::size_t& ibeta_low  )
{
  const auto& alphaGrid = data.alphaGrid();
  nc_assert(alphaGrid.size()>1);
  nc_assert(nc_is_grid(alphaGrid));
  nc_assert_always(alphaGrid.size()<std::numeric_limits<std::uint16_t>::max());

  ibeta_low = 0;
  out_alpharanges.clear();
  const double agrid_front = alphaGrid.front();
  const double agrid_back = alphaGrid.back();
  //At each beta grid point, we search for grid points marking the front and
  //back of the active grid range. For (large) gains of efficiency, we start off
  //the search for these alpha grid points at a given beta grid point, based on
  //their values at the previous beta grid point. This works well because the
  //kinematic boundaries are smooth curves.
  const auto itBegin = alphaGrid.begin();
  const auto itLast = std::prev(alphaGrid.end());
  auto itLow = alphaGrid.begin();
  auto itUpp = std::prev(alphaGrid.end());

  for (auto&& beta : enumerate(data.betaGrid())) {
    double alow(-1.0),aupp(-2.0);
    if ( beta.val > -ekin_div_kT ) {
      auto alims = getAlphaLimits( ekin_div_kT, beta.val );
      alow = alims.first;
      aupp = alims.second;
    }
    if ( agrid_back <= alow || agrid_front >= aupp || aupp < alow ) {
      //No kinematically accessible alpha grid ranges at this beta point (or
      //energy is so ultra low that numerical imprecision led to aupp=alow) .
      if ( out_alpharanges.empty() ) {
        //Still didn't encounter any beta point with accessible alpha ranges, so
        //simply increment ibeta_low.
        ibeta_low = beta.idx + 1;
      } else {
        //Already encountered at least one beta point with accessible alpha
        //range, so must insert an empty alpha range at this point:
        nc_assert( alphaGrid.size()
                   < static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()));
        auto na = static_cast<std::uint16_t>(alphaGrid.size());
        out_alpharanges.emplace_back(na,na);
      }
      continue;
    }

    //Ok, there is a accessible alpha grid range at this beta point, time to
    //find it.

    //Move itLow down or up as needed:
    while ( *itLow > alow && itLow > itBegin )
      --itLow;
    while ( itLow < itLast && *std::next(itLow) <= alow )
      ++itLow;
    //Move itUpp down or up as needed:
    if ( itUpp < itLow )
      itUpp = itLow;
    while ( *itUpp < aupp && itUpp < itLast )
      ++itUpp;
    while ( itUpp > itBegin && *std::prev(itUpp) >= aupp )
      --itUpp;

    //Register result. Notice that in case aupp==alow (for ultra ultra cold
    //neutrons, due to numerical issues) and that value is exactly on an alpha
    //grid point, we could have an empty range (itLow=itUpp). That is OK, we
    //simply register an empty range:
    nc_assert( itUpp < alphaGrid.end() );
    nc_assert( (aupp==alow) ? (itLow <= itUpp) : (itLow < itUpp) );
    out_alpharanges.emplace_back( static_cast<std::uint16_t>(std::distance(itBegin,itLow)),
                                  static_cast<std::uint16_t>(std::distance(itBegin,itUpp)) );

  }

}

NC::SABUtils::TailedBreakdown NC::SABUtils::createTailedBreakdown( const NC::Span<const double>& alphaGrid,
                                                                   const NC::Span<const double>& sab,
                                                                   const NC::Span<const double>& logsab,
                                                                   const NC::Span<const double>& alphaIntegrals_cumul,
                                                                   double alpha_low, double alpha_upp,
                                                                   const unsigned aidx_low, const unsigned aidx_upp )
{
  nc_assert( alpha_low <= alpha_upp );
  nc_assert( aidx_low <= aidx_upp );
  nc_assert( aidx_upp < alphaGrid.size() );

  //Constrain ranges to grid (xs outside is modelled as 0):

  alpha_low = ncclamp(alpha_low,alphaGrid.front(),alphaGrid.back());
  alpha_upp = ncclamp(alpha_upp,alphaGrid.front(),alphaGrid.back());

  TailedBreakdown tb;
  if ( aidx_low == aidx_upp || alpha_low == alpha_upp ) {
    //vanishing!
    return tb;
  }

  nc_assert( aidx_upp+1==alphaGrid.size() || alpha_upp <= alphaGrid[aidx_upp] );
  nc_assert( aidx_low==0 || alpha_low >= alphaGrid[aidx_low] );
  nc_assert( aidx_low+1==alphaGrid.size() || alpha_low < alphaGrid[aidx_low+1] );
  nc_assert( aidx_upp==0 || alpha_upp > alphaGrid[aidx_upp-1] );

  auto interpSVal = [&alphaGrid,&sab,&logsab](const std::size_t alphaidx_lowedge, const double alpha)
                    {
                      nc_assert( alphaidx_lowedge + 1 < (unsigned)alphaGrid.size() );
                      const double alpha0(alphaGrid[alphaidx_lowedge]);
                      const double alpha1(alphaGrid[alphaidx_lowedge+1]);
                      nc_assert( valueInInterval(alpha0,alpha1,alpha) );
                      auto i0(alphaidx_lowedge), i1(alphaidx_lowedge+1);
                      return interpolate_loglin_fallbacklinlin_fast(alpha0,sab[i0],alpha1,sab[i1],alpha,
                                                                    logsab[i0],logsab[i1]);
                    };
  auto setTailPoint = [&interpSVal](TailedBreakdown::TailPoint& tp,unsigned aidx,double alpha)
                      {
                        tp.alpha = alpha;
                        tp.sval = interpSVal(aidx,alpha);
                        tp.logsval = tp.sval > 0.0 ? std::log(tp.sval) : -kInfinity;
                      };

  //Enough setting up, time to analyse the tails and how they fall wrt the grid:
  if ( aidx_low + 1 == aidx_upp ) {
    //Special "narrow" case. Only a single alpha bin is touched by the
    //non-expanded region.  Integrate entire alpha-range in one go:
    tb.narrow = true;
    setTailPoint(tb.front,aidx_low,alpha_low);
    setTailPoint(tb.back,aidx_low,alpha_upp);
    tb.xs_front = integrateAlphaInterval_fast( tb.front.alpha, tb.front.sval,
                                               tb.back.alpha, tb.back.sval,
                                               tb.front.logsval, tb.back.logsval );
    return tb;
  }

  //Non-narrow case, treat "front" and "back" tails separately.
  tb.imiddle_low = aidx_low;
  tb.imiddle_upp = aidx_upp;

  //Front (not there if alpha_low is outside the grid range):
  if ( alpha_low >= alphaGrid[aidx_low] ) {
    nc_assert( alpha_low <= alphaGrid[aidx_low + 1] );
    setTailPoint(tb.front,aidx_low,alpha_low);
    tb.xs_front = integrateAlphaInterval_fast( tb.front.alpha, tb.front.sval,
                                               alphaGrid[aidx_low+1], sab[aidx_low+1],
                                               tb.front.logsval, logsab[aidx_low+1] );
    ++tb.imiddle_low;
  }
  //Back (not there if alpha_upp is outside the grid range):
  if ( alpha_upp <= alphaGrid[aidx_upp] ) {
    nc_assert( aidx_upp != 0 );
    nc_assert( alpha_upp >= alphaGrid[aidx_upp-1] );
    setTailPoint(tb.back,aidx_upp-1,alpha_upp);
    tb.xs_back = integrateAlphaInterval_fast( alphaGrid[aidx_upp-1], sab[aidx_upp-1],
                                              tb.back.alpha, tb.back.sval,
                                              logsab[aidx_upp-1], tb.back.logsval );
    --tb.imiddle_upp;
  }
  tb.xs_middle = ( tb.imiddle_upp > tb.imiddle_low ?
                   alphaIntegrals_cumul[tb.imiddle_upp] - alphaIntegrals_cumul[tb.imiddle_low]
                   : 0.0 );
  return tb;
}
