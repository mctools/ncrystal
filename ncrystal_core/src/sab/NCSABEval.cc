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

#include "NCrystal/internal/sab/NCSABEval.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;
namespace NCSD = NCrystal::SABUtils::detail_sce;

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    namespace detail_sce {

      template< class TInt >
      inline double extract_logS( const SCE_SLogVals& data, TInt i ) ncnoexceptndebug { nc_assert(i<4); return data.logS[i]; }

      template< class TInt >
      inline double extract_logS( const SCE_Empty&, TInt ) ncnoexceptndebug { return -1.0; }//never called

      template <class TData>
      inline double evalAtAlpha0 ( const TData& data, double beta )
      {
        nc_assert( valueInInterval( data.beta, beta ) );
        return interpolate<TData::betaInterpType>(data.beta.first, data.S[0], data.beta.second, data.S[2], beta );
      }

      template <class TData>
      inline double evalAtAlpha1 ( const TData& data, double beta )
      {
        nc_assert( valueInInterval( data.beta, beta ) );
        return interpolate<TData::betaInterpType>(data.beta.first, data.S[1], data.beta.second, data.S[3], beta );
      }

      template <class TData>
      inline double evalAtBeta0 ( const TData& data, double alpha )
      {
        nc_assert( valueInInterval( data.alpha, alpha ) );
        nc_assert( TData::interpOrder == SABInterpolationOrder::ALPHA_FIRST );
        return ( ( TData::alphaInterpType == InterpolationScheme::LOGLIN && TData::hasLogValues )
                 ? interpolate_loglin_fallbacklinlin_fast(data.alpha.first, data.S[0], data.alpha.second, data.S[1], alpha, extract_logS(data,0), extract_logS(data,1) )
                 : interpolate<TData::alphaInterpType>(data.alpha.first, data.S[0], data.alpha.second, data.S[1], alpha ) );
      }

      template <class TData>
      inline double evalAtBeta1 ( const TData& data, double alpha )
      {
        nc_assert( valueInInterval( data.alpha, alpha ) );
        nc_assert( TData::interpOrder == SABInterpolationOrder::ALPHA_FIRST );
        return ( ( TData::alphaInterpType == InterpolationScheme::LOGLIN && TData::hasLogValues )
                 ? interpolate_loglin_fallbacklinlin_fast(data.alpha.first, data.S[2], data.alpha.second, data.S[3], alpha, extract_logS(data,2), extract_logS(data,3) )
                 : interpolate<TData::alphaInterpType>(data.alpha.first, data.S[2], data.alpha.second, data.S[3], alpha ) );
      }

      template <class Func>
      inline double SCE_integrateRombergFlex(const Func& f, double a, double b, Optional<double> splitpoint = NullOpt)
      {
        //TODO: tune next values further? Also, consider swapping romberg17<->romberg33 elsewhere in this file.
        constexpr double prec = 1e-10;
        constexpr unsigned minlvl = 5;
        constexpr unsigned maxlvl = 8;
        if ( !splitpoint.has_value() || splitpoint.value() <=a || splitpoint.value() >= b ) {
          return integrateRombergFlex(f,a,b,prec,minlvl,maxlvl);
        } else {
          nc_assert(splitpoint.value()>=a&&splitpoint.value()<=b);
          return ( integrateRombergFlex(f,a,splitpoint.value(),prec,minlvl,maxlvl)
                   + integrateRombergFlex(f,splitpoint.value(),b,prec,minlvl,maxlvl) );

        }
      }

      template <class TData>
      void intCrossedKB( StableSum& sum, const TData& data, double ekin_div_kT )
      {
        //Crosses kinematic borders (KB) and has already been trimmed (with
        //trimToKB). We need to:
        // 1) Find the complete integral over the entire region
        // 2) Subtract the integral above beta+(alpha)
        // 3) Subtract the integral below beta-(alpha)

        //////////////////////////////////////////////////////////
        // 1) Find the complete integral over the entire region //
        //////////////////////////////////////////////////////////
        data.integral(sum);

        /////////////////////////////////////////////////
        // 2) Subtract the integral above beta+(alpha) //
        /////////////////////////////////////////////////

        const double betaplus_at_a0 = getBetaPlus(ekin_div_kT,data.alpha.first);
        if ( data.beta.second > betaplus_at_a0 ) {
          nc_assert(valueInInterval(data.beta,ncmax(betaplus_at_a0,data.beta.first)));
          nc_assert(valueInInterval(data.beta,data.beta.second));

          sum.add( - integrateRomberg17([ekin_div_kT,&data]( double beta )
          {
            nc_assert( valueInInterval(data.beta,beta) || beta<data.beta.second+0.001*(data.beta.second-data.beta.first) );
            beta = ncmin( data.beta.second, beta );//NCRomberg::evalFuncMany might evaluate slightly above upper limit
            const double alow = data.alpha.first;
            const double ahigh = ncmin(getAlphaMinus(ekin_div_kT,beta),data.alpha.second);
            const double dA = ahigh-alow;
            if ( !(dA>0.0) )
              return 0.0;
            const double Slow = evalAtAlpha0(data,beta);
            const double Shigh = data.eval(ahigh, beta );
            if ( TData::alphaInterpType == InterpolationScheme::LOGLIN )
              return integrateAlphaInterval( alow, Slow, ahigh, Shigh);
            else
              return (ahigh-alow)*0.5*(Slow+Shigh);
          },ncmax(betaplus_at_a0,data.beta.first),data.beta.second));
        };

        /////////////////////////////////////////////////
        // 3) Subtract the integral below beta-(alpha) //
        /////////////////////////////////////////////////


        //First the region at alpha<E/kT:
        const double betaminus_at_a0 = getBetaMinus(ekin_div_kT,data.alpha.first);
        if ( data.alpha.first < ekin_div_kT && data.beta.first < betaminus_at_a0 ) {
          Optional<double> splitpt;
          const double ahighmax = ncmin(data.alpha.second,ekin_div_kT);
#if 1
          const double b2 = ncmin(0.0,ncmin(data.beta.second,betaminus_at_a0));
          if ( ahighmax > 0.8*ekin_div_kT && data.beta.first<= -0.9*ekin_div_kT ) {
            splitpt = data.beta.first + 0.01*(b2-data.beta.first);
          }
          sum.add( - SCE_integrateRombergFlex([ekin_div_kT,ahighmax,&data]( double beta )
          {
            beta = ncmin( data.beta.second, beta );//NCRomberg::evalFuncMany might evaluate slightly above upper limit
            const double alow = data.alpha.first;
            const double ahigh = ncmin(getAlphaMinus(ekin_div_kT,beta),ahighmax);
            const double dA = ahigh-alow;
            if ( !(dA>0.0) )
              return 0.0;
            const double Slow = evalAtAlpha0(data,beta);
            const double Shigh = data.eval(ahigh, beta );
            if ( TData::alphaInterpType == InterpolationScheme::LOGLIN )
              return integrateAlphaInterval( alow, Slow, ahigh, Shigh);
            else
              return (ahigh-alow)*0.5*(Slow+Shigh);
          },
              data.beta.first,b2,splitpt));
#else
          sum.add( - integrateRomberg17([ekin_div_kT,ahighmax,&data]( double beta )
          {
            beta = ncmin( data.beta.second, beta );//NCRomberg::evalFuncMany might evaluate slightly above upper limit
            const double alow = data.alpha.first;
            const double ahigh = ncmin(getAlphaMinus(ekin_div_kT,beta),ahighmax);
            const double dA = ahigh-alow;
            if ( !(dA>0.0) )
              return 0.0;
            const double Slow = evalAtAlpha0(data,beta);
            const double Shigh = data.eval( ahigh, beta );
            if ( TData::alphaInterpType == InterpolationScheme::LOGLIN )
              return integrateAlphaInterval( alow, Slow, ahigh, Shigh);
            else
              return (ahigh-alow)*0.5*(Slow+Shigh);
          },data.beta.first,ncmin(0.0,ncmin(data.beta.second,betaminus_at_a0))));
#endif
        }
        //Finally the region at alpha>E/kT:
#if 0
        if ( data.alpha.second > ekin_div_kT ) {
          const double alowmin = ncmax(data.alpha.first,ekin_div_kT);
          sum.add( - SCE_integrateRombergFlex([ekin_div_kT,alowmin,&data]( double beta )
          {
            beta = ncmin( data.beta.second, beta );//NCRomberg::evalFuncMany might evaluate slightly above upper limit
            const double alow = ncmax(getAlphaPlus(ekin_div_kT,beta),alowmin);
            const double ahigh = data.alpha.second;
            const double dA = ahigh-alow;
            if ( !(dA>0.0) )
              return 0.0;
            const double Slow = data.eval( alow, beta );
            const double Shigh = evalAtAlpha1( data, beta );
            if ( TData::alphaInterpType == InterpolationScheme::LOGLIN )
              return integrateAlphaInterval( alow, Slow, ahigh, Shigh);
            else
              return (ahigh-alow)*0.5*(Slow+Shigh);
          },data.beta.first,ncclamp(getBetaMinus(ekin_div_kT,data.alpha.second),data.beta.first,data.beta.second)));
        }
#else
        if ( data.alpha.second > ekin_div_kT ) {
          auto gen_integrand = [ekin_div_kT,&data]( PairDD clamp_alpha )
          {
            const double ahigh = ncclamp( data.alpha.second, clamp_alpha );

            return [ekin_div_kT,ahigh,clamp_alpha,&data]( double beta )
            {
              const double alow = ncclamp( getAlphaPlus(ekin_div_kT,beta), clamp_alpha );
              const double dA = ahigh-alow;
              if ( !(dA>0.0) )
                return 0.0;
              const double Slow = data.eval( alow, beta );
              const double Shigh = evalAtAlpha1( data, beta );
              if ( TData::alphaInterpType == InterpolationScheme::LOGLIN )
                return integrateAlphaInterval( alow, Slow, ahigh, Shigh);
              else
                return (ahigh-alow)*0.5*(Slow+Shigh);
            };
          };
          const double alowmin = ncmax(data.alpha.first,ekin_div_kT);
          const double b0 = data.beta.first;
          const double b2 = ncclamp( getBetaMinus(ekin_div_kT,data.alpha.second), data.beta );

          if ( b0 < (-0.9)*ekin_div_kT ) {

            const double b1 = b0 + 0.01*(b2-b0);//ncmin(0.01,0.01*(b2-b0));
            const double ap_at_b1 = ncclamp(getAlphaPlus(ekin_div_kT,b1),{alowmin,data.alpha.second});
            sum.add( -SCE_integrateRombergFlex(gen_integrand({alowmin,data.alpha.second}),b0,b1) );
            sum.add( -SCE_integrateRombergFlex(gen_integrand({ap_at_b1,data.alpha.second}),b1,b2) );
          } else {
            sum.add( -integrateRomberg17(gen_integrand({alowmin,data.alpha.second}),b0,b2) );
          }
        }
#endif

      }
      template <class TData>
      Optional<TData> trimToKB( const TData& data_orig, double ekin_div_kT, BetaLimits_t blim_at_a0, BetaLimits_t blim_at_a1 )
      {
        //Todo: Just "inplace_splitAtBeta etc. and avoid this Optional + dataptr stuff!

        Optional<TData> res;
        const TData* data = &data_orig;
        auto doSplitBeta = [&res,&data](double beta, IgnorePart ip)
        {
          if ( beta <= data->beta.first ) {
            if ( ip == IgnorePart::Below )
              return;
            beta = std::nextafter ( data->beta.first, data->beta.second );
          } else if ( beta >= data->beta.second ) {
            if ( ip == IgnorePart::Above )
              return;
            beta = std::nextafter ( data->beta.second, data->beta.first );
          }
          res = data->splitAtBeta( ncclamp(beta,data->beta), ip );
          data = &res.value();
        };
        auto doSplitAlpha = [&res,&data](double alpha, IgnorePart ip)
        {
          if ( alpha <= data->alpha.first ) {
            if ( ip == IgnorePart::Below )
              return;
            alpha = std::nextafter ( data->alpha.first, data->alpha.second );
          } else if ( alpha >= data->alpha.second ) {
            if ( ip == IgnorePart::Above )
              return;
            alpha = std::nextafter ( data->alpha.second, data->alpha.first );
          }
          res = data->splitAtAlpha( ncclamp(alpha,data->alpha), ip );
          data = &res.value();
        };


        //Trim excess edge towards positive beta:
        if ( data->beta.second > blim_at_a1.second )
          doSplitBeta( blim_at_a1.second, IgnorePart::Above );

        //Trim excess edge towards negative beta:
        const double bminus_lower_reach = ( data->alpha.first >= ekin_div_kT
                                            ? blim_at_a0.first
                                            : ( data->alpha.second <= ekin_div_kT
                                                ? blim_at_a1.first
                                                : -ekin_div_kT ) );

        if ( data->beta.first < bminus_lower_reach )
          doSplitBeta( bminus_lower_reach, IgnorePart::Below );
        //Trim excess edge towards positive alpha:
        if ( data->beta.second < blim_at_a1.first ) {
          nc_assert(data->alpha.second>=ekin_div_kT);
          doSplitAlpha(getAlphaPlus(ekin_div_kT,data->beta.second),IgnorePart::Above);//NB: After calling doSplitAlpha we can no longer trust blim_at_a0/blim_at_a1 [SO DONT USE THEM!!!]
        }
        blim_at_a1 = { -1,-1 };//invalidate since we might have called doSplitAlpha(...,Above)!!!

        //Trim excess edge towards negative alpha:
        if ( data->beta.first > 0 && data->beta.first > blim_at_a0.second ) {
          //handle case above beta=0:
          doSplitAlpha(getAlphaMinus(ekin_div_kT,data->beta.first),IgnorePart::Below);
        } else if ( data->beta.second < 0 && data->alpha.first < ekin_div_kT && data->beta.second < blim_at_a0.first ) {
          //handle case below beta=0:
          doSplitAlpha(getAlphaMinus(ekin_div_kT,data->beta.second),IgnorePart::Below);
        }
        blim_at_a0 = { -1,-1 };//invalidate since we might have called doSplitAlpha(...,Below)!!!
        //Trim excess edge towards positive alpha:
        if ( data->alpha.second > ekin_div_kT && data->beta.second < getBetaMinus(ekin_div_kT,data->alpha.second) )
          doSplitAlpha(getAlphaPlus(ekin_div_kT,data->beta.second),IgnorePart::Above);

        return res;
      }
   }
  }
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline typename NCS::SABCellEval<AIT,IO>::SCE_Data NCS::SABCellEval<AIT,IO>::SCE_Data::splitAtBeta( double betaval, IgnorePart ip ) const
{
  nc_assert( valueInInterval( this->beta, betaval ) );
  SCE_Data copy = *this;
  if ( ip == IgnorePart::Above ) {
    copy.beta.second = betaval;
    copy.S[2] = evalAtAlpha0( *this, betaval );
    copy.S[3] = evalAtAlpha1( *this, betaval );
    copy.update_logS( 2, copy.S[2] );
    copy.update_logS( 3, copy.S[3] );
  } else {
    copy.beta.first = betaval;
    copy.S[0] = evalAtAlpha0( *this, betaval );
    copy.S[1] = evalAtAlpha1( *this, betaval );
    copy.update_logS( 0, copy.S[0] );
    copy.update_logS( 1, copy.S[1] );
  }
  copy.validate();
  return copy;
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline typename NCS::SABCellEval<AIT,IO>::SCE_Data NCS::SABCellEval<AIT,IO>::SCE_Data::splitAtAlpha( double alphaval, IgnorePart ip ) const
{
  nc_assert( valueInInterval( this->alpha, alphaval ) );
  SCE_Data copy = *this;
  if ( ip == IgnorePart::Above ) {
    copy.alpha.second = alphaval;
    copy.S[1] = evalAtBeta0( *this, alphaval );
    copy.S[3] = evalAtBeta1( *this, alphaval );
    copy.update_logS( 1, copy.S[1] );
    copy.update_logS( 3, copy.S[3] );
  } else {
    copy.alpha.first = alphaval;
    copy.S[0] = evalAtBeta0( *this, alphaval );
    copy.S[2] = evalAtBeta1( *this, alphaval );
    copy.update_logS( 0, copy.S[0] );
    copy.update_logS( 2, copy.S[2] );
  }
  copy.validate();
  return copy;
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline NC::PairDD NCS::SABCellEval<AIT,IO>::SCE_Data::sample( RNG& rng ) const
{
  const double a0 = this->alpha.first;
  const double b0 = this->beta.first;
  const double dA = this->alpha.second - a0;
  const double dB = this->beta.second - b0;
  const double Smaxval = this->sMax();
  while ( true ) {
    const double a = a0 + dA*rng();
    const double b = b0 + dB*rng();
    if ( this->eval(a,b) >= Smaxval * rng() )//"=" ok to avoid inf loop when Smax=0
      return { a, b };
  }
}

// template <class TData>
//  NCSD::sample( const TData& data, RNG& rng )
// {
// }

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline NCS::KBStatus NCS::SABCellEval<AIT,IO>::SCE_Data::kbStatus( double ekin_div_kT ) const
{
  const double ee = ekin_div_kT;
  const double b0 = this->beta.first;
  const double b1 = this->beta.second;
  const double a0 = this->alpha.first;
  const double a1 = this->alpha.second;
  //First check if we are completely outside:
  if ( b1 <= -ee )
    return KBStatus::FullyOutside;
  const double b0_m_a1 = b0 - a1;
  if ( b0_m_a1 >= 0.0 ) {
    //Fully above betaplus, if:
    //  b0 >= a1 + 2*sqrt(ee*a1)  <=>
    //  b0-a1 >= 2*sqrt(ee*a1)  <=> (since both sides non-negative)
    //  (b0-a1)^2 >= 4*ee*a1
    if ( ncsquare(b0_m_a1) >= 4.0*ee*a1 )
      return KBStatus::FullyOutside;
  } else if ( a0 >= b1 ) {
    if ( b1 <= -ee )
      return KBStatus::FullyOutside;
    //Fully below betaminus, if (ai is a1 if a1<e and a0 if a0>=e):
    //  b1 <= ai - 2*sqrt(ee*ai)  <=>
    //  ai-b1 >= 2*sqrt(ee*ai)  <=> (since both sides non-negative)
    //  (ai-b1)^2 >= 4*ee*ai
    if ( a1 <= ee ) {
      if ( ncsquare(a1-b1) >= 4*ee*a1 )
        return KBStatus::FullyOutside;
    } else if ( a0 >= ee ) {
      if ( ncsquare(a0-b1) >= 4*ee*a0 )
        return KBStatus::FullyOutside;
    } else {
      return KBStatus::Crossing;
    }
  }
  //Ok, we are not fully outside.
  // auto blim_at_a1 = getBetaLimits( ekin_div_kT, a1 );
  // if ( b0 >= blim_at_a1.second )
  //   return KBStatus::FullyOutside;
  // if ( a1 <= ee && b1 <= blim_at_a1.first )
  //   return KBStatus::FullyOutside;
  // if ( a0 >= ee && b1 <= blim_at_a0.first )
    // return KBStatus::FullyOutside;

  //Next check if we are completely inside:
  auto blim_at_a0 = getBetaLimits( ekin_div_kT, a0 );
  if ( b0 >= -ee && b1 <= blim_at_a0.second ) {
    bool crossBMinusLeftSide = ( a0 < ee && b0 < blim_at_a0.first );
    bool crossBMinusRightSide = ( a1 > ee && b0 < getBetaMinus(ee,a1) );
    if ( !crossBMinusLeftSide && !crossBMinusRightSide )
      return KBStatus::FullyInside;
  }
  return KBStatus::Crossing;
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline double NCS::SABCellEval<AIT,IO>::SCE_Data::sOverlayWKB( double ekin_div_kT ) const
{
  const double ee = ekin_div_kT;
  const double b0 = this->beta.first;
  const double b1 = this->beta.second;
  const double a0 = this->alpha.first;
  const double a1 = this->alpha.second;
  //First check if we are completely outside:
  if ( b1 <= -ee )
    return 0.0;
  auto blim_at_a1 = getBetaLimits( ekin_div_kT, a1 );
  if ( b0 >= blim_at_a1.second )
    return 0.0;
  if ( a1 <= ee && b1 <= blim_at_a1.first )
    return 0.0;
  auto blim_at_a0 = getBetaLimits( ekin_div_kT, a0 );
  if ( a0 >= ee && b1 <= blim_at_a0.first )
    return 0.0;
  //Next check if we are completely inside:
  if ( b0 >= -ee && b1 <= blim_at_a0.second ) {
    bool crossBMinusLeftSide = ( a0 < ee && b0 < blim_at_a0.first );
    bool crossBMinusRightSide = ( a1 > ee && b0 < blim_at_a1.first );
    if ( !crossBMinusLeftSide && !crossBMinusRightSide ) {
      //yes, completely inside!
      return this->sMax();
    }
  }//  else {

  //Crosses kinematic borders, strip edges:
  const auto trimmed_data = trimToKB( *this, ekin_div_kT, blim_at_a0, blim_at_a1 );

  //For efficiency we just use the values at the corners as the overlay values: (todo: revisit??)
  return ( trimmed_data.has_value() ? trimmed_data.value() : *this ).sMax();
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline void NCS::SABCellEval<AIT,IO>::SCE_Data::integralWKB( StableSum& sum, double ekin_div_kT ) const
{
  const double ee = ekin_div_kT;
  const double b0 = this->beta.first;
  const double b1 = this->beta.second;
  const double a0 = this->alpha.first;
  const double a1 = this->alpha.second;
  //First check if we are completely outside:
  if ( b1 <= -ee )
    return;
  auto blim_at_a1 = getBetaLimits( ekin_div_kT, a1 );
  if ( b0 >= blim_at_a1.second )
    return;
  if ( a1 <= ee && b1 <= blim_at_a1.first )
    return;
  auto blim_at_a0 = getBetaLimits( ekin_div_kT, a0 );
  if ( a0 >= ee && b1 <= blim_at_a0.first )
    return;
  //Next check if we are completely inside:
  if ( b0 >= -ee && b1 <= blim_at_a0.second ) {
    bool crossBMinusLeftSide = ( a0 < ee && b0 < blim_at_a0.first );
    bool crossBMinusRightSide = ( a1 > ee && b0 < blim_at_a1.first );
    if ( !crossBMinusLeftSide && !crossBMinusRightSide ) {
      this->integral(sum);//yes, completely inside!
      return;
    }
  }//  else {

  //Crosses kinematic borders, needs more careful treatment:
  const auto trimmed_data = trimToKB( *this, ekin_div_kT, blim_at_a0, blim_at_a1 );

  intCrossedKB( sum,
                ( trimmed_data.has_value() ? trimmed_data.value() : *this ),
                ekin_div_kT );

  //todo: just call kbStatus? And get blim_at_a0/blim_at_a1 back from there if needed???
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline void NCS::SABCellEval<AIT,IO>::SCE_Data::validate() const {
#ifndef NDEBUG
  nc_assert( !ncisnanorinf(this->alpha.first) );
  nc_assert( !ncisnanorinf(this->alpha.second) );
  nc_assert( !ncisnanorinf(this->beta.first) );
  nc_assert( !ncisnanorinf(this->beta.second) );
  nc_assert( this->alpha.first < this->alpha.second );
  nc_assert( this->beta.first < this->beta.second );
  for ( auto i : ncrange(4) ) {
    nc_assert( !ncisnanorinf(this->S[i]) );
    nc_assert( !this->hasLogValues || !ncisnan(extract_logS(*this,i)) );
    nc_assert( this->S[i]>=0.0 );
  }
#endif
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline double NCS::SABCellEval<AIT,IO>::SCE_Data::integral() const
{
  static_assert( interpOrder == SABInterpolationOrder::ALPHA_FIRST,"not implemented" );
  static_assert( betaInterpType == InterpolationScheme::LINLIN,"not implemented" );
  if ( alphaInterpType == InterpolationScheme::LOGLIN ) {
    static_assert( ( alphaInterpType != InterpolationScheme::LOGLIN || SCE_Data::hasLogValues ), "not implemented" );
    const double alphaIntegralAtBeta0 = integrateAlphaInterval_fast(this->alpha.first, this->S[0], this->alpha.second, this->S[1],
                                                                    extract_logS(*this,0),extract_logS(*this,1));
    const double alphaIntegralAtBeta1 = integrateAlphaInterval_fast(this->alpha.first, this->S[2], this->alpha.second, this->S[3],
                                                                    extract_logS(*this,2),extract_logS(*this,3));
    const double k = 0.5*(this->beta.second-this->beta.first);
    return k * ( alphaIntegralAtBeta0 + alphaIntegralAtBeta1 );
  } else {
    //double lin-lin: area times average height at the four corners:
    const double k = 0.25*(this->alpha.second-this->alpha.first)*(this->beta.second-this->beta.first);
    return k*(this->S[0]+this->S[1]+this->S[2]+this->S[3]);
  }
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline void NCS::SABCellEval<AIT,IO>::SCE_Data::integral( StableSum& sum ) const
{
  static_assert( interpOrder == SABInterpolationOrder::ALPHA_FIRST,"not implemented" );
  static_assert( betaInterpType == InterpolationScheme::LINLIN,"not implemented" );
  if ( alphaInterpType == InterpolationScheme::LOGLIN ) {
    static_assert( ( alphaInterpType != InterpolationScheme::LOGLIN || SCE_Data::hasLogValues ), "not implemented" );
    const double alphaIntegralAtBeta0 = integrateAlphaInterval_fast(this->alpha.first, this->S[0], this->alpha.second, this->S[1],
                                                                    extract_logS(*this,0),extract_logS(*this,1));
    const double alphaIntegralAtBeta1 = integrateAlphaInterval_fast(this->alpha.first, this->S[2], this->alpha.second, this->S[3],
                                                                    extract_logS(*this,2),extract_logS(*this,3));
    const double k = 0.5*(this->beta.second-this->beta.first);
    sum.add( k * alphaIntegralAtBeta0 );
    sum.add( k * alphaIntegralAtBeta1 );
  } else {
    //double lin-lin: area times average height at the four corners:
    const double k = 0.25*(this->alpha.second-this->alpha.first)*(this->beta.second-this->beta.first);
    sum.add( k * this->S[0] );
    sum.add( k * this->S[1] );
    sum.add( k * this->S[2] );
    sum.add( k * this->S[3] );
  }
}


// template <class TData>
// void NCSD::integral(const TData& data, StableSum& sum)
// {
//   static_assert( TData::interpOrder == SABInterpolationOrder::ALPHA_FIRST,"not implemented" );
//   static_assert( TData::betaInterpType == InterpolationScheme::LINLIN,"not implemented" );
//   if ( TData::alphaInterpType == InterpolationScheme::LOGLIN ) {
//     static_assert( ( TData::alphaInterpType != InterpolationScheme::LOGLIN || TData::hasLogValues), "not implemented" );
//     const double alphaIntegralAtBeta0 = integrateAlphaInterval_fast(data.alpha.first, data.S[0], data.alpha.second, data.S[1],
//                                                                     extract_logS(data,0),extract_logS(data,1));
//     const double alphaIntegralAtBeta1 = integrateAlphaInterval_fast(data.alpha.first, data.S[2], data.alpha.second, data.S[3],
//                                                                     extract_logS(data,2),extract_logS(data,3));
//     const double k = 0.5*(data.beta.second-data.beta.first);
//     sum.add( k * alphaIntegralAtBeta0 );
//     sum.add( k * alphaIntegralAtBeta1 );
//   } else {
//     //double lin-lin: area times average height at the four corners:
//     const double k = 0.25*(data.alpha.second-data.alpha.first)*(data.beta.second-data.beta.first);
//     sum.add( k * data.S[0] );
//     sum.add( k * data.S[1] );
//     sum.add( k * data.S[2] );
//     sum.add( k * data.S[3] );
//   }
// }

// template <class TData>
// double NCSD::integral(const TData& data)
// {
//   static_assert( TData::interpOrder == SABInterpolationOrder::ALPHA_FIRST,"not implemented" );
//   static_assert( TData::betaInterpType == InterpolationScheme::LINLIN,"not implemented" );
//   if ( TData::alphaInterpType == InterpolationScheme::LOGLIN ) {
//     static_assert( ( TData::alphaInterpType != InterpolationScheme::LOGLIN || TData::hasLogValues), "not implemented" );
//     const double alphaIntegralAtBeta0 = integrateAlphaInterval_fast(data.alpha.first, data.S[0], data.alpha.second, data.S[1],
//                                                                     extract_logS(data,0),extract_logS(data,1));
//     const double alphaIntegralAtBeta1 = integrateAlphaInterval_fast(data.alpha.first, data.S[2], data.alpha.second, data.S[3],
//                                                                     extract_logS(data,2),extract_logS(data,3));
//     const double k = 0.5*(data.beta.second-data.beta.first);
//     return k * ( alphaIntegralAtBeta0 + alphaIntegralAtBeta1 );
//   } else {
//     //double lin-lin: area times average height at the four corners:
//     const double k = 0.25*(data.alpha.second-data.alpha.first)*(data.beta.second-data.beta.first);
//     return k*(data.S[0]+data.S[1]+data.S[2]+data.S[3]);
//   }
// }


template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline double NCS::SABCellEval<AIT,IO>::SCE_Data::eval( double alphaval, double betaval ) const
{
  if ( interpOrder == SABInterpolationOrder::ALPHA_FIRST ) {
    //First figure out values at (beta0,alpha) and (beta1,alpha):
    const double s0 = evalAtBeta0( *this,  alphaval );
    const double s1 = evalAtBeta1( *this,  alphaval );
    //Then interpolate to actual beta:
    return interpolate<betaInterpType>( this->beta.first, s0, this->beta.second, s1, betaval );
  } else {
    //First figure out values at (beta,alpha0) and (beta,alpha1):
    const double s0 = evalAtAlpha0( *this,  betaval );
    const double s1 = evalAtAlpha1( *this,  betaval );
    //Then interpolate to actual alpha:
    return interpolate<alphaInterpType>( this->alpha.first, s0, this->alpha.second, s1, alphaval );
  }
}

NCS::CellIndex NCS::getCellIndex(const SABData& data, double alpha, double beta )
{
  auto& agrid = data.alphaGrid();
  auto itAB = agrid.begin();
  auto itAE = agrid.end();
  auto itA = std::lower_bound( itAB, itAE, alpha );
  if ( itA == itAE || ( itA == itAB && alpha < *itAB ) )
    return { NullOpt };
  auto& bgrid = data.betaGrid();
  auto itBB = bgrid.begin();
  auto itBE = bgrid.end();
  auto itB = std::lower_bound( itBB, itBE, beta );
  if ( itB == itBE || ( itB == itBB && beta < *itBB) )
    return { NullOpt };
  auto da = std::distance( itAB, itA );
  auto db = std::distance( itBB, itB );
  return CellIndex{ static_cast<CellIndex::index_t>( da ? --da : da),
                    static_cast<CellIndex::index_t>( db ? --db : db) };
}

template<NCS::InterpolationScheme AIT, NCS::SABInterpolationOrder IO>
inline typename NCS::SABEval<AIT,IO>::celleval_t NCS::SABEval<AIT,IO>::getCell( CellIndex idx ) const
{
  nc_assert( idx.isValid() );
  auto& sab = m_sab->sab();
  auto& agrid = m_sab->alphaGrid();
  auto& bgrid = m_sab->betaGrid();
  const auto nalpha = agrid.size();
  const auto itA = std::next(agrid.begin(),idx.ia());
  nc_assert( std::next(itA) < agrid.end());
  const auto itB = std::next(bgrid.begin(),idx.ib());
  nc_assert( std::next(itB) < bgrid.end());
  const auto itSB = sab.begin();
  const auto itS_a0b0 = std::next( itSB, idx.ib()*nalpha+idx.ia() );
  const auto itS_a0b1 = std::next( itS_a0b0, nalpha );
  const auto itS_a1b0 = std::next( itS_a0b0 );
  const auto itS_a1b1 = std::next( itS_a0b1 );
  const double svals[4] = { *itS_a0b0, *itS_a1b0, *itS_a0b1, *itS_a1b1 };
  return celleval_t{ PairDD{ *itA, *std::next(itA) }, PairDD{ *itB, *std::next(itB) }, svals };
}

//Explicit instantiations (reducing the amount of templated code in the header
//file). We might add a LINLIN specialisation as well:
template class NCS::SABCellEval<NCS::InterpolationScheme::LOGLIN,NCS::SABInterpolationOrder::ALPHA_FIRST>;
template class NCS::SABEval<NCS::InterpolationScheme::LOGLIN,NCS::SABInterpolationOrder::ALPHA_FIRST>;

//todo: do we really need the SABInterpolationOrder???
