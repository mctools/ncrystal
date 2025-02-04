#ifndef NCrystal_SABUtils_hh
#define NCrystal_SABUtils_hh

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

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/sab/NCScatKnlData.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    //Some input kernels are not properly filling the entire grid on which they
    //are specified. The following function detects and removes such parts of
    //the grid:
    ScatKnlData trimZeroEdgesFromKernel(ScatKnlData&&);

    //Transform ScatKnlData to the standard unscaled and asymmetric
    //S(alpha,beta) format. Internally this calls validateScatKnlData(...) on
    //the input and, in debug builds only, also on the output. Thus, creating
    //SABData objects via this function alleviates the need for littering code
    //elsewhere with validate calls. The function also invokes
    //trimZeroEdgesFromKernel if needed, but will trigger a WARNING message in
    //case anything is trimmed.
    SABData transformKernelToStdFormat(ScatKnlData&&);

    //If beta grid is defined as [0,b1,b2,..,bn] it is assumed that it is a
    //space-saving shorthand for [-bn,...,-b2,-b1,0,b1,b2,..,bn]. Expand, using
    //S(alpha,-beta)=S(alpha,beta):
    void expandBetaAndSABToAllBetas( Span<const double> halfbetagrid,
                                     Span<const double> alphagrid,
                                     Span<const double> sab_for_halfbetagrid,
                                     VectD& complete_betagrid,
                                     VectD& complete_sab );

    //Access S(alpha|beta_i) slice (const and non-const version):
    Span<const double> sliceSABAtBetaIdx_const( Span<const double> sab, std::size_t nalpha, std::size_t beta_idx);
    Span<double> sliceSABAtBetaIdx( Span<double> sab, std::size_t nalpha, std::size_t beta_idx);

    //Or just calculate the index:
    std::size_t calcSABIdx( std::size_t nalpha, std::size_t alpha_idx, std::size_t beta_idx);

    //interpolate "loglin" (linear in log(f), fallback to linear when undefined)
    double interpolate_loglin_fallbacklinlin(double a, double fa, double b, double fb, double x);
    double interpolate_loglin_fallbacklinlin_fast(double a, double fa, double b, double fb, double x, double logfa, double logfb);
    double interpolate_linear(double a, double fa, double b, double fb, double x);

    //Templated selection based on enum:
    enum class InterpolationScheme { LOGLIN, LINLIN };
    template<InterpolationScheme scheme>
    double interpolate(double a, double fa, double b, double fb, double x);

    //Functions for sampling:
    double sampleLogLinDist(double a, double fa, double b, double fb, double rand);
    double sampleLogLinDist_fast(double a, double fa, double b, double fb, double rand, double logfa, double logfb);

    //The next two functions can be used to evaluate a SAB integral over a bin
    //of alpha-values, assuming a linear relationship in log(s).  The analytical
    //solution when s1,s2 are both nonzero is (a2-a1)*(s2-s1)/log(s2/s1), but
    //the evaluation will take performance and numerical precision issues into
    //account. In case either s1 or s2 is zero, the function reverts to a simple
    //trapezoidal integral.

    //Version without pre-calculated logs. Returns full (better than 1e-15) precision results.
    double integrateAlphaInterval(double a1,double s1, double a2 , double s2 );

    //Faster version with pre-calculated logs. Returns slightly reduced (better than 1e-14) precision results.
    double integrateAlphaInterval_fast(double a1,double s1, double a2 , double s2 , double logs1, double logs2);

    //Find the grid cells touched by the kinematically accessible region for
    //ekin_div_kt = ekin/kT. The ibeta_low index will indicate the lowest
    //beta-bin, which will span [ibeta_low,ibeta_low+1], and the alpha-values
    //indices in the out_alpharanges vector will donate the upper and lower
    //edges of the touched grid cells at the corresponding beta indices
    //(starting from ibeta_low of course, no need to lead with a number of empty
    //entries). For some grids it is possible that the allowed alpha-range will
    //fall above alpha_max at very high beta values, or that there are other
    //reasons (such as numerically alpha+=alpha- for extremely low energies) for
    //empty cells. Such are marked by having the indices set to invalid values
    //(alphaGrid.size()).
    void activeGridCells( const SABData& data,
                          double ekin_div_kT,
                          std::vector<std::pair<uint16_t,uint16_t>>& out_alpharanges,
                          std::size_t& ibeta_low );

    //Expand kinematic alpha-limits at each beta point outwards to nearest
    //points in the alpha-grid.  Thus, this function appears similar to the
    //activeGridCells function, but concerns 1D ranges at each beta-grid point,
    //rather than the range of 2D cells. Of course, ranges are limited to the
    //extent of the alpha grid. Empty ranges are marked by indices both equal to
    //alphaGrid.size() in the out_alpharanges vector.
    void activeGridRanges( const SABData& data,
                           double ekin_div_kT,
                           std::vector<std::pair<uint16_t,uint16_t>>& out_alpharanges,
                           std::size_t& ibeta_low );

    //Break down how XS contributions falls over whole grid cells, as well as
    //front and back tails covering partial cells.
    struct TailedBreakdown {
      double xs_front=0.0, xs_middle=0.0, xs_back=0.0;
      unsigned imiddle_low=0, imiddle_upp=0;
      struct TailPoint {
        double alpha=0.0, sval=0.0, logsval=0.0;
      } front, back;
      bool narrow = false;//if range inside single grid bin
    };
    TailedBreakdown createTailedBreakdown( const Span<const double>& alphaGrid,
                                           const Span<const double>& sab,
                                           const Span<const double>& logsab,
                                           const Span<const double>& alphaIntegrals_cumul,
                                           double alpha_low, double alpha_upp,
                                           const unsigned aidx_low, const unsigned aidx_upp );

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::Span<const double> NCrystal::SABUtils::sliceSABAtBetaIdx_const( NCrystal::Span<const double> sab, std::size_t nalpha, std::size_t beta_idx)
{
  nc_assert( (beta_idx+1) * nalpha <= static_cast<std::size_t>(sab.size()) );
  return {sab.begin() + beta_idx*nalpha, sab.begin() + (beta_idx+1)*nalpha};
}

inline NCrystal::Span<double> NCrystal::SABUtils::sliceSABAtBetaIdx( NCrystal::Span<double> sab, std::size_t nalpha, std::size_t beta_idx)
{
  nc_assert( (beta_idx+1) * nalpha <= static_cast<std::size_t>(sab.size()) );
  return {sab.begin() + beta_idx*nalpha, sab.begin() + (beta_idx+1)*nalpha};
}

inline std::size_t NCrystal::SABUtils::calcSABIdx( std::size_t nalpha, std::size_t alpha_idx, std::size_t beta_idx)
{
  return beta_idx * nalpha + alpha_idx;
}

inline double NCrystal::SABUtils::interpolate_loglin_fallbacklinlin(double a, double fa, double b, double fb, double x)
{
  nc_assert ( fa>=0.0 && fb >= 0.0 );
  nc_assert ( x >= a && x <= b );
  const double midpoint = 0.5 * ( b + a );
  const bool linlin_mode = ( fa*fb == 0.0 );
  if ( x < midpoint ) {
    //choose form most numerically stable for x near a
    const double r = (x-a) / (b-a);
    return ( linlin_mode ? ( fa + (fb-fa)*r ) : fa * std::pow(fb/fa,r) );
  } else {
    //choose form most numerically stable for x near b
    const double s = (b-x) / (b-a);
    return ( linlin_mode ? ( fb + (fa-fb)*s ) : fb * std::pow(fa/fb,s) );
  }
}

inline double NCrystal::SABUtils::interpolate_loglin_fallbacklinlin_fast(double a, double fa, double b, double fb, double x, double logfa, double logfb)
{
  nc_assert ( fa>=0.0 && fb >= 0.0 );
  nc_assert ( x >= a && x <= b );
  nc_assert ( !ncisnan(logfa) && !ncisnan(logfb) );
  const double midpoint = 0.5 * ( b + a );
  const bool linlin_mode = ( fa*fb == 0.0 );
  if ( x < midpoint ) {
    //choose form most numerically stable for x near a
    const double r = (x-a) / (b-a);
    return ( linlin_mode ? ( fa + (fb-fa)*r ) : std::exp(logfa+(logfb-logfa)*r) );
  } else {
    //choose form most numerically stable for x near b
    const double s = (b-x) / (b-a);
    return ( linlin_mode ? ( fb + (fa-fb)*s ) : std::exp(logfb+(logfa-logfb)*s) );
  }
}

inline double NCrystal::SABUtils::integrateAlphaInterval(double a1,double s1, double a2 , double s2 )
{
  nc_assert( s1 >= 0.0 && s2 >= 0.0 );
  nc_assert( a2 > a1 );

  //With ds=s2-s1, ps=s1+s2, da=a2-a1, y=ds/ps, one gets 2*s2=ps+ds,
  //2*s1=ps-ds and thus the analytical result can be expressed as:
  //  da*ps*y/log( (ps+ds) / (ps-ds) ) = da*ps* [ y/log((1+y)/(1-y)) ]
  //This is an even function in y. As s1,s2>0, |y|<=1. When |y| is
  //sufficiently small, the expression can be evaluated via a Taylor expansion:
  //
  // y/log((1+y)/(1-y)) = [1/2]-[1/6]y^2-[2/45]y^4-[22/945]y^6-[214/14175]y^8-[5098/467775]y^10-[5359534/638512875]y^12
  //
  //The Taylor expansion is not only faster, as it avoids a log call, but it
  //is also numerically safer when s1~=s2.
  const double da = a2-a1;
  const double ps = s1+s2;
  const double ds = s2-s1;
  if ( ncabs(ds) < 0.10*ps ) {
    //both s1 and s2 are nonzero and s1 ~= s2 so we can use the Taylor expansion
    nc_assert( ps > 0. );
    const double y = ds / ps;
    const double c1           = 0.166666666666666666666666666666666666666666667   ; // = 1/6
    const double c2           = 0.0444444444444444444444444444444444444444444444  ; // = 2/45
    const double c3           = 0.0232804232804232804232804232804232804232804233  ; // = 22/945
    const double c4           = 0.0150970017636684303350970017636684303350970018  ; // = 214/14175
    const double c5           = 0.0108984020095131206242317353428464539575650687  ; // = 5098/467775
    const double c6           = 0.00839377592816746255899695053134206573360012513 ; // = 5359534/638512875
    const double ysq = y*y;
    return da*ps*(0.5-ysq*(c1+ysq*(c2+ysq*(c3+ysq*(c4+ysq*(c5+ysq*c6))))));
  }
  //Evaluate via analytical expression or fall-back to trapezoidal integration:
  return ncmin(s1,s2)<1e-300 ? 0.5*da*ps : da*ds/std::log(s2/s1);
}

inline double NCrystal::SABUtils::integrateAlphaInterval_fast(double a1,double s1, double a2 , double s2 , double logs1, double logs2)
{
  //Similar implementation as integrateAlphaInterval, but due to the
  //precalculated logs, the Taylor expansion will be used more rarely and
  //with fewer terms (it is not needed for speed reasons, just numerical
  //stability). We also reorder, to put the most likely code first:
  nc_assert(!ncisnan(logs1)&&!ncisnan(logs2));
  nc_assert( s1 >= 0.0 && s2 >= 0.0 );
  nc_assert( a2 > a1 );
  const double da = a2-a1;
  const double ps = s1+s2;
  const double ds = s2-s1;
  if ( ncmin(s1,s2) < 1e-300 )//1e-300, because we want to stay clear of log(s)=-inf issues.
    return 0.5*da*ps;//trapezoidal integration
  if ( ncabs(ds) > 0.006*ps ) {
    nc_assert_always(!ncisinf(logs1)&&!ncisinf(logs2));//should have been caught by checks above
    return da*ds/(logs2-logs1);
  }
  nc_assert( ps > 0. );
  const double y = ds / ps;
  const double c1           = 0.166666666666666666666666666666666666666666667   ; // = 1/6
  const double c2           = 0.0444444444444444444444444444444444444444444444  ; // = 2/45
  const double c3           = 0.0232804232804232804232804232804232804232804233  ; // = 22/945
  const double ysq = y*y;
  return da*ps*(0.5-ysq*(c1+ysq*(c2+ysq*c3)));
}

inline double NCrystal::SABUtils::sampleLogLinDist(double a, double fa, double b, double fb, double rand)
{
  //(NB: Code duplicated below!)
  nc_assert(b>a);
  nc_assert(fa>=0.);
  nc_assert(fb>=0.);
  double df = fb-fa;
  if ( fa*fb*df != 0.0 ) {
    //usual, non-degenerate case:
    double a_sub_b = a-b;
    double logfa_fb = std::log(fb/fa);
    if ( a_sub_b * logfa_fb != 0.0 )
      return a_sub_b*std::log(fa*std::exp(a*logfa_fb/a_sub_b)/(fa+rand*df))/logfa_fb;
    df = 0.0;
  }
  if (!df)
    return a + rand*(b-a);//fa=fb, select uniformly in [a,b]
  //exactly one of fa and fb is 0:
  nc_assert(fa||fb);
  double x = (b-a)*std::sqrt(rand);
  return (fa ? b-x : a + x);
}

inline double NCrystal::SABUtils::sampleLogLinDist_fast(double a, double fa, double b, double fb, double rand, double logfa, double logfb)
{
  //(NB: Code duplicated above!)
  //Slightly faster. Although, still log+exp remains so not sure that it is really worth the extra caching overhead.
  nc_assert(b>a);
  nc_assert(fa>=0.);
  nc_assert(fb>=0.);
  double df = fb-fa;
  if ( fa*fb*df != 0.0 ) {
    //usual, non-degenerate case:
    double a_sub_b = a-b;
    double logfa_fb = logfb-logfa;
    if ( a_sub_b * logfa_fb != 0.0 )
      return a_sub_b*std::log(fa*std::exp(a*logfa_fb/a_sub_b)/(fa+rand*df))/logfa_fb;
    df = 0.0;
  }
  if (!df)
    return a + rand*(b-a);//fa=fb, select uniformly in [a,b]
  //exactly one of fa and fb is 0:
  nc_assert(fa||fb);
  double x = (b-a)*std::sqrt(rand);
  return (fa ? b-x : a + x);
}

inline double NCrystal::SABUtils::interpolate_linear(double a, double fa, double b, double fb, double x)
{
  nc_assert( b > a );
  nc_assert( !ncisnanorinf(x) );
  nc_assert( x >= a && x <= b );
  double r = ( x - a ) / ( b - a );
  return (1-r) * fa + r * fb;
}

template<NCrystal::SABUtils::InterpolationScheme scheme>
inline double NCrystal::SABUtils::interpolate(double a, double fa, double b, double fb, double x) {
  return ( scheme == InterpolationScheme::LOGLIN
           ? interpolate_loglin_fallbacklinlin(a, fa, b, fb, x)
           : interpolate_linear(a, fa, b, fb, x) );
}

#endif
