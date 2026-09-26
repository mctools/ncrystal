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

#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSKnlGrid.hh"
#include "NCrystal/internal/vdos/NCVDOSExpand.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"//fixme?
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCMsg.hh"//fixme?
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"

namespace NC=NCrystal;
namespace NCV=NCrystal::VDOS;

namespace NCRYSTAL_NAMESPACE {
  namespace VDOS {
    namespace detail {
      static inline double stirlingsSeriesSum9thOrder(double);
      static VectD fillSABFromVDOS( const VDOSGn&, const double, const VDOSLux&,
                                    const VectD&, const VectD&,
                                    ScaleGnContributionFct );

      namespace {
        static bool s_verbose = ncgetenv_bool("DEBUG_PHONON");
      }

      static VectD setupLegacyAlphaGrid( double, double, unsigned );
      static VectD setupLegacyBetaGrid( const VDOSGn&, double,
                                        VDOSLux, unsigned );
    }
  }
}

inline double NCV::detail::stirlingsSeriesSum9thOrder(double inv_n)
{
  //Calculate and return the sum_{k=0}^9{ gk/x^k }, needed to estimate
  //Stirling's series (see eq. 5.11.3 and 5.11.4 in
  //https://dlmf.nist.gov/5.11): Thus, pass in inv_n = 1/n and get back
  //1+1/(12n)+1/(288n^2)+... (up to terms with 1/n^9) from Stirling's
  //series, The result must be multiplied by sqrt(2pi*n)*(n/e)^n to provide
  //an estimate of the faculty of n ("n!"). Factors ci are taken from:
  //https://oeis.org/A001163 and https://oeis.org/A001164
  //
  //For arguments n>=9 this can be used to evaluate n! to a relative
  //precision of O(1e-13) or better.
  constexpr double c1 = 1./12.;
  constexpr double c2 = 1./288.;
  constexpr double c3 = -139./51840.;
  constexpr double c4 = -571./2488320.;
  constexpr double c5 = 163879./209018880.;
  constexpr double c6 = 5246819./75246796800.;
  constexpr double c7 = -534703531./902961561600.;
  constexpr double c8 = -4483131259./86684309913600.;
  constexpr double c9 = 432261921612371./514904800886784000.;
  //Horner's method, explicit std::fma per step (reproducible, and this
  //multiply-add pattern silently contracts under -mfma otherwise; see
  //doc/devel_fma_attribute.md):
  double p = c9;
  p = std::fma( inv_n, p, c8 );
  p = std::fma( inv_n, p, c7 );
  p = std::fma( inv_n, p, c6 );
  p = std::fma( inv_n, p, c5 );
  p = std::fma( inv_n, p, c4 );
  p = std::fma( inv_n, p, c3 );
  p = std::fma( inv_n, p, c2 );
  p = std::fma( inv_n, p, c1 );
  return std::fma( inv_n, p, 1.0 );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //sab[i] += aFact[i]*c, for i in [0,n) with explicit std::fma. Used in the
    //hottest loop below:
    NCRYSTAL_FMADISPATCH_ATTR
    void vdosScatKnlAccumulate( double* sab, const double* aFact,
                                double c, std::size_t n )
    {
      for ( std::size_t i = 0; i < n; ++i )
        sab[i] = std::fma( aFact[i], c, sab[i] );
    }
  }
}

NC::VectD
NCV::detail::fillSABFromVDOS( const VDOSGn& Gn_asym,
                              const double alpha2x,
                              const VDOSLux& vdoslux,
                              const VectD& alphaGrid,
                              const VectD& betaGrid,
                              ScaleGnContributionFct scaleGnContribFct )
{
  // Evaluate S(alpha,beta) from Sjolander's II.28, recast to alpha/beta
  // and excluding sigma*kT/4E from the definition of S.

  const unsigned maxOrder = Gn_asym.maxOrder().value();
  const double kT = Gn_asym.kT();
  nc_assert( kT > 0.0 );
  const auto nbeta = betaGrid.size();
  const auto nalpha = alphaGrid.size();
  VectD sab( nalpha * nbeta, 0.0);
  VectD beta_kT( betaGrid );
  for ( auto& e : beta_kT )
    e *= kT;
  nc_assert( nc_is_grid( beta_kT ) );
  VectD workbuf_gneval;
  VectD gneval;
  gneval.reserve(nbeta);


  //Alpha-dependency is contained in f(x,n)=exp(-x)*x^n/n! where
  //x=alpha2x*alpha.  For n<stirling_threshold, f(x,n) are calculated
  //directly and recursively, using f(x,n)=f(x,n-1)*(x/n). For higher n,
  //this becomes numerically unstable and we instead use Stirling's series
  //to rewrite the factor of "n!" and evaluate f(x,n) directly. NB,
  //stirling_threshold should not be changed without validation: raise it
  //too much, and the direct method runs into issues at high n, decrease
  //it too much there will be precision (and speed) issues.
  //
  //The reason for picking 16 is that it seems the stirling formula as
  //implemented here reaches "a few ulps" precision at n=16, and
  //presumably it is better to switch to it as early as possible for
  //numerical stability (although pure computation speed considerations
  //would imply switching as late as possible).
  constexpr unsigned stirling_threshold = 16;//See comment ^^^
  auto x_vals = vectorTrf(alphaGrid,
                          [alpha2x](double alpha){ return alpha*alpha2x; } );
  auto expmhalfx_vals = vectorTrf(x_vals,
                                  [](double x){ return std::exp(-0.5*x); } );
  VectD logx_vals;
  if ( maxOrder >= stirling_threshold )
    logx_vals = vectorTrf(x_vals,[](double x){
      return NC::SABUtils::safeLogOrElse( x, -kInfinity ); } );

  auto fxn_cache = expmhalfx_vals;//We apply half of the exp(-x) factor before
  //the x^n/n! factor is added, and the other
  //half of exp(-x) we add later. This is done
  //to extend range of valid x-values.
  VectD alpha_factors;
  alpha_factors.resize(x_vals.size());

  const double alpha_factor_floor = (vdoslux.isLegacy() ? -1.0 : 1e-100);//fixme tune? Or just 0?

  //Now we loop and fill the S-table. First over phonon order, n, next over
  //beta and finally alpha. We take care to keep as many calculations as
  //possible in the outer loops, while avoiding unnecessarily repeating
  //calculations or utilising enormous memory caches.

  for ( unsigned n = 1; n <= maxOrder; ++n) {
    const double contribScaleFactor = ( scaleGnContribFct
                                        ? scaleGnContribFct(n) : 1.0 );
    nc_assert(contribScaleFactor>=0.0);

    //Prepare for Gn(beta)-evaluations:
    const double invn = 1.0/n;

    //Preparation of f(x)=exp(-x)*x^n/n! is more tricky, due to reasons of
    //efficiency and numerical issues. As explained above, for orders below
    //stirling_threshold, we build up recursively, and above we use
    //Stirling's formula:
    if ( n < stirling_threshold ) {
      //for reasonably low orders we can build up slowly and cheaply (leaving
      //out of fxn_cache a final factor of expmhalfx for numerical stability):
      for (auto x : enumerate(x_vals) )
        vectAt(fxn_cache,x.idx) *= x.val*invn;
      for (auto fxn : enumerate(fxn_cache) )
        vectAt(alpha_factors,fxn.idx) = fxn.val * vectAt(expmhalfx_vals,fxn.idx) * kT;
    } else {
      //Very high order phonons, f(x) is non-zero at very high values of x, but
      //exp(-0.5*x) becomes 0, precluding the direct/cheap evaluation. Instead
      //we evaluate directly, with the help of Stirling's series for the
      //factorial, n!. Everything is suitably rearranged so cancellations happen
      //before the exponential is evaluated.
      const double gn = stirlingsSeriesSum9thOrder(invn);
      const double fact= kT * kInvSqrt2Pi/(std::sqrt(n)*gn);
      if (!fact)
        continue;//nothing can contribute at this order (should not happen?)
      const double logn = std::log(n);
      const double nAsDbl = static_cast<double>(n);
      for (auto x : enumerate(x_vals) ) {
        //Explicit std::fma for the final combination: "n*X-x.val" is
        //exactly the a*b-c shape a compiler may silently fuse under -mfma,
        //feeding straight into exp() below (see docs/devel_fma_attribute.md):
        const double exparg = std::fma( nAsDbl,
                                        vectAt(logx_vals,x.idx) - logn + 1.0,
                                        -x.val );
        nc_assert_always(exparg <= 708.0);
        vectAt(alpha_factors,x.idx) = fact * std::exp(exparg);
        //fixme: break if the alpha_factor just calculated is below
        //alpha_factor_floor and x>n, since the rest will just be zero??
      }
    }

    std::size_t ialphaB = 0;
    std::size_t ialphaE = nalpha;
    while ( ialphaB != ialphaE
            && vectAt(alpha_factors,ialphaB) <= alpha_factor_floor  )
      ++ialphaB;
    while ( ialphaB != ialphaE
            && vectAt(alpha_factors,ialphaE-1) <= alpha_factor_floor  )
      --ialphaE;

    //Evaluate Gn function:
    nc_assert( nc_is_grid( beta_kT ) );
    Gn_asym.evalMany( n, beta_kT, gneval, workbuf_gneval );

    for ( std::size_t ibeta = 0; ibeta < nbeta; ++ibeta ) {
      //fixme: skip earlier if contribScaleFactor is 0
      const double Gn_asym_eval = contribScaleFactor * vectAt(gneval,ibeta);
      if ( !(Gn_asym_eval>0.0) ) {
        continue;
      }
      vdosScatKnlAccumulate( &sab[ibeta*nalpha+ialphaB],
                             alpha_factors.data() + ialphaB,
                             Gn_asym_eval, ialphaE - ialphaB );
    }//beta loop
  }//phonon order loop

  return sab;
}

NC::ScatKnlData
NCV::createScatteringKernel( const VDOSData& vdosdata,
                             VDOSLux vdoslux,
                             Optional<NeutronEnergy> targetEmax_requested,
                             ScaleGnContributionFct scaleGnContributionFct )
{
  //Expand VDOS to Gn functions:
  const auto gnexpn = expandVDOSToGnFcts( vdosdata, vdoslux,
                                          targetEmax_requested);
  //fixmerequest_override_max_order );
  if (detail::s_verbose)
    NCRYSTAL_MSG("VDOS2SK expandVDOSToGnFcts done (maxOrder="
                 <<gnexpn.Gn.maxOrder().value()<<")");
  const auto& Gn_asym = gnexpn.Gn;

  //Setup the alpha/beta grids. The grid-spacing is not even, rather it attempts
  //to best accomodate features of the distributions:

  VectD alphaGrid, betaGrid;

  unsigned nalpha, nbeta;
  std::tie(nalpha, nbeta) = VDOS::gridDimFromLux( vdoslux );

  if ( vdoslux.isLegacy() ) {
    const double betaMin = gnexpn.sabRange.y0()*1.01;
    if ( !( betaMin < 0.0 ) ) {
      NCRYSTAL_THROW2(CalcError,"Beta range does not extend to negative beta. "
                      <<"(increasing vdoslux might help)");
    }
    betaGrid = detail::setupLegacyBetaGrid( gnexpn.Gn, -betaMin,
                                            vdoslux, nbeta );
    nc_assert_always( gnexpn.sabRange.x1()*1.01 > 0.0 );
    alphaGrid = detail::setupLegacyAlphaGrid( gnexpn.alpha2x,
                                              gnexpn.sabRange.x1()*1.01, nalpha );
  } else {
    std::tie( alphaGrid, betaGrid )
      = VDOS::determineAlphaBetaGridFromGn( gnexpn, nalpha, nbeta );
  }

  nc_assert_always(nc_is_grid(alphaGrid));
  nc_assert_always(nc_is_grid(betaGrid));
  if ( !vdoslux.isLegacy() ) {
    nc_assert_always(alphaGrid.front()==0.0);
  }
  if (detail::s_verbose)
    NCRYSTAL_MSG("VDOS2SK alpha/beta grid done (nalpha="
                 <<alphaGrid.size()<<", nbeta="<<betaGrid.size()<<")");

  //All done, now all that remains is to go through the (alpha,beta) pts in the
  //grid and use Sjolander's II.28 equation to calculate S(alpha,beta) there as
  //the sum of individual phonon orders:

  auto sab = detail::fillSABFromVDOS( Gn_asym, gnexpn.alpha2x,
                                      vdoslux, alphaGrid,
                                      betaGrid, scaleGnContributionFct );

  //Type of SABData.suggestedEmax is unfortunately double and not NeutronEnergy:
  double suggestedEmax = gnexpn.suggestedEmax.dbl();

  const auto max_phonon_order = Gn_asym.maxOrder().value();
  if ( scaleGnContributionFct!=nullptr
       && scaleGnContributionFct(max_phonon_order) == 0.0 ) {
    //Caller might have essentially removed the last order(s), so it is unknown
    //how far the kernel can be used.
    suggestedEmax = 0.0;
  }

  if (detail::s_verbose)
    NCRYSTAL_MSG("VDOS2SK created SK with vdos expansion order N="
                 <<max_phonon_order<<", Emax="<<suggestedEmax
                 <<"eV, nalpha="<<alphaGrid.size()
                 << " nbeta="<<betaGrid.size());

  ScatKnlData out;
  out.alphaGrid = std::move(alphaGrid);
  out.betaGrid  = std::move(betaGrid);
  out.sab       = std::move(sab);
  out.temperature = vdosdata.temperature();
  out.boundXS = vdosdata.boundXS();
  out.elementMassAMU = vdosdata.elementMassAMU();
  out.knltype = ScatKnlData::KnlType::SAB;
  out.suggestedEmax = suggestedEmax;
  //prevent beta-thickening code upon conversion to SABData: (fixme: revisit?)
  out.betaGridOptimised = true;
  return out;
}

NC::VectD NCV::detail::setupLegacyAlphaGrid( double alpha2x, double alphaMax,
                                             unsigned npts )
{
  nc_assert(npts>=20);

  //Strategy (region 0 ensures proper binning for upscattering of ultra-low
  //energy neutrons, and are simply merged onto the grid pts from the other 3
  //regions at the end, using the finalise_grid function further down):
  //
  // Region 0 (N=15% of pts): x=linspace from x=1e-3 to x=alpha_upscatmax*alpha2x
  // Region 1 (N=29% of pts): x=1e-50 + linspace[N-1] from x=1e-10 to x=1.0 (including 1.0).
  // Region 2 (N=23% of pts): linspace[N+2] from x=1.0 to x=15.0 (excluding both endpoints).
  // Region 3 (rest, N~=33% of pts): geomspace[N+1] from x=15.0 to final limit (including 15.0).
  //
  //We give special treatment to cases where alphaMax is not above alphaG15Maxx.
  //
  //alpha_upscatmax is roughly 10, but slightly higher or lower depending on
  //number of points.

  nc_assert(alpha2x>0.0);
  const double x2alpha = 1.0 / alpha2x;
  nc_assert(x2alpha>0.0);
  const double alphaMin = x2alpha*1e-50;
  const double alphaMin2 = x2alpha*1e-10;
  const double alphaG1Maxx = x2alpha*1.0;
  const double alphaG15Maxx = x2alpha*15.0;

  unsigned n0 = static_cast<unsigned>(npts*0.15+0.5);
  unsigned n1 = static_cast<unsigned>(npts*0.29+0.5);
  unsigned n2 = static_cast<unsigned>(npts*0.23+0.5);
  unsigned n3 = npts-(n1+n2+n0);
  unsigned npts_123 = n1+n2+n3;
  nc_assert(n1>=3&&n2>=3&&n3>=3&&n0>=3&&n0+n1+n2+n3==npts);

  const double alpha_upscatmax = (n0<10?6.0:(n0>50?14.0:10.0));
  const auto grid_region0 = linspace( ncmin(1e-3,alphaMax*0.01),
                                      ncmin(alpha_upscatmax,alphaMax*0.99),
                                      n0 );

  auto finalise_grid = [&grid_region0,npts](const VectD& grid) -> VectD
  {
    nc_assert_always(grid.size()+grid_region0.size()==npts);

    //We want to merge the grid_region0 points into the
    //grid points. However, the grid_region0 points will be
    //moved a bit to not accidentally coincide with
    //existing points. Thus, for each grid_region0 pt, find
    //the two neighbouring values (looking in *both* grid
    //and grid_region0). Place the point exactly between
    //those two. I.e. just merge+sort all numbers, but make
    //sure we know which are from region0. Those will be
    //placed halfway between neighbours afterwards.
    std::vector<std::pair<double,bool>> all_merged;
    for (auto e: grid)
      all_merged.emplace_back(e,false);
    for (auto e: grid_region0)
      all_merged.emplace_back(e,true);
    std::stable_sort(all_merged.begin(),all_merged.end());
    auto it = std::next(all_merged.begin());
    auto itLast = std::prev(all_merged.end());
    for (;it!=itLast;++it) {
      if ( it->second == false )
        continue;
      it->first = 0.5*(std::prev(it)->first+std::next(it)->first);
    }
    VectD out;
    out.reserve(npts);
    for (const auto& e: all_merged)
      out.push_back(e.first);
    nc_assert_always(nc_is_grid(out));
    nc_assert_always(out.size()==npts);
    return out;
  };

  if ( alphaMax <= alphaMin*100.0 ) {
    //special case #1, insanely small alphaMax
    return finalise_grid(linspace(alphaMax*0.001,alphaMax,npts_123));
  }

  VectD grid;
  grid.reserve(npts_123);
  grid.push_back(alphaMin);
  if ( alphaMax <= alphaG1Maxx*10.0 ) {
    //special case #2, small alphaMax.
    vectorAppend(grid, linspace(alphaMin2,alphaMax,npts_123-1) );
    nc_assert( grid.size() == npts_123 && nc_is_grid(grid) );
    return finalise_grid(grid);
  }

  //Region 1:
  vectorAppend(grid, linspace(alphaMin2,alphaG1Maxx,n1-1) );

  //Region 2:
  if ( alphaMax < 2.0 * alphaG15Maxx) {
    //Special case #3: somewhat small alphaMax => absorb region3 into region 2
    auto ls2 = linspace(alphaG1Maxx,alphaMax,n2+n3+2);
    grid.insert( grid.end(), std::next(ls2.begin()), std::prev(ls2.end()) );
    nc_assert( grid.size() == npts_123 && nc_is_grid(grid) );
    return finalise_grid(grid);
  }
  auto ls2 = linspace(alphaG1Maxx,alphaG15Maxx,n2+2);
  grid.insert( grid.end(), std::next(ls2.begin()), std::prev(ls2.end()) );
  //Region 3:
  vectorAppend(grid, geomspace(alphaG15Maxx,alphaMax,n3) );
  nc_assert( grid.size() == npts_123 && nc_is_grid(grid) );
  return finalise_grid(grid);
}

NC::VectD NCV::detail::setupLegacyBetaGrid( const VDOSGn& Gn,
                                            double betaMax,
                                            VDOSLux vdoslux,
                                            unsigned ntotal )
{
  nc_assert( vdoslux.isLegacy() );
  nc_assert(Gn.maxOrder().value()>=1);
  nc_assert_always(vdoslux.lvl()<=5);
  nc_assert_always(betaMax>0.0);

  //Reach of 1-phonon and 3-phonon spectrums (converted E->beta):
  const double invkT = 1.0/Gn.kT();
  const double G1 = ncabs(Gn.eRange(1).first*invkT);
  nc_assert(G1>0.0);
  double G3 = ncabs(Gn.eRange( std::min<unsigned>(3,Gn.maxOrder().value()) ).first*invkT);//nb: actually =G1 or =G2 in case of maxOrder<3.
  nc_assert( Gn.maxOrder().value()==1 || G3>G1 );
  if (G1==G3)
    G3 = G1*1.0001;


  //Construct betagrid which will run from [-betaMax,D] with
  //0<D<=betaMax. Assuming betaMax is large enough, D will always be large
  //enough to encompass the first phonon order, and energy transfers of 20*kT,
  //whichever is larger (at luxlevel 4 and 5 the 20*kT limit is increased to
  //30*kT and 40*kT respectively).
  //
  //We divide pts in regions:
  //Region0: The single point beta=0.
  //Region1: [-G1,0) (all points flipped and duplicated in (0,G1]).
  //Region2: [-D,-G1) (all points flipped and duplicated in (G1,D]).
  //Region3: [-betaMax,-D) (no duplication of points). Must be odd number of points.
  //
  //always encompass single-phonons, and clip G3/D to betaMax if needed, to
  //ensure no regions will have 0 range:
  betaMax = ncmax(betaMax,G1*1.01);
  G3 = ncmin(G3,betaMax*0.9999);
  const double D = ncmin(betaMax*0.9999,
                         ncmax(G3,std::max<int>(2,vdoslux.lvl()-1)*10.0));
  nc_assert( D>0.0 && G3 > G1 && D>=G3 && D<=betaMax );

  VectD grid;
  grid.reserve(ntotal);

  //First do a bit of pre-analysis for the G1 spectrum in Region1:
  //Here we don't just space the points uniformly, but try to place them to
  //get the single-phonon curve best described. Additionally, we use 20% (but
  //at least 5) of the points to put some very small values around
  //beta=0. This is not needed to describe the actual data, but is used to
  //reduce artefacts in some integration/sampling algorithms (ideally those
  //algs should instead be updated to avoid such artefacts all by
  //themselves...).

  auto rawspec_G1 = Gn.getRawSpectrum(1);
  auto erange_G1 = Gn.eRange(1);
  //Convert E->beta:
  erange_G1.first *= invkT;
  erange_G1.second *= invkT;
  nc_assert(floateq(erange_G1.first,-G1));
  auto evals_G1 = linspace(erange_G1.first,erange_G1.second,rawspec_G1.size());
  nc_assert_always(evals_G1.front()<0.0);
  nc_assert_always(evals_G1.back()>0.0);
  //Ignore positive part and 0 (NB: value representing zero can actually be
  //slightly non-zero due to numerical issues, hence need for epsilon):
  const double epsilon = -0.1*Gn.binWidth(1);
  nc_assert_always(evals_G1.front()<epsilon);
  while ( evals_G1.back()>epsilon )
    evals_G1.pop_back();
  nc_assert(evals_G1.size()>=2);
  rawspec_G1.resize(evals_G1.size());

  //The absolute maximum number of points which can be useful in region 1 is
  //evals_G1.size() plus a bit for fine-grained grid near beta=0:
  const unsigned n1_max = static_cast<unsigned>(evals_G1.size()*1.25+6.5) & ~1;//& ~1 ensures even
  nc_assert(n1_max>=8);

  //Determine how many pts to assign to the different regions, depending on the
  //length of those regions:

  const unsigned n0 = 1;//the beta=0 point

  //Length of beta-scale covered by various regions:
  nc_assert_always(D>G1);
  nc_assert_always(betaMax>D);
  double L1 = 2*G1;
  double L2 = 2*(D-G1);
  double L3 = betaMax-D;
  //Ad-hoc boost of central parts where we expect more variation:
  L1 *= 4;
  L2 *= 2;
  double Lnorm = 1.0/(L1+L2+L3);
  double f1 = ncmax(0.2,L1*Lnorm);
  double f2 = L2*Lnorm;
  if (f1+f2>0.99) {
    double tmp = 0.99/(f1+f2);
    f1 *= tmp;
    f2 *= tmp;
  }

  unsigned n1 = std::min<unsigned>(n1_max,static_cast<unsigned>(f1*ntotal+0.5) / 2);//div 2 due to duplication
  unsigned n2 = static_cast<unsigned>(f2*ntotal+0.5) / 2;//div 2 due to duplication
  n1 = std::max<unsigned>(15,n1);
  n2 = std::max<unsigned>(1,n2);
  unsigned n012;
  do {
    n012 = 2*(n1+n2)+n0;
    if (n012 >= ntotal - 1) {
      if (n2>n1)
        --n2;
      else
        --n1;
    }
    else
      break;
  } while (true);
  nc_assert_always(n1>=10);

  unsigned n3 = ntotal-n012;

  nc_assert(n0==1);
  nc_assert(n1>=10);
  nc_assert(n2>=1);
  nc_assert(n3>=1);
  nc_assert( ntotal == n0 + n3 + 2 * ( n1 + n2) );
  nc_assert_always(n1 <= ntotal );
  nc_assert_always(n2 <= ntotal );
  nc_assert_always(n3 <= ntotal );
  nc_assert_always(D>G1);
  nc_assert_always(betaMax>D);

  //Time to add values to the grid! First Region3:
  {
    auto vals = linspace(-betaMax, -D, n3+1);
    grid.insert(grid.begin(),vals.begin(),std::prev(vals.end()));
    nc_assert(grid.size()==n3);
  }

  const auto idx_R2Start = grid.size();

  //Then Region 2 - negative values.
  {
    auto vals_r2 = linspace(-D,-G1,n2+1);
    vals_r2.resize(vals_r2.size()-1);
    vectorAppend(grid,vals_r2);
  }

  //Now Region 1 - negative values.
  {
    unsigned n1_near0 = std::max<unsigned>(5,static_cast<unsigned>(n1*0.2+0.5));
    unsigned n1_spectrum = n1-n1_near0;
    if (n1_spectrum >= evals_G1.size()) {
      //Great: can keep all points in G1 spectrum:
      n1_spectrum = evals_G1.size();
      n1_near0 = n1 - n1_spectrum;
    } else {
      //Remove least important points from G1 spectrum:
#if 0
      //Old way, just run reducePtsInDistribution. This might leave huge gaps if
      //for instance VDOS is zero over a wide internal region:
      std::tie(evals_G1,rawspec_G1) = reducePtsInDistribution( evals_G1,rawspec_G1, n1_spectrum );
#else
      //New way, avoid leaving huge gaps without pts. Doing so can leave large
      //gaps in the beta points, which when can lead to numerical artifacts when
      //a scattering kernel is later integrated and sampled. If at some point
      //(TODO!) we implement better integration/sampling algorithms, we can
      //hopefully stop doing this - and perhaps also stop needing the n1_near0
      //values.
      unsigned n_for_gaps = ( ( n1_spectrum > 30 && evals_G1.size()-n1_spectrum > 10 )
                              ? std::max<unsigned>(5,static_cast<unsigned>(n1_spectrum*0.1+0.5))
                              : 0 );
      n1_spectrum -= n_for_gaps;
      const double G1binwidth = evals_G1.at(1)-evals_G1.at(0);
      PtReduceCfg prcfg;
      prcfg.equidistant_fraction = 0.0;
      prcfg.tail_floor = 1e-20;
      std::tie(evals_G1,rawspec_G1) = reducePtsInDistribution( evals_G1,
                                                               rawspec_G1,
                                                               n1_spectrum,
                                                               prcfg );
      if ( n_for_gaps > 0 ) {
        struct Gap {
          Gap(double bbb0, double bbb1)
            : b0(bbb0), b1(bbb1)
          {}
          double b0, b1;//edges of gap
          unsigned n = 0;//number of points allocated to fill the gap
          bool operator<( const Gap& o ) const {
            double a = ( b1 - b0 ) / ( n + 1 );
            double b = ( o.b1 - o.b0 ) / ( o.n + 1 );
            if ( floateq(a,b,1e-13,1e-13) )
              return b0 > o.b0;
            return a > b;//reverse sort
          }
        };
        //Find all gaps:
        std::vector<Gap> gaps;
        gaps.reserve(256);
        double dmin = 1.5 * G1binwidth;
        auto it = evals_G1.begin();
        auto itLast = std::prev(evals_G1.end());
        for ( ; it != itLast; ++it ) {
          double dd = *std::next(it)-*it;
          if ( dd > dmin )
            gaps.emplace_back( *it, *std::next(it) );
        }
        //Allocate pts one by one into largest remaining gap:
        while (n_for_gaps>0 && !gaps.empty() ) {
          std::stable_sort(gaps.begin(),gaps.end());
          gaps.front().n += 1;
          --n_for_gaps;
        }
        //Transfer points in gaps to evals_G1:
        for ( const auto& g : gaps ) {
          if ( g.n > 0 ) {
            double bw = (g.b1-g.b0)/(g.n+1.0);
            for ( auto i : ncrange(g.n) ) {
              evals_G1.push_back( g.b0 + (i+1)*bw );
            }
          }
        }
        std::sort(evals_G1.begin(),evals_G1.end());
      }
      if ( n_for_gaps > 0 ) {
        //unused gap filler points, put somewhere else:
        n1_near0 += n_for_gaps;
      }
#endif
    }

    for ( auto e : evals_G1 )
      grid.push_back(e);

    //To reduce artefacts in sampling/integration algs, add points near 0:
    nc_assert( n1_near0 > 0 );
    nc_assert( grid.back() < 0.0 );
    auto v = geomspace(ncmin(1e-50,-0.001*grid.back()),-grid.back()*0.1,n1_near0);
    std::reverse(v.begin(),v.end());
    for (auto e: v)
      grid.push_back( -e );
  }
  const auto idx_R1Back = grid.size()-1;

  //Then Region0, just the beta=0 point:
  grid.push_back( 0.0 );

  //Now replay region1 and region2 in flipped reverse order to generate
  //positive beta values:
  for ( auto i = idx_R1Back; i>=idx_R2Start; --i )
    grid.push_back( -vectAt(grid,i) );

  nc_assert( grid.front()==-betaMax );
  nc_assert( grid.back()==D );
  nc_assert( grid.size() == ntotal );
  nc_assert( nc_is_grid(grid) );
  return grid;
}
