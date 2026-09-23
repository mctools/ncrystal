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

////////////////////////////////////////////////////////////////////////////////
// Testbed for NC::VDOS::estimateGnErange -- a standalone utility, extracted
// from (but NOT YET wired into) VDOSGn::eRange, for estimating the energy
// interval outside of which a Gn spectrum is everywhere below
// relcontriblvl*max(spec). See docs/claude_session_vdos_fma_reprod.md for the
// background: production VDOSGn::eRange snaps the crossing to whichever grid
// point first clears the threshold. This (a) discretises a continuous
// root-finding problem to bin-width resolution, amplifying a tiny (last-ULP,
// cross-platform-varying) input difference into a whole-bin-sized output
// difference whenever the true crossing sits close to a bin edge, and (b) is
// "one-sided": the single point that happens to first clear the threshold is
// trusted as ground truth, with no check on whether that very point's own
// value might be the one affected by such noise.
//
// This app builds synthetic curves (a clean Gaussian-shaped bump, the same
// bump with noise injected right at the naive crossing points, and a hard
// step) plus real Gn spectra captured from two cases that showed up as
// cross-platform-sensitive on actual CI runs (a Li2O-derived Li spectrum at
// order 129, and a synthetic Debye-model spectrum at order 153 -- both
// vdoslux=2000, i.e. relcontriblvl=1e-6), and compares four crossing
// estimators against each other (and, for the synthetic cases, against the
// true analytic answer):
//
//  * snapToGrid:       today's production VDOSGn::eRange behaviour.
//  * linear2pt:        linear interpolation between the immediate bracket.
//  * logLinear2pt:      log-space linear interpolation between the bracket.
//  * estimateGnErange: the new windowed quadratic-log-fit utility
//                      (NCVDOSUtils.cc): high-order Gn spectra approach a
//                      Gaussian shape via the central limit theorem, i.e.
//                      ln(spec) is close to quadratic (not linear) in the
//                      tails, so a windowed quadratic fit was tried after an
//                      initial windowed LINEAR fit turned out to be *worse*
//                      than the plain 2-point logLinear2pt on the synthetic
//                      Gaussian case (curvature bias outweighing the
//                      noise-averaging benefit of widening the window) --
//                      the quadratic fit gets both: essentially exact on the
//                      noise-free Gaussian case (errors at the level of
//                      floating-point round-off), and far less sensitive
//                      than logLinear2pt to noise injected at the bracket
//                      points (a 1e-3 relative perturbation there moves the
//                      quadratic fit's answer by only ~1e-6 relative,
//                      instead of directly showing up at the ~1e-3 level).
//
// Deliberately verbose/exploratory rather than a narrow pass/fail check --
// this is meant to be read, and rerun with variations, while the algorithm in
// estimateGnErange is refined.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;

namespace {

  //////////////////////////////////////////////////////////
  // Reference implementations of the other three methods //
  //////////////////////////////////////////////////////////

  NC::PairDD snapToGrid( NC::Span<const double> egrid,
                         NC::Span<const double> spec, double threshold )
  {
    NC::PairDD r( egrid.front(), egrid.back() );
    for ( std::size_t i = 0; i < spec.size(); ++i ) {
      if ( spec[i] >= threshold ) {
        r.first = egrid[i];
        break;
      }
    }
    for ( std::size_t i = spec.size(); i > 0; --i ) {
      if ( spec[i-1] >= threshold ) {
        r.second = NC::ncmin( r.second, egrid[i-1] );
        break;
      }
    }
    return r;
  }

  double bracketInterp( NC::Span<const double> egrid, NC::Span<const double> spec,
                        std::size_t idxLo, std::size_t idxHi, double threshold,
                        bool logspace )
  {
    const double v0 = spec[idxLo];
    const double v1 = spec[idxHi];
    double t;
    if ( logspace && v0 > 0.0 && v1 > 0.0 )
      t = ( std::log(threshold)-std::log(v0) ) / ( std::log(v1)-std::log(v0) );
    else
      t = (threshold-v0)/(v1-v0);
    return NC::nclerp( egrid[idxLo], egrid[idxHi], NC::ncclamp(t,0.0,1.0) );
  }

  NC::PairDD twoPoint( NC::Span<const double> egrid, NC::Span<const double> spec,
                       double threshold, bool logspace )
  {
    NC::PairDD r( egrid.front(), egrid.back() );
    for ( std::size_t i = 0; i < spec.size(); ++i ) {
      if ( spec[i] >= threshold ) {
        r.first = ( i==0 ? egrid.front()
                   : bracketInterp(egrid,spec,i-1,i,threshold,logspace) );
        break;
      }
    }
    for ( std::size_t i = spec.size(); i > 0; --i ) {
      if ( spec[i-1] >= threshold ) {
        const double x = ( i==spec.size() ? egrid[i-1]
                          : bracketInterp(egrid,spec,i-1,i,threshold,logspace) );
        r.second = NC::ncmin( r.second, x );
        break;
      }
    }
    return r;
  }

  //////////////////////
  // Reporting helper //
  //////////////////////

  void report( const char* label,
              NC::Span<const double> egrid, NC::Span<const double> spec,
              double relcontriblvl,
              const NC::PairDD* truth = nullptr )
  {
    const double spec_max = *std::max_element(spec.begin(),spec.end());
    const double threshold = relcontriblvl*spec_max;
    std::cout << "=== " << label << " (npts=" << spec.size()
              << ", relcontriblvl=" << NC::fmtg(relcontriblvl) << ") ===" << std::endl;
    auto pr = [&]( const char* name, NC::PairDD r )
    {
      std::cout << "  " << name << ": [" << NC::fmt(r.first,"%.10g")
                << ", " << NC::fmt(r.second,"%.10g") << "]";
      if ( truth ) {
        const double elo = r.first - truth->first;
        const double ehi = r.second - truth->second;
        std::cout << "  (err_lo=" << NC::fmt(elo,"%.3g")
                  << ", err_hi=" << NC::fmt(ehi,"%.3g") << ")";
      }
      std::cout << std::endl;
    };
    pr( "snapToGrid      ", snapToGrid(egrid,spec,threshold) );
    pr( "linear2pt       ", twoPoint(egrid,spec,threshold,false) );
    pr( "logLinear2pt    ", twoPoint(egrid,spec,threshold,true) );
    pr( "estimateGnErange", NC::VDOS::estimateGnErange(egrid,spec,relcontriblvl) );
    if ( truth )
      std::cout << "  truth           : [" << NC::fmt(truth->first,"%.10g")
                << ", " << NC::fmt(truth->second,"%.10g") << "]" << std::endl;
    std::cout << std::endl;
  }

  ///////////////////////////////////////////////
  // Synthetic case 1+2: clean/noisy Gaussian.  //
  ///////////////////////////////////////////////

  //Gaussian bump spec(x) = exp(-k*(x-mu)^2), sampled on an equidistant grid
  //of npts points spanning [mu-halfwidth,mu+halfwidth]. Truth (continuous
  //crossing of relcontriblvl*spec_max=relcontriblvl, since spec_max=1 at
  //x=mu) is exactly mu +- sqrt(ln(1/relcontriblvl)/k), modelling the
  //Gaussian-like shape real high-order Gn spectra approach (central limit
  //theorem acting on repeated self-convolution):
  void gaussianCase( const char* label, double k, double mu, double halfwidth,
                     std::size_t npts, double relcontriblvl,
                     double noiseAmplRel = 0.0 )
  {
    NC::VectD egrid(npts), spec(npts);
    const double binwidth = 2.0*halfwidth/(npts-1);
    for ( auto i : NC::ncrange(npts) ) {
      const double x = mu - halfwidth + i*binwidth;
      egrid[i] = x;
      spec[i] = std::exp( -k*(x-mu)*(x-mu) );
    }
    if ( noiseAmplRel != 0.0 ) {
      //Perturb specifically the points nearest the naive crossing on both
      //sides, mimicking a last-ULP-level fluctuation landing exactly on the
      //single point the old snap-to-grid/2pt methods would trust as their
      //sole anchor -- this is the "one-sidedness" concern:
      const double threshold = relcontriblvl;//spec_max=1
      for ( auto i : NC::ncrange(std::size_t(1),npts) ) {
        const bool crossesUp = spec[i-1] < threshold && spec[i] >= threshold;
        const bool crossesDown = spec[i-1] >= threshold && spec[i] < threshold;
        if ( crossesUp || crossesDown ) {
          spec[i-1] *= ( 1.0 + noiseAmplRel );
          spec[i] *= ( 1.0 - noiseAmplRel );
        }
      }
    }
    //Truth must be relative to the actual (discretely sampled) spec_max, not
    //the continuous peak value of 1.0 -- these differ whenever the grid does
    //not happen to include a point exactly at mu, which would otherwise
    //contaminate the error comparison below with a discretisation artifact
    //unrelated to which crossing-estimator is more accurate:
    const double spec_max = *std::max_element(spec.begin(),spec.end());
    const double d = std::sqrt( std::log( 1.0/(relcontriblvl*spec_max) )/k );
    NC::PairDD truth( mu-d, mu+d );
    report( label, egrid, spec, relcontriblvl, &truth );
  }

  ///////////////////////////////////
  // Synthetic case 3: a hard step //
  ///////////////////////////////////

  //spec(x) = 1 for x in [mu-w,mu+w], 0 outside: no interpolation scheme can
  //improve on the true crossing here (it coincides with a grid point by
  //construction), but this checks that none of the methods misbehave
  //(overshoot outside the window, nan, etc) when confronted with a genuine
  //discontinuity rather than a smooth tail:
  void stepCase( double mu, double w, std::size_t npts, double relcontriblvl )
  {
    NC::VectD egrid(npts), spec(npts);
    const double halfwidth = 2.0*w;
    const double binwidth = 2.0*halfwidth/(npts-1);
    for ( auto i : NC::ncrange(npts) ) {
      const double x = mu - halfwidth + i*binwidth;
      egrid[i] = x;
      spec[i] = ( NC::ncabs(x-mu) <= w ) ? 1.0 : 0.0;
    }
    report( "hard step (discontinuous, robustness check)",
            egrid, spec, relcontriblvl );
  }

  ///////////////////////////////////////////////////////////////
  // Real data, captured via temporary debug instrumentation in //
  // VDOSGn::eRange while investigating the sabxs/nctool CI     //
  // divergences documented in                                  //
  // docs/claude_session_vdos_fma_reprod.md.                    //
  ///////////////////////////////////////////////////////////////

  //Debye-model VDOS (debyeTemp=410K, as in app_fmavdos's testExpansionGrids),
  //VDOSLux(2000) => relcontriblvl=1e-6, order=153 (=maxOrder for this case):
  void debyeOrder153Case()
  {
    const double lower = -2.8995131203490909;
    const double binwidth = 0.021638157614545454;
    const NC::VectD spec = {
      3.5285831638648959e-13, 7.9416557395114186e-13, 1.7616423881341544e-12,
      3.8535173829849065e-12, 8.3130665343855631e-12, 1.7690057548626921e-11,
      3.7138413843954209e-11, 7.6932143671248009e-11, 1.5727079816181896e-10,
      3.1732752472983032e-10, 6.3204303227586635e-10, 1.2428655777294966e-09,
      2.4132245839823095e-09, 4.6272549350747501e-09, 8.7630593331003398e-09,
      1.639260436272431e-08, 3.0293671816193269e-08, 5.5311858116976203e-08,
      9.9792139760772752e-08, 1.7792393440103246e-07, 3.1353021974917828e-07,
      5.4610651540732562e-07, 9.4031435330186411e-07, 1.60070135957299e-06,
      2.694212615586963e-06, 4.4841433595257316e-06, 7.3806393942868516e-06,
      1.2014763459072707e-05, 1.9345597850927932e-05, 3.0812919117973241e-05,
      4.8551705751815892e-05, 7.5689121329201379e-05, 0.00011674957341874873,
      0.0001781987930855037, 0.00026916327608871672, 0.00040236635979588295,
      0.00059532595530774812, 0.00087186056344022158, 0.0012639485083900427,
      0.0018139789748601025, 0.0025774209792698744, 0.003625916409142115,
      0.0050507744862018054, 0.0069668066224631347, 0.0095163924911285704,
      0.012873611017277105, 0.01724820590251613, 0.022889087636759339,
      0.030087007636656948, 0.039175981575036593, 0.050532995765671961,
      0.06457551108032486, 0.081756291869027622, 0.10255514055491548,
      0.12746721803513283, 0.15698777892859239, 0.1915933483538089,
      0.23171960782628281, 0.27773653130281006, 0.32992160234025386,
      0.38843222899600005, 0.45327873011230807, 0.52429946883131384,
      0.60113983090130474, 0.68323676397653088, 0.76981049288754322,
      0.85986479622426193, 0.95219687316918444, 1.0454173592716456,
      1.1379804900430215, 1.2282237963953515, 1.3144160885298022,
      1.394811892074874, 1.4677099898955579, 1.5315133391269804,
      1.5847874116907077, 1.6263139721490338, 1.6551374690954923,
      1.670601569354927, 1.6723738863069157, 1.6604576087853797,
      1.6351894783179104, 1.5972243360850582, 1.5475072108704018,
      1.4872345921650256, 1.4178070825457314, 1.3407760156183326,
      1.2577868390439659, 1.1705220901795268, 1.0806466430574726,
      0.98975760129216728, 0.89934078393764594, 0.81073523923224211,
      0.72510666663609902, 0.64343007211197034, 0.56648146279255385,
      0.49483793597390952, 0.4288851563010837, 0.36883095730788124,
      0.31472365314921691, 0.26647359906596729, 0.22387658360134088,
      0.18663775557302195, 0.15439496503057465, 0.12674060961018338,
      0.10324130629190324, 0.083454936158531889, 0.06694482205918019,
      0.053290985467928265, 0.042098582474756188, 0.033003736526340789,
      0.025677067156876489, 0.019825261927261625, 0.015191057337780528,
      0.011551988904443325, 0.0087182466599474047, 0.0065299357493555938,
      0.0048539977512721145, 0.0035810013307342255, 0.00262196438798126,
      0.0019053266418256065, 0.0013741533234781578, 0.00098361831818445634,
      0.00069878899386986585, 0.00049271492018100117, 0.00034480818475462945,
      0.00023949332222448103, 0.00016509915287753028, 0.00011296223845621545,
      7.6711410463931823e-05, 5.170421104229645e-05, 3.4588520003554188e-05,
      2.2965654051122935e-05, 1.5134459047519326e-05, 9.899118289663627e-06,
      6.4263998885247204e-06, 4.1407642578271976e-06, 2.6481007848426415e-06,
      1.6808512754949705e-06, 1.0589225106143327e-06, 6.621224902628531e-07,
      4.1091380125808288e-07, 2.5310464637482963e-07, 1.5473376424129942e-07,
      9.388669205016139e-08, 5.6539915605942196e-08, 3.3793754635371761e-08,
      2.0046830161559381e-08, 1.1802677886406865e-08, 6.8966529691565971e-09,
      3.999598201576595e-09, 2.3020326485236863e-09, 1.3149838373954515e-09,
      7.4548638421514508e-10, 4.1943592926081002e-10, 2.3420376771729576e-10,
      1.2978378303610078e-10, 7.1374007102099835e-11, 3.8953576988022259e-11,
      2.1097728815440002e-11, 1.1339667041341752e-11, 6.0483158719457259e-12,
      3.2012842168437777e-12, 1.6814112515985862e-12, 8.7621725469248586e-13,
      4.5313957497504441e-13, 2.3248425706281867e-13,
    };
    NC::VectD egrid(spec.size());
    for ( auto i : NC::ncrange(spec.size()) )
      egrid[i] = NC::VDOS::equidistantGridPoint(lower,binwidth,i);
    report( "real: Debye VDOS (T=410K), order=153, vdoslux=2000",
            egrid, spec, 1e-6 );
  }

  //Li2O_sg225_LithiumOxide.ncmat, Li element, temp=293.15K,
  //VDOSLux(2000) => relcontriblvl=1e-6, order=129 (=maxOrder for Li here):
  void li2oOrder129Case()
  {
    const double lower = -6.737992804182225;
    const double binwidth = 0.037853892158327108;
    const NC::VectD spec = {
      2.4018326197233627e-13, 5.847091923795665e-13, 1.3985335023138133e-12,
      3.287551586393825e-12, 7.596400478881977e-12, 1.7258200862047844e-11,
      3.8559030717202e-11, 8.47406107970998e-11, 1.8322458087583387e-10,
      3.8984279065390105e-10, 8.163842450527712e-10, 1.6829970141712412e-09,
      3.41615842389853e-09, 6.82871988226768e-09, 1.3445169911653304e-08,
      2.607933420058272e-08, 4.984319441666212e-08, 9.387908367336454e-08,
      1.7428435116762332e-07, 3.189672173223916e-07, 5.755749961979074e-07,
      1.0242234790262413e-06, 1.7975949519039726e-06, 3.1121411982908906e-06,
      5.315716593613001e-06, 8.959068821075785e-06, 1.4901317214939004e-05,
      2.4462912365000082e-05, 3.964360593032082e-05, 6.342759378283182e-05,
      0.00010020303324737984, 0.00015632821563798445, 0.00024088123157143182,
      0.0003666331939276956, 0.0005512859273664137, 0.0008190122092935482,
      0.0012023287032246702, 0.0017443171596581363, 0.00250118689701343,
      0.0035451400040914117, 0.004967459760020927, 0.0068816930338309315,
      0.009426740733254944, 0.012769610023051513, 0.017107522910621523,
      0.022669024299717278, 0.02971369636459355, 0.03853007340092262,
      0.04943137035602322, 0.06274869609882694, 0.07882152397280144,
      0.09798533869817813, 0.12055656727732275, 0.14681512416549122,
      0.17698514428986048, 0.21121472336576252, 0.24955571142298513,
      0.29194478852064343, 0.3381871675276844, 0.387944296653267,
      0.44072685872415535, 0.49589417755123383, 0.5526608463628436,
      0.6101110021019235, 0.6672202057075243, 0.722884384676335,
      0.775954789069422, 0.8252774475854637, 0.8697352272901724,
      0.9082903348874516, 0.9400249759513069, 0.9641779263240422,
      0.9801749682318384, 0.9876514897638029, 0.9864660143480362,
      0.9767039802632257, 0.9586716854680145, 0.9328809036263005,
      0.9000252181618456, 0.8609495731110959, 0.8166148720955443,
      0.7680596509053499, 0.7163608984397459, 0.6625960108264516,
      0.6078076510820789, 0.5529729768967365, 0.49897832287329263,
      0.446600014076817, 0.3964915775813423, 0.3491772369139407,
      0.3050512444104519, 0.26438234514586095, 0.2273224825021407,
      0.19391875160972946, 0.16412757851972604, 0.13783014065447943,
      0.11484813496705472, 0.0949591295554225, 0.07791088709949039,
      0.06343421021411948, 0.051254017333390925, 0.04109850328483777,
      0.03270636432305377, 0.02583216892929663, 0.0202500315724535,
      0.01575579744366068, 0.012167974142593291, 0.009327654738372113,
      0.007097669485752819, 0.005361184844767866, 0.004019942240453038,
      0.0029922986845395006, 0.0022111998512400702, 0.0016221856890777466,
      0.0011815007447430535, 0.0008543570629823144, 0.0006133772910623087,
      0.00043722953688567727, 0.0003094533895800244, 0.00021746789730275166,
      0.0001517466861688891, 0.00010514223345749383, 7.234002615491207e-05,
      4.9423438097211494e-05, 3.353120906414468e-05, 2.2591047821786635e-05,
      1.5114826546658785e-05, 1.0042876466822723e-05, 6.626886978987173e-06,
      4.342758016959359e-06, 2.826404139961336e-06, 1.826935439416417e-06,
      1.172843343117397e-06, 7.478113449978027e-07, 4.735725374732798e-07,
      2.9787248272133846e-07, 1.8609319443206175e-07, 1.1547645930055734e-07,
      7.117471589422507e-08, 4.357462045230393e-08, 2.649861429232971e-08,
      1.600664483874134e-08, 9.60441181497638e-09, 5.7245336130949956e-09,
      3.3893286667507067e-09, 1.9934100899941806e-09, 1.1646460583422373e-09,
      6.759435042001941e-10, 3.8971800517563627e-10, 2.232124441278612e-10,
      1.27004272419174e-10, 7.178828636653034e-11, 4.031145331279735e-11,
      2.2487619041829344e-11, 1.2462436797248393e-11, 6.861325871258659e-12,
      3.7528431107566266e-12, 2.039225600339778e-12, 1.1008248179772373e-12,
      5.903264456279235e-13, 3.1450496474651e-13, 1.6645611076888047e-13,
    };
    NC::VectD egrid(spec.size());
    for ( auto i : NC::ncrange(spec.size()) )
      egrid[i] = NC::VDOS::equidistantGridPoint(lower,binwidth,i);
    report( "real: Li2O (Li), order=129, vdoslux=2000",
            egrid, spec, 1e-6 );
  }
}

int main() {
  std::cout << "--- Synthetic: clean Gaussian bump, varying resolution ---" << std::endl;
  gaussianCase( "gaussian, coarse (24 pts)", 5.0, 0.0, 3.0, 24, 1e-6 );
  gaussianCase( "gaussian, medium (60 pts)", 5.0, 0.0, 3.0, 60, 1e-6 );
  gaussianCase( "gaussian, fine (400 pts)",  5.0, 0.0, 3.0, 400, 1e-6 );

  std::cout << "--- Synthetic: same coarse Gaussian, with noise injected"
               " right at the naive crossing points (the \"one-sidedness\""
               " concern) ---" << std::endl;
  gaussianCase( "gaussian, coarse, +-1e-6 rel noise at crossing", 5.0, 0.0, 3.0,
               24, 1e-6, 1e-6 );
  gaussianCase( "gaussian, coarse, +-1e-3 rel noise at crossing", 5.0, 0.0, 3.0,
               24, 1e-6, 1e-3 );

  std::cout << "--- Synthetic: hard step (discontinuity robustness) ---" << std::endl;
  stepCase( 0.0, 1.0, 40, 1e-6 );

  std::cout << "--- Real data captured from CI-sensitive cases ---" << std::endl;
  debyeOrder153Case();
  li2oOrder129Case();

  return 0;
}
