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

#include "NCrystal/internal/sanshardsphere/NCSANSSphScat.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <sstream>
namespace NC = NCrystal;

/////////////////////////////////////////////////////////////////////////////////////////
//
// Theory:
//
// Starting from [eq.1] from NCSANSUtils.hh:
//
//   sigma_sans(k) = C_sans * (Vp/k^2) * integral_0^2k{Q*P(Q)*S(Q)}dQ
//
// We define this diluted hard-sphere model by putting S(Q)=1 (dilute approx.) and
// P(Q) = [ 3*(sin(QR)-QR*cos(QR))/(QR)^3 ]^2, R being the sphere radius.
//
// Now with x=QR:
//
//      integral_0^2k{Q*P(Q)*S(Q)}dQ = (1/R)*integral_0^2k{QR*([ 3*(sin(QR)-QR*cos(QR))/(QR)^3 ]^2)}dQ
//                                   = (1/R^2)*integral_0^2kR{x*([ 3*(sin(x)-x*cos(x))/(x)^3 ]^2)}dx
//                                   = (9/R^2)*integral_0^2kR{[sin(x)-x*cos(x)]^2/x^5}dx
// so we can write:
//
//  sigma_sans(k) = C_sans*Vsphere*(9/((Rk)^2)*integral_0^2Rk{[sin(x)-x*cos(x)]^2/x^5}dx
//
// To sample a scattering, we must sample a value of x=QR according to the
// integrand above. To get the total cross section evaluated, we can carry out
// the integration analytically (done using SageMath for convenience):
//
//   integral_0^2Rk{[sin(x)-x*cos(x)]^2/x^5}dx = (32*(Rk)^4 - 8*(Rk)^2 + 4*Rk*sin(4*Rk) + cos(4*Rk) - 1)/ ( 128 * (Rk)^4 )
//
// Thus we get :
//
//  sigma_sans(k) = C_sans*(3pi/32) * R^3 *{ (Rk)^(-6)*[ 32*(Rk)^4 - 8*(Rk)^2 + 4*Rk*sin(4*Rk) + cos(4*Rk) - 1 ] }
//
// For low values of Rk, the this should be evaluated by a Taylor expansion for numerical stability.
//
//////////////////////////////////////

NC::SANSSphereScatter::SANSSphereScatter( SANSScaleFactor sfact_, sphere_radius r_ )
  : m_r(r_.value),
    m_scale([]( SANSScaleFactor sfact, sphere_radius rr )
    {
      double r = rr.value;
      if ( !(r>0) || !(r<1e9) )
        NCRYSTAL_THROW2(BadInput,"SANSSphereScatter radius value invalid or out of range: "<< r <<" Aa");
      constexpr double kkk = kPi * (3./32.);
      return sfact.dbl() * kkk * nccube(r);
    }(sfact_,r_))
{
}

NC::CrossSect NC::SANSSphereScatter::crossSectionIsotropic(CachePtr&, NeutronEnergy ekin ) const
{
  //Result is m_scale multiplied with rk^(-6)*[ 32*(rk)^4 - 8*(rk)^2 + 4*rk*sin(4*rk) + cos(4*rk) - 1 ].
  const double ksq = ekin2ksq( ekin.dbl() );
  const double rk = m_r * std::sqrt(ksq);//Unit is Aa * Aa^-1 = 1
  const double rksq = rk * rk;
  if ( rksq > 0.9 ) {
    if ( rksq > 1e20 )
      return CrossSect{ 0.0 };//guard against infinities
    const double fourrk=4.0*rk;
    auto sincos_4rk = sincos_fast(fourrk);
    const double rk4 = rksq * rksq;
    //Using StableSum in case some large terms cancel out:
    StableSum sum;
    sum.add(32.0*rk4);
    sum.add( - 8.0*rksq);
    sum.add(fourrk*sincos_4rk.sin);
    sum.add(sincos_4rk.cos);
    sum.add(-1.0);
    return CrossSect{ m_scale * sum.sum() / ( rksq * rk4 ) };
    //    return CrossSect{ m_scale * ( 32.0*rk4 - 8.0*rksq + fourrk*sin4rk + cos4rk - 1.0 ) / (rksq * rk4) };
  } else {
    constexpr double c0 = 28.44444444444444444444444;// 256/9
    constexpr double c1 = -11.37777777777777777777778;// -512/45
    constexpr double c2 = 2.600634920634920634920635;// 4096/1575
    constexpr double c3 = -0.3852792475014697236919459;// -16384/42525
    constexpr double c4 = 0.04002901272742542583812425;// 131072/3274425
    constexpr double c5 = -0.00307915482518657121831725;// -131072/42567525
    constexpr double c6 = 0.0001824684340851301462706519;// 1048576/5746615875
    constexpr double c7 = -0.00000858674983930024217744244;// -4194304/488462349375
    constexpr double c8 = 0.0000003286794196861336718638255;// 33554432/102088631019375
    constexpr double c9 = -1.043426729162329117028017e-8;// -67108864/6431583754220625
    constexpr double c10 = 2.791777201772118038870949e-10;// 536870912/1923043542511966875
    constexpr double c11 = -6.381205032621984088847884e-12;// -2147483648/336532619939594203125
    constexpr double c12 = 1.260484944715453647179829e-13;// 17179869184/136295711075535652265625
    constexpr double c13 = -2.173249904681816633068671e-15;// -8589934592/3952575621190533915703125
    const double y = rksq;
    return CrossSect{ m_scale * (c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*(c8+y*(c9+y*(c10+y*(c11+y*(c12+y*c13))))))))))))) };
  }
}

namespace NCRYSTAL_NAMESPACE {

  namespace {
    double hardSphereQRDensityFct( double x )
    {
      double xsq = x*x;
      if ( xsq < 0.8 ) {
        //evaluate with Taylor expansion:
        const double y = xsq;
        constexpr double c10 = 1.144019857775170297458039e-17;// 2/174822140228360625
        constexpr double c9 = -1.554826988521708722454335e-15;// -2/1286316750844125
        constexpr double c8 = 1.763173804983617691263216e-13;// 2/11343181224375
        constexpr double c7 = -1.637792556629227099884499e-11;// -8/488462349375
        constexpr double c6 = 1.218108213992987655539096e-9;// 1/820945125
        constexpr double c5 = -7.047626095245142864190483e-8;// -1/14189175
        constexpr double c4 = 0.000003053971307939561907815876;// 2/654885
        constexpr double c3 = -0.00009406231628453850676072898;// -4/42525
        constexpr double c2 = 0.001904761904761904761904762;// 1/525
        constexpr double c1 = -0.02222222222222222222222222;// -1/45
        constexpr double c0 = 0.1111111111111111111111111;// 1/9
        return x * (c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*(c8+y*(c9+y*c10))))))))));
      } else {
        auto sincos_x = sincos_fast(x);
        return ncsquare(sincos_x.sin-x*sincos_x.cos) / ( xsq*xsq*x );
      }
    }

    double sampleHardSphereQRUnbounded( RNG& rng ) {
      //Sample x in [0,infinity] according to hardSphereQRDensityFct(x) using an
      //overlay function which is a constant 0.105 below 4, and 1.05/x**3 above
      //4.
      constexpr double xthr = 4.0;
      constexpr double hflat = 0.105;
      constexpr double ktail = 1.05;
      constexpr double overlay_integral_below4 = hflat * xthr;
      constexpr double overlay_integral_above4 = ktail / 32.0;
      constexpr double overlay_probflat = overlay_integral_below4 / ( overlay_integral_below4 + overlay_integral_above4 );
      while (true) {
        double xgen;
        double overlayofx;
        if ( rng.generate() < overlay_probflat ) {
          xgen = rng.generate() * xthr;
          overlayofx = hflat;
        } else {
          xgen = xthr / std::sqrt(rng.generate());
          overlayofx = ktail / (xgen*xgen*xgen);
        }
        if ( overlayofx*rng.generate() < hardSphereQRDensityFct(xgen) )
          return xgen;
      }
    }

    double sampleHardSphereQR( RNG& rng, double qrmax ) {
      //According to notes above, we must sample x=Q*R over [0,xmax=2*k*R]
      //according to the density function f(x)=[sin(x)-x*cos(x)]^2/x^5.

      const double xmax = qrmax;

      //Strategy: We note that the first peak around xmax~=2 dominates
      //heavily. The second peak is located at xmax~=6 and the second and
      //higher peaks can be overlay by a function const/x^3. Thus, if xmax is
      //less than 4, we simply generate values in [0,xmax] and reject them
      //using a flat overlay function. If xmax is higher, we generate xmax in
      //[0,infinity] using a two-part overlay function (flat when x<4, ~1/x^3
      //for x>4). In case the generated x value is beyond xmax, we try
      //again. Since the majority of x values will be generated in the initial
      //peak, the need to throw a value away and try again should be rare.

      constexpr double xthr = 4.0;
      constexpr double xfirstpeak = 1.525526411927935;
      if ( xmax <= xthr ) {
        //below the first peak f(x) is strictly rising, so we pick the overlay
        //height accordingly (just using 0.105 everywhere would give terrible
        //acceptance rates for very long wavelength neutrons):
        const double hflat = ( xmax < xfirstpeak
                               ? 1.001*hardSphereQRDensityFct(xmax)
                               : 0.105 );
        while (true) {
          double x = rng.generate() * xmax;
          if ( hflat * rng.generate() <= hardSphereQRDensityFct(x) )
            return x;
        }
      } else {
        while (true) {
          double x = sampleHardSphereQRUnbounded(rng);
          if ( x <= xmax )
            return x;
        }
      }
    }
  }
}

NC::ScatterOutcomeIsotropic NC::SANSSphereScatter::sampleScatterIsotropic(CachePtr&, RNG& rng, NeutronEnergy ekin ) const
{
  //Implement as elastic scattering, i.e. sample Q over 0..2k according to the hard-sphere form factor.
  auto ksq = ekin2ksq( ekin.dbl() );
  if ( !(ksq>0.0) )
    return { ekin, CosineScatAngle{1.0} };//Do nothing

  const double qrmax = 2 * std::sqrt(ksq) * m_r;
  const double qr = sampleHardSphereQR( rng, qrmax );
  const double qsq = ncsquare( qr/m_r );
  double muval = ncclamp( 1.0 - qsq / ( 2.0 * ksq ), -1.0, 1.0 );
  return { ekin, CosineScatAngle{muval} };
}

NC::Optional<std::string> NC::SANSSphereScatter::specificJSONDescription() const
{
  std::ostringstream ss;
  CachePtr dummy;
  auto xs_at10Aa = crossSectionIsotropic(dummy, NeutronWavelength{10.0} );
  {
    std::ostringstream tmp;
    tmp << "radius="<<m_r<<"Aa;xs@10Aa="<<xs_at10Aa;
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "radius", m_r );
  streamJSONDictEntry( ss, "scale", m_scale  );
  streamJSONDictEntry( ss, "xsAt10Aa", xs_at10Aa.dbl(), JSONDictPos::LAST  );
  return ss.str();
}

struct NC::SANSSphereScatter::detail_direct_t
{
};

NC::SANSSphereScatter::SANSSphereScatter( detail_direct_t, double r, double scale)
  : m_r(r),
    m_scale(scale)
{
}

std::shared_ptr<NC::ProcImpl::Process> NC::SANSSphereScatter::createMerged( const Process& oraw,
                                                                            double scale_self,
                                                                            double scale_other ) const
{
  auto optr = dynamic_cast<const SANSSphereScatter*>(&oraw);
  return ( ( optr && m_r == optr->m_r )
           ? std::make_shared<SANSSphereScatter>( detail_direct_t(), m_r,
                                                  scale_self * m_scale + scale_other * optr->m_scale )
           : nullptr );
}
