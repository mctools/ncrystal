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

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCStableDbl.hh"
#include <iostream>

namespace NC=NCrystal;

#define REQUIRE(x) nc_assert_always(x)
#define REQUIREFLTEQ(x, y) nc_assert_always(::NC::floateq((x), (y)))

namespace {
  template <class Func>
  void testIntegration(const Func& f, double a, double b,
                       double expected_integral,
                       std::string descr,
                       bool expect_exact_simpsons ) {
    std::cout<<"\nTest integral of f(x)="<<descr<<" over ["<<a<<", "<<b<<"]:"<<std::endl;
    std::cout<<"  Exact                : "<<expected_integral<<std::endl;
    auto rd = [expected_integral]( double val ) { return (val/expected_integral-1.0); };
    auto rd_expect_exact = [&rd]( double val )
    {
      if ( rd(val) < 4e-15 )
        return 0.0;//cheating, but works for the purpose of stabilising the unit test
      return rd(val);
    };

    auto rd_simpsons = [&rd,&rd_expect_exact,expect_exact_simpsons]( double v )
    {
      return ( expect_exact_simpsons ? rd_expect_exact(v) : rd(v) );
    };
    double val;
    val = NC::integrateTrapezoidal(f,a,b,5);
    std::cout<<"  Trapezoidal(n=5)     : "<<val<<" (relative deviation: "<<rd(val)<<")"<<std::endl;
    val = NC::integrateTrapezoidal(f,a,b,10000);
    std::cout<<"  Trapezoidal(n=10000) : "<<val<<" (relative deviation: "<<rd(val)<<")"<<std::endl;
    val = NC::integrateSimpsons(f,a,b,2);
    std::cout<<"  Simpsons(n=2)        : "<<val<<" (relative deviation: "<<rd_simpsons(val)<<")"<<std::endl;
    val = NC::integrateSimpsons(f,a,b,10000);
    //rd_expect_exact always needed in the next line, to stabilise tests:
    std::cout<<"  Simpsons(n=10000)    : "<<val<<" (relative deviation: "<<rd_expect_exact(val)<<")"<<std::endl;

    val = NC::integrateTrapezoidal(f,a,b,17);
    std::cout<<"  Trapezoidal(n=17)     : "<<val<<" (relative deviation: "<<rd(val)<<")"<<std::endl;
    val = NC::integrateSimpsons(f,a,b,17);
    std::cout<<"  Simpsons(n=17)        : "<<val<<" (relative deviation: "<<rd_simpsons(val)<<")"<<std::endl;
    val = NC::integrateRomberg17(f,a,b);
    std::cout<<"  Romberg(n=17)        : "<<val<<" (relative deviation: "<<rd(val)<<")"<<std::endl;


    val = NC::integrateTrapezoidal(f,a,b,33);
    std::cout<<"  Trapezoidal(n=33)     : "<<val<<" (relative deviation: "<<rd(val)<<")"<<std::endl;
    val = NC::integrateSimpsons(f,a,b,33);
    std::cout<<"  Simpsons(n=33)        : "<<val<<" (relative deviation: "<<rd_simpsons(val)<<")"<<std::endl;
    val = NC::integrateRomberg33(f,a,b);
    std::cout<<"  Romberg(n=33)        : "<<val<<" (relative deviation: "<<rd(val)<<")"<<std::endl;
  }


  // Tests powspace for endpoint, size, monotonicity, and shape.
  // Covers ascending, descending, constant, linear, and nonlinear ranges.
  // Checks p values below, equal to, and above one.
  // Also checks small and large point counts and representative values.
  // Invalid inputs are not called because their behavior is not specified.
  // Requires the project VectD, powspace, floateq, and assertion utilities.

  static void checkPowspace(double a, double b, unsigned n, double p)
  {
    using NC::fmtg;
    std::cout<<"  -> powspace("<<fmtg(a)<<","<<fmtg(b)<<", "<<n<<","<<fmtg(p)
             <<")"<<std::endl;

    using namespace NCrystal;
    if ( !(a>0.0) || !(b>a) || !(p>0) || !(n>=2)
         || !std::isfinite(a)
         || !std::isfinite(b)
         || !std::isfinite(p) ) {
      std::cout<<"    -> invalid powspace input -> trigger LogicError in dbg blds"
               <<std::endl;
#ifndef NDEBUG //powspace checks for invalid input happens in debug builds only.
      bool threw = false;
      try {
        powspace(a, b, n, p);
      } catch (const NC::Error::LogicError&) {
        threw = true;
      }
      REQUIRE(threw);
#endif
      return;
    }
    const VectD x = powspace(a, b, n, p);
    REQUIRE(x.size() == n);
    REQUIREFLTEQ(x.front(), a);
    REQUIREFLTEQ(x.back(), b);

    for (unsigned i = 0; i < n; ++i) {
      const double t = static_cast<double>(i) /
        static_cast<double>(n - 1);
      const double want = a + (b - a) * std::pow(t, p);
      REQUIREFLTEQ(vectAt(x, i), want);
      REQUIRE(std::isfinite(vectAt(x, i)));
    }

    if (b >= a) {
      for (unsigned i = 1; i < n; ++i)
        REQUIRE(vectAt(x, i) >= vectAt(x, i - 1));
    } else {
      for (unsigned i = 1; i < n; ++i)
        REQUIRE(vectAt(x, i) <= vectAt(x, i - 1));
    }
  }

  void testPowspace()
  {
    using NC::VectD;
    using NC::vectAt;
    using NC::powspace;

    constexpr auto nan = std::numeric_limits<double>::quiet_NaN();

    checkPowspace(0.5, 1.0, 2, 0.25);
    checkPowspace(0.5, 1.0, 2, 0.0);
    checkPowspace(0.5, 1.0, 2, -0.1);

    checkPowspace(0.5, 1.0, 1, 0.25);
    checkPowspace(0.5, 1.0, 0, 0.25);

    checkPowspace(0.5, 1.0, 2, nan);
    checkPowspace(nan, 1.0, 2, 0.25);
    checkPowspace(0.5, nan, 2, 0.25);

    checkPowspace(0.5, 1.0, 2, NC::kInfinity);
    checkPowspace(NC::kInfinity, 1.0, 2, 0.25);
    checkPowspace(0.5, NC::kInfinity, 2, 0.25);

    checkPowspace(0.0, 1.0, 2, 0.25);
    checkPowspace(0.0, 1.0, 2, 1.0);
    checkPowspace(0.0, 1.0, 2, 100.0);

    checkPowspace(0.0, 1.0, 3, 0.5);
    checkPowspace(0.0, 1.0, 3, 1.0);
    checkPowspace(0.0, 1.0, 3, 2.0);

    checkPowspace(1e-10, 1.0, 2, 0.25);
    checkPowspace(1e-10, 1.0, 2, 1.0);
    checkPowspace(1e-10, 1.0, 2, 100.0);

    checkPowspace(1e-10, 1.0, 3, 0.5);
    checkPowspace(1e-10, 1.0, 3, 1.0);
    checkPowspace(1e-10, 1.0, 3, 2.0);

    checkPowspace(2.0, 10.0, 5, 2.0);
    checkPowspace(10.0, 2.0, 5, 2.0);
    checkPowspace(-7.0, 13.0, 17, 0.25);
    checkPowspace(-7.0, 13.0, 17, 1.0);
    checkPowspace(-7.0, 13.0, 17, 4.0);
    checkPowspace(3.25, 3.25, 101, 0.1);
    checkPowspace(3.25, 3.25, 101, 10.0);
    checkPowspace(-1.0e6, 1.0e6, 1001, 3.0);

    {
      const VectD x = powspace(2.0, 10.0, 5, 2.0);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 2.5);
      REQUIREFLTEQ(vectAt(x, 2), 4.0);
      REQUIREFLTEQ(vectAt(x, 3), 6.5);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 0.4);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 6.59479341998814);
      REQUIREFLTEQ(vectAt(x, 2), 8.06286626604159);
      REQUIREFLTEQ(vectAt(x, 3), 9.13040983191);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 0.5);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 6.0);
      REQUIREFLTEQ(vectAt(x, 2), 7.65685424949238);
      REQUIREFLTEQ(vectAt(x, 3), 8.92820323027551);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 1.0);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 4.0);
      REQUIREFLTEQ(vectAt(x, 2), 6.0);
      REQUIREFLTEQ(vectAt(x, 3), 8.0);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 1.5);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 3.0);
      REQUIREFLTEQ(vectAt(x, 2), 4.82842712474619);
      REQUIREFLTEQ(vectAt(x, 3), 7.19615242270663);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 2.0);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 2.5);
      REQUIREFLTEQ(vectAt(x, 2), 4.0);
      REQUIREFLTEQ(vectAt(x, 3), 6.5);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 2.1);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 2.43527528164806);
      REQUIREFLTEQ(vectAt(x, 2), 3.86606598307359);
      REQUIREFLTEQ(vectAt(x, 3), 6.37238676007);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 2.5);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 2.25);
      REQUIREFLTEQ(vectAt(x, 2), 3.41421356237310);
      REQUIREFLTEQ(vectAt(x, 3), 5.89711431702997);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 3.0);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 2.125);
      REQUIREFLTEQ(vectAt(x, 2), 3.0);
      REQUIREFLTEQ(vectAt(x, 3), 5.375);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 3.5);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 2.0625);
      REQUIREFLTEQ(vectAt(x, 2), 2.70710678118655);
      REQUIREFLTEQ(vectAt(x, 3), 4.92283573777248);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 4.0);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 2.03125);
      REQUIREFLTEQ(vectAt(x, 2), 2.5);
      REQUIREFLTEQ(vectAt(x, 3), 4.53125);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }

    {
      const VectD x = powspace(2.0, 10.0, 5, 8.0);
      REQUIREFLTEQ(vectAt(x, 0), 2.0);
      REQUIREFLTEQ(vectAt(x, 1), 2.0001220703125);
      REQUIREFLTEQ(vectAt(x, 2), 2.03125);
      REQUIREFLTEQ(vectAt(x, 3), 2.8009033203125);
      REQUIREFLTEQ(vectAt(x, 4), 10.0);
    }
  }

  void testFindRoot()
  {

    using NC::findRoot2;
    using NC::floateq;

    // Linear root.
    REQUIREFLTEQ(findRoot2([](double x) {
      return x - 2.5;
    }, 0.0, 5.0), 2.5);

    // Linear root with a negative result.
    REQUIREFLTEQ(findRoot2([](double x) {
      return x + 3.0;
    }, -5.0, 0.0), -3.0);

    // Positive quadratic root.
    REQUIREFLTEQ(findRoot2([](double x) {
      return x * x - 4.0;
    }, 0.0, 3.0), 2.0);

    // Negative quadratic root.
    REQUIREFLTEQ(findRoot2([](double x) {
      return x * x - 4.0;
    }, -3.0, 0.0), -2.0);

    // Exponential root.
    REQUIREFLTEQ(findRoot2([](double x) {
      return std::exp(x) - 3.0;
    }, 0.0, 2.0), std::log(3.0));

    // Trigonometric root.
    REQUIREFLTEQ(findRoot2([](double x) {
      return std::sin(x);
    }, 2.0, 4.0), std::acos(-1.0));

    // Root near the lower end of the interval.
    REQUIRE(floateq(findRoot2([](double x) {
      return x - 1.0e-8;
    }, -1.0, 1.0), 1.0e-8, 1.0e-10, 1.0e-12));

    // Root near the upper end of the interval.
    REQUIRE(floateq(findRoot2([](double x) {
      return x - (1.0 - 1.0e-8);
    }, 0.0, 1.0), 1.0 - 1.0e-8,
        1.0e-10, 1.0e-12));

    // Tighter requested tolerance and a nonlinear root.
    REQUIRE(floateq(findRoot2([](double x) {
      return x * x * x - 2.0;
    }, 0.0, 2.0, 1.0e-14),
        std::pow(2.0, 1.0 / 3.0), 1.0e-12, 1.0e-14));

    // A valid bracket for a non-monotonic function.
    REQUIREFLTEQ(findRoot2([](double x) {
      return (x - 0.25) * (x - 1.75);
    }, 0.0, 1.0), 0.25);

    // Invalid: the interval bounds are reversed.
    {
      bool threw = false;
      try {
        findRoot2([](double x) {
          return x - 1.0;
        }, 2.0, 0.0);
      } catch (const NC::Error::CalcError&) {
        threw = true;
      }
      REQUIRE(threw);
    }

    // Invalid: there is no sign change.
    {
      bool threw = false;
      try {
        findRoot2([](double x) {
          return x + 1.0;
        }, 0.0, 1.0);
      } catch (const NC::Error::CalcError&) {
        threw = true;
      }
      REQUIRE(threw);
    }

    // Both endpoints are positive, so this bracket must be rejected.
    // The buggy endpoint check instead accepts it and finds an interior
    // root (this test would have failed with code in NCrystal <= 4.4.6):
    {
      bool threw = false;
      try {
        findRoot2([](double x) {
          return (x - 0.1) * (x - 0.2);
        }, 0.0, 1.0);
      } catch (const NC::Error::CalcError&) {
        threw = true;
      }
      REQUIRE(threw);
    }

  }

  void testStableDbl() {
    using NC::StableDbl;
    const double u = std::numeric_limits<double>::epsilon();

    {
      StableDbl x;
      REQUIRE(x.value() == 0.0);
      x = 1.5;
      REQUIRE(x.value() == 1.5);
      StableDbl y(x);
      REQUIRE(y.value() == 1.5);
      y = 2.0;
      REQUIRE(y.value() == 2.0);
      y = x;
      REQUIRE(y.value() == 1.5);
      y += 1e-100;
      REQUIRE(y.value() == 1.5);
      y -= 1.5;
      REQUIRE(y.value() == 1e-100);
    }

    {
      StableDbl x(3.0);
      REQUIRE(x.value() == 3.0);
    }

    {
      auto x = StableDbl::fromState({3.0, u});
      REQUIRE(x.value() == 3.0);
    }

    {
      auto x = StableDbl::fromState({1.0, u});
      StableDbl y(-1.0);
      StableDbl z = x + y;
      REQUIRE(z.value() == u);
    }

    {
      auto x = StableDbl::fromState({1.0,u});
      StableDbl y(-1.0);
      StableDbl z = x - y;
      REQUIRE(z.value() == 2.0 + u);
    }

    {
      auto x = StableDbl::fromState({1.0,u});
      StableDbl z = -x;
      REQUIRE(z.value() == -1.0-u);
      REQUIRE((-z).value() == x.value());
    }

    {
      StableDbl x(1.5);
      StableDbl y(2.0);
      REQUIRE((x * y).value() == 3.0);
    }

    {
      StableDbl x(1.5);
      StableDbl y(2.0);
      REQUIRE((x * 2.0).value() == 3.0);
      REQUIRE((2.0 * x).value() == 3.0);
    }

    {
      auto x = StableDbl::fromState({1.0,u});
      auto y = StableDbl::fromState({1.0,u});
      StableDbl z = x * y;
      REQUIRE(z.value() == 1.0 + 2.0 * u);
    }

    // The exact result is u^2. Ordinary double arithmetic
    // rounds (1+u)^2 to 1+2u before the subtraction.
    {
      StableDbl a(1.0 + u);
      StableDbl b(1.0 + 2.0 * u);
      StableDbl z = a * a - b;
      REQUIRE(z.value() == u * u);
      REQUIRE(z.value() != 0.0);
    }

    {
      StableDbl x(1.25);
      StableDbl y(2.5);
      x += y;
      REQUIRE(x.value() == 3.75);
      x -= y;
      REQUIRE(x.value() == 1.25);
      x *= y;
      REQUIRE(x.value() == 3.125);
      x *= 2.0;
      REQUIRE(x.value() == 6.25);
    }

    {
      auto x = StableDbl::fromState({1.0,u});
      REQUIRE( x.state().first == 1.0 );
      REQUIRE( x.state().second == u );
      auto y = StableDbl::fromState({1.0,2.0*u});
      StableDbl z = (x + y) - y;
      REQUIRE(z.value() == x.value());
    }

    {
      StableDbl x(0.0);
      StableDbl y(-0.0);
      StableDbl z = x * y;
      REQUIRE(z.value() == 0.0);
    }
  }

}



int main() {

  testIntegration([](double x){return x*x+0.1;}, -1,2, 3.3, "x^2+1/10", true);
  testIntegration([](double x){return std::exp(-x*x)*x*x*x;}, 0,2,0.5-2.5*std::exp(-4.0), "exp(-x^2)*x^3 ", false);

  for( auto val : {0.5, 1.0, 2.0, 0.0, NC::kInfinity} ) {
    std::cout<<"NC::ekin2wl("<<val<<") = "<<NC::ekin2wl(val)<<std::endl;
    std::cout<<"NC::ekin2wlsq("<<val<<") = "<<NC::ekin2wlsq(val)<<std::endl;
    std::cout<<"NC::ekin2wlsqinv("<<val<<") = "<<NC::ekin2wlsqinv(val)<<std::endl;
    std::cout<<"NC::constexpr_ekin2wl("<<val<<") = "<<NC::constexpr_ekin2wl(val)<<std::endl;
    std::cout<<"NC::wl2ekin("<<val<<") = "<<NC::wl2ekin(val)<<std::endl;
    std::cout<<"NC::ekin2ksq("<<val<<") = "<<NC::ekin2ksq(val)<<std::endl;
    std::cout<<"NC::ekin2k("<<val<<") = "<<NC::ekin2k(val)<<std::endl;
    std::cout<<"NC::constexpr_ekin2k("<<val<<") = "<<NC::constexpr_ekin2k(val)<<std::endl;
    std::cout<<"NC::wl2k("<<val<<") = "<<NC::wl2k(val)<<std::endl;
    std::cout<<"NC::wl2ksq("<<val<<") = "<<NC::wl2ksq(val)<<std::endl;

    std::cout<<"NC::wlsq2ekin("<<val<<") = "<<NC::wlsq2ekin(val)<<std::endl;
    std::cout<<"NC::ksq2ekin("<<val<<") = "<<NC::ksq2ekin(val)<<std::endl;
    std::cout<<"NC::k2ekin("<<val<<") = "<<NC::k2ekin(val)<<std::endl;
    std::cout<<"NC::k2wl("<<val<<") = "<<NC::k2wl(val)<<std::endl;

  }

  std::cout<<"linspace test: ";
  for ( auto val : NC::linspace(0.1,2.0,20) )
    std::cout<<", "<<NC::fmt(val,"%.15g");
  std::cout<<std::endl;


  std::cout<<"fmt tests:"<<std::endl;
  NC::VectD fmtvals = {
    std::numeric_limits<double>::infinity(),
    -std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::denorm_min(),
    2*std::numeric_limits<double>::denorm_min(),
    0.5*std::numeric_limits<double>::denorm_min(),
    0.0,
    -0.0,
    -2*std::numeric_limits<double>::denorm_min(),
    -0.5*std::numeric_limits<double>::denorm_min(),
    std::nextafter(std::numeric_limits<double>::min(),9999),
    std::numeric_limits<double>::min(),
    std::nextafter(std::numeric_limits<double>::min(),0),
    -std::nextafter(std::numeric_limits<double>::min(),9999),
    -std::numeric_limits<double>::min(),
    -std::nextafter(std::numeric_limits<double>::min(),0),
    3.63069e-313,
    3.630691234567891234567e-313,
    1.234567891234567891232e-320,
    -3.63069e-313,
    -3.630691234567891234567e-313,
    -1.234567891234567891232e-320
  };
  for ( auto& e :
          NC::geomspace( std::numeric_limits<double>::denorm_min() * (1e100*0.1),
                         std::numeric_limits<double>::min() * (1e100*10),
                         107 ) ) {//107 is a prime
    fmtvals.push_back( e*1e-100 );
  }
  for ( auto val : fmtvals )
    std::cout<<"  fmt("<<val<<",\"%.6g\") = >>"
             <<std::flush<<NC::fmt(val,"%.6g")<<"<<"<<std::endl;
  for ( auto val : fmtvals ) {
    std::cout<<"  test fmt("<<val<<")"<<std::endl;
    std::ostringstream ss;
    ss<<NC::fmt(val);
    std::string tmp = ss.str();
    double backconv = NC::str2dbl(tmp);
    nc_assert_always( (NC::ncisnan(val)&&NC::ncisnan(backconv))
                      || val==backconv);
  }
  std::cout<<"powspace testing start..."<<std::endl;
  testPowspace();
  std::cout<<"powspace testing done."<<std::endl;
  std::cout<<"root finding testing start..."<<std::endl;
  testFindRoot();
  std::cout<<"root finding testing done."<<std::endl;
  std::cout<<"StableDbl testing start..."<<std::endl;
  testStableDbl();
  std::cout<<"StableDbl testing done."<<std::endl;
  return 0;
}
