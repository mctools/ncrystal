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
// Verifies that NC::stable_exp/stable_expm1/stable_tanh/stable_sinh/         //
// stable_log                                                                //
// (NCMath.hh/.cc) return the exactly correctly-rounded double result for a   //
// broad set of representative and edge-case arguments.                      //
//                                                                            //
// The reference values below were generated externally with Python's        //
// mpmath library at 60 decimal digits of precision, then rounded to the     //
// nearest double and printed with 20 significant decimal digits -- enough   //
// for the C++ compiler to parse them back to bit-for-bit the same double    //
// mpmath rounded to (a double needs at most 17 significant decimal digits   //
// to round-trip uniquely). Deliberately plain decimal literals throughout   //
// (no hex floats, no ldexp), so this doubles as its own reproducibility     //
// test: if these literals do not parse to exactly the intended doubles on   //
// some platform, this test itself would fail first.                        //
//                                                                            //
// Checks the stable_ functions are within a small, fixed ULP tolerance of   //
// the reference (NOT bit-exact: composing individually correctly-rounded    //
// steps does not itself guarantee an exactly correctly-rounded final        //
// result -- see the comment on the stable_ functions themselves in          //
// NCMath.hh). Deliberately does NOT report how often the platform's plain   //
// std:: call itself disagrees with the reference: that count is expected    //
// to be platform-dependent (that is the whole premise here), so printing it //
// would make this test's own golden log non-reproducible.                  //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  //Maximum ULP distance we require the stable_ functions to be within of
  //the correctly-rounded reference (see the rationale on the stable_
  //functions themselves in NCMath.hh: a handful of ULP, not necessarily
  //exact, is what a single Newton-Raphson refinement step can guarantee
  //in practice). Comfortably above what was observed when this test was
  //written (at most 3 ULP, across exp/expm1/tanh/sinh):
  constexpr unsigned kMaxUlpTol = 8;

  //True if b is reachable from a in at most maxUlps steps of
  //std::nextafter (so also true for a==b, incl. +-0.0). Simple and
  //obviously correct, rather than relying on a signed-integer ULP-key
  //trick, since maxUlps is always small here:
  bool withinUlps( double a, double b, unsigned maxUlps )
  {
    if ( a == b )
      return true;
    if ( NC::ncisnan(a) || NC::ncisnan(b) )
      return false;
    double v = a;
    for ( unsigned i = 0; i < maxUlps; ++i ) {
      v = std::nextafter( v, b );
      if ( v == b )
        return true;
    }
    return false;
  }

  struct TestCase { double x, y; };

  template<class StableFct>
  void checkFct( const char* name, StableFct stablefct,
                 const TestCase* cases, std::size_t n )
  {
    unsigned nMismatch(0);
    for ( std::size_t i = 0; i < n; ++i ) {
      const double x = cases[i].x;
      const double yref = cases[i].y;
      const double ystable = stablefct(x);
      if ( !withinUlps( ystable, yref, kMaxUlpTol ) ) {
        ++nMismatch;
        std::cout << "ERROR: " << name << "(" << NC::fmt(x,"%.17g")
                  << ") = " << NC::fmt(ystable,"%.17g")
                  << " not within " << kMaxUlpTol << " ULP of expected "
                  << NC::fmt(yref,"%.17g") << std::endl;
      }
    }
    std::cout << name << ": " << n << " cases checked" << std::endl;
    REQUIRE( nMismatch == 0 );
  }

}

int main()
{
  //Reference (x,y) pairs generated with mpmath (60 decimal digits), y
  //rounded to the nearest double and printed with 20 significant digits:

  static const TestCase exp_cases[] = {
    { -700, 9.8596765437597707718e-305 },
    { -50, 1.9287498479639178206e-22 },
    { -20, 2.0611536224385578699e-9 },
    { -10, 0.000045399929762484854173 },
    { -5, 0.0067379469990854670008 },
    { -2, 0.13533528323661270232 },
    { -1, 0.36787944117144233402 },
    { -0.5, 0.60653065971263342426 },
    { -0.1, 0.90483741803595962860 },
    { -0.01, 0.99004983374916810668 },
    { -0.001, 0.99900049983337502191 },
    { -0.00000001, 0.99999999000000006077 },
    { 0, 1.0000000000000000000 },
    { 0.00000001, 1.0000000099999999392 },
    { 0.001, 1.0010005001667083846 },
    { 0.01, 1.0100501670841679491 },
    { 0.1, 1.1051709180756477124 },
    { 0.5, 1.6487212707001281942 },
    { 1, 2.7182818284590450908 },
    { 2, 7.3890560989306504069 },
    { 5, 148.41315910257659993 },
    { 10, 22026.465794806717895 },
    { 20, 485165195.40979027748 },
    { 50, 5.1847055285870720451e+21 },
    { 100, 2.6881171418161356094e+43 },
    { 300, 1.9424263952412558252e+130 },
    { 700, 1.0142320547350044928e+304 },
    { 1.234567890123456, 3.4368930843460052671 },
    { -3.7182818284590452, 0.024275641750774682587 },
    { 12.9, 400312.19132988248020 },
    { -0.070000000001, 0.93239381990501579711 },
  };

  static const TestCase expm1_cases[] = {
    { -20, -0.99999999793884641885 },
    { -10, -0.99995460007023750926 },
    { -5, -0.99326205300091452433 },
    { -2, -0.86466471676338729768 },
    { -1, -0.63212055882855766598 },
    { -0.5, -0.39346934028736657574 },
    { -0.1, -0.095162581964040426907 },
    { -0.05, -0.048770575499285991061 },
    { -0.01, -0.0099501662508319470984 },
    { -0.001, -0.00099950016662500823332 },
    { -0.00001, -9.9999500001666662992e-6 },
    { -0.0000001, -9.9999995000000162626e-8 },
    { -0.00000000001, -9.9999999999500002043e-12 },
    { 0, 0.0 },
    { 0.00000000001, 1.0000000000050000201e-11 },
    { 0.0000001, 1.0000000500000017243e-7 },
    { 0.00001, 0.000010000050000166666518 },
    { 0.001, 0.0010005001667083416629 },
    { 0.01, 0.010050167084168056680 },
    { 0.05, 0.051271096376024040409 },
    { 0.1, 0.10517091807564762918 },
    { 0.5, 0.64872127070012819416 },
    { 1, 1.7182818284590453128 },
    { 2, 6.3890560989306504069 },
    { 5, 147.41315910257659993 },
    { 10, 22025.465794806717895 },
    { 20, 485165194.40979027748 },
    { 50, 5.1847055285870720451e+21 },
    { 3.7182818284590452, 40.193555674716122894 },
    { -8.9061234567890123, -0.99986444369335825133 },
  };

  static const TestCase tanh_cases[] = {
    { -20, -1.0000000000000000000 },
    { -10, -0.99999999587769272669 },
    { -5, -0.99990920426259510823 },
    { -2, -0.96402758007581690336 },
    { -1, -0.76159415595576485103 },
    { -0.5, -0.46211715726000973659 },
    { -0.1, -0.099667994624955819072 },
    { -0.07, -0.069885890316428986302 },
    { -0.05, -0.049958374957879969624 },
    { -0.01, -0.0099996666799994603225 },
    { -0.001, -0.00099999966666680000296 },
    { 0, 0.0 },
    { 0.001, 0.00099999966666680000296 },
    { 0.01, 0.0099996666799994603225 },
    { 0.05, 0.049958374957879969624 },
    { 0.07, 0.069885890316428986302 },
    { 0.1, 0.099667994624955819072 },
    { 0.5, 0.46211715726000973659 },
    { 1, 0.76159415595576485103 },
    { 2, 0.96402758007581690336 },
    { 5, 0.99990920426259510823 },
    { 10, 0.99999999587769272669 },
    { 20, 1.0000000000000000000 },
    { 3.7182818284590452, 0.99882208059107291387 },
    { -6.9061234567890123, -0.99999799346406137435 },
  };

  static const TestCase sinh_cases[] = {
    { -20, -242582597.70489513874 },
    { -10, -11013.232874703393463 },
    { -5, -74.203210577788752289 },
    { -2, -3.6268604078470185748 },
    { -1, -1.1752011936438013784 },
    { -0.5, -0.52109530549374738495 },
    { -0.1, -0.10016675001984402804 },
    { -0.07, -0.070057180674134120202 },
    { -0.05, -0.050020835937655015735 },
    { -0.01, -0.010000166667500002757 },
    { -0.001, -0.0010000001666666750565 },
    { 0, 0.0 },
    { 0.001, 0.0010000001666666750565 },
    { 0.01, 0.010000166667500002757 },
    { 0.05, 0.050020835937655015735 },
    { 0.07, 0.070057180674134120202 },
    { 0.1, 0.10016675001984402804 },
    { 0.5, 0.52109530549374738495 },
    { 1, 1.1752011936438013784 },
    { 2, 3.6268604078470185748 },
    { 5, 74.203210577788752289 },
    { 10, 11013.232874703393463 },
    { 20, 242582597.70489513874 },
    { 3.7182818284590452, 20.584640016482673275 },
    { -6.9061234567890123, -499.18425343581833431 },
  };

  static const TestCase log_cases[] = {
    { 1e-308, -709.19620864216608425 },
    { 1e-300, -690.77552789821368151 },
    { 1e-50, -115.12925464970228973 },
    { 1e-20, -46.051701859880914469 },
    { 1e-10, -23.025850929940457235 },
    { 1e-05, -11.512925464970228617 },
    { 0.01, -4.6051701859880909140 },
    { 0.1, -2.3025850929940454570 },
    { 0.5, -0.69314718055994528623 },
    { 0.9, -0.10536051565782628137 },
    { 0.99, -0.010050335853501450589 },
    { 0.999, -0.0010005003335835343636 },
    { 0.9999999, -1.0000000494736473776e-7 },
    { 1, 0.0 },
    { 1.0000001, 9.9999995058387043653e-8 },
    { 1.001, 0.00099950033308342321053 },
    { 1.01, 0.0099503308531680920157 },
    { 1.1, 0.095310179804324934860 },
    { 1.5, 0.40546510810816438486 },
    { 2, 0.69314718055994528623 },
    { 5, 1.6094379124341002818 },
    { 10, 2.3025850929940459011 },
    { 100, 4.6051701859880918022 },
    { 100000, 11.512925464970228617 },
    { 10000000000, 23.025850929940457235 },
    { 1e20, 46.051701859880914469 },
    { 1e50, 115.12925464970228973 },
    { 1e100, 230.25850929940457945 },
    { 1e200, 460.51701859880915890 },
    { 1e300, 690.77552789821368151 },
    { 1.7976931348623157e308, 709.78271289338397310 },
    { 2.718281828459045, 1.0000000000000000000 },
    { 1.234567890123456, 0.21072102231565192909 },
    { 3.7182818284590452, 1.3132616875182228089 },
    { 5e-324, -744.44007192138121809 },
    { 2.2250738585072014e-308, -708.39641853226407875 },
  };

  checkFct( "stable_exp", NC::stable_exp,
           exp_cases, sizeof(exp_cases)/sizeof(exp_cases[0]) );
  checkFct( "stable_expm1", NC::stable_expm1,
           expm1_cases, sizeof(expm1_cases)/sizeof(expm1_cases[0]) );
  checkFct( "stable_tanh", NC::stable_tanh,
           tanh_cases, sizeof(tanh_cases)/sizeof(tanh_cases[0]) );
  checkFct( "stable_sinh", NC::stable_sinh,
           sinh_cases, sizeof(sinh_cases)/sizeof(sinh_cases[0]) );
  checkFct( "stable_log", NC::stable_log,
           log_cases, sizeof(log_cases)/sizeof(log_cases[0]) );

  std::cout << "All checks passed" << std::endl;
  return 0;
}
