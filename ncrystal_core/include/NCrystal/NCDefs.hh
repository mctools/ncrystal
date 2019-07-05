#ifndef NCrystal_Defs_hh
#define NCrystal_Defs_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

//Definitions, includes, utility functions, etc., that should be available
//everywhere - including to NCrystal users.

#if __cplusplus >= 201103L
#  include <cstdint>
#else
#  include <stdint.h>
#endif
#include <cstddef>//for std::size_t
#include <limits>
#include <cmath>
#include <string>

#ifndef ncrystal_api_h
#  include "NCrystal/ncapi.h"
#endif
#ifndef NCrystal_NCException_hh
#  include "NCrystal/NCException.hh"
#endif
#ifndef NCrystal_NCMem_hh
#  include "NCrystal/NCMem.hh"
#endif

namespace NCrystal {

  //Utility functions for converting between neutron wavelength [Aa] and kinetic
  //energy [eV], and for providing infinity:
  NCRYSTAL_API double wl2ekin( double wl );     //cost: 1 branch + 1 division
  NCRYSTAL_API double ekin2wl( double ekin );   //cost: 1 branch + 1 division + 1 sqrt
  NCRYSTAL_API double ekin2wlsq( double ekin ); //cost: 1 branch + 1 division
  NCRYSTAL_API double ekin2wlsqinv( double ekin ); //cost: 1 multiplication
  const double infinity = std::numeric_limits<double>::infinity();

  //Math constants (avoid M_PI etc. for portability reasons). Note that these
  //are const and not static const on purpose (see
  //https://stackoverflow.com/questions/12042549/constants-only-header-file-c).
  const double kInf         = std::numeric_limits<double>::infinity()           ; // = infinity
  const double kSqrt2       = 1.41421356237309504880168872420969807856967188    ; // = sqrt(2)
  const double kPi          = 3.1415926535897932384626433832795028841971694     ; // = pi
  const double k2Pi         = 6.2831853071795864769252867665590057683943388     ; // = 2pi
  const double k4Pi         = 12.5663706143591729538505735331180115367886776    ; // = 4pi
  const double kPiHalf      = 1.5707963267948966192313216916397514420985847     ; // = pi/2
  const double kPiSq        = 9.86960440108935861883449099987615113531369941    ; // = pi^2
  const double k4PiSq       = 39.4784176043574344753379639995046045412547976    ; // = 4*pi^2
  const double kInv4PiSq    = 0.0253302959105844428609698658024319097260896937  ; // = 1/(4*pi^2)
  const double kInvPiSq     = 0.101321183642337771443879463209727638904358775   ; // = 1/pi^2
  const double kSqrtPi      = 1.77245385090551602729816748334114518279754946    ; // = sqrt(pi)
  const double kSqrt2Pi     = 2.50662827463100050241576528481104525300698674    ; // = sqrt(2pi)
  const double kInvSqrt2    = 0.707106781186547524400844362104849039284835938   ; // = 1/sqrt(2)
  const double kInvPi       = 0.318309886183790671537767526745028724068919291   ; // = 1/pi
  const double kInv2Pi      = 0.159154943091895335768883763372514362034459646   ; // = 1/(2pi)
  const double kInv4Pi      = 0.0795774715459476678844418816862571810172298229  ; // = 1/(4pi)
  const double kInvSqrt2Pi  = 0.398942280401432677939946059934381868475858631   ; // = 1/sqrt(2pi)
  const double kInvSqrtPi   = 0.564189583547756286948079451560772585844050629   ; // = 1/sqrt(pi)
  const double kSigma2FWHM  = 2.35482004503094938202313865291939927549477138    ; // = sqrt(8ln(2))
  const double kFWHM2Sigma  = 0.424660900144009521360751417051444809857570547   ; // = 1/sqrt(8ln(2))
  const double kDeg         = 0.0174532925199432957692369076848861271344287189  ; // = pi/180
  const double kArcMin      = 0.000290888208665721596153948461414768785573811981; // = pi/(180*60)
  const double kArcSec      = 0.00000484813681109535993589914102357947975956353302; // = pi/(180*3600)
  const double kToDeg       = 57.2957795130823208767981548141051703324054725    ; // = 180/pi
  const double kToArcMin    = 3437.74677078493925260788928884631021994432835    ; // = 180*60/pi
  const double kToArcSec    = 206264.806247096355156473357330778613196659701    ; // = 180*3600/pi

  //Use to mark variables as used and thus to avoid compiler warnings. This
  //should only be used in exceptional cases to workaround compiler issues, or
  //when a variable is otherwise used only in an assert:
  template<class T> inline void markused( const T& ) { }

  //Generic RNG interface:
  class NCRYSTAL_API RandomBase : public RCBase {
  public:
    virtual double generate() = 0;//generate numbers uniformly in [0,1[
  protected:
    RandomBase(){}
    virtual ~RandomBase();
  };

}

//Technically constants like M_PI from cmath/math.h are not dictated by the
//standards, and they are thus absent on some platforms (like windows with
//visual studio). For portability we either add them all (when
//NCRYSTAL_NO_CMATH_CONSTANTS is not defined) or remove them all.
#ifndef NCRYSTAL_NO_CMATH_CONSTANTS
#  ifndef M_E
#    define M_E        2.71828182845904523536  // e
#  endif
#  ifndef M_LOG2E
#    define M_LOG2E    1.44269504088896340736  //  log2(e)
#  endif
#  ifndef M_LOG10E
#    define M_LOG10E   0.434294481903251827651 //  log10(e)
#  endif
#  ifndef M_LN2
#    define M_LN2      0.693147180559945309417 //  ln(2)
#  endif
#  ifndef M_LN10
#    define M_LN10     2.30258509299404568402  //  ln(10)
#  endif
#  ifndef M_PI
#    define M_PI       3.14159265358979323846  //  pi
#  endif
#  ifndef M_PI_2
#    define M_PI_2     1.57079632679489661923  //  pi/2
#  endif
#  ifndef M_PI_4
#    define M_PI_4     0.785398163397448309616 //  pi/4
#  endif
#  ifndef M_1_PI
#    define M_1_PI     0.318309886183790671538 //  1/pi
#  endif
#  ifndef M_2_PI
#    define M_2_PI     0.636619772367581343076 //  2/pi
#  endif
#  ifndef M_2_SQRTPI
#    define M_2_SQRTPI 1.12837916709551257390  //  2/sqrt(pi)
#  endif
#  ifndef M_SQRT2
#    define M_SQRT2    1.41421356237309504880  //  sqrt(2)
#  endif
#  ifndef M_SQRT1_2
#    define M_SQRT1_2  0.707106781186547524401 //  1/sqrt(2)
#  endif
#else
#  ifdef M_E
#    undef M_E
#  endif
#  ifdef M_LOG2E
#    undef M_LOG2E
#  endif
#  ifdef M_LOG10E
#    undef M_LOG10E
#  endif
#  ifdef M_LN2
#    undef M_LN2
#  endif
#  ifdef M_LN10
#    undef M_LN10
#  endif
#  ifdef M_PI
#    undef M_PI
#  endif
#  ifdef M_PI_2
#    undef M_PI_2
#  endif
#  ifdef M_PI_4
#    undef M_PI_4
#  endif
#  ifdef M_1_PI
#    undef M_1_PI
#  endif
#  ifdef M_2_PI
#    undef M_2_PI
#  endif
#  ifdef M_2_SQRTPI
#    undef M_2_SQRTPI
#  endif
#  ifdef M_SQRT2
#    undef M_SQRT2
#  endif
#  ifdef M_SQRT1_2
#    undef M_SQRT1_2
#  endif
#endif


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {

  //The constant 8.1804... in the functions wl2ekin and ekin2wl is based on the
  //equation "h^2 * c^2 / (2.0*m)", using CODATA Internationally recommended
  //2014 values of the fundamental physical constants
  //(http://physics.nist.gov/cuu/Constants/Table/allascii.txt):
  //
  //  h = 4.135667662e-15 [Ev*s] <- Planck constant
  //  c = 299792458.0e10 [Aa/s]  <- speed of light in vacuum
  //  m = 939.5654133e6 [eV]     <- neutron mass energy equivalent
  //
  //  h^2 * c^2 / (2.0*m) = 0.081804209605330899

  inline double wl2ekin( double wl)
  {
    //Aangstrom to eV
    double wl2 = wl*wl;
    return wl2 ? ( 0.081804209605330899 / wl2 )  : std::numeric_limits<double>::infinity();
  }

  inline double ekin2wl( double ekin)
  {
    //eV to Aangstrom
    return ekin ? std::sqrt( 0.081804209605330899 / ekin ) : std::numeric_limits<double>::infinity();
  }

  inline double ekin2wlsq( double ekin)
  {
    //eV to Aangstrom^2
    return ekin ? 0.081804209605330899 / ekin : std::numeric_limits<double>::infinity();
  }

  inline double ekin2wlsqinv( double ekin)
  {
    //eV to 1/Aangstrom^2
    return ekin * 12.22430978582345950656;//constant is 1/0.081804209605330899
  }

  //Some obscure compilers like to complain about unused constants defined
  //through function calls rather than literal constants. Workaround is to mark
  //them as used in a global dummy function:
  inline void dummy_markused_global_constants() { markused(kInf); markused(infinity); }


}

#endif
