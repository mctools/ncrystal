#ifndef NCrystal_Defs_hh
#define NCrystal_Defs_hh

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

//Definitions, includes, utility functions, etc., that should be available
//everywhere - including to NCrystal users.

#include <cstdint>
#include <cstddef>//for std::size_t
#include <limits>
#include <cmath>
#include <utility>//std::move, std::forward
#include <tuple>//std::tuple, std::tie
//These we always include - simply because they otherwise have to be included in
//so many places:
#include <string>
#include <vector>
#include <map>
#include <set>
#include <mutex>
#include <atomic>

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
  //energy [eV]:
  NCRYSTAL_API constexpr double wl2ekin( double wl );     //cost: 1 branch + 1 division + 1 mult
  NCRYSTAL_API constexpr double ekin2wl( double ekin );   //cost: 1 branch + 1 division + 1 sqrt
  NCRYSTAL_API constexpr double ekin2wlsq( double ekin ); //cost: 1 branch + 1 division
  NCRYSTAL_API constexpr double ekin2wlsqinv( double ekin ); //cost: 1 multiplication
  NCRYSTAL_API constexpr double wlsq2ekin( double wl );   //cost: 1 branch + 1 division

  //Math constants (avoid M_PI etc. for portability reasons).
  constexpr double kInfinity    = std::numeric_limits<double>::infinity()           ; // = infinity
  constexpr double kSqrt2       = 1.41421356237309504880168872420969807856967188    ; // = sqrt(2)
  constexpr double kPi          = 3.1415926535897932384626433832795028841971694     ; // = pi
  constexpr double k2Pi         = 6.2831853071795864769252867665590057683943388     ; // = 2pi
  constexpr double k4Pi         = 12.5663706143591729538505735331180115367886776    ; // = 4pi
  constexpr double kPiHalf      = 1.5707963267948966192313216916397514420985847     ; // = pi/2
  constexpr double kPiSq        = 9.86960440108935861883449099987615113531369941    ; // = pi^2
  constexpr double k4PiSq       = 39.4784176043574344753379639995046045412547976    ; // = 4*pi^2
  constexpr double kInv4PiSq    = 0.0253302959105844428609698658024319097260896937  ; // = 1/(4*pi^2)
  constexpr double kInvPiSq     = 0.101321183642337771443879463209727638904358775   ; // = 1/pi^2
  constexpr double kSqrtPi      = 1.77245385090551602729816748334114518279754946    ; // = sqrt(pi)
  constexpr double kSqrt2Pi     = 2.50662827463100050241576528481104525300698674    ; // = sqrt(2pi)
  constexpr double kInvSqrt2    = 0.707106781186547524400844362104849039284835938   ; // = 1/sqrt(2)
  constexpr double kInvPi       = 0.318309886183790671537767526745028724068919291   ; // = 1/pi
  constexpr double kInv2Pi      = 0.159154943091895335768883763372514362034459646   ; // = 1/(2pi)
  constexpr double kInv4Pi      = 0.0795774715459476678844418816862571810172298229  ; // = 1/(4pi)
  constexpr double kInvSqrt2Pi  = 0.398942280401432677939946059934381868475858631   ; // = 1/sqrt(2pi)
  constexpr double kInvSqrtPi   = 0.564189583547756286948079451560772585844050629   ; // = 1/sqrt(pi)
  constexpr double kSigma2FWHM  = 2.35482004503094938202313865291939927549477138    ; // = sqrt(8ln(2))
  constexpr double kFWHM2Sigma  = 0.424660900144009521360751417051444809857570547   ; // = 1/sqrt(8ln(2))
  constexpr double kDeg         = 0.0174532925199432957692369076848861271344287189  ; // = pi/180
  constexpr double kArcMin      = 0.000290888208665721596153948461414768785573811981; // = pi/(180*60)
  constexpr double kArcSec      = 0.00000484813681109535993589914102357947975956353302; // = pi/(180*3600)
  constexpr double kToDeg       = 57.2957795130823208767981548141051703324054725    ; // = 180/pi
  constexpr double kToArcMin    = 3437.74677078493925260788928884631021994432835    ; // = 180*60/pi
  constexpr double kToArcSec    = 206264.806247096355156473357330778613196659701    ; // = 180*3600/pi
  constexpr double kE           = 2.71828182845904523536028747135266249775724709    ; // = Euler's number, e
  constexpr double kInvE        = 0.367879441171442321595523770161460867445811131   ; // = 1/e

  //C++14 provides string_literals, allowing "hello"s as a shorthand for
  //std::string("hello",5). As C++11 does not support this, we implement our
  //own, "hello"_s variant.
  inline std::string operator "" _s(const char* c, std::size_t n) { return std::string(c,n); }

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

  //Disable copy/move semantics, by private inheritance from this class:
  class NCRYSTAL_API NoCopyMove {
  protected:
    constexpr NoCopyMove() = default;
    ~NoCopyMove() = default;
    NoCopyMove( const NoCopyMove& ) = delete;
    NoCopyMove& operator=( const NoCopyMove& ) = delete;
    NoCopyMove( NoCopyMove&& ) = delete;
    NoCopyMove& operator=( NoCopyMove&& ) = delete;
  };

  //Disable copy semantics, retain move semantics:
  class NCRYSTAL_API MoveOnly {
  protected:
    constexpr MoveOnly() = default;
    ~MoveOnly() = default;
    MoveOnly( const MoveOnly& ) = delete;
    MoveOnly& operator=( const MoveOnly& ) = delete;
    MoveOnly( MoveOnly&& ) = default;
    MoveOnly& operator=( MoveOnly&& ) = default;
  };

  struct NCRYSTAL_API UniqueIDValue {
    //type-safe unique id holder.
    uint64_t value;
#if __cplusplus >= 202002L
    auto operator<=>(const UniqueIDValue&) const = default;
#else
    bool operator<(UniqueIDValue const& o) const { return value < o.value; }
    bool operator==(UniqueIDValue const& o) const { return value == o.value; }
#endif
  };

  class NCRYSTAL_API UniqueID : private MoveOnly {
  public:
    //Lock-free (on most platforms) and MT safe (on all platforms) unique ID.
    //Move-only to avoid two objects with same ID (another option would be to
    //generate a new ID for the copied-to object).

    UniqueIDValue getUniqueID() const  { return {m_uid}; }

    UniqueID();
    ~UniqueID() = default;
    UniqueID( const UniqueID& ) = delete;
    UniqueID& operator=( const UniqueID& ) = delete;
    UniqueID( UniqueID&& ) = default;
    UniqueID& operator=( UniqueID&& ) = default;
  private:
    uint64_t m_uid;
  };

  //Pimpl idiom helper (move-only, automatic lifetime mgmt, const-correctness, flexible constructor):
  template<typename T>
  class NCRYSTAL_API Pimpl : private MoveOnly {
  private:
    std::unique_ptr<T> m_ptr;
  public:
    Pimpl() : m_ptr(std::make_unique<T>()) {}
    template<typename ...Args> Pimpl( Args&& ...args ) : m_ptr(std::make_unique<T>(std::forward<Args>(args)... )) {}
    ~Pimpl() = default;
    Pimpl( Pimpl&& ) = default;
    Pimpl& operator=( Pimpl&& ) = default;
    const T* operator->() const { return m_ptr.get(); }
    T* operator->() { return m_ptr.get(); }
    const T& operator*() const  { return *m_ptr.get(); }
    T& operator*() { return *m_ptr.get(); }
  };

  //A few typedefs for very common types:
  typedef std::vector<double> VectD;
  typedef std::vector<std::string> VectS;
  typedef std::pair<double,double> PairDD;

  //Structs which can be used in interfaces accepting cross-section values, to
  //make sure one does not accidentally mix up bound and free cross sections.
  struct NCRYSTAL_API SigmaBound { double val; };
  struct NCRYSTAL_API SigmaFree  { double val; };

}

//For inserting code only in DEBUG builds:
#ifndef NDEBUG
#  define NCRYSTAL_DEBUGONLY(x) x
#else
#  define NCRYSTAL_DEBUGONLY(x) do {} while(0)
#endif


//Technically constants like M_PI from cmath/math.h are not dictated by the
//standards, and they are thus absent on some platforms (like windows with
//visual studio). For portability we remove them all.
#ifdef M_E
#  undef M_E
#endif
#ifdef M_LOG2E
#  undef M_LOG2E
#endif
#ifdef M_LOG10E
#  undef M_LOG10E
#endif
#ifdef M_LN2
#  undef M_LN2
#endif
#ifdef M_LN10
#  undef M_LN10
#endif
#ifdef M_PI
#  undef M_PI
#endif
#ifdef M_PI_2
#  undef M_PI_2
#endif
#ifdef M_PI_4
#  undef M_PI_4
#endif
#ifdef M_1_PI
#  undef M_1_PI
#endif
#ifdef M_2_PI
#  undef M_2_PI
#endif
#ifdef M_2_SQRTPI
#  undef M_2_SQRTPI
#endif
#ifdef M_SQRT2
#  undef M_SQRT2
#endif
#ifdef M_SQRT1_2
#  undef M_SQRT1_2
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

  inline constexpr double wl2ekin( double wl)
  {
    //angstrom to eV
    return wlsq2ekin( wl * wl );
  }

  inline constexpr double ekin2wl( double ekin)
  {
    //eV to angstrom
    return ekin ? std::sqrt( 0.081804209605330899 / ekin ) : kInfinity;
  }

  inline constexpr double wlsq2ekin( double wlsq )
  {
    //angstrom^2 to eV
    return (wlsq ? ( 0.081804209605330899 / wlsq )  : kInfinity);
  }

  inline constexpr double ekin2wlsq( double ekin)
  {
    //eV to angstrom^2
    return ekin ? 0.081804209605330899 / ekin : kInfinity;
  }

  inline constexpr double ekin2wlsqinv( double ekin)
  {
    //eV to 1/angstrom^2
    return ekin * 12.22430978582345950656;//constant is 1/0.081804209605330899
  }

  //Some obscure compilers like to complain about unused constants defined
  //through function calls rather than literal constants. Workaround is to mark
  //them as used in a global dummy function:
  inline void dummy_markused_global_constants() { markused(kInfinity);  }


}

#endif
