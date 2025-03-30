#ifndef NCrystal_Defs_hh
#define NCrystal_Defs_hh

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

//Definitions, includes, utility functions, etc., that should be available
//everywhere - including to NCrystal users.

#include <cstdint>
#include <limits>
#include <cmath>
#include <utility>//std::move, std::forward
#include <tuple>//std::tuple, std::tie
//These we always include - simply because they otherwise have to be included in
//so many places:
#include <type_traits>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <map>
#include <set>
#include <mutex>
#include <atomic>
#include <ostream>
#include <cstring>
#if nc_cplusplus >= 202002L
#  include <compare>
#endif

#ifndef ncrystal_api_h
#  include "NCrystal/ncapi.h"
#endif
#ifndef NCrystal_NCException_hh
#  include "NCrystal/core/NCException.hh"
#endif
#ifndef NCrystal_NCMem_hh
#  include "NCrystal/core/NCMem.hh"
#endif

namespace NCRYSTAL_NAMESPACE {

  //Untyped utility functions for converting between neutron wavelength [Aa],
  //kinetic energy [eV], and wavenumber k=2pi/lambda [1/Aa]:
  NCRYSTAL_API constexpr double wl2ekin( double wl ) noexcept; //cost: 1 branch + 1 division + 1 mult
  NCRYSTAL_API double ekin2wl( double ekin ) noexcept; //cost: 1 branch + 1 division + 1 sqrt
  NCRYSTAL_API constexpr double constexpr_ekin2wl( double ekin ) noexcept;//expensive, use only compiletime
  NCRYSTAL_API constexpr double ekin2wlsq( double ekin ) noexcept; //cost: 1 branch + 1 division
  NCRYSTAL_API constexpr double ekin2wlsqinv( double ekin ) noexcept; //cost: 1 multiplication
  NCRYSTAL_API constexpr double wlsq2ekin( double wl ) noexcept; //cost: 1 branch + 1 division
  NCRYSTAL_API constexpr double ekin2ksq( double ekin ) noexcept; //cost: 1 multiplication
  NCRYSTAL_API double ekin2k( double ekin ) noexcept; //cost: 1 multiplication + 1 sqrt
  NCRYSTAL_API constexpr double constexpr_ekin2k( double ekin ) noexcept;//expensive, use only compiletime
  NCRYSTAL_API constexpr double ksq2ekin( double ksq ) noexcept; //cost: 1 multiplication
  NCRYSTAL_API constexpr double k2ekin( double k ) noexcept; //cost: 2 multiplications
  NCRYSTAL_API constexpr double wl2k( double wl ) noexcept; //cost: 1 branch + 1 division
  NCRYSTAL_API constexpr double wl2ksq( double wl ) noexcept; //cost: 1 branch + 1 division + 1 multiplication
  NCRYSTAL_API constexpr double k2wl( double wl ) noexcept; //cost: 1 branch + 1 division

  //Physics constants (more are in internal NCMath.hh header):
  constexpr double constant_boltzmann = 8.6173303e-5;  // eV/K
  constexpr double const_neutron_mass_amu = 1.00866491588; // [amu]
  constexpr double const_inv_neutron_mass_amu = 1.0/const_neutron_mass_amu; // [amu]
  constexpr double constant_dalton2kg =  1.660539040e-27; // amu to kg (source: NIST/CODATA 2018)

  //Various constexpr functions (not efficient for runtime usage!!):
  NCRYSTAL_API constexpr double constexpr_sqrt(double);
  NCRYSTAL_API constexpr double constexpr_abs(double);
  template <typename TVal>
  constexpr TVal ncconstexpr_max( TVal, TVal );
  template <typename TVal, typename ... Args>
  constexpr TVal ncconstexpr_max( TVal, TVal, Args ... );
  template<class TInt>
  constexpr TInt ncconstexpr_gcd( TInt , TInt );
  template<class TInt>
  constexpr TInt ncconstexpr_lcm( TInt, TInt );
  template <typename TInt, typename ... Args>
  constexpr TInt ncconstexpr_lcm( TInt, TInt, Args ... );
  template<class TInt>
  constexpr TInt ncconstexpr_roundToNextMultipleOf( TInt a, TInt b );//round a to next mult of b (a,b>0)

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
  constexpr unsigned supported_ncmat_format_version_min = 1;
  constexpr unsigned supported_ncmat_format_version_max = 7;

  //C++14 provides string_literals, allowing "hello"s as a shorthand for
  //std::string("hello",5). As C++11 does not support this, we implement our
  //own, "hello"_s variant.
  inline std::string operator ""_s(const char* c, std::size_t n) { return std::string(c,n); }

  //Use to mark variables as used and thus to avoid compiler warnings. This
  //should only be used in exceptional cases to workaround compiler issues, or
  //when a variable is otherwise used only in an assert:
  template<class T> inline void markused( const T& ) { }

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
    MoveOnly() = default;
    ~MoveOnly() = default;
    MoveOnly( const MoveOnly& ) = delete;
    MoveOnly& operator=( const MoveOnly& ) = delete;
    MoveOnly( MoveOnly&& ) = default;
    MoveOnly& operator=( MoveOnly&& ) = default;
  };

  class NCRYSTAL_API RNG : private MoveOnly  {
  public:
    //Random number stream base class with interfaces for generating
    //numbers. Actual implementations of specific streams should derive from RNG
    //which adds additional features concerning for instance state manipulation
    //and jumping. Code merely needing random numbers and not advanced stream
    //control should simply use the RNG interface to make this clear.
    virtual ~RNG();

    //Generate random number uniformly in interval (0.0,1.0].
    double operator()();
    double generate();

    //Generate integer uniformly in { 0, 1, ..., N-1 }:
    std::uint32_t generateInt( std::uint32_t N );
    std::uint64_t generateInt64( std::uint64_t N );

    //Coin flips (50% true, 50% false) and completely randomised bit patterns:
    virtual bool coinflip() = 0;
    virtual std::uint64_t generate64RndmBits() = 0;
    virtual std::uint32_t generate32RndmBits() = 0;

    //Possibly a more efficient version for getting many numbers at once:
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    virtual void generateMany( std::size_t n, double* tgt );
    virtual void generateRandomBits( std::size_t nbytes, std::uint8_t* data );
#endif

  protected:
    virtual double actualGenerate() = 0;//uniformly in (0,1]
  };

  struct NCRYSTAL_API UniqueIDValue {
    //type-safe unique id holder.
    std::uint64_t value;
#if nc_cplusplus >= 202002L
    auto operator<=>(const UniqueIDValue&) const = default;
#else
    constexpr bool operator<(UniqueIDValue const& o) const noexcept { return value < o.value; }
    constexpr bool operator==(UniqueIDValue const& o) const noexcept { return value == o.value; }
    constexpr bool operator!=(UniqueIDValue const& o) const noexcept { return value != o.value; }
#endif
  };

  class NCRYSTAL_API UniqueID : private MoveOnly {
  public:
    //Lock-free (on most platforms) and MT safe (on all platforms) unique ID.
    //Move-only to avoid two objects with same ID (another option would be to
    //generate a new ID for the copied-to object).

    ncconstexpr17 UniqueIDValue getUniqueID() const noexcept { return {m_uid}; }

    UniqueID();
    ~UniqueID() = default;
    UniqueID( const UniqueID& ) = delete;
    UniqueID& operator=( const UniqueID& ) = delete;
    UniqueID( UniqueID&& ) = default;
    UniqueID& operator=( UniqueID&& ) = default;
  private:
    std::uint64_t m_uid;
  };

  ////////////////////////////////////////////////////////////////////////////
  //Very simple optional class for C++11, similar to std::optional from C++17,
  //but with less features. It has copy semantics exactly when the wrapped type
  //has (implemented via a second boolean template parameter and a class
  //specialisation inlined elsewhere, so it works on all four major compilers):

  struct NCRYSTAL_API NullOptType {};
  constexpr const NullOptType NullOpt;

  template<class T, bool = (std::is_copy_constructible<T>::value
                            &&std::is_copy_assignable<T>::value)>
  class Optional {

    static_assert(std::is_nothrow_destructible<T>::value,
                  "Optional can only keep objects with noexcept destructors");
  public:

    static constexpr bool has_copy_semantics = false;
    using value_type = T;

    //Construct with value:
    ncconstexpr17 Optional( T&& ) noexcept(std::is_nothrow_move_constructible<T>::value);

    //Assign/construct value-less:
    constexpr Optional() noexcept;
    constexpr Optional( NullOptType ) noexcept;
    ncconstexpr17 Optional& operator=( NullOptType ) noexcept;

    //Clear value:
    void reset() noexcept;

    //Set value:
    Optional& operator=( T&& o ) noexcept(std::is_nothrow_move_constructible<T>::value);

    template<class TOther, class U = TOther,
             typename = typename std::enable_if<std::is_rvalue_reference<TOther>::value
                                                && std::is_constructible<T,U>::value
                                                && !std::is_base_of<Optional<T,false>,
                                typename std::remove_cv<typename std::remove_reference<U>::type>::type>::value>::type>
    Optional& operator=( TOther&& to ) noexcept(std::is_nothrow_constructible<T,TOther>::value)
    {
      reset();
      new(&m_value)T(std::forward<TOther>(to));
      m_hasValue = true;
      return *this;

    }

    template<typename... Args>
    void emplace( Args&& ... );

    //Access value:
    constexpr bool has_value() const noexcept { return m_hasValue; }
    ncconstexpr17 T& value() ncnoexceptndebug { nc_assert(m_hasValue); return m_value; }
    ncconstexpr17 const T& value() const ncnoexceptndebug { nc_assert(m_hasValue); return m_value; }

    template<class U>
    constexpr T value_or(U&& u) const;

    //move/assign/destruct (moved-from Optional is always without value). Note
    //that if T has copy-semantics, this class will get them as well!:
    Optional( Optional&& ) noexcept(std::is_nothrow_move_constructible<T>::value);
    Optional& operator=( Optional&& ) noexcept(std::is_nothrow_move_constructible<T>::value);
    ~Optional() noexcept;

    //Argh, operator= overloads occasionally gets swallowed by TOther. Provide
    //set methods as workaround:
    Optional& set( const Optional& o );
    Optional& set( Optional&& o );//leaves o without value

    bool operator<( const Optional& o) const;
    bool operator==( const Optional& o) const;
    bool operator!=( const Optional& o) const;

  private:
    union { char m_dummy; T m_value; };
    bool m_hasValue;
    friend class Optional<T,true>;
  };

  //Pimpl idiom helper (move-only, automatic lifetime mgmt, const-correctness, flexible constructor):
  template<typename T>
  class Pimpl final : private MoveOnly {
  private:
    T * m_ptr = nullptr;
  public:
    Pimpl() : m_ptr(new T) {}
    ~Pimpl() { delete m_ptr; }
    template<typename ...Args>
    Pimpl( Args&& ...args ) : m_ptr(new T(std::forward<Args>(args)... )) {}
    Pimpl( Pimpl&& o ) noexcept(std::is_nothrow_destructible<T>::value)
    {
      T* tmp = nullptr; std::swap(tmp,m_ptr); std::swap(m_ptr,o.m_ptr); delete tmp;
    }
    Pimpl& operator=( Pimpl&& o ) noexcept(std::is_nothrow_destructible<T>::value)
    {
      T* tmp = nullptr; std::swap(tmp,m_ptr); std::swap(m_ptr,o.m_ptr); delete tmp; return *this;
    }
    const T* operator->() const { return m_ptr; }
    T* operator->() { return m_ptr; }
    const T& operator*() const  { return *m_ptr; }
    T& operator*() { return *m_ptr; }
  };

  template<class T>
  class COWPimpl final {
  public:

    // Helper class implementing both pimpl and copy-on-write mechanics. This is
    // supposed to be both multi-thread safe and light-weight, but is best
    // suited for classes which potentially keep a large internal structure
    // which is modified seldomly and copied often.

    //Arguments of COWPimpl constructor will be passed along to constructor of T:
    template<typename ...Args> COWPimpl( Args&& ... );

    //Access internal T object through * or -> dereferencing. For modifications,
    //one must first retrieve a modification object through modify() (which
    //takes care of the MT-safe detaching as needed).

    class Modifier;
    ncnodiscard17 Modifier modify() { return Modifier{*this,true}; }
    constexpr const T* operator->() const noexcept { return &m_data->t; }
    constexpr const T& operator*() const noexcept { return m_data->t; }

    //If completely sure there is only one thread accessing the instance
    //(e.g. in the constructor of the class having a COWPimpl data member),
    //modifications can be done without locking first:
    ncnodiscard17 Modifier modifyWithoutLocking() { return Modifier{*this,false}; }

    //Fully supports copy and move operations, through cheap multi-thread safe
    //shallow copying (the overhead is a mutex lock).
    ~COWPimpl();
    COWPimpl( const COWPimpl& );
    COWPimpl& operator=( const COWPimpl& );
    COWPimpl( COWPimpl&& );
    COWPimpl& operator=( COWPimpl&& );

  private:
    struct Data;
    Data * m_data = nullptr;
    void releaseData();
  public:
    class Modifier final : private MoveOnly {
      Data * m_data = nullptr;
      std::mutex * m_mtx = nullptr;
    public:
      Modifier( COWPimpl& c, bool lock );
      ~Modifier();
      Modifier(Modifier&&o);
      Modifier& operator=(Modifier&&o);
      void reset();//Early release of lock and any refs to COWPimpl object (careful!)
      ncconstexpr17 T* operator->() noexcept { return &m_data->t; }
      ncconstexpr17 T& operator*() noexcept { return m_data->t; }
    };
  };

  //A few typedefs for very common types:
  typedef std::vector<double> VectD;
  typedef std::vector<std::string> VectS;
  typedef std::pair<double,double> PairDD;
  typedef std::pair<std::string,std::string> PairSS;

  //Structs which can be used in interfaces accepting cross-section values, to
  //make sure one does not accidentally mix up bound and free cross sections.

  //Convert uniformly randomised std::uint64_t (i.e. 64 independently randomised
  //bits) to a double precision floating point uniformly distributed over
  //(0,1]. This will map uniformly from 1.0 at input 0x0 to epsilon~=5.42e-20 at
  //std::uint64_max. The least significant bits in the input integer will also be
  //least significant for the output value, and due to precision issues the
  //lowest bits will only affect the result when the generated values are low
  //(e.g. the lowest 3 bits will only matter when the generated value is below
  //0.004):
  NCRYSTAL_API ncconstexpr17 double randUInt64ToFP01( std::uint64_t );

  struct NCRYSTAL_API no_init_t {};
  constexpr no_init_t no_init = no_init_t{};

  template <class TValue, std::size_t N>
  class FixedVector {
    //Base class for N-vectors (such as NeutronDirection or the internal Vector
    //class, which are implemented via the StronglyTypedFixedVector class and
    //CRTP, see below). This common base class provides some generic methods for
    //accessing data, printing, and interoperability with various types. This
    //interoperability includes std::array's and C arrays.
  public:
    using stdarray_type = std::array<TValue, N>;
    using carray_type = TValue[N];
    using value_type = TValue;
    using size_type = typename stdarray_type::size_type;
    using iterator = typename stdarray_type::iterator;
    using const_iterator = typename stdarray_type::const_iterator;
    using reverse_iterator = typename stdarray_type::reverse_iterator;
    using const_reverse_iterator = typename stdarray_type::const_reverse_iterator;

    static constexpr size_type size() noexcept { return N; }
    ncconstexpr17 TValue* data() noexcept { return m_data.data(); }
    constexpr const TValue* data() const noexcept { return m_data.data(); }

    //Default construct (will zero-initialise):
    explicit ncconstexpr17 FixedVector() noexcept;

    //Constructor which will NOT zero-initialise:
    explicit ncconstexpr17 FixedVector( no_init_t ) noexcept;

    //For convenience we provide special constructors and setters for 2-vectors
    //and 3-vectors (using one of these with a wrong length vector type will
    //result in a compile-time error):
    template<decltype(N) U = N, typename = typename std::enable_if<U==2>::type>
    ncconstexpr17 FixedVector(double xx, double yy) noexcept { m_data = { xx, yy }; }
    template<decltype(N) U = N, typename = typename std::enable_if<U==3>::type>
    ncconstexpr17 FixedVector(double xx, double yy, double zz) noexcept { m_data = { xx, yy, zz }; }
    template<decltype(N) U = N, typename = typename std::enable_if<U==2>::type>
    ncconstexpr17 void set(double xx, double yy) noexcept { m_data = { xx, yy }; }
    template<decltype(N) U = N, typename = typename std::enable_if<U==3>::type>
    ncconstexpr17 void set(double xx, double yy, double zz) noexcept { m_data = { xx, yy, zz }; }

    //Interoperability with std::array:
    explicit constexpr FixedVector(const stdarray_type& o) noexcept;
    explicit ncconstexpr17 operator stdarray_type&() noexcept { return m_data; }
    explicit constexpr operator const stdarray_type&() const noexcept { return m_data; }
    ncconstexpr17 stdarray_type& array() noexcept { return m_data; }
    constexpr const stdarray_type& array() const noexcept { return m_data; }

    //Interoperability with C arrays:
    explicit constexpr FixedVector(const carray_type& xyz) noexcept;
    explicit ncconstexpr17 operator carray_type&() noexcept { return rawArray(); }
    explicit ncconstexpr17 operator const carray_type&() const noexcept { return rawArray(); }
    ncconstexpr17 carray_type& rawArray() noexcept;
    ncconstexpr17 const carray_type& rawArray() const noexcept;
    ncconstexpr17 void applyTo( carray_type& ) const noexcept;

    //Access:
    ncconstexpr17 value_type operator[]( size_type i ) const ncnoexceptndebug { nc_assert(i<N); return m_data[i]; }
    ncconstexpr17 value_type& operator[]( size_type i ) ncnoexceptndebug { nc_assert(i<N); return m_data[i]; }
    ncconstexpr17 value_type at( size_type i ) const  { return m_data.at(i); }
    ncconstexpr17 value_type& at( size_type i ) { return m_data.at(i); }

    //Comparisons:
#if nc_cplusplus >= 202002L
    auto operator<=>(const FixedVector&) const = default;
#else
    //Can't be constexpr pre-C++20 due to lack of constexpr on std::array:
    bool operator<( const FixedVector& o ) const noexcept { return m_data < o.m_data; }
    bool operator>( const FixedVector& o ) const noexcept { return m_data > o.m_data; }
    bool operator<=( const FixedVector& o ) const noexcept { return m_data <= o.m_data; }
    bool operator>=( const FixedVector& o ) const noexcept { return m_data >= o.m_data; }
    bool operator==( const FixedVector& o ) const noexcept { return m_data == o.m_data; }
    bool operator!=( const FixedVector& o ) const noexcept { return m_data != o.m_data; }
#endif

    //Support standard iteration:
    ncconstexpr17 const_iterator begin() const noexcept { return m_data.cbegin(); }
    ncconstexpr17 const_iterator end() const noexcept { return m_data.cend(); }
    ncconstexpr17 const_iterator cbegin() const noexcept { return m_data.cbegin(); }
    ncconstexpr17 const_iterator cend() const noexcept { return m_data.cend(); }
    ncconstexpr17 iterator begin() noexcept { return m_data.begin(); }
    ncconstexpr17 iterator end() noexcept { return m_data.end(); }
    ncconstexpr17 const_reverse_iterator rbegin() const noexcept { return m_data.crbegin(); }
    ncconstexpr17 const_reverse_iterator rend() const noexcept { return m_data.crend(); }
    ncconstexpr17 const_reverse_iterator crbegin() const noexcept { return m_data.crbegin(); }
    ncconstexpr17 const_reverse_iterator crend() const noexcept { return m_data.crend(); }
    ncconstexpr17 reverse_iterator rbegin() noexcept { return m_data.rbegin(); }
    ncconstexpr17 reverse_iterator rend() noexcept { return m_data.rend(); }

  protected:
    stdarray_type m_data;
  };

  //Functions which do not care about the type of vector and which do not want
  //to be strongly typed (for instance vector-math functions), can simply accept
  //references to the FixedVector type as arguments. We typedef the most
  //important one for convenience:
  using ThreeVector = FixedVector<double,3>;

  //Strongly typed vectors (such as NeutronDirection or CrystalAxis) should be
  //implemented via the following CRTP base class. In addition to the
  //FixedVector interfaces, it also provides a mechanism for explitly (not
  //implicitly!) interpreting a vector as a diffent kind of strongly typed
  //vector:

  template <class Derived, class TValue, std::size_t N>
  class StronglyTypedFixedVector : public FixedVector<TValue,N> {

  public:
    using FixedVector<TValue,N>::FixedVector;

    //Interoperability with other compatible classes (those also implemented via
    //StronglyTypedFixedVector and not adding neither data members nor vtables):
    template<class TOther> ncconstexpr17 TOther& as() noexcept;
    template<class TOther> ncconstexpr17 const TOther& as() const noexcept;

  private:
    template<class TOther>
    static ncconstexpr17 void ensureCompatible() noexcept;
  };

  template <class TValue, std::size_t N>
  ncconstexpr17 std::ostream& operator<<(std::ostream& os,
                                         const FixedVector<TValue,N>& dir);


  //std::as_const is only available in C++17 and std::add_const_t only in C++14:
  template< class T >
  using nc_add_const_t = typename std::add_const<T>::type;
  template <class T>
  constexpr nc_add_const_t<T>& nc_as_const(T& t) noexcept { return t; }

  //Substitute for std::map::try_emplace which only available in C++17. For
  //simplicity we only implement the version which copies the key:
  template <class TMap, class ... Args>
  std::pair<typename TMap::iterator, bool> nc_map_try_emplace( TMap&, const typename TMap::key_type&, Args&&... );

  //For overriding existing keys (like m[k]=v but without default constructed value):
  template <class TMap, class ... Args>
  void nc_map_force_emplace( TMap&, const typename TMap::key_type&, Args&&... );

  //Convenience:
  using voidfct_t = std::function<void()>;

  //isOneOf: To test an argument against multiple values, e.g. write
  //isOneOf(a,"foo","bar","foobar) instead of
  //a=="foo"||a=="bar"||a=="foobar". The search is linear from left to right, so
  //is not meant for very large number of test cases:

  template <class T1>
  inline constexpr bool isOneOf(T1) {
    return false;
  }

  template <class T1, class T2, class... Ts>
  inline constexpr bool isOneOf(T1 needle, T2 haystack0, Ts... haystack_rest) {
    //linear search (so best for very small searches)
    return needle == haystack0 || isOneOf(needle, haystack_rest...);
  }

  template<class T>
  class ValRange {
  public:
    // Small helper class which allows python-like iteration over
    // integer range (with ncrange helper functions below):
    //
    //   for ( auto i : ncrange(i) )
    //         std::cout<<i<<std::endl;
    //
    using value_type = T;
    static_assert(std::is_integral<T>::value,
                  "ValRange and ncrange only works with integral types");
    class Iterator {
      value_type m_v;
      constexpr Iterator(value_type v) noexcept : m_v(v) {}
      friend class ValRange;
    public:
      constexpr const value_type* operator->() const noexcept { return &m_v; }
      constexpr const value_type& operator*() const noexcept { return m_v; }
      ncconstexpr17 Iterator operator++() noexcept { ++m_v; return Iterator(m_v); }
      ncconstexpr17 Iterator operator++(int) noexcept { Iterator r(m_v); ++m_v; return r; }
      constexpr bool operator==(const Iterator& o) const noexcept { return m_v == o.m_v; }
      constexpr bool operator!=(const Iterator& o) const noexcept { return m_v != o.m_v; }
      constexpr bool operator<(const Iterator& o) const noexcept { return m_v < o.m_v; }
    };
    constexpr ValRange(value_type n) noexcept : m_begin(0), m_end(n) {}
    constexpr ValRange(value_type l, value_type n) noexcept : m_begin(l), m_end(n)
    {
#if nc_cplusplus >= 201703L
      assert( l<=n );
#endif
    }
    constexpr Iterator begin() const noexcept { return m_begin; }
    constexpr Iterator end() const noexcept { return m_end; }
  private:
    Iterator m_begin, m_end;
  };

  template<class T> inline ValRange<T> ncrange(T n) { return ValRange<T>(n); }
  template<class T> inline ValRange<T> ncrange(T l, T n) { return ValRange<T>(l,n); }

  //Type-safe access to underlying int value of enum:
  template<class TEnum>
  constexpr typename std::underlying_type<TEnum>::type enumAsInt(TEnum const value) noexcept;

}

//We perform all mutex locking/unlocking with the following defines. This
//allows for easy dead-lock debugging by simply compiling with
//NCRYSTAL_DEBUG_LOCKS defined. Unique variable names are generated based on
//__LINE__, so do not put more than one of these statements on the same line.

//#define NCRYSTAL_DEBUG_LOCKS

#ifdef NCRYSTAL_LOCK_GUARD
#  undef NCRYSTAL_LOCK_GUARD
#endif
#ifdef NCRYSTAL_LOCK_MUTEX
#  undef NCRYSTAL_LOCK_MUTEX
#endif
#ifdef NCRYSTAL_UNLOCK_MUTEX
#  undef NCRYSTAL_UNLOCK_MUTEX
#endif
#ifdef NCRYSTAL_DEBUG_LOCKS
//NB: Not using NCrystal message infrastructure for this (it is anyway purely
//for developers):
#  define NCRYSTAL_LOCK_GUARD(MtxVariable) ::NCrystal::LockGuard ncrystal_join(nc_lock_guard_instance_line,__LINE__)(MtxVariable,__FILE__,__LINE__)
#  define NCRYSTAL_LOCK_MUTEX(MtxVariable) do { std::cout<<"NCrystal::Will lock mtx "<<(void*)&MtxVariable<<" ("<<__FILE__<<" : "<<__LINE__<<")"<<std::endl; MtxVariable.lock(); } while (0)
#  define NCRYSTAL_UNLOCK_MUTEX(MtxVariable) do { std::cout<<"NCrystal::Will unlock mtx "<<(void*)&MtxVariable<<" ("<<__FILE__<<" : "<<__LINE__<<")"<<std::endl; MtxVariable.unlock(); } while (0)
#  include <iostream>
#else
#  define NCRYSTAL_LOCK_GUARD(MtxVariable) ::NCrystal::LockGuard ncrystal_join(nc_lock_guard_instance_line,__LINE__)(MtxVariable)
#  define NCRYSTAL_LOCK_MUTEX(MtxVariable) do { MtxVariable.lock(); } while (0)
#  define NCRYSTAL_UNLOCK_MUTEX(MtxVariable) do { MtxVariable.unlock(); } while (0)
#endif

namespace NCRYSTAL_NAMESPACE {
  class NCRYSTAL_API LockGuard {
    using mutex_t = std::mutex;
    std::lock_guard<mutex_t> m_lg;
#ifdef NCRYSTAL_DEBUG_LOCKS
    std::string m_file; unsigned m_lineno; void * m_mtxaddr;
  public:
    LockGuard(mutex_t& mtx,const char * file, unsigned lineno)
      : m_lg([this,&mtx,&file,&lineno]()->mutex_t& {
        std::cout<<"NCrystal::LockGuard("<<(void*)this<<") Will lock mtx "<<(void*)&mtx<<" ("<<file<<" : "<<lineno<<")"<<std::endl; return mtx; }()),
        m_file(file), m_lineno(lineno), m_mtxaddr(&mtx) {}
    ~LockGuard() { std::cout<<"NCrystal::LockGuard("<<(void*)this<<") Will unlock mtx "<<m_mtxaddr<<" ("<<m_file<<" : "<<m_lineno<<")"<<std::endl; }
#else
  public:
    LockGuard(mutex_t& mtx) : m_lg(mtx) {}
#endif
  };
}

//For inserting code only in DEBUG builds:
#ifdef NCRYSTAL_DEBUGONLY
#  undef NCRYSTAL_DEBUGONLY
#endif
#ifndef NDEBUG
#  define NCRYSTAL_DEBUGONLY(x) x
#else
#  define NCRYSTAL_DEBUGONLY(x) do {} while(0)
#endif

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

namespace NCRYSTAL_NAMESPACE {

  inline constexpr double detail_sqrtNR(double x, double xc, double xp)
  {
    return xc == xp ? xc : detail_sqrtNR(x, 0.5 * (xc + x / xc), xc);
  }

  inline constexpr double constexpr_abs(double x)
  {
    return x < 0 ? -x : x;
  }

  inline constexpr bool constexpr_isinf(double x)
  {
    return constexpr_abs(x) >= std::numeric_limits<double>::infinity();
  }

  inline constexpr double constexpr_sqrt(double x)
  {
    //TODO: Mark consteval in c++20.
    return constexpr_isinf( x ) ? x : detail_sqrtNR(x, x, 0.);
  }

  template <typename TVal>
  inline constexpr TVal ncconstexpr_max(TVal a,TVal b)
  {
    return a > b ? a : b;
  }

  template <typename TVal, typename ... Args>
  inline constexpr TVal ncconstexpr_max(TVal a, TVal b, Args ... args)
  {
    return ncconstexpr_max<TVal>(ncconstexpr_max<TVal>(a,b),args...);
  }

  template<class TInt>
  inline constexpr TInt ncconstexpr_gcd( TInt a, TInt b )
  {
    return b ? ncconstexpr_gcd( b, a % b ) : a;
  }

  template<class TInt>
  inline constexpr TInt ncconstexpr_lcm( TInt a, TInt b )
  {
    return a * b / ncconstexpr_gcd(a,b);
  }

  template <typename TInt, typename ... Args>
  inline constexpr TInt ncconstexpr_lcm(TInt a, TInt b, Args ... args)
  {
    return ncconstexpr_lcm<TInt>(ncconstexpr_lcm<TInt>(a,b),args...);
  }

  template<class TInt>
  inline constexpr TInt ncconstexpr_roundToNextMultipleOf( TInt a, TInt b )
  {
    return ( a%b ? a + b - (a%b) : (a?a:b) );
  }

  template<class TEnum>
  inline constexpr typename std::underlying_type<TEnum>::type enumAsInt( TEnum value ) noexcept
  {
    return static_cast<typename std::underlying_type<TEnum>::type>(value);
  }

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

  inline constexpr double wl2ekin( double wl) noexcept
  {
    //angstrom to eV
    return wlsq2ekin( wl * wl );
  }

  inline double ekin2wl( double ekin) noexcept
  {
    //NB: std::sqrt is NOT constexpr!
    //eV to angstrom
    return ekin ? std::sqrt( 0.081804209605330899 / ekin ) : kInfinity;
  }

  inline constexpr double constexpr_ekin2wl( double ekin) noexcept
  {
    //eV to angstrom
    return ekin ? constexpr_sqrt( 0.081804209605330899 / ekin ) : kInfinity;
  }

  inline constexpr double wlsq2ekin( double wlsq ) noexcept
  {
    //angstrom^2 to eV
    return (wlsq ? ( 0.081804209605330899 / wlsq )  : kInfinity);
  }

  inline constexpr double ekin2wlsq( double ekin) noexcept
  {
    //eV to angstrom^2
    return ekin ? 0.081804209605330899 / ekin : kInfinity;
  }

  inline constexpr double ekin2wlsqinv( double ekin) noexcept
  {
    //eV to 1/angstrom^2
    return ekin * 12.22430978582345950656;//constant is 1/0.081804209605330899
  }

  namespace detail {
    constexpr double const_ekin2ksq_factor = k4PiSq * ekin2wlsqinv(1.0);
    constexpr double const_ksq2ekin_factor = 1.0 / const_ekin2ksq_factor;
    constexpr double constexpr_square_hlpr( double x ) noexcept { return x*x; }
  }

  inline constexpr double ekin2ksq( double ekin ) noexcept
  {
    return ekin * detail::const_ekin2ksq_factor;
  }

  inline double ekin2k( double ekin ) noexcept
  {
    //NB: std::sqrt is NOT constexpr!
    return std::sqrt( ekin * detail::const_ekin2ksq_factor );
  }

  inline constexpr double constexpr_ekin2k( double ekin ) noexcept
  {
    return constexpr_sqrt( ekin * detail::const_ekin2ksq_factor );
  }

  inline constexpr double ksq2ekin( double ksq ) noexcept
  {
    return detail::const_ksq2ekin_factor * ksq;
  }

  inline constexpr double k2ekin( double k ) noexcept
  {
    return detail::const_ksq2ekin_factor * ( k * k );
  }

  inline constexpr double wl2k( double wl ) noexcept
  {
    return wl ? k2Pi / wl : kInfinity;
  }

  inline constexpr double k2wl( double k ) noexcept
  {
    return k ? k2Pi / k : kInfinity;
  }

  inline constexpr double wl2ksq( double wl ) noexcept
  {
    return detail::constexpr_square_hlpr( wl2k(wl) );
  }

  //Some obscure compilers like to complain about unused constants defined
  //through function calls rather than literal constants. Workaround is to mark
  //them as used in a global dummy function:
  inline void dummy_markused_global_constants() { markused(kInfinity);  }

  template <class TValue, std::size_t N>
  inline ncconstexpr17 FixedVector<TValue,N>::FixedVector() noexcept
  {
    m_data.fill(0.0);
  }

  template <class TValue, std::size_t N>
  inline ncconstexpr17 FixedVector<TValue,N>::FixedVector( no_init_t ) noexcept
  {
  }

  template <class TValue, std::size_t N>
  inline constexpr FixedVector<TValue,N>::FixedVector(const stdarray_type& o) noexcept
    : m_data{o}
  {
  }

  template <class TValue, std::size_t N>
  inline constexpr FixedVector<TValue,N>::FixedVector(const carray_type& xyz) noexcept
    : m_data{xyz[0],xyz[1],xyz[2]}
  {
  }

  template <class TValue, std::size_t N>
  inline ncconstexpr17 void FixedVector<TValue,N>::applyTo( carray_type& a ) const noexcept
  {
    static_assert(sizeof(a)==N*sizeof(TValue),"");
    std::copy(std::begin(m_data), std::end(m_data), std::begin(a));
  }

  template <class Derived, class TValue, std::size_t N>
  template<class TOther>
  inline ncconstexpr17 void StronglyTypedFixedVector<Derived,TValue,N>::ensureCompatible() noexcept
  {
    //Verify our class and TOther are safe for the reinterpret_casts used in the as<..>() methods.
    static_assert( N == TOther::size(), "incompatible types (incompatible size)" );
    static_assert( std::is_same<TValue,typename TOther::value_type>::value,
                   "incompatible types (incompatible value type)" );
    //value_type
    using OurExpectedCRTPBase = StronglyTypedFixedVector<Derived,TValue,N>;
    using OtherExpectedCRTPBase = StronglyTypedFixedVector<TOther,TValue,N>;
    static_assert(std::is_base_of<OurExpectedCRTPBase, Derived>::value,
                  "incompatible types (source class not correctly CRTP derived)");
    static_assert(std::is_base_of<OtherExpectedCRTPBase, TOther>::value,
                  "incompatible types (target class not derived from StronglyTypedFixedVector or not correctly CRTP derived)");
    //Both this and the other class must be compatible with a simple
    //implementation. We want to make sure that neither class added vtables or
    //data members, which could make the reinterpret casts problematic.  The
    //checks might in principle fail to catch some wrong usage on machines with
    //extreme alignments, but in practice this is not really an issue.
    class Simple : public StronglyTypedFixedVector<Simple,TValue,N> {};
    static_assert(sizeof(Derived)==sizeof(Simple), "incompatible types (source class adds data-members or virtual functions)");
    static_assert(sizeof(TOther)==sizeof(Simple), "incompatible types (target class adds data-members or virtual functions)");
  }

  template <class Derived, class TValue, std::size_t N>
  template<class TOther>
  ncconstexpr17 TOther& StronglyTypedFixedVector<Derived,TValue,N>::as() noexcept
  {
    //Due to the common base FixedVector<TValue,N> class keeping all data, and
    //the various checks in ensureCompatible, this should be safe. There is a
    //small loophole in the standard for making the cast UB, but unless
    //compilers will go out of their way to make life miserable for us for no
    //reason at all, we should be safe. In any case, compilation or the very
    //first runtime usage should hopefully break down horribly and visibly if we
    //ever encounter such a compiler.
    ensureCompatible<TOther>();
    return *reinterpret_cast<TOther*>(this);
  }

  template <class Derived, class TValue, std::size_t N>
  template<class TOther>
  ncconstexpr17 const TOther& StronglyTypedFixedVector<Derived,TValue,N>::as() const noexcept
  {
    return ensureCompatible<TOther>(), *reinterpret_cast<const TOther*>(this);
  }

  template <class TValue, std::size_t N>
  inline ncconstexpr17 std::ostream& operator<< (std::ostream& os, const FixedVector<TValue,N>& dir)
  {
    if ( N==0 )
      return os << "{}";
    os << "{ " << dir[0];
    for (std::size_t i = 1; i < N; ++i)
      os << ", " << dir[i];
    return os << " }";
  }

  template <class TValue, std::size_t N>
  inline ncconstexpr17 typename FixedVector<TValue,N>::carray_type& FixedVector<TValue,N>::rawArray() noexcept
  {
    return *reinterpret_cast<carray_type*>(m_data.data());
  }

  template <class TValue, std::size_t N>
  inline ncconstexpr17 const typename FixedVector<TValue,N>::carray_type& FixedVector<TValue,N>::rawArray() const noexcept
  {
    return *reinterpret_cast<const carray_type*>(m_data.data());
  }

  template <typename T>
  class Optional<T, true> final : public Optional<T, false>
  {
    //Specialisation adding copy semantics.
    using base_t = Optional<T,false>;
    base_t& asBase() noexcept { return *static_cast<Optional<T,false>*>(this); }
    const base_t& asBase() const noexcept { return *static_cast<const Optional<T,false>*>(this); }
  public:
    static constexpr bool has_copy_semantics = true;
    static_assert(std::is_copy_constructible<T>::value,"");
    static_assert(std::is_copy_assignable<T>::value,"");

    using Optional<T, false>::Optional;

    Optional& operator=( Optional&& o ) noexcept(std::is_nothrow_move_constructible<T>::value)
    {
      asBase() = std::move(o);
      return *this;
    }

    ncconstexpr17 Optional& operator=( NullOptType ) noexcept
    {
      this->reset();
      return *this;
    }

    Optional& operator=( const T& o ) noexcept(std::is_nothrow_copy_constructible<T>::value)
    {
      this->reset();
      new(&asBase().m_value)T(o);
      asBase().m_hasValue = true;
      return *this;
    }

    template<class TOther, class U = TOther,
             typename = typename std::enable_if<
               std::is_constructible<T,U>::value &&
               !std::is_base_of<Optional<T,false>,
                                typename std::remove_cv<typename std::remove_reference<U>::type>::type>::value>::type>
    Optional& operator=( TOther&& to ) noexcept(std::is_nothrow_constructible<T,TOther>::value)
    {
      this->reset();
      new(&asBase().m_value)T(std::forward<TOther>(to));
      asBase().m_hasValue = true;
      return *this;
    }

    ncconstexpr17 Optional( const T& t ) noexcept(std::is_nothrow_copy_constructible<T>::value)
    {
      asBase().m_hasValue = true;
      new(&(asBase().m_value))T(t);
    }

    Optional( const Optional& o ) noexcept(std::is_nothrow_copy_constructible<T>::value)
      : Optional<T, false>()
      {
        if ( o.m_hasValue ) {
          new(&this->m_value) T(o.m_value);
          this->m_hasValue = true;
        } else {
          this->m_dummy = 0;
          this->m_hasValue = false;
        }
      }

    Optional& operator=( const Optional& o ) noexcept(std::is_nothrow_copy_constructible<T>::value)
    {
      if ( &o == this )
        return *this;
      this->reset();
      if ( o.m_hasValue ) {
        new(&this->m_value) T(o.m_value);
        this->m_hasValue = true;
      }
      return *this;
    }

    //Re-expose potentially hidden methods or those with changed return types:
    Optional& set( const Optional& o )
    {
      this->asBase().set(o.asBase());
      return *this;
    }

    Optional& set( Optional&& o )
    {
      this->asBase().set(std::move(o.asBase()));
      return *this;
    }

    constexpr Optional() noexcept : Optional<T,false>() {};
    ~Optional() noexcept {}
    Optional( Optional&& o ) noexcept(std::is_nothrow_move_constructible<T>::value)
      : Optional<T,false>(std::move(o.asBase()))
    {
    }

  };

  //NB: using placement new in most places to add new values to m_value, this
  //ensures that no existing object is assumed to exist in m_value which would
  //be UB.

  template<class T,bool bcopy>
  inline ncconstexpr17 Optional<T,bcopy>::Optional(T&& t) noexcept(std::is_nothrow_move_constructible<T>::value)
    : m_hasValue(true) { new(&m_value)T(std::move(t)); }

  template<class T,bool bcopy>
  inline Optional<T,bcopy>::~Optional() noexcept
  {
    if (m_hasValue)
      m_value.~T();
  }

  template<class T,bool bcopy>
  inline ncconstexpr17 Optional<T,bcopy>& Optional<T,bcopy>::operator=( NullOptType ) noexcept
  {
    reset();
    return *this;
  }

  template<class T,bool bcopy>
  inline Optional<T,bcopy>& Optional<T,bcopy>::operator=( T&& o ) noexcept(std::is_nothrow_move_constructible<T>::value)
  {
    reset();
    new(&m_value)T(std::move(o));
    m_hasValue = true;
    return *this;
  }

  template<class T,bool bcopy>
  inline Optional<T,bcopy>::Optional( Optional&& o ) noexcept(std::is_nothrow_move_constructible<T>::value)
  {
    if ( o.m_hasValue ) {
      new(&m_value)T(std::move(o.m_value));
      m_hasValue = true;
      o.reset();
    } else {
      m_dummy = 0;
      m_hasValue = false;
    }
  }

  template<class T,bool bcopy>
  inline Optional<T,bcopy>& Optional<T,bcopy>::set( const Optional& o )
  {
    if ( &o == this )
      return *this;
    reset();
    if ( o.m_hasValue ) {
      new(&m_value) T(o.m_value);
      m_hasValue = true;
    }
    return *this;
  }

  template<class T,bool bcopy>
  inline Optional<T,bcopy>& Optional<T,bcopy>::set( Optional&& o )
  {
    if ( &o == this )
      return *this;
    reset();
    if ( o.m_hasValue ) {
      new(&m_value)T(std::move(o.m_value));
      m_hasValue = true;
      o.reset();
    }
    return *this;
  }

  template<class T,bool bcopy>
  inline Optional<T,bcopy>& Optional<T,bcopy>::operator=( Optional&& o ) noexcept(std::is_nothrow_move_constructible<T>::value)
  {
    if ( &o == this )
      return *this;
    reset();
    if ( o.m_hasValue ) {
      new(&m_value)T(std::move(o.m_value));
      m_hasValue = true;
      o.reset();
    }
    return *this;
  }

  template<class T,bool bcopy>
  template<typename... Args>
  inline void Optional<T,bcopy>::emplace( Args&& ...args ) {
    reset();
    new(&m_value) T(std::forward<Args>(args)...);//wanted to use curly braces, but got unexpected narrowing errors with Optional<std::string>(5,'c')
    m_hasValue = true;
  }

  template<class T,bool bcopy>
  inline void Optional<T,bcopy>::reset() noexcept {
    if (m_hasValue) {
      m_value.~T();
      m_hasValue = false;
      m_dummy = 0;
    }
  }

  template<class T,bool bcopy>
  template<class U>
  inline constexpr T Optional<T,bcopy>::value_or(U&& u) const
  {
    return m_hasValue ? m_value : std::forward<U>(u);
  }

  template<class T,bool bcopy>
  inline constexpr Optional<T,bcopy>::Optional() noexcept : m_dummy(0), m_hasValue(false) {}

  template<class T,bool bcopy>
  inline constexpr Optional<T,bcopy>::Optional( NullOptType ) noexcept : m_dummy(0), m_hasValue(false) {}

  template<class T,bool bcopy>
  inline bool Optional<T,bcopy>::operator<( const Optional& o) const
  {
    if ( has_value() && o.has_value() )
      return value() < o.value();
    //std::optional: lhs is considered less than rhs if, and only if, rhs contains a value and lhs does not.
    return !has_value() && o.has_value();
  }

  template<class T,bool bcopy>
  inline bool Optional<T,bcopy>::operator==( const Optional& o) const
  {
    if ( has_value() && o.has_value() )
      return value() == o.value();
    return has_value() == o.has_value();
  }

  template<class T,bool bcopy>
  inline bool Optional<T,bcopy>::operator!=( const Optional& o) const
  {
    if ( has_value() && o.has_value() )
      return value() != o.value();
    return has_value() != o.has_value();
  }


  inline double RNG::operator()()
  {
    return generate();
  }

  inline double RNG::generate()
  {
    double r = actualGenerate();
#ifndef NDEBUG
    if ( ! ( r > 0.0 && r <= 1.0 ) )
      NCRYSTAL_THROW2(CalcError,"Random number stream generated number "<<r<<" which is outside (0.0,1.0]");
#endif
    return r;
  }

  inline std::uint32_t RNG::generateInt( std::uint32_t N )
  {
    constexpr std::uint32_t nmax = std::numeric_limits<std::uint32_t>::max();
    const std::uint32_t lim = nmax - nmax % N;//remove bias
    do {
      std::uint32_t r = generate32RndmBits();
      if ( r < lim )
        return r % N;
    } while(true);
  }

  inline std::uint64_t RNG::generateInt64( std::uint64_t N )
  {
    constexpr std::uint64_t nmax = std::numeric_limits<std::uint64_t>::max();
    const std::uint64_t lim = nmax - nmax % N;//remove bias
    do {
      std::uint64_t r = generate64RndmBits();
      if ( r < lim )
        return r % N;
    } while(true);
  }

  inline ncconstexpr17 double randUInt64ToFP01( std::uint64_t x )
  {
    //A note about the implementation: xoroshiro authors recommend: "(x >> 11) *
    //0x1.0p-53".  This selects all k*(2^53) for k=0..(2^53-1) with equal
    //probability. Here we invert r->1.0-r to get numbers in the interval
    //(0,1]. The lowest possible value selected by this is 1.1e-16 and the next
    //one is 2.2e-16. To smooth this out a bit further, and also reach values
    //closer to 0, we use the remaining 11 bits to pick a value uniformly inside
    //each "1.1e-16-width bin". The lowest achievable value is then 2^-53 -
    //0x7FF*2^-64 ~= 5.42e-20. Even though the 3 lowest bits in xoroshiro
    //generated 64bit integers supposedly have some statistical issues, this
    //should be acceptable (whenever the primary term generates numbers greater
    //than 0.004 the lowest three bits have no effect at all).
#if nc_cplusplus >= 201703L
    const double r1 = ( x >> 11 ) * 0x1.0p-53;
    const double r2 = ( x & 0x7FF ) * 0x1.0p-64;
#else
    const double r1 = ( x >> 11 ) * 1.1102230246251565404236316680908203125e-16;
    const double r2 = ( x & 0x7FF ) * 5.42101086242752217003726400434970855712890625e-20;
#endif
    return ( 1.0 - r1 ) - r2;
  }

  template<class T>
  struct COWPimpl<T>::Data : private NoCopyMove {
    template<typename ...Args> Data( Args&& ...args ) : t(std::forward<Args>(args)... ) {}
    T t;
    std::mutex mtx;
    std::uint_fast64_t refcount = 1;
  };

  template<class T>
  inline COWPimpl<T>::Modifier::Modifier(Modifier&&o)
  {
    std::swap(m_data,o.m_data);
    std::swap(m_mtx,o.m_mtx);
  }

  template<class T>
  inline typename COWPimpl<T>::Modifier& COWPimpl<T>::Modifier::operator=(Modifier&&o)
  {
    reset();
    std::swap(m_data,o.m_data);
    std::swap(m_mtx,o.m_mtx);
  }

  template<class T>
  inline void COWPimpl<T>::Modifier::reset()
  {
    if (m_mtx) {
      NCRYSTAL_UNLOCK_MUTEX((*m_mtx));
    }
    m_mtx = nullptr;
    m_data = nullptr;
  }

  template<class T>
  inline COWPimpl<T>::Modifier::~Modifier()
  {
    reset();
  }

  template<class T>
  inline COWPimpl<T>::Modifier::Modifier( COWPimpl& c, bool lock )
    : m_data( c.m_data )
  {
    nc_assert( m_data != nullptr );
    if (!lock)
      return;
    NCRYSTAL_LOCK_MUTEX(m_data->mtx);
    if ( m_data->refcount > 1 ) {
      //Detach:
      auto newdata = new Data( m_data->t );
      --( m_data->refcount );
      NCRYSTAL_UNLOCK_MUTEX(m_data->mtx);
      c.m_data = m_data = newdata;
      NCRYSTAL_LOCK_MUTEX(m_data->mtx);
    }
    m_mtx = & m_data->mtx;//the mutex we have locked
  }

  template<class T>
  inline void COWPimpl<T>::releaseData()
  {
    if (!m_data)
      return;

    //Check refcount while holding lock, but release lock before deleting
    //m_impl!
    Data * data_to_delete = nullptr;
    {
      NCRYSTAL_LOCK_GUARD(m_data->mtx);
      if ( m_data->refcount==1 ) {
        //we are the only object referring to m_data, and should remain so even
        //after releasing the lock, since we are the only remaining source of
        //m_data.
        std::swap( data_to_delete, m_data );
      } else {
        --( m_data->refcount );
      }
    }
    if (data_to_delete)
      delete data_to_delete;
  }

  template<class T>
  inline COWPimpl<T>::~COWPimpl()
  {
    releaseData();
  }

  template<class T>
  inline COWPimpl<T>::COWPimpl( const COWPimpl& o )
  {
    if ( o.m_data ) {
      NCRYSTAL_LOCK_GUARD(o.m_data->mtx);
      m_data = o.m_data;
      ++( m_data->refcount );
    }
  }

  template<class T>
  template<typename ...Args>
  inline COWPimpl<T>::COWPimpl( Args&& ...args )
    : m_data(new Data(std::forward<Args>(args)... ))
  {
  }

  template<class T>
  inline COWPimpl<T>& COWPimpl<T>::operator=( const COWPimpl& o )
  {
    if ( m_data == o.m_data )
      return *this;
    releaseData();
    if ( !o.m_data)
      return *this;//caller assigned moved-from object!
    NCRYSTAL_LOCK_GUARD(o.m_data->mtx);
    m_data = o.m_data;
    ++( m_data->refcount );
    return *this;
  }

  template<class T>
  COWPimpl<T>::COWPimpl( COWPimpl&& o )
  {
    std::swap(m_data,o.m_data);
  }

  template<class T>
  COWPimpl<T>& COWPimpl<T>::operator=( COWPimpl&& o )
  {
    if ( m_data == o.m_data )
      return *this;
    releaseData();
    std::swap(m_data,o.m_data);
    return *this;
  }

  template <class TMap, class ... Args>
  inline std::pair<typename TMap::iterator, bool> nc_map_try_emplace( TMap& themap, const typename TMap::key_type& key, Args&&... args )
  {
#if nc_cplusplus >= 201703L
    return themap.try_emplace(key, std::forward<Args>(args)... );
#else
    //slower workaround
    auto it = themap.find(key);
    if ( it != themap.end() )
      return std::pair<typename TMap::iterator, bool>( it, false );
    auto res = themap.emplace( typename TMap::value_type( key, typename TMap::mapped_type(std::forward<Args>(args)...) ) );
    nc_assert(res.second==true);
    return res;
#endif
  }

  template <class TMap, class ... Args>
  inline void nc_map_force_emplace( TMap& themap, const typename TMap::key_type& key, Args&&... args )
  {
    auto res = nc_map_try_emplace( themap, key, std::forward<Args>(args)... );
    if ( !res.second ) {
      //was not inserted, must override existing (note that try_emplace
      //crucially guarantees that args were not moved from already!):
      res.first->second = typename TMap::mapped_type( std::forward<Args>(args)... );
    }
  }
}

#ifndef NCrystal_Fmt_hh
#  include "NCrystal/core/NCFmt.hh"
#endif

#endif
