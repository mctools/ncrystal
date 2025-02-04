#ifndef NCrystal_Types_hh
#define NCrystal_Types_hh

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

#ifndef NCrystal_SmallVector_hh
#  include "NCrystal/core/NCSmallVector.hh"
#endif
#ifndef NCrystal_ImmutBuf_hh
#  include "NCrystal/core/NCImmutBuf.hh"
#endif
#ifndef NCrystal_Variant_hh
#  include "NCrystal/core/NCVariant.hh"
#endif
#ifndef NCrystal_Fmt_hh
#  include "NCrystal/core/NCFmt.hh"
#endif

//////////////////////////////////////////////
//                                          //
// Common types used in various interfaces  //
//                                          //
//////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  //Base class used for MT-safe caching. Client code can provide polymorphic
  //CachePtr objects and pass them into various interfaces. This allows client
  //code to control the lifetime and number of caches (usually one per thread),
  //while giving NCrystal full control of what to actually put inside the
  //caches.
  class NCRYSTAL_API CacheBase : private NoCopyMove {
  public:
    CacheBase() = default;
    virtual ~CacheBase() = default;

    //Only required function allows callers to invalidate the cache (while
    //callers can otherwise just reset the CachePtr, it would require a new
    //malloc on next usage. This function should ideally be called through the
    //free-standing invalidateCache function:
    virtual void invalidateCache() = 0;
  };
  using CachePtr = std::unique_ptr<CacheBase>;

  //CachePtr Invalidate Using this function
  NCRYSTAL_API inline void invalidateCache( CachePtr& cp )
  {
    if ( cp != nullptr )
      cp->invalidateCache();
  }

  //Types representing various physical quantities. Technically they simply wrap
  //a floating point number, but encapsulating them like this ensures strong
  //type-safety in interfaces. It also means that users can selectively pass
  //either NeutronEnergy or NeutronWavelength objects, in all interfaces
  //accepting NeutronEnergy parameters. For simplicity and efficiency we
  //implement it via the CRTP technique.

  struct DoValidate_t {};
  constexpr DoValidate_t DoValidate = DoValidate_t{};

  template <class Derived, class TValue = double>
  class EncapsulatedValue {
    //CRTP base class for encapsulated values. Just as efficient as using raw
    //TValue objects, but with strong type safety.
  public:
    using value_type = TValue;

    //NB: For gcc 4.7 @ centos7, some constructors are defined in-class rather
    //than below with the other inlines!

    //Default constructor zero-initialises value:
    constexpr EncapsulatedValue() noexcept : m_value() {}

    //Construct from value (marked explicit to ensure strong type safety!):
    explicit constexpr EncapsulatedValue(TValue v) noexcept : m_value(v) {}

    //Constructors which will do the same but also invoke Derived::validate() fct:
    EncapsulatedValue( DoValidate_t, EncapsulatedValue );
    explicit EncapsulatedValue( DoValidate_t, TValue );

    //Accessing the wrapped value always require invocation of a named
    //method. Casts are forbidden (for strong type safety):
    operator TValue() const = delete;
    ncconstexpr17 TValue& get() noexcept;
    constexpr const TValue& get() const noexcept;
    ncconstexpr17 void set( TValue ) noexcept;
    ncconstexpr17 void set( const Derived& ) noexcept;

    //"dbl" is alias for "get" (can only be used when TValue is a double)::
    template<class U = TValue, typename = typename std::enable_if<std::is_same<U,double>::value>::type>
    ncconstexpr17 double& dbl() noexcept { return m_value; }
    template<class U = TValue, typename = typename std::enable_if<std::is_same<U,double>::value>::type>
    constexpr const double& dbl() const noexcept { return m_value; }

    //Can be used in boolean expression (usually this evaluates as "true" if value is non-zero).
    template<class U = TValue, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
    explicit constexpr operator bool() const noexcept { return bool(m_value); }

    //Supports comparisons:
#if nc_cplusplus >= 202002L
    auto operator<=>( const EncapsulatedValue& ) const = default;
#else
    constexpr bool operator<( const EncapsulatedValue& o ) const noexcept { return m_value < o.m_value; }
    constexpr bool operator>( const EncapsulatedValue& o ) const noexcept { return m_value > o.m_value; }
    constexpr bool operator<=( const EncapsulatedValue& o ) const noexcept { return m_value <= o.m_value; }
    constexpr bool operator>=( const EncapsulatedValue& o ) const noexcept { return m_value >= o.m_value; }
    constexpr bool operator==( const EncapsulatedValue& o ) const noexcept { return m_value == o.m_value; }
    constexpr bool operator!=( const EncapsulatedValue& o ) const noexcept { return m_value != o.m_value; }
#endif
    void stream( std::ostream& ) const;
  protected:
    TValue m_value;
  };

  //Values are printed with suitable trailing unit:
  template <class Derived, class TValue>
  std::ostream& operator<<(std::ostream&, const EncapsulatedValue<Derived,TValue>& );

  class NeutronEnergy;

  class NCRYSTAL_API NeutronWavelength final : public EncapsulatedValue<NeutronWavelength> {
    //Neutron wavelength in angstrom
  public:
    using EncapsulatedValue::EncapsulatedValue;
    NeutronWavelength() = default;//for VSCode
    static constexpr const char * unit() noexcept { return "Aa"; }
    //Automatic conversions from/to energy (not constexpr, since std::sqrt is not constexpr):
    NeutronWavelength(NeutronEnergy) noexcept;
    void validate() const;
    ncnodiscard17 constexpr NeutronEnergy energy() const noexcept;
    ncnodiscard17 constexpr double k() const noexcept;//Corresponding wavenumber (k) in units of 1/Aa
    ncnodiscard17 constexpr double ksq() const noexcept;//Corresponding squared wavenumber (k^2) in units of 1/Aa^2
  };

  class NCRYSTAL_API NeutronEnergy final : public EncapsulatedValue<NeutronEnergy> {
  public:
    //Neutron kinetic energy in electronvolt
    using EncapsulatedValue::EncapsulatedValue;
    NeutronEnergy() = default;//For VSCode
    static constexpr const char * unit() noexcept { return "eV"; }
    //Automatic conversions from/to wavelength:
    constexpr NeutronEnergy(NeutronWavelength) noexcept;
    void validate() const;
    ncnodiscard17 NeutronWavelength wavelength() const noexcept;
    ncnodiscard17 double k() const noexcept;//Corresponding wavenumber (k) in units of 1/Aa
    ncnodiscard17 constexpr double ksq() const noexcept;//Corresponding squared wavenumber (k^2) in units of 1/Aa^2
    //Constexpr versions (do not use in runtime!):
    ncnodiscard17 constexpr NeutronWavelength constexpr_wavelength() const noexcept;
    ncnodiscard17 constexpr double constexpr_k() const noexcept;
  };

  class NCRYSTAL_API CosineScatAngle final : public EncapsulatedValue<CosineScatAngle> {
  public:
    //Cosine of scattering angle (aka "mu")
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
    void validate() const;
  };

  class NCRYSTAL_API CrossSect final : public EncapsulatedValue<CrossSect> {
  public:
    //Cross section value (usually calculated per-atom for a particular neutron
    //state and material).
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "barn"; }
    void validate() const;
  };

  class NCRYSTAL_API Temperature final : public EncapsulatedValue<Temperature> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "K"; }
    constexpr double kT() const noexcept;
    void validate() const;

    struct temp_range{ double first, second; };//std::pair not allowed in constexpr
    static constexpr temp_range allowed_range = temp_range{0.001, 1000000.0};
  };

  class NCRYSTAL_API DebyeTemperature final : public EncapsulatedValue<DebyeTemperature> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "K"; }
    constexpr double kT() const noexcept;
    void validate() const;
  };

  class NCRYSTAL_API AtomMass final : public EncapsulatedValue<AtomMass> {
  public:
    //Mass of nuclei or atoms (in atomic mass units, a.k.a. Daltons).
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "u"; }
    constexpr double relativeToNeutronMass() const noexcept;
    void validate() const;
  };

  class SigmaFree;
  class NCRYSTAL_API SigmaBound final : public EncapsulatedValue<SigmaBound> {
  public:
    //Bound cross section value (usually a material-specific constant).
    using EncapsulatedValue::EncapsulatedValue;
    SigmaBound() = default;//for VSCode
    static constexpr const char * unit() noexcept { return "barn"; }
    constexpr SigmaBound( SigmaFree, AtomMass ) noexcept;
    constexpr SigmaFree free( AtomMass am ) const noexcept;
    void validate() const;
  };

  class NCRYSTAL_API SigmaFree final : public EncapsulatedValue<SigmaFree> {
  public:
    //Free-limit cross section value (usually a material-specific constant).
    using EncapsulatedValue::EncapsulatedValue;
    constexpr SigmaFree() noexcept : EncapsulatedValue{0.0} {}//line here for VSCode and gcc4.8 @ centos7
    static constexpr const char * unit() noexcept { return "barn"; }
    constexpr SigmaFree( SigmaBound, AtomMass ) noexcept;
    constexpr SigmaBound bound( AtomMass am ) const noexcept;
    void validate() const;
  };

  class NCRYSTAL_API SigmaAbsorption final : public EncapsulatedValue<SigmaAbsorption> {
  public:
    //Absorption cross section @ 2200m/s (usually a material-specific constant).
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "barn"; }
    void validate() const;
  };

  class NumberDensity;
  class NCRYSTAL_API Density final : public EncapsulatedValue<Density> {
  public:
    //Material density in g/cm3.
    using EncapsulatedValue::EncapsulatedValue;
    Density() = default;//for VSCode
    static constexpr const char * unit() noexcept { return "g/cm3"; }
    Density( NumberDensity, AtomMass averageAtomMass );
    void validate() const;
  };

  class NCRYSTAL_API NumberDensity final : public EncapsulatedValue<NumberDensity> {
  public:
    //Material number density in atoms per Aa^3.
    using EncapsulatedValue::EncapsulatedValue;
    NumberDensity() = default;//for VSCode
    static constexpr const char * unit() noexcept { return "atoms/Aa^3"; }
    NumberDensity( Density, AtomMass averageAtomMass );
    void validate() const;
  };

  struct NCRYSTAL_API DensityState {
    //Keep information about density overrides: Either an absolute override of
    //density the numberdensity value, or a scaling factor. The default state is
    //a scaling factor of 1.0, which is the same as no override.
    //NB: NCCfgManip.hh assumes underlying type is unsigned!
    enum class Type : unsigned { DENSITY, NUMBERDENSITY, SCALEFACTOR };
    Type type = Type::SCALEFACTOR;
    double value = 1.0;
    bool operator==( const DensityState& ) const;
    void validate() const;
    //For C++11:
    DensityState() = default;
    DensityState( Type tt, double vv ) : type(tt), value(vv) {}
  };
  NCRYSTAL_API std::ostream& operator<<( std::ostream&, const DensityState& );

  class NCRYSTAL_API Pressure final : public EncapsulatedValue<Pressure> {
  public:
    //Pressure in pascal (Pa).
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "Pa"; }
    static constexpr Pressure one_atm() noexcept { return Pressure{ 101325.0 }; }
    static constexpr double one_atm_raw() noexcept { return 101325.0; }
    void validate() const;
  };

  class SLDContrast;

  class NCRYSTAL_API ScatLenDensity final : public EncapsulatedValue<ScatLenDensity> {
  public:
    //Scattering length density (usually used for SANS physics). Note that this might be negative.
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "10^-6/Aa^2"; }
    void stream( std::ostream& ) const;
    void validate() const;
    ncnodiscard17 constexpr SLDContrast contrast(ScatLenDensity other) const noexcept;
  };

  class NCRYSTAL_API SLDContrast final : public EncapsulatedValue<SLDContrast> {
  public:
    //Scattering length density contrast for SANS physics, (also known as
    //"delta-rho"). This is an absolute difference between two scattering length
    //densities.
    using EncapsulatedValue::EncapsulatedValue;
    SLDContrast() = default;//for VSCode
    static constexpr const char * unit() noexcept { return "10^-6/Aa^2"; }
    void validate() const;
    constexpr SLDContrast(ScatLenDensity rho1, ScatLenDensity rho2) noexcept;
    ncnodiscard17 constexpr double valuePerAa2() const noexcept;//Value in units of 1/Aa^2
    ncnodiscard17 constexpr double valuePerCM2() const noexcept;//Value in units of 1/cm^2
  };

  class MosaicitySigma;

  class NCRYSTAL_API MosaicityFWHM final : public EncapsulatedValue<MosaicityFWHM> {
  public:
    //Mosaic spread (FWHM) in single crystals.
    using EncapsulatedValue::EncapsulatedValue;
    MosaicityFWHM() = default;//for VSCode
    constexpr MosaicityFWHM(MosaicitySigma) noexcept;
    ncnodiscard17 constexpr MosaicitySigma sigma() const noexcept;
    static constexpr const char * unit() noexcept { return "radians"; }
    void validate() const;
  };

  class NCRYSTAL_API MosaicitySigma final : public EncapsulatedValue<MosaicitySigma> {
  public:
    //Mosaic spread (std dev, i.e. Gaussian sigma parameter) in single crystals.
    using EncapsulatedValue::EncapsulatedValue;
    MosaicitySigma() = default;//for VSCode
    constexpr MosaicitySigma(MosaicityFWHM) noexcept;
    ncnodiscard17 constexpr MosaicityFWHM fwhm() const noexcept;
    static constexpr const char * unit() noexcept { return "radians"; }
    void validate() const;
  };

  class NCRYSTAL_API NeutronDirection final : public StronglyTypedFixedVector<NeutronDirection,double,3> {
    //Neutron direction vector (usually it should be normalised, but this is not
    //imposed by the class itself).
  public:
    using StronglyTypedFixedVector::StronglyTypedFixedVector;
  };

  class NCRYSTAL_API Length final : public EncapsulatedValue<Length> {
  public:
    //Length (in meters) used for macroscopic length scales like geometries and
    //particle positions. Although, this is the same type of unit as
    //NeutronWavelength, we keep the two uses in separate types with separate
    //default units (but we provide explicit conversion options):
    using EncapsulatedValue::EncapsulatedValue;
    Length() = default;//for VSCode
    static constexpr const char * unit() noexcept { return "m"; }
    void validate() const;
    //Automatic conversions from/to NeutronWavelengthwavelength:
    explicit constexpr Length(NeutronWavelength wl) noexcept;
    ncnodiscard17 constexpr NeutronWavelength as_wavelength() const noexcept;
    //Related numeric definitions of the chosen numerical values:
    static constexpr double kilometer = 1000.0;
    static constexpr double meter = 1.0;
    static constexpr double centimeter = 0.01;
    static constexpr double millimeter = 0.001;
    static constexpr double micrometer = 1e-6;
    static constexpr double nanometer = 1e-9;
    static constexpr double angstrom = 1e-10;
    static constexpr double cm = centimeter;
    static constexpr double mm = millimeter;
  };

  class NCRYSTAL_API LabAxis final : public StronglyTypedFixedVector<LabAxis,double,3> {
    //Axis or direction in lab-space.
  public:
    using StronglyTypedFixedVector::StronglyTypedFixedVector;
  };

  class NCRYSTAL_API CrystalAxis final : public StronglyTypedFixedVector<CrystalAxis,double,3> {
    //Axis or direction in crystal-space (i.e. in the coordinate system of the
    //unit cell of a single crystal).
  public:
    using StronglyTypedFixedVector::StronglyTypedFixedVector;
  };

  struct NCRYSTAL_API HKL {
    //Small structure containing hkl point. Comparisons use an inverse
    //lexicographic sort, to ensure the first point in a sorted list always have
    //a standard form where the first non-zero digit is positive.
    int h, k, l;
    HKL( no_init_t ) noexcept {}
    constexpr HKL() noexcept : h(0), k(0), l(0) {}
    constexpr HKL(int h_, int k_, int l_) noexcept;
    bool operator<(const HKL&o) const noexcept;
    bool operator==(const HKL&o) const noexcept;
    constexpr HKL flipped() const noexcept { return HKL{-h,-k,-l}; }
  };

  class NCRYSTAL_API HKLPoint final : public StronglyTypedFixedVector<HKLPoint,double,3> {
    //Point in (h,k,l)-space (unlike the HKL structure this is a floating pointer vector).
  public:
    using StronglyTypedFixedVector::StronglyTypedFixedVector;
  };

  struct NCRYSTAL_API OrientDir {
    //Orientation of crystal axis, combining a laboratory direction with either
    //a direction in the unit cell coordinate system or a HKL point (indicating
    //the plane-normal of that HKL plane).
    using CrysDir = Variant<CrystalAxis,HKLPoint,VariantAllowEmpty::No>;
    using LabDir = LabAxis;
    CrysDir crystal;
    LabDir lab;
    bool operator==(const OrientDir&o) const noexcept { return crystal == o.crystal && lab == o.lab; }
  };
  NCRYSTAL_API std::ostream& operator<<( std::ostream&, const OrientDir& );

  class NCRYSTAL_API LCAxis final : public StronglyTypedFixedVector<LCAxis,double,3> {
    //Layered crystal axis (e.g. (0,0,1) in pyrolytic graphite).
  public:
    using StronglyTypedFixedVector::StronglyTypedFixedVector;
  };

  //Composite value types used in some interfaces:
  struct NCRYSTAL_API EnergyDomain {
    NeutronEnergy elow = NeutronEnergy{ 0.0 };
    NeutronEnergy ehigh = NeutronEnergy{ kInfinity };
    constexpr EnergyDomain(NeutronEnergy,NeutronEnergy) ncnoexceptndebug;
    struct no_check_t{};
    constexpr EnergyDomain(no_check_t,NeutronEnergy,NeutronEnergy) noexcept;

    ncnodiscard17 constexpr bool contains(NeutronEnergy) const ncnoexceptndebug;//checks if in [elow,ehigh] (always false if isNull() is true!)
    ncnodiscard17 constexpr bool isNull() const noexcept { return elow.get() > std::numeric_limits<double>::max() || elow==ehigh; }//nb: std::isinf not constexpr
    static constexpr EnergyDomain null() noexcept { return EnergyDomain{ no_check_t{},
                                                                         NeutronEnergy{kInfinity},
                                                                         NeutronEnergy{kInfinity} }; }
  };

  struct NCRYSTAL_API ScatterOutcome {
    NeutronEnergy ekin;
    NeutronDirection direction;
  };

  struct NCRYSTAL_API ScatterOutcomeIsotropic {
    NeutronEnergy ekin;
    CosineScatAngle mu;
    constexpr static ScatterOutcomeIsotropic noScat( NeutronEnergy ef ) noexcept
    {
      return { ef, CosineScatAngle{1.0} };
    }
  };

  NCRYSTAL_API std::ostream& operator<<( std::ostream&, const ScatterOutcome& );
  NCRYSTAL_API std::ostream& operator<<( std::ostream&, const ScatterOutcomeIsotropic& );

  //Index identifying RNG streams:
  class NCRYSTAL_API RNGStreamIndex final : public EncapsulatedValue<RNGStreamIndex,uint64_t> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
    ncconstexpr17 void validate() const noexcept {}
  };

  //Serialised RNG stream state:
  class NCRYSTAL_API RNGStreamState final : public EncapsulatedValue<RNGStreamState,std::string> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
    ncconstexpr17 void validate() const noexcept {}
  };

  //Process properties:
  enum class MaterialType { Anisotropic, Isotropic };
  enum class ProcessType { Absorption, Scatter };

  NCRYSTAL_API std::ostream& operator<<(std::ostream&, MaterialType);
  NCRYSTAL_API std::ostream& operator<<(std::ostream&, ProcessType);

  class NCRYSTAL_API ThreadCount final : public EncapsulatedValue<ThreadCount,uint32_t> {
  public:
    //Thread counts (for specifying number of multiprocessing threads). Values
    //>= 9999 are used to indicate an automatic number of threads.
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
    bool indicatesAutoDetect() const noexcept;
    static ThreadCount auto_detect() noexcept { return ThreadCount{9999}; }
    ncconstexpr17 void validate() const noexcept {}
  };

  class NCRYSTAL_API DataSourceName {
    //Immutable string which is used to indicate the origins of a given data
    //object. This may or may not be a filename (e.g. "myfile.ncmat",
    //"<unknown>", ...), and should in any case ONLY be used to create more
    //meaningful output messages. In no case should an algorithm start to
    //modify its behaviour depending on the content of such names.
  public:
    DataSourceName();//empty string
    DataSourceName( std::string );
    DataSourceName& operator=( std::string );
    bool operator==( const DataSourceName& o ) const { return *m_str == *o.m_str; }
    const std::string& str() const { return m_str; }
  private:
    shared_obj<const std::string> m_str;
  };

  NCRYSTAL_API std::ostream& operator<<( std::ostream&, const DataSourceName& );

  struct NCRYSTAL_API UCNMode {
    static constexpr NeutronEnergy default_threshold() noexcept { return NeutronEnergy{ 300e-9 }; }
    enum class Mode { Refine, Remove, Only };
    Mode mode = Mode::Refine;
    NeutronEnergy threshold = default_threshold();
    bool operator==( const UCNMode& ) const;
    UCNMode( Mode mde = Mode::Refine,
             NeutronEnergy thr = default_threshold() )
    : mode(mde), threshold(thr) {}//constructor needed by clang
  };

  NCRYSTAL_API std::ostream& operator<<( std::ostream&, const UCNMode& );

  //For message passing and output control:
  enum class MsgType : unsigned { Info = 0, Warning = 1, RawOutput = 2 };

  namespace Cfg {

    namespace detail {

      //Internal types related to configuration and needed in public headers for
      //technical reasons.

      namespace varbuf_calc {
        //Typically VarBuf will end up with alignment of 8, object size of 32,
        //and a local buffer of 27. However, for portability and robustness we
        //express the requirements for all supported types and combine
        //appropriately:
        static constexpr auto ValDbl_buf_align = alignof(double);
        static constexpr auto ValDbl_buf_minsize = sizeof(double) + 16;
        static constexpr auto ValInt_buf_align = alignof(int64_t);
        static constexpr auto ValInt_buf_minsize = sizeof(int64_t);
        static constexpr auto ValStr_buf_align = alignof(char);
        static constexpr auto ValStr_buf_minsize = 24;
        static constexpr auto ValBool_buf_align = alignof(char);
        static constexpr auto ValBool_buf_minsize = 1;
        static constexpr auto ValVector_buf_align = alignof(ThreeVector);
        static constexpr auto ValVector_buf_minsize = sizeof(ThreeVector);
        static constexpr auto ValOrientDir_buf_align = alignof(double);
        static constexpr auto ValOrientDir_buf_minsize = 1;//too large, just keep in remote buffer.
        static constexpr auto buf_align = ncconstexpr_lcm<std::size_t>( ValDbl_buf_align, ValInt_buf_align,
                                                                        ValStr_buf_align, ValBool_buf_align,
                                                                        ValVector_buf_align, ValOrientDir_buf_align );
        static constexpr auto buf_minsize = ncconstexpr_max<std::size_t>( ValDbl_buf_minsize, ValInt_buf_minsize,
                                                                          ValStr_buf_minsize, ValBool_buf_minsize,
                                                                          ValVector_buf_minsize, ValOrientDir_buf_minsize );
      }
      enum class VarId : std::uint32_t;//only fwd decl here
      using VarBuf = ImmutableBuffer<varbuf_calc::buf_minsize,varbuf_calc::buf_align,VarId>;
      using VarBufVector = SmallVector_IC<VarBuf,7,SVMode::FASTACCESS>;
    }

    class CfgManip;
    class NCRYSTAL_API CfgData {
      //Class encapsulating a VarBufVector, ensuring that only CfgManip methods
      //can access it. This allows CfgData objects in public headers, leaving a
      //whole bunch of templated code private.
    private:
      friend class CfgManip;
      detail::VarBufVector m_data;
      ncconstexpr17 const detail::VarBufVector& operator()() const noexcept { return m_data; }
      ncconstexpr17 detail::VarBufVector& operator()() noexcept { return m_data; }
    };
  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  inline std::ostream& operator<<(std::ostream& os, MaterialType mt)
  {
    return os << (mt==MaterialType::Isotropic?"Isotropic":"Anisotropic");
  }

  inline std::ostream& operator<<(std::ostream& os, ProcessType pt)
  {
    return os << (pt==ProcessType::Scatter?"Scatter":"Absorption");
  }

  inline std::ostream& operator<<( std::ostream& os, const ScatterOutcome& outcome )
  {
    return os << "ScatterOutcomeIsotropic(E="<<outcome.ekin<<", dir="<<outcome.direction<<')';
  }

  inline std::ostream& operator<<(std::ostream& os, const ScatterOutcomeIsotropic& outcome)
  {
    return os << "ScatterOutcomeIsotropic(E="<<outcome.ekin<<", mu="<<outcome.mu<<')';
  }

  //Next two are defined in-class due to gcc 4.8 @ centos7:

  // template <class Derived, class TValue>
  // inline constexpr EncapsulatedValue<Derived,TValue>::EncapsulatedValue() noexcept
  //   : m_value()
  // {
  // }

  // template <class Derived, class TValue>
  // inline constexpr EncapsulatedValue<Derived,TValue>::EncapsulatedValue(TValue val) noexcept
  //   : m_value(val)
  // {
  // }

  template <class Derived, class TValue>
  inline EncapsulatedValue<Derived,TValue>::EncapsulatedValue( DoValidate_t, EncapsulatedValue o )
    : m_value(o.m_value)
  {
    static_cast<const Derived*>(this)->validate();
  }

  template <class Derived, class TValue>
  inline EncapsulatedValue<Derived,TValue>::EncapsulatedValue( DoValidate_t, TValue val )
    : m_value(val)
  {
    static_cast<const Derived*>(this)->validate();
  }

  template <class Derived, class TValue>
  inline ncconstexpr17 TValue& EncapsulatedValue<Derived,TValue>::get() noexcept
  {
    return m_value;
  }

  template <class Derived, class TValue>
  inline constexpr const TValue& EncapsulatedValue<Derived,TValue>::get() const noexcept
  {
    return m_value;
  }

  template <class Derived, class TValue>
  inline ncconstexpr17 void EncapsulatedValue<Derived,TValue>::set( TValue v ) noexcept
  {
    m_value = v;
  }

  template <class Derived, class TValue>
  inline ncconstexpr17 void EncapsulatedValue<Derived,TValue>::set( const Derived& o ) noexcept
  {
    m_value = o.m_value;
  }

  namespace detail {
    struct NCRYSTAL_API fmthelper{
      static void dofmt( std::ostream& os, const double& val ) { os << fmtg(val); }
      template<class T>
      static void dofmt( std::ostream& os, const T& val ) { os << val; }
    };
  }

  template <class Derived, class TValue>
  inline void EncapsulatedValue<Derived,TValue>::stream(std::ostream& os) const
  {
    detail::fmthelper::dofmt(os,m_value);
    os << Derived::unit();
  }

  template <class Derived, class TValue>
  inline std::ostream& operator<< (std::ostream& os, const EncapsulatedValue<Derived,TValue>& val)
  {
    static_cast<const Derived*>(&val)->stream(os);
    return os;
  }

  inline constexpr NeutronEnergy NeutronWavelength::energy() const noexcept
  {
    return *this;
  }

  inline NeutronWavelength NeutronEnergy::wavelength() const noexcept
  {
    return *this;
  }

  inline constexpr NeutronWavelength NeutronEnergy::constexpr_wavelength() const noexcept
  {
    return NeutronWavelength{ constexpr_ekin2wl( this->get() ) };
  }

  inline double NeutronEnergy::k() const noexcept
  {
    return ekin2k(m_value);
  }

  inline constexpr double NeutronEnergy::constexpr_k() const noexcept
  {
    return constexpr_ekin2k(m_value);
  }

  inline constexpr double NeutronEnergy::ksq() const noexcept
  {
    return ekin2ksq(m_value);
  }

  inline constexpr double NeutronWavelength::k() const noexcept
  {
    return wl2k(m_value);
  }

  inline constexpr double NeutronWavelength::ksq() const noexcept
  {
    return wl2ksq(m_value);
  }

  //Using nc_as_const so it can also be constexpr in C++11 [update: they can't
  //all anyway, since std::sqrt is not constexpr and so ekin2wl is not]:
  inline NeutronWavelength::NeutronWavelength(NeutronEnergy ekin) noexcept
    : EncapsulatedValue(ekin2wl(nc_as_const(ekin).get())) {}
  inline constexpr NeutronEnergy::NeutronEnergy(NeutronWavelength wl) noexcept
    : EncapsulatedValue(wl2ekin(nc_as_const(wl).get())) {}
  inline constexpr MosaicityFWHM::MosaicityFWHM(MosaicitySigma mos) noexcept
    : EncapsulatedValue( nc_as_const(mos).get() * kSigma2FWHM ) {}
  inline constexpr MosaicitySigma::MosaicitySigma(MosaicityFWHM mos) noexcept
    : EncapsulatedValue( nc_as_const(mos).get() * kFWHM2Sigma ) {}

  inline constexpr Length::Length(NeutronWavelength wl) noexcept
    : EncapsulatedValue( nc_as_const(wl).get() * angstrom ) {}

  inline constexpr NeutronWavelength Length::as_wavelength() const noexcept
  {
    static_assert( angstrom == 1e-10, "" );
    return NeutronWavelength{ m_value * 1e10 };
  }

  inline void Length::validate() const
  {
    if ( !(m_value>=0.0) || !std::isfinite(m_value) )
      NCRYSTAL_THROW2(CalcError,"Length::validate() failed. Invalid value:" << *this );
  }

  inline constexpr MosaicitySigma MosaicityFWHM::sigma() const noexcept
  {
    return *this;
  }
  inline constexpr MosaicityFWHM MosaicitySigma::fwhm() const noexcept
  {
    return *this;
  }

  inline constexpr double Temperature::kT() const noexcept
  {
    return constant_boltzmann * m_value;
  }

  inline constexpr double DebyeTemperature::kT() const noexcept
  {
    return constant_boltzmann * m_value;
  }

  namespace detail {
    constexpr double pow2(double x) noexcept { return x*x; }
  }

  inline constexpr SigmaBound::SigmaBound( SigmaFree sf, AtomMass mass ) noexcept
    : EncapsulatedValue<SigmaBound>( detail::pow2( ( const_neutron_mass_amu + nc_as_const(mass).get() ) / nc_as_const(mass).get() ) * nc_as_const(sf).get() )
  {
  }

  inline constexpr SigmaFree::SigmaFree( SigmaBound sb, AtomMass mass ) noexcept
    : EncapsulatedValue<SigmaFree>( detail::pow2( nc_as_const(mass).get() / ( const_neutron_mass_amu + nc_as_const(mass).get() ) ) * nc_as_const(sb).get() )
  {
  }

  inline constexpr SigmaFree SigmaBound::free( AtomMass am ) const noexcept
  {
    return SigmaFree( *this, am );
  }

  inline constexpr SigmaBound SigmaFree::bound( AtomMass am ) const noexcept
  {
    return SigmaBound( *this, am );
  }

  inline Density::Density( NumberDensity nd_, AtomMass averageAtomMass_ )
    : Density( [](NumberDensity nd, AtomMass averageAtomMass)
    {
      nd.validate();
      averageAtomMass.validate();
      constexpr double kkk = 1e27 * constant_dalton2kg;
      double result = kkk * averageAtomMass.get() * nd.get();
      if ( ! (result >= 0) || std::isinf(result) || std::isnan(result) )
        NCRYSTAL_THROW(CalcError,"Problems calculating Density from NumberDensity from NumberDensity and averageAtomMass");
      return Density{ result };
    }(nd_,averageAtomMass_))
  {}

  inline NumberDensity::NumberDensity( Density dens_, AtomMass averageAtomMass_ )
    : NumberDensity([](Density dens, AtomMass averageAtomMass)
    {
      dens.validate();
      averageAtomMass.validate();
      constexpr double kkk = 1e27 * constant_dalton2kg;
      double denom = kkk * averageAtomMass.get();
      if ( ! (denom > 0) || std::isinf(denom) || std::isnan(denom) )
        NCRYSTAL_THROW(CalcError,"Can not calculate NumberDensity from Density when averageAtomMass is vanishing or invalid");
      return NumberDensity{ dens.get() / denom };
    }(dens_,averageAtomMass_))
  {
  }

  inline constexpr double AtomMass::relativeToNeutronMass() const noexcept
  {
    return const_inv_neutron_mass_amu * m_value;
  }

  inline constexpr EnergyDomain::EnergyDomain(NeutronEnergy el,NeutronEnergy eh) ncnoexceptndebug
    : elow(el), ehigh(eh)
  {
#if nc_cplusplus >= 201703L
    nc_assert( elow.dbl() <= ehigh.dbl() );
#endif
  }

  inline constexpr EnergyDomain::EnergyDomain(EnergyDomain::no_check_t, NeutronEnergy el,NeutronEnergy eh) noexcept
    : elow(el), ehigh(eh)
  {
  }

  inline constexpr bool EnergyDomain::contains( NeutronEnergy ekin ) const ncnoexceptndebug
  {
#if nc_cplusplus >= 201703L
    nc_assert( elow.dbl() <= ehigh.dbl() );
#endif
    //nc_as_const to get constexpr .dbl in before C++17:
    return !isNull() && ( nc_as_const(ekin).dbl() >= elow.dbl() ) && ( nc_as_const(ekin).dbl() <= ehigh.dbl() );
    //NB: More efficient version might trigger FPE's:
    //    return ( nc_as_const(ekin).dbl()-elow.dbl() ) * ( nc_as_const(ekin).dbl()-ehigh.dbl() ) <= 0.0;
  }


  inline void AtomMass::validate() const
  {
    if ( ! ( m_value >= 0.0 && m_value < 1e9 )  )
      NCRYSTAL_THROW2(CalcError,"AtomMass::validate() failed. Invalid value:" << *this );
  }

  inline void Pressure::validate() const
  {
    if ( ! ( m_value > 0.0 && m_value < 1e15 ) )
      NCRYSTAL_THROW2(CalcError,"Pressure::validate() failed. Invalid value:" << *this );
  }

  inline void Temperature::validate() const
  {
    if ( ! ( m_value > 0.0 && m_value < 1e99 ) )
      NCRYSTAL_THROW2(CalcError,"Temperature::validate() failed. Invalid value:" << *this );
    if ( ! ( m_value >= allowed_range.first ) )
      NCRYSTAL_THROW2(CalcError,"Temperature::validate() failed for T="<<*this<<" (temperature values below "<<Temperature{allowed_range.first}<< " are not supported).");
    if ( ! ( m_value <= allowed_range.second ) )
      NCRYSTAL_THROW2(CalcError,"Temperature::validate() failed for T="<<*this<<" (temperature values above "<<Temperature{allowed_range.second}<< " are not supported).");
  }

  inline void DebyeTemperature::validate() const
  {
    if ( ! ( m_value > 0.0 && m_value < 1e9 ) )
      NCRYSTAL_THROW2(CalcError,"DebyeTemperature::validate() failed. Invalid value:" << *this );
  }

  inline void CosineScatAngle::validate() const
  {
    if ( ! ( m_value >= -1.0 && m_value <= 1.0 ) )
      NCRYSTAL_THROW2(CalcError,"CosineScatAngle::validate() failed. Invalid value:" << *this );
  }

  inline void CrossSect::validate() const
  {
    if ( ! ( m_value >= 0.0 && std::isfinite(m_value) ) )
      NCRYSTAL_THROW2(CalcError,"CrossSect::validate() failed. Invalid value:" << *this );
  }

  inline void SigmaBound::validate() const
  {
    if ( ! ( m_value >= 0.0 && m_value < 1e9 ) )
      NCRYSTAL_THROW2(CalcError,"SigmaBound::validate() failed. Invalid value:" << *this );
  }

  inline void SigmaFree::validate() const
  {
    if ( ! ( m_value >= 0.0 && m_value < 1e9 ) )
      NCRYSTAL_THROW2(CalcError,"SigmaFree::validate() failed. Invalid value:" << *this );
  }

  inline void SigmaAbsorption::validate() const
  {
    if ( ! ( m_value >= 0.0 && m_value < 1e9 ) )
      NCRYSTAL_THROW2(CalcError,"SigmaAbsorption::validate() failed. Invalid value:" << *this );
  }

  inline void ScatLenDensity::validate() const
  {
    if ( std::isnan( m_value ) || m_value < -1.0e9 || m_value > 1.0e9 )
      NCRYSTAL_THROW2(CalcError,"ScatLenDensity::validate() failed. Invalid value:" << *this );
  }

  inline void ScatLenDensity::stream(std::ostream& os) const
  {
    os << fmtg(get()) << "x" << unit();//Unit starts with number, so need "x" to separate from value.
  }

  inline constexpr SLDContrast ScatLenDensity::contrast(ScatLenDensity other) const noexcept
  {
    return SLDContrast{ *this, other };
  }

  inline constexpr SLDContrast::SLDContrast( ScatLenDensity rho1, ScatLenDensity rho2 ) noexcept
    : SLDContrast( constexpr_abs( nc_as_const(rho1).dbl()-nc_as_const(rho2).dbl() ) )
  {
  }

  inline constexpr double SLDContrast::valuePerAa2() const noexcept
  {
    return dbl() * 1.0e-6;
  }

  inline constexpr double SLDContrast::valuePerCM2() const noexcept
  {
    return dbl() * 1.0e10;
  }

  inline void SLDContrast::validate() const
  {
    if ( ! ( m_value >= 0.0 && m_value < 1e9 ) )
      NCRYSTAL_THROW2(CalcError,"SLDContrast::validate() failed. Invalid value:" << *this );
  }

  inline void Density::validate() const
  {
    if ( ! ( m_value >= 0.0 && m_value < 1e6 ) )
      NCRYSTAL_THROW2(CalcError,"Density::validate() failed. Invalid value:" << *this );
  }

  inline void NumberDensity::validate() const
  {
    if ( ! ( m_value >= 0.0 && m_value < 1e6 ) )
      NCRYSTAL_THROW2(CalcError,"NumberDensity::validate() failed. Invalid value:" << *this );
  }

  inline void MosaicityFWHM::validate() const
  {
    if ( ! ( m_value > 0.0 && m_value <= kPiHalf ) )
      NCRYSTAL_THROW2(CalcError,"MosaicityFWHM::validate() failed. Invalid value:" << *this );
  }

  inline void MosaicitySigma::validate() const
  {
    if ( ! ( m_value > 0.0 && m_value*kSigma2FWHM <= kPiHalf ) )
      NCRYSTAL_THROW2(CalcError,"MosaicitySigma::validate() failed. Invalid value:" << *this );
  }

  inline void NeutronWavelength::validate() const
  {
    if ( ! ( m_value >= 0.0 ) )
      NCRYSTAL_THROW2(CalcError,"NeutronWavelength::validate() failed. Invalid value:" << *this );
  }

  inline void NeutronEnergy::validate() const
  {
    if ( ! ( m_value >= 0.0 ) )
      NCRYSTAL_THROW2(CalcError,"NeutronEnergy::validate() failed. Invalid value:" << *this );
  }

  constexpr inline HKL::HKL(int hh, int kk, int ll) noexcept : h(hh), k(kk), l(ll) {}

  inline bool HKL::operator<(const HKL&o) const noexcept
  {
    return ( h!=o.h ? o.h<h : ( k!=o.k ? o.k<k : o.l<l ) );
  }

  inline bool HKL::operator==(const HKL&o) const noexcept
  {
    return  h==o.h && k==o.k && l==o.l;
  }

  inline std::ostream& operator<<(std::ostream& os, const DataSourceName& dsn)
  {
    return os << dsn.str();
  }

  inline DataSourceName::DataSourceName( std::string dsn )
    : m_str(makeSO<std::string>(std::move(dsn)))
  {
  }

  inline DataSourceName& DataSourceName::operator=( std::string dsn )
  {
    if ( *m_str != dsn )
      m_str = makeSO<std::string>(std::move(dsn));
    return *this;
  }

  inline bool DensityState::operator==( const DensityState& o ) const
  {
    return this->value == o.value && this->type == o.type;
  }

  inline void DensityState::validate() const
  {
    if ( !(value>0.0) || !( value <= 1.0e200) )
      NCRYSTAL_THROW2(BadInput,"Density value invalid or out of bounds: "<<*this);
  }

  inline bool UCNMode::operator==( const UCNMode& o ) const
  {
    return this->mode == o.mode && this->threshold == o.threshold;
  }

  inline bool ThreadCount::indicatesAutoDetect() const noexcept
  {
    return this->get()>=9999;
  }

}

#endif
