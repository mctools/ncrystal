#ifndef NCrystal_Types_hh
#define NCrystal_Types_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/NCDefs.hh"
#include <ostream>

//////////////////////////////////////////////
//                                          //
// Common types used in various interfaces  //
//                                          //
//////////////////////////////////////////////

namespace NCrystal {

  //Base class used to MT-safe caching. Client code can provide polymorphic
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
  class NCRYSTAL_API EncapsulatedValue {
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

    //Accessing the wrapped value always require invocation of a named method or
    //explicit cast (for strong type safety):
    explicit constexpr operator TValue() const noexcept;
    ncconstexpr17 TValue& get() noexcept;
    constexpr const TValue& get() const noexcept;
    ncconstexpr17 void set( TValue ) noexcept;
    ncconstexpr17 void set( const Derived& ) noexcept;

    //"dbl" is alias for "get" (can only be used when TValue is a double)::
    ncconstexpr17 double& dbl() noexcept;
    constexpr const double& dbl() const noexcept;

    //Can be used in boolean expression (usually this evaluates as "true" if value is non-zero).
    explicit constexpr operator bool() const noexcept;

    //Supports comparisons:
#if __cplusplus >= 202002L
    auto operator<=>(const EncapsulatedValue&) const = default;
#else
    constexpr bool operator<( const EncapsulatedValue& o ) const noexcept { return m_value < o.m_value; }
    constexpr bool operator>( const EncapsulatedValue& o ) const noexcept { return m_value > o.m_value; }
    constexpr bool operator<=( const EncapsulatedValue& o ) const noexcept { return m_value <= o.m_value; }
    constexpr bool operator>=( const EncapsulatedValue& o ) const noexcept { return m_value >= o.m_value; }
    constexpr bool operator==( const EncapsulatedValue& o ) const noexcept { return m_value == o.m_value; }
    constexpr bool operator!=( const EncapsulatedValue& o ) const noexcept { return m_value != o.m_value; }
#endif
  protected:
    TValue m_value;
  };

  //Values are printed with suitable trailing unit:
  template <class Derived, class TValue>
  NCRYSTAL_API std::ostream& operator<< (std::ostream& os, const EncapsulatedValue<Derived,TValue>& val) noexcept;

  class NeutronEnergy;

  class NCRYSTAL_API NeutronWavelength final : public EncapsulatedValue<NeutronWavelength> {
    //Neutron wavelength in angstrom
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "Aa"; }
    //Automatic conversions from/to energy:
    constexpr NeutronWavelength(NeutronEnergy) noexcept;
    ncnodiscard17 constexpr NeutronEnergy energy() const noexcept;
    void validate() const;
  };

  class NCRYSTAL_API NeutronEnergy final : public EncapsulatedValue<NeutronEnergy> {
  public:
    //Neutron kinetic energy in electronvolt
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "eV"; }
    //Automatic conversions from/to wavelength:
    constexpr NeutronEnergy(NeutronWavelength) noexcept;
    ncnodiscard17 constexpr NeutronWavelength wavelength() const noexcept;
    void validate() const;
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
    //Cross section value (usually calculated for a particular neutron state and material).
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
    static constexpr const char * unit() noexcept { return "barn"; }
    constexpr SigmaBound( SigmaFree, AtomMass ) noexcept;
    constexpr SigmaFree free( AtomMass am ) const noexcept;
    void validate() const;
  };

  class NCRYSTAL_API SigmaFree final : public EncapsulatedValue<SigmaFree> {
  public:
    //Free-limit cross section value (usually a material-specific constant).
    using EncapsulatedValue::EncapsulatedValue;
    constexpr SigmaFree() noexcept : EncapsulatedValue{0.0} {}//line here for gcc4.8 @ centos7
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

  class NCRYSTAL_API Density final : public EncapsulatedValue<Density> {
  public:
    //Material density in g/cm3.
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "g/cm3"; }
    void validate() const;
  };

  class NCRYSTAL_API NumberDensity final : public EncapsulatedValue<NumberDensity> {
  public:
    //Material number density in atoms per Aa^3.
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "atoms/Aa^3"; }
    void validate() const;
  };

  class MosaicitySigma;

  class NCRYSTAL_API MosaicityFWHM final : public EncapsulatedValue<MosaicityFWHM> {
  public:
    //Mosaic spread (FWHM) in single crystals.
    using EncapsulatedValue::EncapsulatedValue;
    constexpr MosaicityFWHM(MosaicitySigma) noexcept;
    ncnodiscard17 constexpr MosaicitySigma sigma() const noexcept;
    static constexpr const char * unit() noexcept { return "radians"; }
    void validate() const;
  };

  class NCRYSTAL_API MosaicitySigma final : public EncapsulatedValue<MosaicitySigma> {
  public:
    //Mosaic spread (std dev, i.e. Gaussian sigma parameter) in single crystals.
    using EncapsulatedValue::EncapsulatedValue;
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

  class NCRYSTAL_API HKLPoint final : public StronglyTypedFixedVector<HKLPoint,double,3> {
    //Point in (h,k,l)-space.
  public:
    using StronglyTypedFixedVector::StronglyTypedFixedVector;
  };

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
    ncnodiscard17 constexpr bool contains(NeutronEnergy) const ncnoexceptndebug;
  };

  struct NCRYSTAL_API ScatterOutcome {
    NeutronEnergy ekin;
    NeutronDirection direction;
  };

  struct NCRYSTAL_API ScatterOutcomeIsotropic {
    NeutronEnergy ekin;
    CosineScatAngle mu;
  };

  //Index identifying RNG streams:
  class NCRYSTAL_API RNGStreamIndex final : public EncapsulatedValue<RNGStreamIndex,uint64_t> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
  };

  //Serialised RNG stream state:
  class NCRYSTAL_API RNGStreamState final : public EncapsulatedValue<RNGStreamState,std::string> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
  };

  //Process properties:
  enum class MaterialType { Anisotropic, Isotropic };
  enum class ProcessType { Absorption, Scatter };

  NCRYSTAL_API std::ostream& operator<<(std::ostream&, MaterialType);
  NCRYSTAL_API std::ostream& operator<<(std::ostream&, ProcessType);

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {

  inline std::ostream& operator<<(std::ostream& os, MaterialType mt)
  {
    return os << (mt==MaterialType::Isotropic?"Isotropic":"Anisotropic");
  }

  inline std::ostream& operator<<(std::ostream& os, ProcessType pt)
  {
    return os << (pt==ProcessType::Scatter?"Scatter":"Absorption");
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
  inline constexpr EncapsulatedValue<Derived,TValue>::operator TValue() const noexcept
  {
    return m_value;
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

  template <class Derived, class TValue>
  inline ncconstexpr17 double& EncapsulatedValue<Derived,TValue>::dbl() noexcept
  {
    static_assert(std::is_same<TValue,double>::value,
                  ".dbl() method can only be used for types wrapping a double");
    return m_value;
  }

  template <class Derived, class TValue>
  inline constexpr const double& EncapsulatedValue<Derived,TValue>::dbl() const noexcept
  {
    static_assert(std::is_same<TValue,double>::value,
                  ".dbl() method can only be used for types wrapping a double");
    return m_value;
  }

  template <class Derived, class TValue>
  inline constexpr EncapsulatedValue<Derived,TValue>::operator bool() const noexcept
  {
    return bool(m_value);
  }

  template <class Derived, class TValue>
  NCRYSTAL_API inline std::ostream& operator<< (std::ostream& os, const EncapsulatedValue<Derived,TValue>& val) noexcept
  {
    return os << TValue(val) << Derived::unit();
  }

  inline constexpr NeutronEnergy NeutronWavelength::energy() const noexcept
  {
    return *this;
  }

  inline constexpr NeutronWavelength NeutronEnergy::wavelength() const noexcept
  {
    return *this;
  }

  //Using nc_as_const so it can also be constexpr in C++11
  inline constexpr NeutronWavelength::NeutronWavelength(NeutronEnergy ekin) noexcept
    : EncapsulatedValue(ekin2wl(nc_as_const(ekin).get())) {}
  inline constexpr NeutronEnergy::NeutronEnergy(NeutronWavelength wl) noexcept
    : EncapsulatedValue(wl2ekin(nc_as_const(wl).get())) {}
  inline constexpr MosaicityFWHM::MosaicityFWHM(MosaicitySigma m) noexcept
    : EncapsulatedValue( nc_as_const(m).get() * kSigma2FWHM ) {}
  inline constexpr MosaicitySigma::MosaicitySigma(MosaicityFWHM m) noexcept
    : EncapsulatedValue( nc_as_const(m).get() * kFWHM2Sigma ) {}

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

  inline constexpr double AtomMass::relativeToNeutronMass() const noexcept
  {
    return const_inv_neutron_mass_amu * m_value;
  }

  inline constexpr EnergyDomain::EnergyDomain(NeutronEnergy el,NeutronEnergy eh) ncnoexceptndebug
    : elow(el), ehigh(eh)
  {
#  if __cplusplus >= 201703L
    nc_assert( elow.dbl() <= ehigh.dbl() );
#endif
  }

  inline constexpr bool EnergyDomain::contains( NeutronEnergy ekin ) const ncnoexceptndebug
  {
#  if __cplusplus >= 201703L
    nc_assert( elow.dbl() <= ehigh.dbl() );
#endif
    //nc_as_const to get constexpr .dbl in before C++17:
    return ( nc_as_const(ekin).dbl()-elow.dbl() ) * ( nc_as_const(ekin).dbl()-ehigh.dbl() ) <= 0.0;
  }


  inline void AtomMass::validate() const
  {
    if ( ! ( m_value >= 0.0 && m_value < 1e9 )  )
      NCRYSTAL_THROW2(CalcError,"AtomMass::validate() failed. Invalid value:" << *this );
  }

  inline void Temperature::validate() const
  {
    if ( ! ( m_value > 0.0 && m_value < 1e9 ) )
      NCRYSTAL_THROW2(CalcError,"Temperature::validate() failed. Invalid value:" << *this );
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
    if ( ! ( m_value >= 0.0 && m_value < 1e9 ) )
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

}

#endif
