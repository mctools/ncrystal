#ifndef NCrystal_InfoTypes_hh
#define NCrystal_InfoTypes_hh

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

#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/interfaces/NCAtomData.hh"
#include "NCrystal/core/NCSmallVector.hh"

/////////////////////////////////////////////////////////
// Data structures used on Info objects (in NCInfo.hh) //
/////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  class Info;
  class DynamicInfo;

  struct NCRYSTAL_API StructureInfo final {
    unsigned spacegroup = 0;//From 1-230 if provided, 0 if not available
    double lattice_a = 0.0;//angstrom
    double lattice_b = 0.0;//angstrom
    double lattice_c = 0.0;//angstrom
    double alpha = 0.0;//degree
    double beta = 0.0;//degree
    double gamma = 0.0;//degree
    double volume = 0.0;//Aa^3
    unsigned n_atoms = 0;//Number of atoms per unit cell
  };

  enum class NCRYSTAL_API HKLInfoType : unsigned { SymEqvGroup = 0,
                                                   ExplicitHKLs = 1,
                                                   ExplicitNormals = 2,
                                                   Minimal = 3 };
  NCRYSTAL_API std::ostream& operator<< ( std::ostream&, const HKLInfoType&);

  struct NCRYSTAL_API HKLInfo final : private MoveOnly {
    HKL hkl;
    unsigned multiplicity = 0;
    double dspacing = 0.0;//angstrom
    double fsquared = 0.0;//barn

    //Concerning specification of the full list of (h,k,l) entries in the group
    //(which has length "multiplicity"), and the associated plane normals, there
    //are several options. In the most ideal case ("SymEqvGroup") the HKL
    //groupings are according to symmetry-equivalency, and the Info object
    //contains StructureInfo including spacegroup. Thus, nothing more than the
    //hkl values above are needed. This is because the full list of HKL entries
    //in the group can be provided via the symmetry (cf. NCEqRefl.hh), and the
    //normals can then be constructed through knowledge of the unit cell
    //transformations contained in the StructureInfo (cf. NCLatticeUtils.hh).
    //
    //Alternatively, when the groupings are not by symmetry-equivalency, the
    //full list of HKL entries can be available ("ExplicitHKLs"). Again,
    //structure info (but not necessarily the spacegroup) is needed to have
    //access to normals. As another alternative, it is possibly to have an
    //explicit list of normals ("ExplicitNormals"). In that case, normals are
    //always available, whereas construction of (h,k,l) points would rely on
    //StructureInfo. Finally, all information about (h,k,l) points and normals
    //might be missing ("Minimal"), in which case the HKL lists would likely
    //only be useful to describe powder diffraction.
    //
    //The header NCPlaneProvider.hh contains helper classes for accessing HKL
    //values or normals in a manner independent of the exact manner of
    //specification, so client C++ code should rarely have to check the
    //HKLInfoType directly.

    HKLInfoType type() const;
    class NCRYSTAL_API Normal final : public StronglyTypedFixedVector<Normal,double,3> {
    public:
      using StronglyTypedFixedVector::StronglyTypedFixedVector;
    };
    struct ExplicitVals {
      Variant<std::vector<Normal>,std::vector<HKL>> list;//empty implies Type==Minimal
    };
    std::unique_ptr<ExplicitVals> explicitValues;
  };

  //Using small vector since VisualStudio's vector does not support move-only
  //types:
  using HKLList = SmallVector<HKLInfo,1>;

  class NCRYSTAL_API AtomIndex final : public EncapsulatedValue<AtomIndex,unsigned> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
    static AtomIndex createInvalidObject();
    constexpr bool isInvalid() const noexcept;
    ncconstexpr17 void validate() const noexcept {}
    //NB: Since unsigned is only guaranteed to be at least 16bit, the underlying
    //value representing invalid indices might be as low as 65535 (of course, in
    //practice it is usually 4294967295).
  };

  struct NCRYSTAL_API IndexedAtomData {
    //AtomData and associated index. The index is *only* valid in association
    //with a particular Info object. It exists since it is in principle possible
    //to have the same fundamental atom playing more than one role in a given
    //material (for instance, the same atom could have different displacements
    //on different positions in the unit cell).
    AtomDataSP atomDataSP;
    AtomIndex index;

    const AtomData& data() const ncnoexceptndebug;

    //Sort by index (comparison should only be performed with objects associated
    //with the same single-phase Info object). For multi-phase objects,
    //comparison will be done using the atom data UID:
    bool operator<(const IndexedAtomData& o) const ncnoexceptndebug;
    bool operator==(const IndexedAtomData& o) const ncnoexceptndebug;
  };
  NCRYSTAL_API std::ostream& operator<< ( std::ostream&, const IndexedAtomData&);

  class NCRYSTAL_API AtomInfo final : private MoveOnly {
  public:

    //////////////////////////////////////////////////////////////////////////////
    // Information about one kind of atom in a crystal unit cell, sharing both  //
    // atomic composition and dynamic behaviour (reflected in e.g. mean squared //
    // displacement values and associated DynamicInfo object).                  //
    //////////////////////////////////////////////////////////////////////////////

    //Always contains info about atomic compositions:
    ncconstexpr17 const IndexedAtomData& indexedAtomData() const noexcept { return m_iad; }
    ncconstexpr17 const IndexedAtomData& atom() const noexcept { return m_iad; }
    AtomDataSP atomDataSP() const noexcept { return m_iad.atomDataSP; }
    const AtomData& atomData() const ncnoexceptndebug { return m_iad.data(); }

    //Always contains non-empty list of associated unit cell positions:
    class NCRYSTAL_API Pos final : public StronglyTypedFixedVector<Pos,double,3> {
    public:
      using StronglyTypedFixedVector::StronglyTypedFixedVector;
    };
    using AtomPositions = std::vector<Pos>;
    const AtomPositions& unitCellPositions() const noexcept { return m_pos; }
    unsigned numberPerUnitCell() const noexcept { return static_cast<unsigned>( m_pos.size() ); }

    //Mean-square-displacements (a.k.a. "U_iso") in angstrom^2 is optional. Note
    //that this is the displacement projected onto a linear axis, for direct
    //usage in isotropic Debye-Waller factors:
    ncconstexpr17 const Optional<double>& msd() const noexcept { return m_msd; }

    //Debye temperature is optional::
    ncconstexpr17 const Optional<DebyeTemperature>& debyeTemp() const noexcept { return m_dt; }

    //Corresponding DynamicInfo object on the same Info (returns nullptr if not available):
    const DynamicInfo* correspondingDynamicInfo() const { return m_dyninfo; }

    //Constructor:
    AtomInfo( IndexedAtomData, AtomPositions&&, Optional<DebyeTemperature>, Optional<double> msd );
  private:
    IndexedAtomData m_iad;
    Optional<DebyeTemperature> m_dt;
    Optional<double> m_msd;
    std::vector<Pos> m_pos;
    const DynamicInfo* m_dyninfo = nullptr;
  public:
    //Implementation details for builder infrastructure (todo: avoid with proper
    //builder-pattern for AtomInfo objects):
    void detail_setupLink(DynamicInfo*);
    std::vector<Pos>& detail_accessPos() { return m_pos; }
  };

  //NB: VSCode std::vectors does not support move-only objects apparently, so we
  //are using our own SmallVector for AtomInfoList:
  using AtomInfoList = SmallVector<AtomInfo,4>;
  using AtomList = AtomInfoList;//obsolete alias

  class NCRYSTAL_API DynamicInfo : public UniqueID {
  public:
    DynamicInfo(double fraction, IndexedAtomData, Temperature);
    virtual ~DynamicInfo() = 0;//Make abstract
    double fraction() const;
    void changeFraction(double f) { m_fraction = f; }
    Temperature temperature() const;//same as on associated Info object

    const IndexedAtomData& atom() const { return m_atom; }
    AtomDataSP atomDataSP() const { return m_atom.atomDataSP; }
    const AtomData& atomData() const { return m_atom.data(); }

    //Corresponding AtomInfo object on the same Info (nullptr if not available):
    const AtomInfo* correspondingAtomInfo() const { return m_atomInfo; }

  private:
    double m_fraction;
    IndexedAtomData m_atom;
    Temperature m_temperature;
    AtomInfo* m_atomInfo = nullptr;
    friend class AtomInfo;
  };

  //NB: VSCode std::vectors does not support move-only objects apparently, so we
  //are using our own SmallVector for DynamicInfoList:
  typedef SmallVector<std::unique_ptr<DynamicInfo>,4> DynamicInfoList;

  class NCRYSTAL_API DI_Sterile final : public DynamicInfo {
    //Class indicating elements for which inelastic neutron scattering is absent
    //or disabled.
  public:
    using DynamicInfo::DynamicInfo;
    virtual ~DI_Sterile();
  };

  class NCRYSTAL_API DI_FreeGas final : public DynamicInfo {
    //Class indicating elements for which inelastic neutron scattering should be
    //modelled as scattering on a free gas.
  public:
    using DynamicInfo::DynamicInfo;
    virtual ~DI_FreeGas();
  };

  class NCRYSTAL_API DI_ScatKnl : public DynamicInfo {
  public:
    //Base class for dynamic information which can, directly or indirectly,
    //result in a SABData scattering kernel. The class is mostly semantic, as no
    //SABData access interface is provided on this class, as some derived
    //classes (e.g. VDOS) need dedicated algorithms in order to create the
    //SABData object. This class does, however, provide a unified interface for
    //associated data which is needed in order to use the SABData for
    //scattering. Currently this is just the grid of energy points for which SAB
    //integration will perform analysis and caching.
    virtual ~DI_ScatKnl();

    //If source dictated what energy grid to use for caching cross-sections,
    //etc., it can be returned here. It is ok to return a null ptr, leaving the
    //decision entirely to the consuming code. Grids must have at least 3
    //entries, and grids of size 3 actually indicates [emin,emax,npts], where
    //any value can be 0 to leave the choice for the consuming code. Grids of
    //size >=4 must be proper grids.
    typedef std::shared_ptr<const VectD> EGridShPtr;
    virtual EGridShPtr energyGrid() const = 0;
  protected:
    using DynamicInfo::DynamicInfo;
  };

  class NCRYSTAL_API DI_ScatKnlDirect : public DI_ScatKnl {
  public:
    //Pre-calculated scattering kernel which at most needs a conversion to
    //SABData format before it is available. For reasons of efficiency, this
    //conversion is actually not required to be carried out before calling code
    //calls the MT-safe ensureBuildThenReturnSAB().
    using DI_ScatKnl::DI_ScatKnl;
    virtual ~DI_ScatKnlDirect();

    //Use ensureBuildThenReturnSAB to access the scattering kernel:
    std::shared_ptr<const SABData> ensureBuildThenReturnSAB() const;

    //check if SAB is already built:
    bool hasBuiltSAB() const;

  protected:
    //Implement in derived classes to build the completed SABData object (will
    //only be called once and in an MT-safe context, protected by per-object
    //mutex):
    virtual std::shared_ptr<const SABData> buildSAB() const = 0;
  private:
    mutable std::shared_ptr<const SABData> m_sabdata;
    mutable std::mutex m_mutex;
  };

  class NCRYSTAL_API DI_VDOS : public DI_ScatKnl {
  public:
    //For a solid material, a phonon spectrum in the form of a Vibrational
    //Density Of State (VDOS) parameterisation, can be expanded into a full
    //scattering kernel. The calling code is responsible for doing this,
    //including performing choices as to grid layout, expansion order, etc.
    using DI_ScatKnl::DI_ScatKnl;
    virtual ~DI_VDOS();
    virtual const VDOSData& vdosData() const = 0;

    //The above vdosData() function returns regularised VDOS. The following
    //functions provide optional access to the original curves (returns empty
    //vectors if not available):
    virtual const VectD& vdosOrigEgrid() const = 0;
    virtual const VectD& vdosOrigDensity() const = 0;
  };

  class NCRYSTAL_API DI_VDOSDebye final : public DI_ScatKnl {
  public:
    //An idealised VDOS spectrum, based on the Debye Model in which the spectrum
    //rises quadratically with phonon energy below a cutoff value, kT, where T
    //is the Debye temperature (which must be available on the associated Info
    //object).
    DI_VDOSDebye( double fraction,
                  IndexedAtomData,
                  Temperature temperature,
                  DebyeTemperature debyeTemperature );
    virtual ~DI_VDOSDebye();
    DebyeTemperature debyeTemperature() const;
    EGridShPtr energyGrid() const override { return nullptr; }
  private:
    DebyeTemperature m_dt;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline double DynamicInfo::fraction() const { return m_fraction; }
  inline Temperature DynamicInfo::temperature() const { return m_temperature; }
  inline DI_VDOSDebye::DI_VDOSDebye( double fr, IndexedAtomData aaa, Temperature tt, DebyeTemperature dt )
    : DI_ScatKnl(fr,std::move(aaa),tt),m_dt(dt) { nc_assert(m_dt.get()>0.0); }
  inline DebyeTemperature DI_VDOSDebye::debyeTemperature() const { return m_dt; }
  inline const AtomData& IndexedAtomData::data() const ncnoexceptndebug { return atomDataSP; }
  inline bool IndexedAtomData::operator<(const IndexedAtomData& o) const ncnoexceptndebug {
    //Sanity check (same index means same AtomData instance on a single-phase
    //info object, multi-phase objects have invalid indices):
    nc_assert( atomDataSP == o.atomDataSP || index != o.index || ( index.isInvalid() && o.index.isInvalid() ) );
    nc_assert( atomDataSP!=nullptr );
    nc_assert( o.atomDataSP!=nullptr );
    return ( index.isInvalid()
             ? ( *atomDataSP == *o.atomDataSP
                 ? atomDataSP->getUniqueID() < o.atomDataSP->getUniqueID()
                 : *atomDataSP < *o.atomDataSP )
             : index < o.index );
  }
  inline bool IndexedAtomData::operator==(const IndexedAtomData& o) const ncnoexceptndebug {
    nc_assert( atomDataSP == o.atomDataSP || index != o.index || ( index.isInvalid() && o.index.isInvalid() ) );
    nc_assert( atomDataSP!=nullptr );
    nc_assert( o.atomDataSP!=nullptr );
    return index.isInvalid()
      ? atomDataSP->getUniqueID() == o.atomDataSP->getUniqueID()
      : index == o.index;
  }
  inline AtomIndex AtomIndex::createInvalidObject() { return AtomIndex{std::numeric_limits<unsigned>::max()}; }
  inline constexpr bool AtomIndex::isInvalid() const noexcept { return get() == std::numeric_limits<unsigned>::max(); }
  inline std::ostream& operator<< ( std::ostream& os, const IndexedAtomData& ai) {
    os<<"Atom(descr=\""<<ai.data().description(false)<<"\",index="<<ai.index<<")";
    return os;
  }
  inline HKLInfoType HKLInfo::type() const
  {
    if ( !explicitValues )
      return HKLInfoType::SymEqvGroup;//usual case (ncmat,laz,lau)
    auto& evl = explicitValues->list;
    if ( evl.has_value<std::vector<HKL>>() )//ncmat without spacegroup
      return HKLInfoType::ExplicitHKLs;
    if ( evl.has_value<std::vector<Normal>>() )//rare (handcrafted?)
      return HKLInfoType::ExplicitNormals;
    return HKLInfoType::Minimal;
  }

  inline std::ostream& operator<< ( std::ostream& os, const HKLInfoType& ht )
  {
    switch( ht ) {
    case HKLInfoType::SymEqvGroup: return os << "SymEqvGroup";
    case HKLInfoType::ExplicitHKLs: return os << "ExplicitHKLs";
    case HKLInfoType::ExplicitNormals: return os << "ExplicitNormals";
    case HKLInfoType::Minimal: return os << "Minimal";
    };
    nc_assert_always(false);
    return os;
  }


}

#endif
