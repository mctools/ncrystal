#ifndef NCrystal_LegacyProcs_hh
#define NCrystal_LegacyProcs_hh

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

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Legacy interface classes. Provided here for now for backwards               //
// compatibility with code using NCrystal. Everything is implemented in the    //
// NCrystal::Legacy namespace, which is for now injected into the main         //
// NCrystal namespace. This will eventually change!                            //
//                                                                             //
// It is no longer supported to actually implement physics models by extending //
// these legacy interfaces. This is prevented at compile-time using            //
// private-constructor-friend-statement tricks                                 //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/NCRNG.hh"
#include "NCrystal/NCProcImpl.hh"
#include "NCrystal/NCMatInfo.hh"
#include "NCrystal/NCLegacyMem.hh"

namespace NCrystal {

  namespace Legacy {

    class NCRYSTAL_API RandomBase : public RCBase {
      //Legacy RNG interface.
    public:
      virtual double generate() = 0;//generate numbers uniformly in (0,1]
    protected:
      virtual ~RandomBase();
    };

    NCRYSTAL_API RandomBase * defaultRandomGenerator(bool trigger_default = true);
    NCRYSTAL_API RCHolder<RandomBase> defaultRNG(bool trigger_default = true);

    class Process;
    class Absorption;
    class Scatter;
    class AbsorptionImplWrapper;
    class ScatterImplWrapper;
    class ScatterComp;

    class NCRYSTAL_API Info final : public RCBase {
      /////////////////////////////////////////////////////////////////////////////
      // Backwards compatibility wrapper for the MatInfo class. See NCMatInfo.hh //
      // for description of methods and classes. This wrapper will eventually be //
      // phased out and removed.                                                 //
      /////////////////////////////////////////////////////////////////////////////
      shared_obj<const MatInfo> m_wrapped;
    public:
      Info( shared_obj<const MatInfo> mi ) noexcept : m_wrapped(std::move(mi)) {}
      shared_obj<const MatInfo> wrapped() const noexcept { return m_wrapped; }
      UniqueIDValue getUniqueID() const { return m_wrapped->getUniqueID(); }
      bool isCrystalline() const { return m_wrapped->isCrystalline(); }
      bool hasStructureInfo() const { return m_wrapped->hasStructureInfo(); }
      const StructureInfo& getStructureInfo() const { return m_wrapped->getStructureInfo(); }
      double dspacingFromHKL( int h, int k, int l ) const { return m_wrapped->dspacingFromHKL(h,k,l); }
      bool hasDynamicInfo() const { return m_wrapped->hasDynamicInfo(); }
      const DynamicInfoList& getDynamicInfoList() const { return m_wrapped->getDynamicInfoList(); }
      bool hasXSectAbsorption() const { return m_wrapped->hasXSectAbsorption(); }
      double getXSectAbsorption() const { return m_wrapped->getXSectAbsorption().get(); }
      bool hasXSectFree() const { return m_wrapped->hasXSectFree(); }
      double getXSectFree() const { return m_wrapped->getXSectFree().get(); }
      bool providesNonBraggXSects() const { return m_wrapped->providesNonBraggXSects(); }
      double xsectScatNonBragg(double lambda) const { return m_wrapped->xsectScatNonBragg(NeutronWavelength{lambda}).get(); }
      bool hasTemperature() const { return m_wrapped->hasTemperature(); }
      double getTemperature() const { return m_wrapped->getTemperature().get(); }
      bool hasAnyDebyeTemperature() const { return m_wrapped->hasAtomDebyeTemp(); }
      bool hasGlobalDebyeTemperature() const { return false; }
      double getGlobalDebyeTemperature() const { nc_assert_always(false); return 0.0; }
      bool hasPerElementDebyeTemperature() const { return m_wrapped->hasAtomDebyeTemp(); }
      double getDebyeTemperatureByElement(const AtomIndex& idx) const
      {
        nc_assert_always( hasPerElementDebyeTemperature() );
        for ( auto& ai : m_wrapped->getAtomInfos() )
          if ( ai.indexedAtomData().index == idx )
            return ai.debyeTemp().value().get();
        nc_assert_always(false);
        return -1.0;
      }
      bool hasAtomInfo() const { return m_wrapped->hasAtomInfo(); }
      AtomList::const_iterator atomInfoBegin() const { return m_wrapped->atomInfoBegin(); }
      AtomList::const_iterator atomInfoEnd() const { return m_wrapped->atomInfoEnd(); }
      bool hasAtomPositions() const { return m_wrapped->hasAtomPositions(); }
      bool hasAtomMSD() const { return m_wrapped->hasAtomMSD(); }
      bool hasHKLInfo() const { return m_wrapped->hasHKLInfo(); }
      unsigned nHKL() const { return m_wrapped->nHKL(); }
      HKLList::const_iterator hklBegin() const { return m_wrapped->hklBegin(); }
      HKLList::const_iterator hklLast() const { return m_wrapped->hklLast(); }
      HKLList::const_iterator hklEnd() const { return m_wrapped->hklEnd(); }
      double hklDLower() const { return m_wrapped->hklDLower(); }
      double hklDUpper() const { return m_wrapped->hklDUpper(); }
      double hklDMinVal() const { return m_wrapped->hklDMinVal(); }
      double hklDMaxVal() const { return m_wrapped->hklDMaxVal(); }
      bool hasHKLDemiNormals() const { return m_wrapped->hasHKLDemiNormals(); }
      bool hasExpandedHKLInfo() const { return m_wrapped->hasExpandedHKLInfo(); }
      HKLList::const_iterator searchExpandedHKL(short h, short k, short l) const { return m_wrapped->searchExpandedHKL(h,k,l); }
      bool hasDensity() const { return m_wrapped->hasDensity(); }
      double getDensity() const { return m_wrapped->getDensity().dbl(); }
      bool hasNumberDensity() const { return m_wrapped->hasNumberDensity(); }
      double getNumberDensity() const { return m_wrapped->getNumberDensity().dbl(); }
      bool hasComposition() const { return m_wrapped->hasComposition(); }
      using CompositionEntry = MatInfo::CompositionEntry;
      using Composition = MatInfo::Composition;
      const Composition& getComposition() const { return m_wrapped->getComposition(); }
      const std::string& displayLabel(const AtomIndex& ai) const { return m_wrapped->displayLabel(ai); }
      AtomDataSP atomDataSP( const AtomIndex& ai ) const { return m_wrapped->atomDataSP(ai); }
      const AtomData& atomData( const AtomIndex& ai ) const { return m_wrapped->atomData(ai); }
      IndexedAtomData indexedAtomData( const AtomIndex& ai ) const { return m_wrapped->indexedAtomData(ai); }
      using CustomLine = MatInfo::CustomLine;
      using CustomSectionData = MatInfo::CustomSectionData;
      using CustomSectionName = MatInfo::CustomSectionName;
      using CustomData = MatInfo::CustomData;
      const CustomData& getAllCustomSections() const { return m_wrapped->getAllCustomSections(); }
      unsigned countCustomSections(const CustomSectionName& sectionname ) const { return m_wrapped->countCustomSections(sectionname); }
      const CustomSectionData& getCustomSection( const CustomSectionName& name,
                                                 unsigned index=0 ) const { return m_wrapped->getCustomSection(name,index); }
      bool isLocked() const { return m_wrapped->isLocked(); }
    };

    class NCRYSTAL_API CalcBase : public RCBase {
    public:

      //Base class handling RNG, UID, Ref-counting.
      //rather than public.

      const char * getCalcName() const { return m_name.c_str(); }
      void setCalcName(const char * n) { m_name = n; }

      bool isSubCalc(const CalcBase*) const;

      //Can be used to control and access the random number stream (but note that
      //often the setDefaultRandomGenerator(..) in NCRandom.hh will be an easier
      //way to do this globally). Sub-calcs will have their RNG changed as well:
      // void setRandomGenerator(RandomBase* rg);

      //Access current RNG (the first will init and return default RNG if none was
      //set explicitly - this modifies a mutable data member):
      RandomBase* getRNG() const {
        if (!m_legacyrng)
          initLegacyRng();
        nc_assert(!!m_legacyrng);
        return m_legacyrng.obj();
      }

      void setRNG( shared_obj<Modern::RNGStream> rng ) { m_rng = std::move(rng); }

      Modern::RNGStream& getModernRNG() const
      {
        if (!m_rng)
          m_rng = Modern::getDefaultRNGProducer()->produceByIdx(RNGStreamIndex{1000000000000000000ull});
        return *m_rng;
      }

      UniqueIDValue getUniqueID() const  { return m_uid.getUniqueID(); }

    protected:
      //Registering sub-calcs will result in their ref-counts being incremented
      //for the mother calc's lifetime, and also means that future
      //setRandomGenerator calls will propagate to them:
      void registerSubCalc(CalcBase*);
      virtual ~CalcBase();
    private:
      std::vector<CalcBase*> m_subcalcs;
      std::string m_name;
      UniqueID m_uid;
      mutable optional_shared_obj<Modern::RNGStream> m_rng;
      mutable RCHolder<RandomBase> m_legacyrng;
      void initLegacyRng() const;
      //double initDefaultRand() const;
      CalcBase(const char * calculator_type_name);
      friend class Process;
    };

    class NCRYSTAL_API Process : public CalcBase {
    public:

      /////////////////////////////////////////////////////////////////////
      // Base class for calculations of processes in materials.          //
      //                                                                 //
      // Note that for maximum compatibility vector directions are       //
      // specified via double arrays of length 3.                        //
      //                                                                 //
      // Note that the unit for kinetic energy in the calls below is eV  //
      // (electronvolt).                                                 //
      /////////////////////////////////////////////////////////////////////

      //The process cross-section (in barns):
      virtual double crossSection(double ekin, const double (&neutron_direction)[3] ) const = 0;

      //Override if the process is guaranteed to always produce vanishing
      //cross-sections outside a given domain (the special case
      //ekin_low=ekin_high=infinity is used to indicate a somewhat unusual process
      //with vanishing cross-section everywhere):
      virtual void domain(double& ekin_low, double& ekin_high) const { ekin_low = 0.0; ekin_high = kInfinity; }

      //Check if process is oriented. For non-oriented processes, the results do
      //not depend on the incident direction of the neutron, and outcomes such as
      //scatterings are phi-symmetric:
      virtual bool isOriented() const { return true; };

      //For non-oriented processes, callers can either use the method above to
      //access cross-sections (with free choice of reference frame for the
      //direction vectors), or use the following method for convenience:
      virtual double crossSectionNonOriented( double ekin ) const;

      virtual void validate();//call to perform a quick (incomplete) validation
      //that cross sections are vanishing outside
      //domain(..).

      //Check (via domain(..)) if cross section is 0 everywhere:
      bool isNull() const;
    protected:
      virtual ~Process();
      Process(const char * calculator_type_name);
      friend class Absorption;
      friend class Scatter;
    };

    class NCRYSTAL_API Absorption : public Process {
    public:

      /////////////////////////////////////////////////////////////////////
      // Base class for calculations of absorption in materials, for now //
      // adding no additional methods over the ones from the Process     //
      // class.                                                          //
      /////////////////////////////////////////////////////////////////////

    protected:
      virtual ~Absorption();
    private:
      Absorption(const char* calculator_type_name);
      friend class AbsorptionImplWrapper;
    };

    class NCRYSTAL_API Scatter : public Process {
    public:

      /////////////////////////////////////////////////////////////////////
      // Base class for calculations of scattering in materials, adding  //
      // generateScattering methods in addition to the crossSection      //
      // methods from the Process class.                                 //
      //                                                                 //
      // Note that the unit for kinetic energy in the calls below is eV  //
      // (electronvolt) and angles are in radians.                       //
      /////////////////////////////////////////////////////////////////////

      //Assuming a scattering took place, the following method generate an energy
      //transfer and new direction for the neutron.
      virtual void generateScattering( double ekin, const double (&neutron_direction)[3],
                                       double (&resulting_neutron_direction)[3], double& delta_ekin ) const = 0;

      //For non-oriented scatter calculators, callers can either use the methods
      //above to access cross-sections and scatterings (with free choice of
      //reference frame for the direction vectors), or use the following methods
      //for convenience (note that the returned angle is the angle from incoming
      //to outgoing neutron):
      virtual void generateScatteringNonOriented( double ekin,
                                                  double& angle, double& delta_ekin ) const;

    protected:
      virtual ~Scatter();
    private:
      Scatter(const char * calculator_type_name);
      friend class ScatterComp;
      friend class ScatterImplWrapper;
    };

    class NCRYSTAL_API ScatterComp final : public Scatter {
    public:

      /////////////////////////////////////////////////////////////////////
      // Composition class which combines a list of scatter calculators  //
      // into one, more complete picture.                                //
      /////////////////////////////////////////////////////////////////////

      ScatterComp(const char * calculator_type_name = "ScatterComp");

      //Must add at least one component:
      void addComponent(Scatter*, double scale = 1.0 );

      size_t nComponents() const { return m_calcs.size(); }
      const Scatter * component(size_t i) const { return m_calcs.at(i).scatter; }
      Scatter * component(size_t i) { return m_calcs.at(i).scatter; }
      double scale(size_t i) const { return m_calcs.at(i).scale; }

      virtual double crossSection(double ekin, const double (&neutron_direction)[3] ) const;

      virtual void domain(double& ekin_low, double& ekin_high) const {
        ekin_low = m_threshold_lower; ekin_high = m_threshold_upper;
      }

      virtual void generateScattering( double ekin, const double (&neutron_direction)[3],
                                       double (&resulting_neutron_direction)[3], double& delta_ekin ) const;

      virtual bool isOriented() const;

      //Note about exception safety: In case of errors, addComponent(scat,..)
      //might throw exceptions, but in this case it will always ref+unref the
      //passed scat object. Thus placing components directly sc->addComponent(new
      //Something,..) is exception safe RAII.

    protected:
      virtual ~ScatterComp();
      struct Component {
        double threshold_lower;
        double threshold_upper;
        double scale;
        Scatter* scatter;
        bool operator<(const Component&) const;
      };
      std::vector<Component> m_calcs;//thresholds and calcs
      double m_threshold_lower;
      double m_threshold_upper;
      mutable int m_isOriented;
      void checkIsOriented() const;
    };

//Convenience defines:
#ifndef NC_VECTOR_CAST
#  define NC_VECTOR_CAST(v) (reinterpret_cast<double(&)[3]>(v))
#  define NC_CVECTOR_CAST(v) (reinterpret_cast<const double(&)[3]>(v))
#endif

    NCRYSTAL_API RCHolder<const Scatter> wrapModernProcPtrInLegacyScatterClass( ::NCrystal::ProcImpl::ProcPtr );
    NCRYSTAL_API RCHolder<const Absorption> wrapModernProcPtrInLegacyAbsorptionClass( ::NCrystal::ProcImpl::ProcPtr );

  }

}

#endif
