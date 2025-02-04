#ifndef NCrystal_ProcImpl_hh
#define NCrystal_ProcImpl_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/interfaces/NCRNG.hh"

///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
// Interfaces used for the implementation of physics processes. Implementations  //
// should derive from one of the following base classes for scattering or        //
// absorption processes:                                                         //
//                                                                               //
//   ScatterIsotropicMat   (isotropic materials)                                 //
//   ScatterAnisotropicMat (anisotropic materials such as single crystals)       //
//   AbsorptionIsotropicMat (isotropic materials)                                //
//                                                                               //
// Note: We have not yet added an AbsorptionAnisotropicMat class, but it could   //
// of course be done if ever needed.                                             //
//                                                                               //
// Also note that the term "isotropic material" refers to a property of the      //
// material, it does not mean that neutrons will necessarily be scattered        //
// isotropically. Also note that all single crystals are by definition           //
// anisotropic in the sense used here, even those laid out according to an ideal //
// Gaussian mosaicity distribution.                                              //
//                                                                               //
// A few rules should be obeyed in all implementations:                          //
//                                                                               //
// 1) All methods on classes should in general be marked const (exception is     //
//    move-assignment operator).                                                 //
// 2) Never use mutables and statics (exceptions: statics inside functions       //
//    are allowed by experts if *really* careful about MT safety and avoiding    //
//    keeping any info which would amount to "hidden-state".                     //
// 3) Non-mutable data members can be used for process data (e.g. reflection     //
//    plane lists or scattering kernels). Caches for neutron state parameters    //
//    should use the CachePtr objects (see next item).                           //
// 4) If your class is marked "final", you are allowed customised caching        //
//    which can depend on neutron state (e.g. list of reflection planes and      //
//    their contribution to the cross section for a given neutron energy and/or  //
//    direction). Add cache data on a custom subclass of CacheBase, and access   //
//    an instance of it by calling:                                              //
//        auto& cache = accessCache<MyCacheClass>(cachePtr).                     //
//    The reason for doing this is that calling code will take care of managing  //
//    issues such as cache lifetimes and having different caches in each thread  //
//    of a multi-threaded programme.                                             //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace ProcImpl {

    ///////////////////////////////////////////////////////////////////////////////
    //
    // Common base class of all processes.
    //
    class NCRYSTAL_API Process : private MoveOnly {
    public:
      virtual ~Process() = default;

      //Name of process (for debug output, should usually match name of C++
      //class implementing the process):
      virtual const char * name() const noexcept = 0;

      //Unique instance ID:
      UniqueIDValue getUniqueID() const noexcept { return m_uniqueID.getUniqueID(); }

      //Isotropic state of material associated with process:
      virtual MaterialType materialType() const noexcept = 0;

      bool isOriented() const noexcept { return materialType()==MaterialType::Anisotropic; }

      //Does process describe scattering or absorption:
      virtual ProcessType processType() const noexcept = 0;

      //Domain of process, i.e. range of neutron kinetic energies for which
      //process might have non-zero cross section (default implementation
      //returns [0,infinity]):
      virtual EnergyDomain domain() const noexcept;

      //Domain might indicate that this is a null-process, vanishing everywhere:
      bool isNull() const noexcept;

      //In some cases it might be possible, and desirable from an efficiency
      //POV, for a process to be merged with another (usually of its own type)
      //into a single instance, possibly with scaling factors. Classes
      //supporting this can reimplement this method (returns nullptr if merging
      //is not possible, which is why we return std::shared_ptr and not
      //shared_obj):
      virtual std::shared_ptr<Process> createMerged( const Process& other,
                                                     double scale_other,
                                                     double scale_self ) const;

      //The next four functions implement cross section calculations and
      //scattering samplings. Depending on material and process type, some of
      //them might be unavailable and will throw an exception if called:

      //It is always allowed to ask for direction-dependent cross section:
      virtual CrossSect crossSection(CachePtr&, NeutronEnergy, const NeutronDirection& ) const = 0;

      //This can be called only when materialType is Isotropic:
      virtual CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const = 0;

      //This can be called only when processType is Scatter:
      virtual ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const = 0;

      //This can be called only when processType is Scatter and materialType is Isotropic:
      virtual ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const = 0;

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      //Experimental ABI-breaking additions. For now, it is best to use
      //functions in NCABIUtils.hh to access these.

      //Advanced cross section methods more suitable for vectorisation. These do
      //not use the same level of type-safety, and any passed arrays are
      //required to be non-aliasing (i.e. not overlap in memory):
      virtual void evalManyXS( CachePtr&, const double* ekin,
                               const double* ux, const double* uy, const double* uz,
                               std::size_t N, double* out_xs ) const = 0;
      virtual void evalManyXSIsotropic( CachePtr&, const double* ekin, std::size_t N,
                                        double* out_xs ) const = 0;

      //Elastic scattering processes that *never* change the neutrons energy can
      //advertise themselves as such (this is optional but highly recommended).
      virtual bool isPureElasticScatter() const = 0;

      //If you know you need both the cross section and to sample a scattering,
      //it is usually more efficient to request both at once:
      virtual std::pair<CrossSect,ScatterOutcome>
        evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const = 0;

      virtual std::pair<CrossSect,ScatterOutcomeIsotropic>
        evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const = 0;

#endif
      //Callers who would like to reduce allocations during their event loops,
      //can call the following function which is likely (but not 100%
      //guaranteed) to preallocate any required objects in the cache:
      void initCachePtr(CachePtr& cp) const;

      //Summarise meta-data in JSON dictionary:
      std::string jsonDescription() const;

    protected:
      //For unpacking CachePtr arguments into custom classes (should only be
      //used in derived classes marked "final"):
      template<class CacheClass>
      CacheClass& accessCache(CachePtr& cpbase) const;

      //Optionally override to provide JSON dictionary with specific info
      //relevant only for a given model. It is a good idea if such a dictionary
      //contains an entry with the key named "summarystr", providing a short
      //summary of the process:
      virtual Optional<std::string> specificJSONDescription() const { return NullOpt; };

    private:
      UniqueID m_uniqueID;
      //Restrict direct inheritance from this class to framework classes listed
      //as friends here:
      Process() = default;
      //Base classes for actual physics models:
      friend class ScatterIsotropicMat;
      friend class ScatterAnisotropicMat;
      friend class AbsorptionIsotropicMat;
      //Infrastructure classes (final):
      friend class ProcComposition;
      friend class NullProcess;
    };

    //Convenience:
    using ProcPtr = shared_obj<const Process>;
    using OptionalProcPtr = std::shared_ptr<const Process>;

    ///////////////////////////////////////////////////////////////////////////////
    //
    // Base class of all scattering processes in isotropic materials.
    //
    class NCRYSTAL_API ScatterIsotropicMat : public Process {
    public:
      MaterialType materialType() const noexcept final { return MaterialType::Isotropic; }
      ProcessType processType() const noexcept final { return ProcessType::Scatter; }

      //Isotropic material, anisotropic methods are implemented in terms of the isotropic ones:
      CrossSect crossSection(CachePtr&, NeutronEnergy, const NeutronDirection& ) const final;
      ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      bool isPureElasticScatter() const override;
      void evalManyXS( CachePtr&, const double*, const double*, const double*,
                       const double*, std::size_t, double* ) const override;
      void evalManyXSIsotropic( CachePtr&, const double*,
                                std::size_t, double* ) const override;
      std::pair<CrossSect,ScatterOutcome>
        evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;
      std::pair<CrossSect,ScatterOutcomeIsotropic>
        evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;
#endif

      // NB: Several methods are marked as "override" here instead of "final",
      // since some models might be able to do something more efficient than the
      // default implementation).
    };

    ///////////////////////////////////////////////////////////////////////////////
    //
    // Base class of all scattering processes in anisotropic materials, such as
    // single crystals.
    //
    class NCRYSTAL_API ScatterAnisotropicMat : public Process {
    public:
      MaterialType materialType() const noexcept final { return MaterialType::Anisotropic; }
      ProcessType processType() const noexcept final { return ProcessType::Scatter; }

      //Anisotropic material, it is not allowed to call isotropic methods (will
      //throw LogicError if attempted):
      CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const final;
      ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const final;

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      bool isPureElasticScatter() const override;
      void evalManyXS( CachePtr&, const double*, const double*, const double*,
                       const double*, std::size_t, double* ) const override;
      void evalManyXSIsotropic( CachePtr&, const double*,
                                std::size_t, double* ) const final;
      std::pair<CrossSect,ScatterOutcome>
        evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;
      std::pair<CrossSect,ScatterOutcomeIsotropic>
        evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;
#endif


    };

    ///////////////////////////////////////////////////////////////////////////////
    //
    // Base class of all absorption processes in isotropic materials.
    //
    class NCRYSTAL_API AbsorptionIsotropicMat : public Process {
    public:
      MaterialType materialType() const noexcept final { return MaterialType::Isotropic; }
      ProcessType processType() const noexcept final { return ProcessType::Absorption; }

      //Isotropic material, anisotropic xsect is implemented in terms of the isotropic one:
      CrossSect crossSection(CachePtr& cp, NeutronEnergy ekin, const NeutronDirection& ) const final;

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      void evalManyXS( CachePtr& cp, const double* ekin,
                       const double*, const double*, const double*,
                       std::size_t N, double* out_xs ) const override;
      void evalManyXSIsotropic( CachePtr&, const double* ekin, std::size_t N,
                                double* out_xs ) const override;
#endif

      //Absorption process, it is not allowed to call scattering methods (will
      //throw LogicError if attempted):
      ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const final;
      ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const final;
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      bool isPureElasticScatter() const final;
      std::pair<CrossSect,ScatterOutcome>
        evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const final;
      std::pair<CrossSect,ScatterOutcomeIsotropic>
        evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const final;
#endif
    };

    ///////////////////////////////////////////////////////////////////////////////
    //
    // Composition class. This is a technical class which can be used to
    // additively (and incoherently) combine other processes. Cross sectons of
    // all components will simply be added up (possibly after scaling in case of
    // non-unity scale factors). If applicable, scattering events will be
    // sampled by first randomly selecting a component according to its
    // contribution to the cross section.
    //
    // NB: If using scale factors to combine processes within a single material
    // phase, scale factors will often be either 1.0 or represent atomic
    // fractions (e.g. a free-gas process for Al in Al2O3 might have a scale
    // factor of 2/5=0.4). But when using scale factors to combine processes
    // from different phases, it is important to take the atomic number
    // densities of each phase into account, as counted over the entire material
    // (e.g. two phases might both occupy 50% of the whole volume, but if one
    // has 10 times as many atoms inside then the resulting per-atom cross
    // section should of course reflect that).
    //
    // All processes added to a particular ProcComposition instance must all
    // have same processType, but might mix materialType (the material becomes
    // anisotropic when at least one component is anisotropic).
    //

    class NCRYSTAL_API ProcComposition final : public Process {
    public:

      const char * name() const noexcept final { return "ProcComposition"; }
      ProcessType processType() const noexcept final  { return m_processType; }
      MaterialType materialType() const noexcept final { return m_materialType; }

      struct NCRYSTAL_API Component {
        double scale;
        ProcPtr process;
        Component( double sc, ProcPtr pp ) : scale(sc), process(std::move(pp)) {}
        Component( ProcPtr pp ) : scale(1.0), process(std::move(pp)) {}
      };
      using ComponentList = SmallVector<Component,6>;

      ProcComposition( ComponentList components = ComponentList(),
                       ProcessType processType = ProcessType::Scatter );
      const ComponentList& components() const noexcept;
      void addComponent( ProcPtr process, double scale = 1.0 );
      void addComponents( ComponentList components, double scale = 1.0 );

      EnergyDomain domain() const noexcept final { return m_domain; }
      CrossSect crossSection(CachePtr& cacheptr, NeutronEnergy ekin, const NeutronDirection& dir ) const final;
      CrossSect crossSectionIsotropic(CachePtr& cacheptr, NeutronEnergy ekin ) const final;
      ScatterOutcome sampleScatter(CachePtr& cacheptr, RNG& rng, NeutronEnergy ekin, const NeutronDirection& dir ) const final;
      ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr& cacheptr, RNG& rng, NeutronEnergy ekin ) const final;

      //Combine components efficiently. Often this result in their placement in
      //a ProcComposition object, but in case of no null-processes in the list a
      //NullScatter object is returned instead. And if the list contains only a
      //single unscaled component, that component is returned.
      static ProcPtr combine(const ComponentList&, ProcessType processType = ProcessType::Scatter);
      static ProcPtr consumeAndCombine(ComponentList&&, ProcessType processType = ProcessType::Scatter);

      //Combine (with no scaling, i.e. all scale=1.0):
      template<typename... Args>
      static ProcPtr combine(Args &&... args);

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      bool isPureElasticScatter() const override;
      std::pair<CrossSect,ScatterOutcome>
        evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;
      std::pair<CrossSect,ScatterOutcomeIsotropic>
        evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;
      void evalManyXS( CachePtr&, const double* ekin,
                       const double* ux, const double* uy, const double* uz,
                       std::size_t N,
                       double* out_xs ) const override;
      void evalManyXSIsotropic( CachePtr&, const double* ekin, std::size_t N,
                                double* out_xs ) const override;
#endif

    protected:
      Optional<std::string> specificJSONDescription() const override;
    private:
      unsigned m_nHistory = 1;//increment whenever m_components change.
      ComponentList m_components;
      ProcessType m_processType;
      MaterialType m_materialType;
      EnergyDomain m_domain = { NeutronEnergy{0.0}, NeutronEnergy{0.0} };
      class Impl;
      friend class Impl;
    };

    ///////////////////////////////////////////////////////////////////////////////
    //
    // For technical reasons, it might occasionally be convenient to use
    // null-processes:
    //
    class NullProcess : public Process {
    public:
      MaterialType materialType() const noexcept final { return MaterialType::Isotropic; }
      EnergyDomain domain() const noexcept final { return EnergyDomain::null(); }
      CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const final { return CrossSect{0.0}; }
      CrossSect crossSection(CachePtr&, NeutronEnergy, const NeutronDirection& ) const final { return CrossSect{0.0}; }
      ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const final;
      ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const final;
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      void evalManyXS( CachePtr&, const double* ekin,
                       const double* ux, const double* uy, const double* uz,
                       std::size_t N,
                       double* out_xs ) const override;
      void evalManyXSIsotropic( CachePtr&, const double* ekin, std::size_t N,
                                double* out_xs ) const override;
#endif

    private:
      //Only NullScatter/NullAbsorption can inherit from this class:
      NullProcess() = default;
      friend class NullScatter;
      friend class NullAbsorption;
    };

    class NullScatter final : public NullProcess {
    public:
      const char * name() const noexcept override { return "NullScatter"; }
      ProcessType processType() const noexcept override { return ProcessType::Scatter; }
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      bool isPureElasticScatter() const override { return true; }
      std::pair<CrossSect,ScatterOutcome>
        evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;
      std::pair<CrossSect,ScatterOutcomeIsotropic>
        evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;
#endif
    };

    class NullAbsorption final : public NullProcess {
    public:
      const char * name() const noexcept override { return "NullAbsorption"; }
      ProcessType processType() const noexcept override { return ProcessType::Absorption; }
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      bool isPureElasticScatter() const override { return false; }
      std::pair<CrossSect,ScatterOutcome>
        evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;
      std::pair<CrossSect,ScatterOutcomeIsotropic>
        evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;
#endif
    };

    //Global instances (better caching):
    ProcPtr getGlobalNullScatter();
    ProcPtr getGlobalNullAbsorption();
    ProcPtr getGlobalNullProcess(ProcessType);

  }

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace ProcImpl {
    inline EnergyDomain Process::domain() const noexcept { return { NeutronEnergy{0.0}, NeutronEnergy{kInfinity} }; }
    inline bool Process::isNull() const noexcept { return domain().isNull(); }
    inline std::shared_ptr<Process> Process::createMerged( const Process&, double, double ) const { return nullptr; }

    template<class CacheClass>
    inline CacheClass& Process::accessCache(CachePtr& cpbase) const {
      if (!cpbase)
        cpbase = std::make_unique<CacheClass>();
#ifndef NDEBUG
      if ( !dynamic_cast<CacheClass*>(cpbase.get()) ) {
        NCRYSTAL_THROW2( LogicError, "accessCache: type mismatch in class \"" << name()
                         << "\". Was accessCache<..> used with conflicting cache class types?" );
      }
#endif
      return *static_cast<CacheClass*>(cpbase.get());
    }

    inline CrossSect ScatterIsotropicMat::crossSection( CachePtr& cp, NeutronEnergy ekin, const NeutronDirection& ) const
    {
      return crossSectionIsotropic(cp,ekin);
    }

    inline CrossSect AbsorptionIsotropicMat::crossSection( CachePtr& cp, NeutronEnergy ekin, const NeutronDirection& ) const
    {
      return crossSectionIsotropic(cp,ekin);
    }

    inline const ProcComposition::ComponentList& ProcComposition::components() const noexcept
    {
      return m_components;
    }

    template<typename... Args>
    inline ProcPtr ProcComposition::combine(Args &&... args)
    {
      std::array<ProcPtr, sizeof...(args)> proc_ptrs{std::move(args)...};
      static_assert( sizeof...(args) >= 1, "combine(..) needs at least one argument" );
      ComponentList cl;
      cl.reserve_hint(sizeof...(args));
      for(auto&& pp : proc_ptrs)
        cl.push_back({1.0,std::move(pp)});
      auto proctype = cl.front().process->processType();
      return consumeAndCombine( std::move(cl), proctype );
    }

    inline ProcPtr getGlobalNullProcess(ProcessType pt)
    {
      nc_assert( pt == ProcessType::Scatter || pt == ProcessType::Absorption );
      return ( pt == ProcessType::Scatter ? getGlobalNullScatter() : getGlobalNullAbsorption() );
    }
  }
}

#endif
