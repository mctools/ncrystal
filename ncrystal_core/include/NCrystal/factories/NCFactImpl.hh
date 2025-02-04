#ifndef NCrystal_FactImpl_hh
#define NCrystal_FactImpl_hh

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

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/factories/NCFactRequests.hh"
#include "NCrystal/factories/NCFactTypes.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/text/NCTextData.hh"

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Detailed factory infrastructure used to produce text and info objects, as   //
// well as  physics-process objects which leaves the management of RNG streams //
// and CachePtrs to the caller (see NCProcImpl.hh). This is the underlying     //
// form used to represent all physics models in NCrystal, and any plugins      //
// extending NCrystal with new physics capabilities should register factories  //
// here. In addition, factory infrastructure for Info objects is likewise kept //
// here, as is those concerning text data (such as input files). For the       //
// latter, more convenient modification functions are available in the file    //
// NCDataSources.hh.                                                           //
//                                                                             //
// The "unmanaged" physics processes dealt with here are suitable for clients  //
// wishing to deal with caching and RNG streams in the context of advanced     //
// application hooks or multi-threaded programming. Users simply wishing to    //
// access the physics in a single-threaded context can instead use the factory //
// methods in NCFact.hh which provide "managed" physics processes (see         //
// NCProc.hh).                                                                 //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace FactImpl {

    //Use registered factories:
    NCRYSTAL_API shared_obj<const TextData>          createTextData( const TextDataPath& );
    NCRYSTAL_API shared_obj<const Info>              createInfo( const MatCfg& );
    NCRYSTAL_API shared_obj<const ProcImpl::Process> createScatter( const MatCfg& );
    NCRYSTAL_API shared_obj<const ProcImpl::Process> createAbsorption( const MatCfg& );

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // The rest of this file concerns interfaces needed by code extending the //
    // capabilities of NCrystal with new factories. For the case of managing  //
    // sources of input data (e.g. data files) see also NCDataSources.hh.     //
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    //Create objects via derived requests:
    NCRYSTAL_API shared_obj<const Info>              createInfo( const InfoRequest& );
    NCRYSTAL_API shared_obj<const ProcImpl::Process> createScatter( const ScatterRequest& );
    NCRYSTAL_API shared_obj<const ProcImpl::Process> createAbsorption( const AbsorptionRequest& );

    //Same, but resolved via overload resolution:
    NCRYSTAL_API shared_obj<const Info>              create( const InfoRequest& );
    NCRYSTAL_API shared_obj<const ProcImpl::Process> create( const ScatterRequest& );
    NCRYSTAL_API shared_obj<const ProcImpl::Process> create( const AbsorptionRequest& );

    ///////////////////////////
    // Factory base classes: //
    ///////////////////////////

    enum class FactoryType { TextData, Info, Scatter, Absorption  };

    template<FactoryType factoryType,
             class TProduced,
             class TKey,
             class TProdRV = shared_obj<const TProduced> >
    class NCRYSTAL_API Factory : private MoveOnly {
    public:
      using Priority = ::NCrystal::Priority;
      using factory_type = FactoryType;
      using produced_type = TProduced;
      using key_type = TKey;
      virtual const char * name() const noexcept = 0;

      //Analyse cfg and determine if, and if so with what priority, the factory
      //can service the request:
      virtual Priority query( const key_type& ) const = 0;

      //Service the request and produce the requested object (the framework will
      //only call this method after a previous call to query(..) indicated that
      //the production is possible - i.e. did not return Priority::Unable):
      virtual TProdRV produce( const key_type& ) const = 0;

      virtual ~Factory() = default;
    };

    //Text data factories:
    class NCRYSTAL_API TextDataFactory
      : public Factory<FactoryType::TextData,TextDataSource,TextDataPath,TextDataSource>
    {
    public:
      using Factory::Factory;
      //In addition to the usual Factory methods, TextData factories can
      //optionally be browsed, providing list of available entries
      //(cf. listAvailableFiles() fct in NCDataSources.hh):
      struct NCRYSTAL_API BrowseEntry {
        std::string name;
        std::string source;//directory path or description
        Priority priority;//the priority it can be delivered with.
      };
      using BrowseList = std::vector<BrowseEntry>;
      virtual BrowseList browse() const = 0;
    };


    //Info factories operates on TextData objects and the reduced set of
    //parameters available in InfoRequests, and always represent a single phase
    //at a time at the configuration level (the FactImpl infrastructure will
    //take care of combinations for multiple phases). The can of course return
    //multiphased Info objects (e.g. a single mysansspheres.ncmat file might
    //result in a two-phased material being loaded):
    using InfoFactory = Factory<FactoryType::Info,Info,InfoRequest>;

    //Analogously, scatter and absorption factories operates on Info objects and
    //the reduced set of parameters available in ScatterRequests and Absorption
    //requests, respectively. Additionally, they must declare capabilities for
    //handling multi- and single-phase configurations (most plugins with the
    //exceptions of SANS plugins should provide SPOnly factories). Additionally,
    //scatter factories also gets a few convenience methods added, which are
    //useful for a plugin factory wanting to add new physics on top of existing
    //capabilities.

    enum class MultiPhaseCapability { MPOnly, SPOnly, Both };

    class NCRYSTAL_API AbsorptionFactory : public Factory<FactoryType::Absorption,ProcImpl::Process,AbsorptionRequest> {
    public:
      using Factory::Factory;
      using ProcPtr = ProcImpl::ProcPtr;
      using MultiPhaseCapability = FactImpl::MultiPhaseCapability;
      virtual MultiPhaseCapability multiPhaseCapability() const { return MultiPhaseCapability::SPOnly; }
    };

    class NCRYSTAL_API ScatterFactory : public Factory<FactoryType::Scatter,ProcImpl::Process,ScatterRequest> {
    public:
      using Factory::Factory;
      using ProcPtr = ProcImpl::ProcPtr;
      using MultiPhaseCapability = FactImpl::MultiPhaseCapability;

      virtual MultiPhaseCapability multiPhaseCapability() const { return MultiPhaseCapability::SPOnly; }

      //Wraps ::NCrystal::FactImpl::createScatter but excludes our own factory for consideration:
      ProcImpl::ProcPtr globalCreateScatter( const ScatterRequest& ) const;

      //Combine multiple ProcPtr's into one:
      template<typename... Args>
      static ProcImpl::ProcPtr combineProcs(Args &&... args);
    };

    //////////////////////////////////
    // Add to the factory registry: //
    //////////////////////////////////

    NCRYSTAL_API void registerFactory( std::unique_ptr<const TextDataFactory> );
    NCRYSTAL_API void registerFactory( std::unique_ptr<const InfoFactory> );
    NCRYSTAL_API void registerFactory( std::unique_ptr<const ScatterFactory> );
    NCRYSTAL_API void registerFactory( std::unique_ptr<const AbsorptionFactory> );

    /////////////////////////////////
    // Query the factory registry: //
    /////////////////////////////////

    //The following functions will cause initial plugin list to be loaded if not
    //already. Therefore they should NOT be used when registering factories in
    //plugins (use the currentlyHasFactory function below instead if needed):
    NCRYSTAL_API bool hasFactory( FactoryType, const std::string& name );
    NCRYSTAL_API bool hasTextDataFactory( const std::string& name );
    NCRYSTAL_API bool hasInfoFactory( const std::string& name );
    NCRYSTAL_API bool hasScatterFactory( const std::string& name );
    NCRYSTAL_API bool hasAbsorptionFactory( const std::string& name );

    NCRYSTAL_API std::vector<shared_obj<const TextDataFactory>> getTextDataFactoryList();
    NCRYSTAL_API std::vector<shared_obj<const InfoFactory>> getInfoFactoryList();
    NCRYSTAL_API std::vector<shared_obj<const ScatterFactory>> getScatterFactoryList();
    NCRYSTAL_API std::vector<shared_obj<const AbsorptionFactory>> getAbsorptionFactoryList();

    NCRYSTAL_API bool currentlyHasFactory( FactoryType, const std::string& name );

    ////////////////////
    // Miscellaneous: //
    ////////////////////

    //Utility function for guessing data type from data (optionally also from
    //the file name). Returns empty string if unable. We might eventually make
    //this extendable so plugins can provide new type recognition capabilities:
    NCRYSTAL_API std::string guessDataType( const RawStrData&, const std::string& filename = {} );

    //Advanced usage (for NCDataSources.cc, not recommended for other usage):
    NCRYSTAL_API void removeTextDataFactoryIfExists( const std::string& name );

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {
    template<typename... Args>
    inline ProcImpl::ProcPtr ScatterFactory::combineProcs(Args &&... args)
    {
      return ProcImpl::ProcComposition::combine(std::forward<Args>(args)...);
    }
    inline bool hasTextDataFactory( const std::string& name ) { return hasFactory(FactoryType::TextData,name); }
    inline bool hasInfoFactory( const std::string& name ) { return hasFactory(FactoryType::Info,name); }
    inline bool hasScatterFactory( const std::string& name ) { return hasFactory(FactoryType::Scatter,name); }
    inline bool hasAbsorptionFactory( const std::string& name ) { return hasFactory(FactoryType::Absorption,name); }

    inline shared_obj<const Info> create( const InfoRequest& req ) { return createInfo(req); }
    inline shared_obj<const ProcImpl::Process> create( const ScatterRequest& req ) { return createScatter(req); }
    inline shared_obj<const ProcImpl::Process> create( const AbsorptionRequest& req ) { return createAbsorption(req); }


  }
}

#endif
