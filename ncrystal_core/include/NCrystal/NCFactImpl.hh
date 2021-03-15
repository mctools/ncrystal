#ifndef NCrystal_FactImpl_hh
#define NCrystal_FactImpl_hh

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

#include "NCrystal/NCMatCfg.hh"
#include "NCrystal/NCFactTypes.hh"
#include "NCrystal/NCProcImpl.hh"
#include "NCrystal/NCMatInfo.hh"
#include "NCrystal/NCTextData.hh"

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Detailed factory infrastructure used to produce text and info objects, but  //
// also physics-process objects which leaves the management of RNG streams     //
// and CachePtrs to the caller (see NCProcImpl.hh). This is the underlying     //
// form used to represent all physics models in NCrystal, and any plugins      //
// extending NCrystal with new physics capabilities should register factories  //
// here. In addition, factory infrastructure for MatInfo objects is likewise,  //
// kept here, as is those concerning text data (such as input files). For the  //
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

namespace NCrystal {

  namespace FactImpl {

    //Use registered factories:
    NCRYSTAL_API shared_obj<const TextData>          createTextData( const TextDataPath& cfg );
    NCRYSTAL_API shared_obj<const MatInfo>           createInfo( const MatCfg& cfg );
    NCRYSTAL_API shared_obj<const ProcImpl::Process> createScatter( const MatCfg& cfg );
    NCRYSTAL_API shared_obj<const ProcImpl::Process> createAbsorption( const MatCfg& cfg );

    //Disable and enable caching in these factories (default state upon startup
    //is for caching to be enabled, unless the environment variable
    //NCRYSTAL_NOCACHE is set):
    NCRYSTAL_API void setCachingEnabled(bool);
    NCRYSTAL_API bool getCachingEnabled();

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // The rest of this file concerns interfaces needed by code extending the //
    // capabilities of NCrystal with new factories. For the case of managing  //
    // sources of input data (e.g. data files) see also NCDataSources.hh.     //
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


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

    //Info and Absorption factories are simple:
    using InfoFactory = Factory<FactoryType::Info,MatInfo,MatInfoCfg>;
    using AbsorptionFactory = Factory<FactoryType::Absorption,ProcImpl::Process,MatCfg>;

    class NCRYSTAL_API TextDataFactory
      : public Factory<FactoryType::TextData,TextDataSource,TextDataPath,TextDataSource>
    {
    public:
      using Factory::Factory;
      //In addition usual Factory methods, TextData factories can also
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

    //Scatter factories gets a few convenience methods:
    class NCRYSTAL_API ScatterFactory : public Factory<FactoryType::Scatter,ProcImpl::Process,MatCfg> {
    public:
      using Factory::Factory;
      using ProcPtr = ProcImpl::ProcPtr;

      //Wraps ::NCrystal::FactImpl::createScatter but avoids our own factory for consideration:
      ProcImpl::ProcPtr globalCreateScatter( const MatCfg&, bool allowself=false ) const;

      //Scatter factories often need to access MatInfo, so make it easily available:
      static shared_obj<const MatInfo> globalCreateInfo( const MatCfg& );
      static shared_obj<const MatInfo> createInfo( const MatCfg& );

      //Combine multiple ProcPtr's into one:
      template<typename... Args>
      static ProcImpl::ProcPtr combineProcs(Args &&... args);
    };

    /////////////////////////////////////////
    // Access/manipulate factory registry: //
    /////////////////////////////////////////

    enum class RegPolicy { ERROR_IF_EXISTS, OVERRIDE_IF_EXISTS, IGNORE_IF_EXISTS };
    NCRYSTAL_API void registerFactory( std::unique_ptr<const TextDataFactory>, RegPolicy = RegPolicy::OVERRIDE_IF_EXISTS );
    NCRYSTAL_API void registerFactory( std::unique_ptr<const InfoFactory>, RegPolicy = RegPolicy::OVERRIDE_IF_EXISTS );
    NCRYSTAL_API void registerFactory( std::unique_ptr<const ScatterFactory>, RegPolicy = RegPolicy::OVERRIDE_IF_EXISTS );
    NCRYSTAL_API void registerFactory( std::unique_ptr<const AbsorptionFactory>, RegPolicy = RegPolicy::OVERRIDE_IF_EXISTS );

    NCRYSTAL_API bool hasFactory( FactoryType, const std::string& name );
    NCRYSTAL_API bool hasTextDataFactory( const std::string& name );
    NCRYSTAL_API bool hasInfoFactory( const std::string& name );
    NCRYSTAL_API bool hasScatterFactory( const std::string& name );
    NCRYSTAL_API bool hasAbsorptionFactory( const std::string& name );

    NCRYSTAL_API std::vector<shared_obj<const TextDataFactory>> getTextDataFactoryList();
    NCRYSTAL_API std::vector<shared_obj<const InfoFactory>> getInfoFactoryList();
    NCRYSTAL_API std::vector<shared_obj<const ScatterFactory>> getScatterFactoryList();
    NCRYSTAL_API std::vector<shared_obj<const AbsorptionFactory>> getAbsorptionFactoryList();

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

namespace NCrystal {
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
  }
}

#endif
