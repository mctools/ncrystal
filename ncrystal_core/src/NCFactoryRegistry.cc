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

#include "NCrystal/NCFactoryRegistry.hh"
#include "NCrystal/NCFactory.hh"
#include "NCrystal/NCPluginMgmt.hh"
#include "NCrystal/internal/NCString.hh"
#include <iostream>
#include <cstdlib>

namespace NC = NCrystal;

namespace NCrystal {

  namespace {
    static std::atomic<bool> s_debug_factoryreg( ncgetenv_bool("DEBUGFACTORY") );

    struct Registry : public NoCopyMove {
      FactoryList factories;
    };
    Registry& getRegistry()
    {
      static Registry s_reg;
      return s_reg;
    }
    std::mutex& getRegistryMutex()
    {
      static std::mutex s_mtx;
      return s_mtx;
    }
  }
}



bool NC::hasFactory(const char* factoryname)
{
  nc_assert(factoryname);
  const std::string fn(factoryname);
  std::lock_guard<std::mutex> guard(getRegistryMutex());
  for (auto fptr : getRegistry().factories) {
    if ( fn == fptr->getName() )
      return true;
  }
  return false;
}

void NC::registerFactory(std::unique_ptr<NC::FactoryBase> f)
{
  nc_assert_always(!!f);
  Plugins::ensurePluginsLoaded();
  if (s_debug_factoryreg)
    std::cout<<"NCrystal::registerFactory - called with factory named \""<<f->getName()<<"\""<<std::endl;

  //Factory names must be unique, and be simple in nature, since they have to
  //specifiable in NCMatCfg input strings.
  if (contains_any(f->getName()," \"'|><(){}[]")||!isSimpleASCII(f->getName()))
    NCRYSTAL_THROW2(BadInput,"Forbidden characters in factory name: \""<<f->getName()<<"\"");
  if (hasFactory(f->getName()))
    NCRYSTAL_THROW2(BadInput,"Factory named \""<<f->getName()<<"\" already registered");
  if (s_debug_factoryreg)
    std::cout<<"NCrystal::registerFactory - no obvious errors found - adding factory succesfully"<<std::endl;
  std::lock_guard<std::mutex> guard(getRegistryMutex());
  getRegistry().factories.push_back(std::move(f));
}

const NC::FactoryList& NC::getFactories()
{
  Plugins::ensurePluginsLoaded();
  std::lock_guard<std::mutex> guard(getRegistryMutex());
  return getRegistry().factories;//MT TODO: not safe, we should use shared ptrs instead.
}

NC::FactoryBase::~FactoryBase()
{
}

NC::RCHolder<const NC::Info> NC::FactoryBase::createInfo( const MatCfg& ) const
{
  NCRYSTAL_THROW(LogicError,"Unsupported createInfo request");
}

NC::RCHolder<const NC::Scatter> NC::FactoryBase::createScatter( const MatCfg& ) const
{
  NCRYSTAL_THROW(LogicError,"Unsupported createScatter request");
}

NC::RCHolder<const NC::Absorption> NC::FactoryBase::createAbsorption( const MatCfg& ) const
{
  NCRYSTAL_THROW(LogicError,"Unsupported createAbsorption request");
}

#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCScatter.hh"
#include "NCrystal/NCAbsorption.hh"

NC::RCHolder<const NC::Info> NC::FactoryBase::globalCreateInfo( const NC::MatCfg& cfg ) const
{
  return RCHolder<const Info>(::NC::createInfo(cfg));
}

NC::RCHolder<const NC::Scatter> NC::FactoryBase::globalCreateScatter( const NC::MatCfg& cfg, bool allowself ) const
{
  auto cfg2 = cfg.clone();
  if ( ! allowself ) {
    std::string ourname=getName();
    auto factrequest = cfg2.get_scatfactory_parsed();
    factrequest.excluded.insert(ourname);
    if (factrequest.specific==ourname)
      factrequest.specific.clear();
    cfg2.set_scatfactory(factrequest);
  }
  return RCHolder<const Scatter>(::NC::createScatter(cfg2));
}

NC::RCHolder<const NC::Absorption> NC::FactoryBase::globalCreateAbsorption( const NC::MatCfg& cfg, bool allowself ) const
{
  auto cfg2 = cfg.clone();
  if ( ! allowself ) {
    std::string ourname=getName();
    auto factrequest = cfg2.get_absnfactory_parsed();
    factrequest.excluded.insert(ourname);
    if (factrequest.specific==ourname)
      factrequest.specific.clear();
    cfg2.set_absnfactory(factrequest);
  }
  return RCHolder<const Absorption>(::NC::createAbsorption(cfg2));
}

//Temporary workaround function, until we redesign the interfaces:
#include "NCrystal/NCScatterComp.hh"
NC::RCHolder<const NC::Scatter> NC::FactoryBase::combineScatterObjects( NC::RCHolder<const NC::Scatter> sc1,
                                                                        NC::RCHolder<const NC::Scatter> sc2 )
{
  if (sc1->isNull())
    return sc2;
  if (sc2->isNull())
    return sc1;
  auto sc_comp = makeRC<ScatterComp>();
  sc_comp->addComponent(const_cast<Scatter*>(sc1.obj()));
  sc_comp->addComponent(const_cast<Scatter*>(sc2.obj()));
  return sc_comp;
}
