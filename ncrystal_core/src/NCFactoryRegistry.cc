////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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
#include "NCString.hh"
#include "NCrystal/NCDefs.hh"
#include <iostream>
#include <cstdlib>

namespace NCrystal {

  static bool s_debug_factoryreg = (std::getenv("NCRYSTAL_DEBUGFACTORY") ? true : false);

  struct Registry {
    FactoryList factories;
    ~Registry()
    {
      for (size_t i = 0; i < factories.size(); ++i)
        delete factories.at(i);
    }
  };
  static Registry * s_factory_registry = 0;
  Registry* getRegistry()
  {
    if (!s_factory_registry)
      s_factory_registry = new Registry;
    return s_factory_registry;
  }

}

bool NCrystal::hasFactory(const char* factoryname)
{
  if (!s_factory_registry)
    return false;
  nc_assert(factoryname);
  const std::string fn(factoryname);
  const FactoryList& fl = getFactories();
  FactoryList::const_iterator it(fl.begin()), itE(fl.end());
  for (;it!=itE;++it)
    if ( fn == (*it)->getName() )
      return true;
  return false;
}

void NCrystal::registerFactory(FactoryBase*f)
{
  nc_assert_always(f);
  if (s_debug_factoryreg)
    std::cout<<"NCrystal::registerFactory - called with factory named \""<<f->getName()<<"\""<<std::endl;

  //Factory names must be unique, and be simple in nature, since they have to
  //specifiable in NCMatCfg input strings.
  if (hasFactory(f->getName()))
    NCRYSTAL_THROW2(BadInput,"Factory named \""<<f->getName()<<"\" already registered");
  if (contains_any(f->getName()," \"'|><(){}[]")||!isSimpleASCII(f->getName()))
    NCRYSTAL_THROW2(BadInput,"Forbidden characters in factory name: \""<<f->getName()<<"\"");
  if (s_debug_factoryreg)
    std::cout<<"NCrystal::registerFactory - no obvious errors found - adding factory succesfully"<<std::endl;
  getFactories();//trigger inbuilt registry calls
  getRegistry()->factories.push_back(f);
}

//Automatic enablement of .nxs/.laz support is controlled via the
//NCRYSTAL_ENABLE_NXSLAZ macro, and the nxs/laz factories must obviously have
//been built for this to work. Support for .ncmat files is on the other hand
//ALWAYS enabled by default, unless disabled with the macro
//NCRYSTAL_DISABLE_NCMAT. Likewise, the standard scatter/absorption factories
//are always enabled, unless respectively NCRYSTAL_DISABLE_STDSCAT and
//NCRYSTAL_DISABLE_STDABS is defined.


#ifndef NCRYSTAL_DISABLE_STDSCAT
extern "C" void ncrystal_register_stdscat_factory();
#endif
#ifndef NCRYSTAL_DISABLE_STDABS
extern "C" void ncrystal_register_stdabs_factory();
#endif
#ifndef NCRYSTAL_DISABLE_NCMAT
extern "C" void ncrystal_register_ncmat_factory();
#endif
#ifdef NCRYSTAL_ENABLE_NXSLAZ
extern "C" void ncrystal_register_nxslaz_factories();
#endif

NCrystal::FactoryList& NCrystal::getFactories()
{
  static bool first(true);
  if (first) {
    first = false;
#ifndef NCRYSTAL_DISABLE_STDSCAT
    ncrystal_register_stdscat_factory();
#endif
#ifndef NCRYSTAL_DISABLE_STDABS
    ncrystal_register_stdabs_factory();
#endif
#ifndef NCRYSTAL_DISABLE_NCMAT
    ncrystal_register_ncmat_factory();
#endif
#ifdef NCRYSTAL_ENABLE_NXSLAZ
    ncrystal_register_nxslaz_factories();
#endif
  }
  return getRegistry()->factories;
}

void NCrystal::clearFactoryRegistry()
{
  if (s_debug_factoryreg)
    std::cout<<"NCrystal::registerFactory clearFactoryRegistry called - removing all factories and cleaing all Info caches"<<std::endl;
  clearInfoCaches();
  delete s_factory_registry;
  s_factory_registry = 0;
}

NCrystal::FactoryBase::~FactoryBase()
{
}

const NCrystal::Info * NCrystal::FactoryBase::createInfo( const MatCfg& ) const
{
  NCRYSTAL_THROW(LogicError,"Unsupported createInfo request");
}

const NCrystal::Scatter * NCrystal::FactoryBase::createScatter( const MatCfg& ) const
{
  NCRYSTAL_THROW(LogicError,"Unsupported createScatter request");
}

const NCrystal::Absorption * NCrystal::FactoryBase::createAbsorption( const MatCfg& ) const
{
  NCRYSTAL_THROW(LogicError,"Unsupported createAbsorption request");
}
