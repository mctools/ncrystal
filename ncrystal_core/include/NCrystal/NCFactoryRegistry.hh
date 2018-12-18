#ifndef NCrystal_FactoryRegistry_hh
#define NCrystal_FactoryRegistry_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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
#include <vector>

namespace NCrystal {

  //Interface for dynamically registering code capable of transforming MatCfg
  //objects into Info, Scatter or Absorption objects.

  //Note: This interface is only likely to be of interest for NCrystal
  //developers and experts wishing to extend the capabilities of NCrystal with
  //support for new file formats or algorithms.

  class MatCfg;
  class Info;
  class Scatter;
  class Absorption;

  class NCRYSTAL_API FactoryBase {
  public:
    virtual const char * getName() const = 0;
    FactoryBase(){}
    virtual ~FactoryBase();

    //The canCreate methods should return 0 if unable to respond to a request,
    //otherwise they should return a non-zero priority. The available factory
    //returning the highest non-zero priority will subsequently be asked to
    //actually create the object in question:

    virtual int canCreateInfo( const MatCfg& ) const { return 0; }
    virtual int canCreateScatter( const MatCfg& ) const { return 0; }
    virtual int canCreateAbsorption( const MatCfg& ) const { return 0; }

    //Methods actually used to create the desired object, once the factory was
    //selected for the task (the factory should throw a LogicError exception if
    //unable to service the request, since it should only ever be called if the
    //corresponding canCreate method returned a non-zero value):

    virtual const Info * createInfo( const MatCfg& ) const;
    virtual const Scatter * createScatter( const MatCfg& ) const;
    virtual const Absorption * createAbsorption( const MatCfg& ) const;

  };

  //Methods used to actually register factories. The factories will subsequently
  //be owned by the NCrystal factory registry, and if it is required (for
  //valgrind studies for instance) clearFactoryRegistry() can be called.
  NCRYSTAL_API void registerFactory(FactoryBase*);
  typedef std::vector<const FactoryBase*> FactoryList;
  NCRYSTAL_API FactoryList& getFactories();//Access factories
  NCRYSTAL_API void clearFactoryRegistry();//clears factories and info caches
  NCRYSTAL_API bool hasFactory(const char*);

}

#endif
