#ifndef NCrystal_FactoryRegistry_hh
#define NCrystal_FactoryRegistry_hh

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

#include "NCrystal/NCDefs.hh"

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

    virtual RCHolder<const Info> createInfo( const MatCfg& ) const;
    virtual RCHolder<const Scatter> createScatter( const MatCfg& ) const;
    virtual RCHolder<const Absorption> createAbsorption( const MatCfg& ) const;

  protected:

    //TODO: We should rename the above methods or split the factory base into
    //three parts. The name conflicts with the global NCrystal::createInfo
    //etc. functions. For now, we provide these convenience functions which
    //redirects to the global functions (after inserting the current factory
    //into the exclusion list of the relevant xxxfactory cfg parameter, unless
    //allowself=true - MT TODO: support this for info facts as well!):

    RCHolder<const Info> globalCreateInfo( const MatCfg& cfg ) const;
    RCHolder<const Scatter> globalCreateScatter( const MatCfg& cfg, bool allowself=false ) const;
    RCHolder<const Absorption> globalCreateAbsorption( const MatCfg& cfg, bool allowself=false ) const;

    //Combine two Scatter objects into one. This should only be used with
    //freshly created Scatter objects in a factory's createScatter method.
    //(this is due to a design flaw which is planned to be fixed in a future
    //release):
    static RCHolder<const Scatter> combineScatterObjects( RCHolder<const Scatter> sc1,
                                                          RCHolder<const Scatter> sc2 );
  };

  //Methods used to register factories and query available factories:
  NCRYSTAL_API void registerFactory(std::unique_ptr<FactoryBase>);
  typedef std::vector<std::shared_ptr<const FactoryBase>> FactoryList;
  NCRYSTAL_API const FactoryList& getFactories();//Access factories
  NCRYSTAL_API bool hasFactory(const char*);

}

#endif
