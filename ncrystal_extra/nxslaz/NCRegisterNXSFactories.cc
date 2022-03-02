////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2022 NCrystal developers                                   //
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

#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/NCDataSources.hh"
#include "NCrystal/NCMatCfg.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCFactory_NXS.hh"
#include "NCLazLoader.hh"
#include <iostream>

namespace NC = NCrystal;

namespace NCrystal {

  class NXSFactory final : public FactImpl::InfoFactory {
  public:
    const char * name() const noexcept final { return "stdnxs"; }

    Priority query( const FactImpl::InfoRequest& cfg ) const final
    {
      return cfg.getDataType()=="nxs" ? Priority{100} : Priority::Unable;
    }

    shared_obj<const Info> produce( const FactImpl::InfoRequest& cfg ) const final
    {
      nc_assert_always(cfg.getDataType()=="nxs");
      if ( !trim2(cfg.get_atomdb()).empty() )
        std::cout<<"NCrystal WARNING: atomdb parameter is ignored for .nxs files"<<std::endl;
      auto builder = loadNXSCrystal( cfg.textData(),
                                     cfg.get_temp().dbl()==-1.0?Temperature{293.15}:cfg.get_temp(),
                                     cfg.get_dcutoff(),
                                     cfg.get_dcutoffup() );
      builder.dataSourceName = cfg.dataSourceName();
      return buildInfoPtr(std::move(builder));
    }
  };

  class LazFactory final : public FactImpl::InfoFactory {
  public:
    const char * name() const noexcept final { return "stdlaz"; }
    Priority query( const FactImpl::InfoRequest& cfg ) const final
    {
      std::string ext = cfg.getDataType();
      return (ext=="laz"||ext=="lau") ? Priority{100} : Priority::Unable;
    }

    shared_obj<const Info> produce( const FactImpl::InfoRequest& cfg ) const final
    {
      if ( !trim2(cfg.get_atomdb()).empty() )
        std::cout<<"NCrystal WARNING: atomdb parameter is ignored for .laz/.lau files"<<std::endl;
      LazLoader ld (cfg.textData(),
                    cfg.get_dcutoff(),
                    cfg.get_dcutoffup(),
                    cfg.get_temp().dbl()==-1.0?Temperature{293.15}:cfg.get_temp());
      auto builder = ld.read();
      builder.dataSourceName = cfg.dataSourceName();
      return buildInfoPtr(std::move(builder));
    }
  };

}

//Finally, a function which can be used to enable the above factories. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void ncrystal_register_nxslaz_factories()
{
  NC::FactImpl::registerFactory( std::make_unique<NC::NXSFactory>(),
                                 NC::FactImpl::RegPolicy::IGNORE_IF_EXISTS );
  NC::FactImpl::registerFactory( std::make_unique<NC::LazFactory>(),
                                 NC::FactImpl::RegPolicy::IGNORE_IF_EXISTS );
  NC::DataSources::addRecognisedFileExtensions("nxs");
  NC::DataSources::addRecognisedFileExtensions("laz");
  NC::DataSources::addRecognisedFileExtensions("lau");
}
