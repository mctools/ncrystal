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

#include "NCrystal/NCMatInfo.hh"
#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/NCDataSources.hh"
#include "NCrystal/NCMatCfg.hh"
#include "NCFactory_NXS.hh"
#include "NCLazLoader.hh"

namespace NC = NCrystal;

namespace NCrystal {

  class NXSFactory final : public FactImpl::InfoFactory {
  public:
    const char * name() const noexcept final { return "stdnxs"; }

    Priority query( const MatInfoCfg& cfg ) const final
    {
      return cfg.getDataType()=="nxs" ? Priority{100} : Priority::Unable;
    }

    shared_obj<const MatInfo> produce( const MatInfoCfg& cfg ) const final
    {
      nc_assert_always(cfg.getDataType()=="nxs");
      const char * flag_bkgdlikemcstas = "mcstaslikebkgd";
      const char * flag_fixpolyatom = "fixpolyatoms";
      cfg.infofactopt_validate({flag_bkgdlikemcstas,flag_fixpolyatom});
      return makeSO<const MatInfo>(loadNXSCrystal( cfg.textData(),
                                                   cfg.get_temp().dbl()==-1.0?Temperature{293.15}:cfg.get_temp(),
                                                   cfg.get_dcutoff(),
                                                   cfg.get_dcutoffup(),
                                                   cfg.get_infofactopt_flag(flag_bkgdlikemcstas),
                                                   cfg.get_infofactopt_flag(flag_fixpolyatom)
                                                   ));
    }
  };

  class LazFactory final : public FactImpl::InfoFactory {
  public:
    const char * name() const noexcept final { return "stdlaz"; }
    Priority query( const MatInfoCfg& cfg ) const final
    {
      std::string ext = cfg.getDataType();
      return (ext=="laz"||ext=="lau") ? Priority{100} : Priority::Unable;
    }

    shared_obj<const MatInfo> produce( const MatInfoCfg& cfg ) const final
    {
      LazLoader ld (cfg.textData(),
                    cfg.get_dcutoff(),
                    cfg.get_dcutoffup(),
                    cfg.get_temp().dbl()==-1.0?Temperature{293.15}:cfg.get_temp());
      ld.read();
      return ld.getCrystalInfo();
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
