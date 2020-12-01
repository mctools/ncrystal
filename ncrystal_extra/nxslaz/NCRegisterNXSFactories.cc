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

#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCFactoryRegistry.hh"
#include "NCrystal/NCMatCfg.hh"
#include "NCFactory_NXS.hh"
#include "NCLazLoader.hh"

namespace NC = NCrystal;

namespace NCrystal {

  class NXSFactory final : public FactoryBase {
  public:
    const char * getName() const final { return "stdnxs"; }

    int canCreateInfo( const MatCfg& cfg ) const final {
      return cfg.getDataFileExtension()=="nxs" ? 100 : 0;
    }
    RCHolder<const Info> createInfo( const MatCfg& cfg ) const final
    {
      nc_assert_always(canCreateInfo(cfg));
      const char * flag_bkgdlikemcstas = "mcstaslikebkgd";
      const char * flag_fixpolyatom = "fixpolyatoms";
      cfg.infofactopt_validate({flag_bkgdlikemcstas,flag_fixpolyatom});
      return loadNXSCrystal( cfg.getDataFile().c_str(),
                             cfg.get_temp()==-1.0?293.15:cfg.get_temp(),
                             cfg.get_dcutoff(),
                             cfg.get_dcutoffup(),
                             cfg.get_infofactopt_flag(flag_bkgdlikemcstas),
                             cfg.get_infofactopt_flag(flag_fixpolyatom)
                             );
    }
  };

  class LazFactory final : public FactoryBase {
  public:
    const char * getName() const final { return "stdlaz"; }
    int canCreateInfo( const MatCfg& cfg ) const final {
      std::string ext = cfg.getDataFileExtension();
      return (ext=="laz"||ext=="lau") ? 100 : 0;
    }
    RCHolder<const Info> createInfo( const MatCfg& cfg ) const final
    {
      LazLoader ld (cfg.getDataFile().c_str(),
                    cfg.get_dcutoff(),
                    cfg.get_dcutoffup(),
                    cfg.get_temp()==-1.0?293.15:cfg.get_temp());
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
  if (!NC::hasFactory("stdnxs"))
    NC::registerFactory( std::make_unique<NC::NXSFactory>() );
  if (!NC::hasFactory("stdlaz"))
    NC::registerFactory( std::make_unique<NC::LazFactory>() );
}
