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

#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCFactoryRegistry.hh"
#include "NCrystal/NCMatCfg.hh"
#include "NCFactory_NXS.hh"
#include "NCLazLoader.hh"

namespace NCrystal {

  class NXSFactory : public FactoryBase {
  public:
    const char * getName() const { return "stdnxs"; }

    virtual int canCreateInfo( const MatCfg& cfg ) const {
      return cfg.getDataFileExtension()=="nxs" ? 100 : 0;
    }
    virtual const Info * createInfo( const MatCfg& cfg ) const
    {
      nc_assert_always(canCreateInfo(cfg));
      const char * flag_bkgdlikemcstas = "mcstaslikebkgd";
      const char * flag_fixpolyatom = "fixpolyatoms";
#if __cplusplus >= 201103L
      cfg.infofactopt_validate({flag_bkgdlikemcstas,flag_fixpolyatom});
#else
      std::set<std::string> allowed_infofactopts;
      allowed_infofactopts.insert(flag_bkgdlikemcstas);
      allowed_infofactopts.insert(flag_fixpolyatom);
      cfg.infofactopt_validate(allowed_infofactopts);
#endif
      return loadNXSCrystal( cfg.getDataFile().c_str(),
                             cfg.get_temp(),
                             cfg.get_dcutoff(),
                             cfg.get_dcutoffup(),
                             cfg.get_infofactopt_flag(flag_bkgdlikemcstas),
                             cfg.get_infofactopt_flag(flag_fixpolyatom)
                             );
    }
  };

  class LazFactory : public FactoryBase {
  public:
    const char * getName() const { return "stdlaz"; }
    virtual int canCreateInfo( const MatCfg& cfg ) const {
      std::string ext = cfg.getDataFileExtension();
      return (ext=="laz"||ext=="lau") ? 100 : 0;
    }
    virtual const Info * createInfo( const MatCfg& cfg ) const
    {
      const Info *ci(0);
      {
        LazLoader ld (cfg.getDataFile().c_str(),
                      cfg.get_dcutoff(),
                      cfg.get_dcutoffup(),
                      cfg.get_temp());
        ld.read();
        ci = ld.getCrystalInfo();
        ci->ref();//Avoid deletion when ld unrefs.
      }
      ci->unrefNoDelete();
      return ci;
    }
  };

}

//Finally, a function which can be used to enable the above factories. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void ncrystal_register_nxslaz_factories()
{
  if (!NCrystal::hasFactory("stdnxs"))
    NCrystal::registerFactory(new NCrystal::NXSFactory);
  if (!NCrystal::hasFactory("stdlaz"))
    NCrystal::registerFactory(new NCrystal::LazFactory);
}
