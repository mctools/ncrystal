////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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
#include "NCrystal/NCMatCfg.hh"
#include "NCrystal/NCInfo.hh"

#include "NCrystal/NCSCOrientation.hh"
#include "NCrystal/NCPCBragg.hh"
#include "NCrystal/NCSCBragg.hh"
#include "NCrystal/NCSimpleBkgd.hh"
#include "NCrystal/NCScatterComp.hh"

namespace NCrystal {

  class NCStdScatFact : public FactoryBase {
  public:
    const char * getName() const { return "NCrystalStdScatFactory"; }

    virtual int canCreateScatter( const MatCfg& cfg ) const {
      int priority;
      RCHolder<const Info> info;
      std::string bkgd;
      return analyseCfg( cfg, priority, info, bkgd ) ? priority : 0;
    }

    virtual const Scatter * createScatter( const MatCfg& cfg ) const
    {
      int priority;
      RCHolder<const Info> info;
      std::string bkgd;
      if ( !analyseCfg( cfg, priority, info, bkgd ) )
        return 0;
      nc_assert(info);
      Scatter *scatter_bragg(0), *scatter_bkgd(0);
      if (!cfg.get_skipbragg()) {
        if (cfg.isSingleCrystal()) {
          SCOrientation sco = cfg.createSCOrientation();
          scatter_bragg = new SCBragg(info.obj(),sco,cfg.get_mosaicity());
        } else {
          scatter_bragg = new PCBragg(info.obj());
        }
      }
      if (!bkgd.empty()) {
        nc_assert_always(bkgd=="simplethermalising"||bkgd=="simpleelastic");
        scatter_bkgd = new SimpleBkgd(info.obj(),(bkgd=="simplethermalising"));
      }
      if (scatter_bkgd && scatter_bragg) {
        ScatterComp * sc = new ScatterComp("ScatterComp");
        sc->addComponent(scatter_bragg);
        sc->addComponent(scatter_bkgd);
        return sc;
      }
      if (scatter_bragg) return scatter_bragg;
      if (scatter_bkgd) return scatter_bkgd;
      return 0;
    }

  private:
    //Common analysis function shared between canCreateScatter and createScatter
    //methods. If succesful, info object will have already been reffed.
    bool analyseCfg( const MatCfg& cfg, int& priority, RCHolder<const Info>& info, std::string& bkgd ) const {
      info = ::NCrystal::createInfo( cfg );
      if ( !info )
        return false;
      bkgd.clear();
      bool can_do_genscatter = true;
      bool can_do_xsects = true;
      if (!cfg.get_skipbragg()) {
        if ( !info.obj()->hasHKLInfo() )
          return false;//both PCBragg and SCBragg needs this.
        if (!info.obj()->hasStructureInfo()) {
          if (cfg.isSingleCrystal())
            return false;//SCBragg always needs this
          can_do_xsects = false;//PCBragg can generate scatterings, but not xsects, without it.
        }
      }
      if (!cfg.get_braggonly()) {
        //Background as well. For the initial NCrystal release, that means NCSimpleBkgd.
        bkgd = cfg.get_scatterbkgdmodel();
        if (bkgd=="best") {
          //Pick the best available background model. For first NCrystal
          //releasenow it really just means the simple background.
          //todo for nc2: check for phonon dos / S(q,w) availability.
          bkgd = "simpleelastic";
          if (info.obj()->hasTemperature()) {
            bkgd = "simplethermalising";
          }
        }
        if (bkgd=="simpleelastic"||bkgd=="simplethermalising") {
          if (!info.obj()->providesNonBraggXSects())
            can_do_xsects = false;
          if (bkgd=="simplethermalising"&&!info.obj()->hasTemperature())
            can_do_genscatter = false;
        } else {
          bkgd.clear();
          return false;//unknown or unsupported choice
        }
      }
      priority = (can_do_xsects?50:0) + (can_do_genscatter?50:0);
      if (!priority)
        return 0;
      return priority;
    }

  };
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void ncrystal_register_stdscat_factory()
{
  if (!NCrystal::hasFactory("NCrystalStdScatFactory"))
    NCrystal::registerFactory(new NCrystal::NCStdScatFact);
}
