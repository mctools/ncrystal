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
#include "NCrystal/NCMatCfg.hh"

#include "NCPlaneProvider.hh"

#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCSCOrientation.hh"
#include "NCrystal/NCScatterComp.hh"

#include "NCrystal/NCPCBragg.hh"
#include "NCrystal/NCSCBragg.hh"
#include "NCrystal/NCLCBragg.hh"
#include "NCrystal/NCBkgdPhonDebye.hh"
#include "NCrystal/NCBkgdExtCurve.hh"

namespace NCrystal {

  class PlaneProviderWCutOff : public PlaneProvider {
  public:

    //Wraps plane provider and provides only those planes above or below a given
    //threshold. Will assume ownership of wrapped plane provider; If
    //select_above is true, only planes with dspacing >= cutoff will be
    //returned, if above is false, only planes with dspacing<cutoff will be returned;
    PlaneProviderWCutOff(double dcut, PlaneProvider* pp)
      : PlaneProvider(), m_pp(pp), m_dcut(dcut) { nc_assert(pp); pp->prepareLoop(); }
    virtual ~PlaneProviderWCutOff() { delete m_pp; }

    virtual bool getNextPlane(double& dspacing, double& fsq, Vector& demi_normal) {
      while (m_pp->getNextPlane(dspacing,fsq,demi_normal)) {
        if ( dspacing>=m_dcut ) {
          return true;
        } else {
          fsq*=2;//getNextPlane provides demi-normals, e.g. only half of the normals.

          if (m_withheldPlanes.empty()||m_withheldPlanes.back().first!=dspacing) {
#if __cplusplus >= 201103L
            m_withheldPlanes.emplace_back(dspacing,fsq);
#else
            m_withheldPlanes.push_back(std::make_pair(dspacing,fsq));
#endif
          } else {
            m_withheldPlanes.back().second += fsq;
          }
        }
      }
      return false;
    }

    virtual void prepareLoop() { m_pp->prepareLoop(); m_withheldPlanes.clear(); }
    virtual bool canProvide() const { return m_pp->canProvide(); }

    std::vector<std::pair<double,double> >& planesWithheldInLastLoop() { return m_withheldPlanes; };

  private:
    PlaneProvider* m_pp;
    double m_dcut;
    std::vector<std::pair<double,double> > m_withheldPlanes;
  };

  class NCStdScatFact : public FactoryBase {
  public:
    const char * getName() const { return "stdscat"; }

    virtual int canCreateScatter( const MatCfg& cfg ) const {
      RCHolder<const Info> info;
      std::string bkgd;
      return analyseCfg( cfg, info, bkgd ) ? 100 : 0;
    }

    virtual const Scatter * createScatter( const MatCfg& cfg ) const
    {
      //Analyse and extract info:
      RCHolder<const Info> info;
      std::string bkgd;
      if ( !analyseCfg( cfg, info, bkgd ) )
        return 0;
      nc_assert(info);

      //Collect components on ScatterComp:
      RCHolder<ScatterComp> sc(new ScatterComp);

      //Bragg component:
      if (cfg.get_bragg()) {
        if (cfg.isSingleCrystal()) {
          //TODO for NC2: factory function somewhere for this, so can be easily created directly in test-code wo matcfg?
          UniquePtr<PlaneProvider> sc_pp(createStdPlaneProvider(info.obj()));
          PlaneProviderWCutOff* ppwcutoff(0);
          nc_assert(info.obj()->hasHKLInfo());
          if ( cfg.get_sccutoff() && cfg.get_sccutoff() > info.obj()->hklDMinVal() ) {
            //Improve efficieny by treating planes with dspacing less than
            //sccutoff as having isotropic mosaicity distribution.
            sc_pp = ( ppwcutoff = new PlaneProviderWCutOff(cfg.get_sccutoff(),sc_pp.release()) );
          }
          SCOrientation sco = cfg.createSCOrientation();
          if (cfg.isLayeredCrystal()) {
            double lcdir[3];
            cfg.get_lcaxis(lcdir);
            sc.obj()->addComponent(new LCBragg(info.obj(), sco, cfg.get_mos(), lcdir, cfg.get_lcmode(),
                                               0,sc_pp.obj(),cfg.get_mosprec(),0.0));
          } else {
            sc.obj()->addComponent(new SCBragg(info.obj(),sco,cfg.get_mos(),0.0,sc_pp.obj(),cfg.get_mosprec(),0.));
          }
          if ( ppwcutoff && !ppwcutoff->planesWithheldInLastLoop().empty() ) {
            nc_assert_always(info.obj()->hasStructureInfo());
            sc.obj()->addComponent(new PCBragg(info.obj()->getStructureInfo(),ppwcutoff->planesWithheldInLastLoop()));
          }
        } else {
          sc.obj()->addComponent(new PCBragg(info.obj()));
          //NB: Layered polycrystals get same treatment as unlayered
          //polycrystals in our current modelling.
        }
      }

      //Allowed bkgd options depends on main bkgd mode:
      std::set<std::string> allowed_bkgdopts;

      //Common options for cross-section-curve-only bkgd modes:
      bool bkgdcurve_dothermalise = false;
      if ( bkgd == "external" || bkgd=="phonondebye" ) {
        allowed_bkgdopts.insert("elastic");
        allowed_bkgdopts.insert("thermalise");
        bool opt_elastic = cfg.get_bkgdopt_flag("elastic");
        bool opt_therm = cfg.get_bkgdopt_flag("thermalise");
        if (opt_elastic&&opt_therm)
          NCRYSTAL_THROW2(BadInput,"Can not specify both elastic and thermalise flags to bkgd="<<bkgd);
        //Use flags and temp availability to determine mode:
        bkgdcurve_dothermalise = info.obj()->hasTemperature();//default depends on temp availability
        if (opt_elastic) bkgdcurve_dothermalise = false;//user overrides to elastic
        else if (opt_therm) bkgdcurve_dothermalise = true;//user overrides to thermalising
      }
      //Common options for modes needing phonzeroinco flags:
      bool opt_no_pzi(false), opt_only_pzi(false);
      if ( bkgd == "phonondebye"
           ) {
        opt_no_pzi = cfg.get_bkgdopt_flag("no_phonzeroinco");
        opt_only_pzi = cfg.get_bkgdopt_flag("only_phonzeroinco");
        allowed_bkgdopts.insert("no_phonzeroinco");
        allowed_bkgdopts.insert("only_phonzeroinco");
        if (opt_no_pzi&&opt_only_pzi)
          NCRYSTAL_THROW2(BadInput,"Can not specify both no_phonzeroinco and only_phonzeroinco flags to bkgd="<<bkgd);
      }


      if ( bkgd == "external" ) {
        sc.obj()->addComponent(new BkgdExtCurve(info.obj(),bkgdcurve_dothermalise));
      } else if ( bkgd == "phonondebye" ) {
        allowed_bkgdopts.insert("nphonon");
        allowed_bkgdopts.insert("no_extrap");
        sc.obj()->addComponent(new BkgdPhonDebye( info.obj(),
                                                  bkgdcurve_dothermalise,
                                                  cfg.get_bkgdopt_int("nphonon",0),
                                                  !opt_no_pzi, opt_only_pzi,
                                                  cfg.get_bkgdopt_flag("no_extrap") ) );
      } else {
        nc_assert_always(bkgd=="none");
      }

      //Complain if user specified options not supported by bkgd mode in question:
      cfg.bkgdopt_validate(allowed_bkgdopts);

      //Wrap it up and return:
      if (sc.obj()->nComponents()==0) {
        //No components available, represent this with a NullScatter instead:
        return new NullScatter;
      } else if ( sc.obj()->nComponents()==1 && sc.obj()->scale(0) == 1.0 ) {
        //Single component with unit scale - no need to wrap in ScatterComp:
        RCHolder<Scatter> comp(sc.obj()->component(0));
        sc.clear();
        return comp.releaseNoDelete();
      } else {
        //Usual case, return ScatterComp:
        return sc.releaseNoDelete();
      }
    }

  private:
    //Common analysis function shared between canCreateScatter and createScatter
    //methods. If succesful, info object will have already been reffed.
    bool analyseCfg( const MatCfg& cfg, RCHolder<const Info>& info, std::string& bkgd ) const {
      info = ::NCrystal::createInfo( cfg );
      if ( !info )
        return false;
      bkgd.clear();

      //Bragg
      if (cfg.get_bragg()) {
        if ( !info.obj()->hasHKLInfo() || !info.obj()->hasStructureInfo() ) {
          //in principle PCBragg could generate scatterings (but not xsects) without
          //structure info, but we keeps things simple.
          return false;
        }
      }


      //Bkgd
      bkgd = cfg.get_bkgd_name();

      if (bkgd=="none")
        return true;

      if (bkgd=="best") {
        //Automatically select bkgd model. This will choose the most realistic
        //model available:
        if (false) {
        } else if (BkgdPhonDebye::hasSufficientInfo(info.obj())) {
          bkgd = "phonondebye";
        } else if (info.obj()->providesNonBraggXSects()) {
          bkgd="external";
        } else {
          bkgd="none";
        }
        return true;
      }

      //User specified bkgd explicitly:
      if (bkgd=="phonondebye"||bkgd=="external") {
        bool can_do = ( bkgd=="phonondebye"
                        ? BkgdPhonDebye::hasSufficientInfo(info.obj())
                        : info.obj()->providesNonBraggXSects() );
        if (!can_do)
          return false;
        if (!info.obj()->hasTemperature()&&cfg.get_bkgdopt_flag("thermalise"))
          return false;
        return true;
      }
      return false;//bkgd mode which is not support by this factory
    }

  };
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void ncrystal_register_stdscat_factory()
{
  if (!NCrystal::hasFactory("stdscat"))
    NCrystal::registerFactory(new NCrystal::NCStdScatFact);
}
