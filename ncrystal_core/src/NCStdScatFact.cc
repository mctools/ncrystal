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
#include "NCrystal/NCMatCfg.hh"

#include "NCPlaneProvider.hh"
#include "NCDynInfoUtils.hh"
#include "NCNeutronSCL.hh"

#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCSCOrientation.hh"
#include "NCrystal/NCScatterComp.hh"

#include "NCrystal/NCPCBragg.hh"
#include "NCrystal/NCSCBragg.hh"
#include "NCrystal/NCLCBragg.hh"
#include "NCrystal/NCBkgdExtCurve.hh"
#include "NCrystal/NCFreeGas.hh"
#include "NCrystal/NCElIncScatter.hh"
#include "NCrystal/NCSABScatter.hh"
#include "NCSABFactory.hh"

///////////////////////////////////////////////////////////////////////////////////
// The standard Scatter factory. Will accept with priority 100 any request for   //
// which the inelas parameter has a known value ("auto","external","dyninfo",    //
// etc.). Thus, custom non-standard Scatter factories will have to use a         //
// non-standard value of the inelas parameter in order to get precedence (for    //
// when the standard Scatter factory should still be accessible for some         //
// configurations), or merely return a priority below 100 (for when the standard //
// Scatter factory should be completely replaced).                               //
///////////////////////////////////////////////////////////////////////////////////

namespace NCrystal {

  class PlaneProviderWCutOff : public PlaneProvider {
  public:

    //Wraps plane provider and provides only those planes above or below a given
    //threshold. Will assume ownership of wrapped plane provider; If
    //select_above is true, only planes with dspacing >= cutoff will be
    //returned, if above is false, only planes with dspacing<cutoff will be returned;
    PlaneProviderWCutOff(double dcut, std::unique_ptr<PlaneProvider> pp)
      : PlaneProvider(), m_pp(std::move(pp)), m_dcut(dcut) { nc_assert(m_pp); m_pp->prepareLoop(); }
    virtual ~PlaneProviderWCutOff() {}

    virtual bool getNextPlane(double& dspacing, double& fsq, Vector& demi_normal) {
      while (m_pp->getNextPlane(dspacing,fsq,demi_normal)) {
        if ( dspacing>=m_dcut ) {
          return true;
        } else {
          fsq*=2;//getNextPlane provides demi-normals, e.g. only half of the normals.
          if (m_withheldPlanes.empty()||m_withheldPlanes.back().first!=dspacing)
            m_withheldPlanes.emplace_back(dspacing,fsq);
          else
            m_withheldPlanes.back().second += fsq;
        }
      }
      return false;
    }

    virtual void prepareLoop() { m_pp->prepareLoop(); m_withheldPlanes.clear(); }
    virtual bool canProvide() const { return m_pp->canProvide(); }

    std::vector<PairDD >& planesWithheldInLastLoop() { return m_withheldPlanes; };

  private:
    std::unique_ptr<PlaneProvider> m_pp;
    double m_dcut;
    std::vector<PairDD > m_withheldPlanes;
  };

  class NCStdScatFact : public FactoryBase {
  public:
    const char * getName() const { return "stdscat"; }

    virtual int canCreateScatter( const MatCfg& cfg ) const {

      RCHolder<const Info> info;
      std::string inelas;
      return analyseCfg( cfg, info, inelas ) ? 100 : 0;
    }

    virtual const Scatter * createScatter( const MatCfg& cfg ) const
    {
      //Analyse and extract info:
      RCHolder<const Info> info_holder;
      std::string inelas;
      if ( !analyseCfg( cfg, info_holder, inelas ) )
        return 0;
      nc_assert_always(info_holder.obj());
      nc_assert_always(isOneOf(inelas,"none","external","dyninfo","vdosdebye","freegas"));
      const Info& info = *info_holder.obj();

      //Collect components on ScatterComp:
      RCHolder<ScatterComp> sc_holder(new ScatterComp);
      ScatterComp& sc = *sc_holder.obj();

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Incoherent-elastic component:
      if ( cfg.get_incoh_elas() && info.isCrystalline() ) {
        const bool has_msd  = info.hasAtomMSD() || ( info.hasTemperature() && info.hasAnyDebyeTemperature() );
        if ( has_msd )
          sc.addComponent( new ElIncScatter( &info ) );
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Coherent-elastic (Bragg) component:
      if ( cfg.get_coh_elas() && info.isCrystalline() && info.hasHKLInfo() ) {
        if (cfg.isSingleCrystal()) {
          //TODO for NC2: factory function somewhere for this, so can be easily created directly in test-code wo matcfg?
          UniquePtr<PlaneProvider> sc_pp(createStdPlaneProvider(&info));
          PlaneProviderWCutOff* ppwcutoff(0);
          nc_assert(info.hasHKLInfo());
          if ( cfg.get_sccutoff() && cfg.get_sccutoff() > info.hklDMinVal() ) {
            //Improve efficiency by treating planes with dspacing less than
            //sccutoff as having isotropic mosaicity distribution.
            sc_pp = ( ppwcutoff = new PlaneProviderWCutOff(cfg.get_sccutoff(),std::unique_ptr<PlaneProvider>(sc_pp.release())) );
          }
          SCOrientation sco = cfg.createSCOrientation();
          if (cfg.isLayeredCrystal()) {
            double lcdir[3];
            cfg.get_lcaxis(lcdir);
            sc.addComponent(new LCBragg(&info, sco, cfg.get_mos(), lcdir, cfg.get_lcmode(),
                                        0,sc_pp.obj(),cfg.get_mosprec(),0.0));
          } else {
            sc.addComponent(new SCBragg(&info,sco,cfg.get_mos(),0.0,sc_pp.obj(),cfg.get_mosprec(),0.));
          }
          if ( ppwcutoff && !ppwcutoff->planesWithheldInLastLoop().empty() ) {
            nc_assert_always(info.hasStructureInfo());
            sc.addComponent(new PCBragg(info.getStructureInfo(),ppwcutoff->planesWithheldInLastLoop()));
          }
        } else {
          sc.addComponent(new PCBragg(&info));
          //NB: Layered polycrystals get same treatment as unlayered
          //polycrystals in our current modelling.
        }
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Inelastic components:
      if ( inelas == "none" ) {

        //do not add anything.

      } else if ( inelas == "external" ) {

        if ( !info.providesNonBraggXSects() )
          NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires input source which provides direct"
                          " parameterisation of (non-Bragg) scattering cross sections (try e.g. inelas=auto instead)");

        sc.addComponent( new BkgdExtCurve( &info ) );

      } else {

        nc_assert_always( isOneOf(inelas,"dyninfo","freegas", "vdosdebye" ) );

        if ( !info.hasComposition() )
          NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires specification of material composition");

        if ( !info.hasTemperature() )
          NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires specification of material temperature");

        if ( inelas == "dyninfo" ) {

          if ( !info.hasDynamicInfo() )
            NCRYSTAL_THROW(BadInput,"inelas=dyninfo does not work for input without specific dynamic information. It is possible that"
                           " other modes might work (try e.g. inelas=auto instead).");

          for (auto& di : info.getDynamicInfoList()) {
            const DI_ScatKnl* di_scatknl = dynamic_cast<const DI_ScatKnl*>(di.get());
            if (di_scatknl) {
              sc.addComponent( new SABScatter( *di_scatknl, cfg.get_vdoslux() ), di->fraction() );
            } else if (dynamic_cast<const DI_Sterile*>(di.get())) {
              continue;//just skip past sterile components
            } else if (dynamic_cast<const DI_FreeGas*>(di.get())) {
              sc.addComponent( new FreeGas( info.getTemperature(), di->elementName() ), di->fraction() );
            } else {
              NCRYSTAL_THROW(LogicError,"Unsupported DynamicInfo entry encountered.");
            }
          }

        } else if ( inelas=="freegas" ) {

          for ( auto& e : info.getComposition() )
            sc.addComponent( new FreeGas( info.getTemperature(), e.first), e.second );

        } else {

          nc_assert_always(inelas=="vdosdebye");

          if ( !info.hasAnyDebyeTemperature() )
            NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires specification of Debye temperature");
          if ( !info.isCrystalline() || !info.hasAtomInfo() )
            NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires crystalline material with atomic information");

          unsigned ntot = 0.0;
          for (auto it = info.atomInfoBegin(); it!= info.atomInfoEnd(); ++it)
            ntot += it->number_per_unit_cell;
          for (auto it = info.atomInfoBegin(); it!= info.atomInfoEnd(); ++it) {
            const double debyeTemp = it->debye_temp > 0.0 ? it->debye_temp : info.getGlobalDebyeTemperature();
            const double elementMassAMU = NeutronSCL::instance()->getAtomicMass(it->element_name);
            const auto boundXS = SigmaBound{NeutronSCL::instance()->getBoundXS(it->element_name)};
            auto sabdata =  extractSABDataFromVDOSDebyeModel( debyeTemp,
                                                              info.getTemperature(),
                                                              boundXS,
                                                              elementMassAMU,
                                                              cfg.get_vdoslux() );
            auto scathelper = SAB::createScatterHelperWithCache( std::move(sabdata) );
            sc.addComponent( new SABScatter( std::move(scathelper) ), it->number_per_unit_cell*1.0/ntot );
          }
        }
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Wrap it up and return:
      if (sc.nComponents()==0) {
        //No components available, represent this with a NullScatter instead:
        return new NullScatter;
      } else if ( sc.nComponents()==1 && sc.scale(0) == 1.0 ) {
        //Single component with unit scale - no need to wrap in ScatterComp:
        RCHolder<Scatter> comp(sc.component(0));
        sc_holder.clear();
        return comp.releaseNoDelete();
      } else {
        //Usual case, return ScatterComp:
        return sc_holder.releaseNoDelete();
      }
    }

  private:
    //Common analysis function shared between canCreateScatter and createScatter
    //methods.
    bool analyseCfg( const MatCfg& cfg, RCHolder<const Info>& info_holder, std::string& inelas ) const {

      info_holder = ::NCrystal::createInfo( cfg );//NB: ::NCrystal prefix is important to call correct function.
      if ( !info_holder )
        return false;
      const Info& info = *info_holder.obj();

      inelas = cfg.get_inelas();

      if ( inelas == "none" )
        return true;

      if ( isOneOf(inelas,"external","dyninfo","vdosdebye","freegas" ) ) {
        return true;
      }

      if (inelas!="auto")
        return false;//inelas mode not supported by this standard factory

      //Automatic selection requested.  Try to select in a way which respects
      //the behaviour outlined in ncmat_doc.md as well as respecting non-ncmat
      //sources of info (i.e. laz/lau files gives no inelas modelling and .nxs
      //files use external xs curves (and isotropic/deltaE=0 scattering even for
      //inelastic components).

      const bool has_temperature_and_composition = info.hasTemperature() && info.hasComposition();
      if ( info.providesNonBraggXSects() ) {
        inelas = "external";//.nxs files end up here
        return true;
      }
      if ( info.hasDynamicInfo() ) {
        nc_assert(has_temperature_and_composition);//guaranteed by NCInfo validation
        inelas = "dyninfo";//.ncmat files end up here
        return true;
      }
      if ( has_temperature_and_composition ) {
        inelas = info.hasAnyDebyeTemperature() ? "vdosdebye" : "freegas";
        return true;
      }

      inelas = "none";//.laz/.lau files end up here
      return true;
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
