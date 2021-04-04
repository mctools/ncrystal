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

#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/internal/NCPlaneProvider.hh"
#include "NCrystal/internal/NCDynInfoUtils.hh"
#include "NCrystal/NCProcImpl.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCSCOrientation.hh"
#include "NCrystal/internal/NCPCBragg.hh"
#include "NCrystal/internal/NCSCBragg.hh"
#include "NCrystal/internal/NCLCBragg.hh"
#include "NCrystal/internal/NCBkgdExtCurve.hh"
#include "NCrystal/internal/NCFreeGas.hh"
#include "NCrystal/internal/NCElIncScatter.hh"
#include "NCrystal/internal/NCSABScatter.hh"
#include "NCrystal/internal/NCSABFactory.hh"
#include "NCrystal/internal/NCString.hh"//for safe_str2dbl

namespace NC = NCrystal;

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

    bool hasPlanesWithheldInLastLoop() const { return !m_withheldPlanes.empty(); };

    PCBragg::VectDFM&& consumePlanesWithheldInLastLoop() { return std::move(m_withheldPlanes); };

  private:
    std::unique_ptr<PlaneProvider> m_pp;
    double m_dcut;
    PCBragg::VectDFM m_withheldPlanes;
  };

  class StdScatFact : public FactImpl::ScatterFactory {
  public:
    const char * name() const noexcept final { return "stdscat"; }

    Priority query( const MatCfg& cfg ) const final
    {
      return analyseCfg( cfg ).able ? Priority{100} : Priority{Priority::Unable};
    }

    ProcImpl::ProcPtr produce( const MatCfg& cfg ) const final
    {
      //Analyse and extract info:
      auto ana = analyseCfg( cfg );
      if ( !ana.able )
        NCRYSTAL_THROW2(LogicError,"Could not create ProcImpl process for cfg=\""<<cfg
                        <<"\" (this factory-produce method should not have been called)");

      const Info& info = ana.info;
      const auto& inelas = ana.inelas;

      nc_assert_always(isOneOf(inelas,"none","external","dyninfo","vdosdebye","freegas"));

      //Collect components:
      ProcImpl::ProcComposition::ComponentList components;

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Incoherent-elastic component:
      if ( cfg.get_incoh_elas() && info.isCrystalline() ) {
        const bool has_msd = info.hasAtomMSD() || ( info.hasTemperature() && info.hasDebyeTemperature() );
        if ( has_msd )
          components.push_back({1.0,makeSO<ElIncScatter>(ElIncScatter::msd_from_atominfo_t(),info)});
      }

      if ( info.countCustomSections( "SPECIALINCOHELAS" ) > 0 ) {
        //Special "secret" test code for testing purposes (before we perhaps
        //decide to make it an official feature). Add @CUSTOM_SPECIALINCOHELAS
        //section to enable incoherent-elastic treatment in non-crystals.
        if ( info.countCustomSections( "SPECIALINCOHELAS" ) != 1 )
          NCRYSTAL_THROW(BadInput,"Only one CUSTOM_SPECIALINCOHELAS section is allowed in input.");
        bool has_atominfo_msd = info.hasAtomMSD() || ( info.hasTemperature() && info.hasDebyeTemperature() );
        bool has_dyninfovdos_msd = false;
        if ( !has_atominfo_msd ) {
          for ( auto& di : info.getDynamicInfoList() ) {
            if ( dynamic_cast<const DI_VDOS*>(di.get()) ) {
              //Allow if even a single VDOS is present (it might for instance be
              //e.g. H-in-polyethylene where the almost irrelevant C atom is
              //simply modelled as free gas).
              has_dyninfovdos_msd = true;
              break;
            }
          }
        }
        if (!has_atominfo_msd && !has_dyninfovdos_msd)
          NCRYSTAL_THROW(BadInput,"Presence of CUSTOM_SPECIALINCOHELAS section requires ability to"
                         " determine atomic mean squared displacements for all atoms.");
        if ( info.isCrystalline() )
          NCRYSTAL_THROW(BadInput,"A @CUSTOM_SPECIALINCOHELAS section can currently only be used when there is otherwise no unit cell information in the input.");
        if ( cfg.get_incoh_elas() ) {
          //Decode:
          //
          //Syntax is single line with factor (>0.0,<=10.0), followed optionally by the keyword "include_sigma_coh"
          const auto& cd = info.getCustomSection("SPECIALINCOHELAS");//vector<VectS>
          if ( cd.size() != 1 )
            NCRYSTAL_THROW(BadInput,"CUSTOM_SPECIALINCOHELAS section should only have 1 line.");
          const auto& line = cd.at(0);
          nc_assert_always(!line.empty());
          if ( line.empty() || line.size() > 2  )
            NCRYSTAL_THROW(BadInput,"Bad syntax in CUSTOM_SPECIALINCOHELAS section (too many entries).");
          double scale_factor;
          if (!safe_str2dbl( line.at(0), scale_factor ) || !(scale_factor>0.0) || !(scale_factor<=10.0) )
            NCRYSTAL_THROW(BadInput,"Invalid scale factor in CUSTOM_SPECIALINCOHELAS section. Value should be in interval (0.0,10.0].");
          bool include_sigma_coh(false);
          if ( line.size()==2 )  {
            if ( line.at(1)!="include_sigma_coh" )
              NCRYSTAL_THROW2(BadInput,"Invalid keyword in CUSTOM_SPECIALINCOHELAS section. Syntax requires either"
                              " just the scale factor, or a scale factor followed by the keyword \"include_sigma_coh\".");
            include_sigma_coh = true;
          }
          if (has_atominfo_msd)
            components.push_back({1.0,makeSO<ElIncScatter>(ElIncScatter::msd_from_atominfo_t(),info,scale_factor, include_sigma_coh)});
          else
            components.push_back({1.0,makeSO<ElIncScatter>(ElIncScatter::msd_from_dyninfo_t(),info,scale_factor, include_sigma_coh)});
        }
      }


      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Coherent-elastic (Bragg) component:
      if ( cfg.get_coh_elas() && info.isCrystalline() && info.hasHKLInfo() ) {
        if (cfg.isSingleCrystal()) {
          //TODO: factory function somewhere for this, so can be easily created directly in test-code wo matcfg?
          auto sc_pp = createStdPlaneProvider( ana.info );
          PlaneProviderWCutOff* ppwcutoff(nullptr);
          nc_assert(info.hasHKLInfo());
          if ( cfg.get_sccutoff() && cfg.get_sccutoff() > info.hklDMinVal() ) {
            //Improve efficiency by treating planes with dspacing less than
            //sccutoff as having isotropic mosaicity distribution.
            auto tmp = std::make_unique<PlaneProviderWCutOff>(cfg.get_sccutoff(),std::move(sc_pp));
            ppwcutoff = tmp.get();
            sc_pp = std::move(tmp);
            nc_assert( sc_pp!=nullptr && (void*)sc_pp.get()==(void*)ppwcutoff );
          }
          SCOrientation sco = cfg.createSCOrientation();
          if (cfg.isLayeredCrystal()) {
            components.push_back({1.0,makeSO<LCBragg>( info, sco, cfg.get_mos(), cfg.get_lcaxis(), cfg.get_lcmode(),
                                                       0,sc_pp.get(),cfg.get_mosprec(),0.0 )});
          } else {
            components.push_back({1.0,makeSO<SCBragg>( info, sco,cfg.get_mos(),0.0,
                                                       sc_pp.get(),cfg.get_mosprec(),0.)});


          }
          if ( ppwcutoff && ppwcutoff->hasPlanesWithheldInLastLoop() ) {
            nc_assert_always(info.hasStructureInfo());
            components.push_back({1.0,makeSO<PCBragg>(info.getStructureInfo(),ppwcutoff->consumePlanesWithheldInLastLoop())});
          }
        } else {
          components.push_back({1.0,makeSO<PCBragg>(info)});
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

        components.push_back({1.0,makeSO<BkgdExtCurve>(ana.info)});

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
              components.push_back({di->fraction(),makeSO<SABScatter>(*di_scatknl, cfg.get_vdoslux())});
            } else if (dynamic_cast<const DI_Sterile*>(di.get())) {
              continue;//just skip past sterile components
            } else if (dynamic_cast<const DI_FreeGas*>(di.get())) {
              components.push_back({di->fraction(),makeSO<FreeGas>(info.getTemperature(), di->atomData())});
            } else {
              NCRYSTAL_THROW(LogicError,"Unsupported DynamicInfo entry encountered.");
            }
          }

        } else if ( inelas=="freegas" ) {

          for ( auto& e : info.getComposition() ) {
            components.push_back({e.fraction,makeSO<FreeGas>(info.getTemperature(), e.atom.data())});
          }

        } else {

          nc_assert_always(inelas=="vdosdebye");

          if ( !info.hasDebyeTemperature() )
            NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires specification of Debye temperature");
          if ( !info.isCrystalline() || !info.hasAtomInfo() )
            NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires crystalline material with atomic information");

          unsigned ntot = 0.0;
          for (auto it = info.atomInfoBegin(); it!= info.atomInfoEnd(); ++it)
            ntot += it->numberPerUnitCell();
          for (auto it = info.atomInfoBegin(); it!= info.atomInfoEnd(); ++it) {
            nc_assert_always( it->debyeTemp().has_value() );
            auto sabdata =  extractSABDataFromVDOSDebyeModel( it->debyeTemp().value(),
                                                              info.getTemperature(),
                                                              it->atomData().scatteringXS(),
                                                              it->atomData().averageMassAMU(),
                                                              cfg.get_vdoslux() );
            auto scathelper = SAB::createScatterHelperWithCache( std::move(sabdata) );
            components.push_back({it->numberPerUnitCell()*1.0/ntot,makeSO<SABScatter>(std::move(scathelper))});

          }
        }
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Wrap it up and return:
      return ProcImpl::ProcComposition::consumeAndCombine( std::move(components), ProcessType::Scatter );
    }

  private:
    //Common analysis function shared between canCreateScatter and createScatter
    //methods.
    struct CfgAnalysis {
      CfgAnalysis(shared_obj<const Info> ii) : info(std::move(ii)) {}
      bool able = true;
      shared_obj<const Info> info;
      std::string inelas;
    };
    CfgAnalysis analyseCfg( const MatCfg& cfg ) const {
      CfgAnalysis result(FactImpl::createInfo( cfg ));
      const Info& info = result.info;

      result.inelas = cfg.get_inelas();

      if ( result.inelas == "none" )
        return result;

      if ( isOneOf(result.inelas,"external","dyninfo","vdosdebye","freegas" ) ) {
        return result;
      }

      if (result.inelas!="auto") {
        //inelas mode not supported by this standard factory
        result.able = false;
        return result;
      }

      //Automatic selection requested.  Try to select in a way which respects
      //the behaviour outlined in ncmat_doc.md as well as respecting non-ncmat
      //sources of info (i.e. laz/lau files gives no inelas modelling and .nxs
      //files use external xs curves (and isotropic/deltaE=0 scattering even for
      //inelastic components).

      const bool has_temperature_and_composition = info.hasTemperature() && info.hasComposition();
      if ( info.providesNonBraggXSects() ) {
        result.inelas = "external";//.nxs files end up here
        return result;
      }
      if ( info.hasDynamicInfo() ) {
        nc_assert(has_temperature_and_composition);//guaranteed by NCInfo validation
        result.inelas = "dyninfo";//.ncmat files end up here
        return result;
      }
      if ( has_temperature_and_composition ) {
        result.inelas = info.hasDebyeTemperature() ? "vdosdebye" : "freegas";
        return result;
      }

      result.inelas = "none";//.laz/.lau files end up here
      return result;
    }

  };
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void ncrystal_register_stdscat_factory()
{
  if (!NC::FactImpl::hasScatterFactory("stdscat"))
    NC::FactImpl::registerFactory(std::make_unique<NC::StdScatFact>());
}
