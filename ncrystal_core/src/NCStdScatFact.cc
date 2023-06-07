////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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
#include "NCrystal/internal/NCSABUCN.hh"
#include "NCrystal/internal/NCString.hh"

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

    Optional<Plane> getNextPlane() override {
      //double& dspacing, double& fsq, Vector& demi_normal
      Optional<Plane> res;
      while ( ( res = m_pp->getNextPlane() ).has_value() ) {
        if ( res.value().dspacing>=m_dcut ) {
          return res;
        } else {
          const double fsq = res.value().fsq * 2;//getNextPlane provides demi-normals, e.g. only half of the normals.
          if (m_withheldPlanes.empty()||m_withheldPlanes.back().first!=res.value().dspacing)
            m_withheldPlanes.emplace_back(res.value().dspacing,fsq);
          else
            m_withheldPlanes.back().second += fsq;
        }
      }
      return NullOpt;
    }

    void prepareLoop() override { m_pp->prepareLoop(); m_withheldPlanes.clear(); }
    bool canProvide() const override { return m_pp->canProvide(); }
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

    Priority query( const FactImpl::ScatterRequest& cfg ) const final
    {
      return analyseCfg( cfg ).able ? Priority{100} : Priority{Priority::Unable};
    }

    ProcImpl::ProcPtr produce( const FactImpl::ScatterRequest& cfg ) const final
    {
      //Analyse and extract info:
      auto ana = analyseCfg( cfg );
      if ( !ana.able )
        NCRYSTAL_THROW2(LogicError,"Could not create ProcImpl process for cfg=\""<<cfg
                        <<"\" (this factory-produce method should not have been called)");

      const Info& info = cfg.info();
      const auto& inelas = ana.inelas;
      const auto ucnmode = cfg.get_ucnmode();
      const auto vdoslux = cfg.get_vdoslux();

      nc_assert_always(isOneOf(inelas,"0","external","dyninfo","vdosdebye","freegas"));

      //Unofficial hacks in @CUSTOM_UNOFFICIALHACKS section for special hacks
      //that are needed for various tests, preliminary support of new materials,
      //or special plugins.
      auto getUnofficialHack = [&info]( const std::string& keyword ) -> Optional<VectS>
      {
        auto n = info.countCustomSections( "UNOFFICIALHACKS" );
        if ( n == 0 )
          return NullOpt;
        if ( n > 1 )
          NCRYSTAL_THROW(BadInput,"Only one CUSTOM_UNOFFICIALHACKS section is allowed in input.");
        for ( auto& line : info.getCustomSection("UNOFFICIALHACKS") )
          if ( line.at(0) == keyword )
            return VectS(std::next(line.begin()),line.end());
        return NullOpt;
      };

      //Collect components:
      ProcImpl::ProcComposition::ComponentList components;

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Crystals: Incoherent-elastic component:
      if ( cfg.get_incoh_elas() && info.isCrystalline() ) {
        if ( ElIncScatter::hasSufficientInfo(info) )
          components.push_back({1.0,makeSO<ElIncScatter>(info)});
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Crystals: Coherent-elastic (Bragg) component:
      if ( cfg.get_coh_elas() && info.isCrystalline() ) {
        nc_assert(info.hasHKLInfo());
        if (cfg.isSingleCrystal()) {
          //TODO: factory function somewhere for this, so can be easily created directly in test-code?
          auto sc_pp = createStdPlaneProvider( cfg.infoPtr() );
          PlaneProviderWCutOff* ppwcutoff(nullptr);
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
      //Amorphous (non-crystal) solids: Incoherent and coherent elastic (via incoh. approximation) component
      if ( !info.isCrystalline() && info.stateOfMatter() == Info::StateOfMatter::Solid ) {
        bool add_inc = cfg.get_incoh_elas();
        bool add_coh = cfg.get_coh_elas();
        //Unofficial hack disabling the coherent part (for instance so plugins
        //can create their own coherent elastic based on S(Q)):
        if ( getUnofficialHack("no_cohelas_via_incohapprox_for_amorphous_solids").has_value() )
          add_coh = false;
        if ( add_inc || add_coh ) {
          ElIncScatterCfg elinc_cfg;
          elinc_cfg.use_sigma_incoherent = add_inc;
          elinc_cfg.use_sigma_coherent = add_coh;
          if ( ElIncScatter::hasSufficientInfo(info, elinc_cfg) )
            components.push_back({1.0,makeSO<ElIncScatter>(info,elinc_cfg)});
        }
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Inelastic components:
      if ( inelas == "0" ) {

        //do not add anything.

      } else if ( inelas == "external" ) {

        if ( ucnmode.has_value() )
          NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode is not compatible with any ucnmode)");

        if ( !info.providesNonBraggXSects() )
          NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires input source which provides direct"
                          " parameterisation of (non-Bragg) scattering cross sections (try e.g. inelas=auto instead)");

        components.push_back({1.0,makeSO<BkgdExtCurve>(cfg.infoPtr())});

      } else {
        nc_assert_always( isOneOf(inelas,"dyninfo","freegas", "vdosdebye" ) );

        if ( !info.hasTemperature() )
          NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires specification of material temperature");

        if ( inelas == "dyninfo" ) {

          if ( !info.hasDynamicInfo() )
            NCRYSTAL_THROW(BadInput,"inelas=dyninfo does not work for input without specific dynamic information. It is possible that"
                           " other modes might work (try e.g. inelas=auto instead).");

          uint32_t vdos2sabExcludeFlag = 0;
          auto specialIgnoreContribs = getUnofficialHack("vdos2sab_ignorecontrib");
          if ( specialIgnoreContribs.has_value() ) {
            //Parse syntax:  vdos2sab_ignorecontrib low [high] [coherent|incoherent]
            auto& l = specialIgnoreContribs.value();
            nc_assert_always(l.size()>=1);
            unsigned mode = 3;
            if ( isOneOf(l.back(),"coherent","incoherent") ) {
              mode = l.back()=="coherent" ? 1 : 2;
              l.pop_back();
            }
            nc_assert_always(l.size()>=1);
            unsigned low = str2int(l.at(0));
            unsigned high = l.size()==2 ? str2int(l.at(1)) : low;
            nc_assert_always(low>=1&&high>=low);
            nc_assert(high<=9999);
            vdos2sabExcludeFlag = mode + 4*low + 40000*high;
          }

          for (auto& di : info.getDynamicInfoList()) {
            const DI_ScatKnl* di_scatknl = dynamic_cast<const DI_ScatKnl*>(di.get());
            if (di_scatknl) {
              auto sabdata = NC::extractSABDataFromDynInfo( di_scatknl, vdoslux, true/*use cache*/, vdos2sabExcludeFlag );
              if ( sabdata->boundXS() ) {
                auto sab_scatter = makeSO<SABScatter>( sabdata, di_scatknl->energyGrid() );
                if ( ucnmode.has_value() ) {
                  auto scUCN = UCN::UCNScatter::createWithCache( sabdata, ucnmode.value().threshold );
                  nc_assert(isOneOf(ucnmode.value().mode,UCNMode::Mode::Refine,UCNMode::Mode::Remove,UCNMode::Mode::Only));
                  if ( scUCN->isNull() ) {
                    //Just the normal process, the UCN process is apparently null.
                    if ( isOneOf( ucnmode.value().mode, UCNMode::Mode::Refine,  UCNMode::Mode::Remove ) )
                      components.push_back({di->fraction(),sab_scatter});
                  } else {
                    if ( isOneOf( ucnmode.value().mode, UCNMode::Mode::Refine,  UCNMode::Mode::Remove ) )
                      components.push_back({ di->fraction(), makeSO<UCN::ExcludeUCNScatter>( sab_scatter, scUCN ) });
                    if ( isOneOf( ucnmode.value().mode, UCNMode::Mode::Refine,  UCNMode::Mode::Only ) )
                      components.push_back({ di->fraction(), scUCN });
                  }
                } else {
                  components.push_back({di->fraction(),sab_scatter});
                }
              }
            } else if (dynamic_cast<const DI_Sterile*>(di.get())) {
              continue;//just skip past sterile components
            } else if (dynamic_cast<const DI_FreeGas*>(di.get())) {
              components.push_back({di->fraction(),makeSO<FreeGas>(info.getTemperature(), di->atomData())});
              if ( ucnmode.has_value() )
                NCRYSTAL_THROW2(BadInput,"components with freegas dynamics are currently not supported when using any ucnmode)");
            } else {
              NCRYSTAL_THROW(LogicError,"Unsupported DynamicInfo entry encountered.");
            }
          }

        } else if ( inelas=="freegas" ) {
          if ( ucnmode.has_value() )
            NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode is not compatible with any ucnmode)");
          for ( auto& e : info.getComposition() ) {
            components.push_back({e.fraction,makeSO<FreeGas>(info.getTemperature(), e.atom.data())});
          }

        } else {
          nc_assert_always(inelas=="vdosdebye");
          if ( ucnmode.has_value() )
            NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode is not compatible with any ucnmode)");
          if ( !info.hasDebyeTemperature() )
            NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires specification of Debye temperature");//TODO: This should be allowed also for elements with actual VDOS
          if ( !info.isCrystalline() || !info.hasAtomInfo() )
            NCRYSTAL_THROW2(BadInput,"inelas="<<inelas<<" mode requires crystalline material with atomic information");

          unsigned ntot = 0.0;
          for ( auto& ai : info.getAtomInfos() )
            ntot += ai.numberPerUnitCell();
          for ( auto& ai : info.getAtomInfos() ) {
            nc_assert_always( ai.debyeTemp().has_value() );
            auto sabdata =  extractSABDataFromVDOSDebyeModel( ai.debyeTemp().value(),
                                                              info.getTemperature(),
                                                              ai.atomData().scatteringXS(),
                                                              ai.atomData().averageMassAMU(),
                                                              vdoslux );
            auto scathelper = SAB::createScatterHelperWithCache( std::move(sabdata) );
            components.push_back({ai.numberPerUnitCell()*1.0/ntot,makeSO<SABScatter>(std::move(scathelper))});

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
      bool able = true;
      std::string inelas;
    };
    CfgAnalysis analyseCfg( const FactImpl::ScatterRequest& cfg ) const {
      CfgAnalysis result;
      const Info& info = cfg.info();

      result.inelas = cfg.get_inelas();

      if ( result.inelas == "0" )
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
      //sources of info (i.e. .nxs files might use external xs curves).

      if ( info.providesNonBraggXSects() ) {
        result.inelas = "external";//.nxs files end up here
        return result;
      }
      if ( info.hasDynamicInfo() ) {
        nc_assert( info.hasTemperature() );//guaranteed by InfoBuilder validation
        result.inelas = "dyninfo";//.ncmat files end up here
        return result;
      }
      if ( info.hasTemperature() ) {
        //No direct information provided about material dynamics available. If
        //Debye temperatures are available and the material could be a solid
        //(e.g. state of matter is solid or unknown, not explicitly
        //gas/liquid/...), we can in principle revert to a free-gas
        //model. However, a free-gas model accounts for the entire cross
        //section, and should usually not be combined with e.g. Bragg
        //diffraction. So we only add the free-gas modelling if there is no HKL
        //information available.
        const bool could_be_solid = ( info.stateOfMatter() == Info::StateOfMatter::Solid
                                      || info.stateOfMatter() == Info::StateOfMatter::Unknown );
        if ( could_be_solid && info.hasDebyeTemperature() ) {
          result.inelas = "vdosdebye";
        } else {
          result.inelas = info.hasHKLInfo() ? "0" : "freegas";
        }
        return result;
      }

      result.inelas = "0";
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
