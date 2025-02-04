////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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

#include "NCrystal/dump/NCDump.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"
#include <sstream>

void NCrystal::dump( std::ostream& os, const Info&c, DumpVerbosity verbosityVal )
{

  const auto verbose = enumAsInt( verbosityVal );
  nc_assert( verbose <= 2 );
  const bool isSinglePhase = c.isSinglePhase();

  //Display-labels - falling back to the AtomData description for multiphase:
  auto safeDisplayLabelImpl = [](const Info& info, const IndexedAtomData& ai)
  {
    if (ai.index.isInvalid())
      return ai.data().description(false);
    return info.displayLabel(ai.index);
  };
  auto safeDisplayLabel = [&safeDisplayLabelImpl,&c](const IndexedAtomData& ai)
  {
    return safeDisplayLabelImpl(c,ai);
  };

  //Figure out max display label width for column alignment:
  unsigned longestDisplayLabel(0);
  for ( auto& ce: c.getComposition() ) {
    longestDisplayLabel = std::max<unsigned>(longestDisplayLabel,
                                             static_cast<unsigned>(safeDisplayLabel(ce.atom).size()));
  }
  nc_assert_always(longestDisplayLabel<1000);

  const char * hr = "----------------------------------------------------------------------------------------------------\n";
  os << hr <<"------------------------------------   NCrystal Material Info   ------------------------------------\n";

  if ( !c.getDataSourceName().str().empty() && c.getDataSourceName().str() != "<unknown>" ) {
    os<<hr;
    os<<"Data source: "<<c.getDataSourceName()<<'\n';
  }

  if (!Cfg::CfgManip::empty(c.getCfgData())) {
    os<<hr;
    os<<"Process cfg hints: ";
    Cfg::CfgManip::stream( c.getCfgData(), os );
    os<<'\n';
  }

  //Densities and composition-related information:
  os<<hr;
  const auto averageMass = c.getAverageAtomMass();
  os<<"Density : "<<fmtg(c.getDensity().dbl())<<" g/cm3, "<<fmtg(c.getNumberDensity().dbl())<<" atoms/Aa^3\n";

  {//Composition
    os<<hr;
    os<<"Composition (by mole):";
    for (auto& e : c.getComposition()) {
      auto lbl = safeDisplayLabel(e.atom);
      os<<' '<<fmtg(e.fraction*100.0)<<"% "<<lbl;
    }
    os<<'\n'<<hr;
    os<<"Composition (by mass):";
    for (auto& e : c.getComposition()) {
      auto lbl = safeDisplayLabel(e.atom);
      os<<' '<<fmtg(e.fraction*100.0 * e.atom.data().averageMassAMU().dbl()/averageMass.dbl())
        <<"% "<<lbl;
    }
    os<<'\n'<<hr;
    os<<"Atom data:\n";
    for (auto& e : c.getComposition()) {
      auto lbl = safeDisplayLabel(e.atom);
      nc_assert(longestDisplayLabel>0);
      os<<"   "<<str_ljust(lbl,longestDisplayLabel)<<" = "<<*e.atom.atomDataSP<<'\n';
    }
  }
  os<<hr;
  os<<"Averaged quantities:\n";
  os<<"   Atomic mass               : "<<fmtg(averageMass.dbl())<<"u\n";//
  os<<"   Absorption XS at 2200m/s  : "<<fmtg(c.getXSectAbsorption().dbl())<<" barn\n";
  os<<"   Free scattering XS        : "<<fmtg(c.getXSectFree().dbl())<<" barn\n";
  auto sld = c.getSLD();
  os<<"   Scattering length density : "<<fmtg(c.getSLD().dbl())<<' '<<decltype(sld)::unit()<<'\n';

  if (c.hasTemperature()) {
    os<<hr;
    os<<"Temperature : "<<fmtg(c.getTemperature().dbl())<<" kelvin\n";
  }

  if ( c.stateOfMatter() != Info::StateOfMatter::Unknown ) {
    os<<hr;
    const char * subtype="";
    if ( c.stateOfMatter() == Info::StateOfMatter::Solid )
      subtype = isSinglePhase ? ( c.isCrystalline() ? " (crystalline)" : " (amorphous)" ) : "";
    os<<"State of matter: "<<Info::toString(c.stateOfMatter())<<subtype<<'\n';
  }

  if ( !isSinglePhase ) {
    ///////////////////
    //  Multi phase  //
    ///////////////////

    auto& phases = c.getPhases();
    os<<hr<<"Multi-phase material with "<<phases.size()<<" phases:\n";
    for ( auto& phidx : ncrange( phases.size() ) ) {
      const auto& ph = phases.at(phidx);
      const auto& volfrac = ph.first;
      const auto& phinfo = *ph.second;
      const double mole_frac = volfrac*phinfo.getNumberDensity().dbl()/c.getNumberDensity().dbl();
      os<<"  "<<phidx<<") "<< fmt(100.*volfrac,"%7g") <<"% (volume) "
        <<fmt(100.0*mole_frac,"%7g")<<"% (mole) "
        <<fmt(100.*mole_frac*phinfo.getAverageAtomMass().dbl()/averageMass.dbl(),"%7g")<<"% (mass): ";
      bool firstph(true);
      for ( auto& phc : phinfo.getComposition() ) {
        if (firstph)
          firstph=false;
        else
          os<<" / ";
        os<<safeDisplayLabelImpl(phinfo,phc.atom);
      }
      if ( phinfo.isMultiPhase() ) {
        os<<"  {"<<phinfo.getPhases().size()<<" sub-phases}";
      }

      os<<'\n';
    }
  } else {
    ////////////////////
    //  Single phase  //
    ////////////////////

    if ( c.hasStructureInfo() ) {
      const StructureInfo& si = c.getStructureInfo();
      os<<hr;
      if (si.spacegroup!=0)
        os<<"Space group number      : "<<si.spacegroup<<'\n';
      os<<"Lattice spacings   [Aa] : "<<fmtg(si.lattice_a)<<' '<<fmtg(si.lattice_b)<<' '<<fmtg(si.lattice_c)<<'\n';
      os<<"Lattice angles    [deg] : "<<fmtg(si.alpha)<<' '<<fmtg(si.beta)<<' '<<fmtg(si.gamma)<<'\n';
      os<<"Unit cell volume [Aa^3] : "<<fmtg(si.volume)<<'\n';
      os<<"Atoms / unit cell       : "<<si.n_atoms<<'\n';
    };

    if ( c.hasAtomInfo() ) {
      os<<hr;
      unsigned ntot = 0;
      for ( auto& ai : c.getAtomInfos() )
        ntot += ai.numberPerUnitCell();
      os<<"Atoms in unit cell (total "<<ntot<<"):\n";
      for ( auto& ai : c.getAtomInfos() ) {
        auto lbl = safeDisplayLabel( ai.atom() );
        nc_assert(longestDisplayLabel>0);
        os << "     "<< ai.numberPerUnitCell() <<' '
           << str_ljust(lbl,longestDisplayLabel)
           <<" atom"<<(ai.numberPerUnitCell()==1?"":"s");
        if ( ai.debyeTemp().has_value() || ai.msd().has_value() ) {
          os <<" [";
          if ( ai.debyeTemp().has_value() ) {
            os <<"T_Debye="<<ai.debyeTemp().value();
            if ( ai.msd().has_value() )
              os << ", ";
          }
          if ( ai.msd().has_value() )
            os <<"MSD="<<ai.msd().value()<<"Aa^2";
          os << "]";
        }
        os<<'\n';
      }
      {
        os<<hr;
        os<<"Atomic coordinates:\n";
        auto& atomlist = c.getAtomInfos();
        if ( !verbose && ntot>30) {
          os<<"  (suppressed due to their large number, increase verbosity to show)\n";
        } else {
          auto fmtg_frac_str = [](double value)->std::string
          {
            std::ostringstream ss;
            ss << fmtg_frac(value);
            return ss.str();
          };
          for ( auto& ai : atomlist ) {
            auto lbl = safeDisplayLabel(ai.atom());
            for ( auto& pos : ai.unitCellPositions() ) {
              nc_assert(longestDisplayLabel>0);
              os<<"     "<<str_ljust(lbl,longestDisplayLabel)<<"   "
                <<str_rjust(fmtg_frac_str(pos[0]),10)<<"   "
                <<str_rjust(fmtg_frac_str(pos[1]),10)<<"   "
                <<str_rjust(fmtg_frac_str(pos[2]),10)<<'\n';
            }
          }
        }
      }
    }

    if ( c.hasDynamicInfo()) {
      os<<hr;
      for (auto& di: c.getDynamicInfoList()) {
        auto lbl = safeDisplayLabel(di->atom());
        os<<"Dynamic info for "<<lbl<<" ("<<fmtg(di->fraction()*100.0)<<"%):\n";
        auto di_knl = dynamic_cast<const DI_ScatKnl*>(di.get());
        if (di_knl) {
          auto di_skd = dynamic_cast<const DI_ScatKnlDirect*>(di_knl);
          auto di_vdos = dynamic_cast<const DI_VDOS*>(di_knl);
          auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di_knl);
          os<<"   type: S(alpha,beta)"
            <<(di_vdos?" [from VDOS]":(di_vdosdebye?" [from VDOSDebye]":""))
            <<'\n';
          auto sp_egrid = di_knl->energyGrid();
          bool eg_is_auto = !sp_egrid || (sp_egrid->size()==3 && sp_egrid->at(0)==0.0 && sp_egrid->at(1)==0.0 && sp_egrid->at(2)==0.0);
          if (!eg_is_auto) {
            const auto& eg = *sp_egrid;
            //"Grids must have at least 3 entries, and grids of size 3 actually
            //indicates [emin,emax,npts], where any value can be 0 to leave the
            //choice for the consuming code. Grids of size >=4 must be proper
            //grids."
            if ( eg.size()==3 ) {
              double eg_emin = eg.at(0);
              double eg_emax = eg.at(1);
              double eg_npts = eg.at(2);
              os<<"   egrid: ";
              if ( eg_emin == 0.0 )
                os<<"auto";
              else
                os<<fmtg(eg_emin)<<" eV";
              os<<" -> ";
              if ( eg_emax == 0.0 )
                os<<"auto";
              else
                os<<fmtg(eg_emax)<<" eV";
              if ( eg_npts > 0)
                os<<" ("<<fmtg(eg_npts)<<" points)";
              os<<'\n';
            } else {
              nc_assert(eg.size()>=4);
              os<<"   egrid: "<<fmtg(eg.front())<<" eV -> "<<fmtg(eg.back())<<" eV ("<<eg.size()<<" points)\n";
            }
          }

          if (di_skd) {
            const auto& sabdata = *(di_skd->ensureBuildThenReturnSAB());
            const auto& ag = sabdata.alphaGrid();
            const auto& bg = sabdata.betaGrid();
            const auto& sab = sabdata.sab();
            os<<"   alpha-grid   : "<<fmtg(ag.front())<<" -> "<<fmtg(ag.back())<<" ("<<ag.size()<<" points)\n";
            os<<"   beta-grid    : "<<fmtg(bg.front())<<" -> "<<fmtg(bg.back())<<" ("<<bg.size()<<" points)\n";
            os<<"   S(alpha,beta): "<<sab.size()<<" points, S_max = "<<fmtg(*std::max_element(sab.begin(),sab.end()))<<'\n';
          }
          if (di_vdos) {
            os<<"   VDOS Source: "<<di_vdos->vdosData().vdos_density().size()<<" points\n";
            os<<"   VDOS E_max: "<<fmtg(di_vdos->vdosData().vdos_egrid().second*1000.0)<<" meV\n";
          } else if (di_vdosdebye) {
            os<<"   VDOS E_max: "<<fmtg(di_vdosdebye->debyeTemperature().kT()*1000.0)<<" meV\n";
          }
        } else if (dynamic_cast<const DI_Sterile*>(di.get())) {
          os<<"   type: sterile\n";
        } else if (dynamic_cast<const DI_FreeGas*>(di.get())) {
          os<<"   type: freegas\n";
        } else {
          nc_assert_always(false);
        }
      }
    }

    if ( c.providesNonBraggXSects() ) {
      os<<hr;
      os<<"Provides non-bragg cross section calculations. A few values are:\n";
      os<<"   lambda[Aa]  sigma_scat[barn]\n";
      double ll[] = {0.5, 1.0, 1.798, 2.5, 5, 10, 20 };
      for (unsigned i = 0; i < sizeof(ll)/sizeof(ll[0]); ++i)
        os<<fmt(ll[i],"%13g")<<' '
          <<fmt(c.xsectScatNonBragg(NeutronWavelength{ll[i]}).dbl(),"%17g")<<'\n';
    }

    {
      const auto& customsections = c.getAllCustomSections();
      if ( !customsections.empty() ) {
        os<<hr;
        os<<"Custom data sections:\n";
        for (const auto& e: customsections) {
          os<<"  "<<e.first<<":\n";
          for (const auto& line : e.second) {
            os<<"    ";
            for (std::size_t i = 0; i < line.size(); ++i ) {
              os<<line.at(i);
              if ( i+1 != line.size() )
                os<<' ';
            }
            os<<'\n';
          }
        }
      }
    }

    if ( c.hasHKLInfo() ) {
      //Figure out how many planes to show, and try to get them available fast.
      enum class HideResult { Show, Hide, HideRest };
      std::function<HideResult(std::size_t)> hide_entry = [](std::size_t i) { return i < 10 ? HideResult::Show : HideResult::HideRest; };
      std::size_t nhklplanes_max = 10;
      const HKLList * thelist = nullptr;
      Optional<HKLList> localHKLList;//for lifetime extension
      bool some_hidden = false;
      if ( verbose || c.hklListIsFullyInitialised() ) {
        //Need full list or full list is already available anyway:
        thelist = &c.hklList();
        if (verbose) {
          nhklplanes_max = thelist->size();
          if ( verbose == 1 && nhklplanes_max > 50 ) {
            //hide all but initial 30 and last 10
            hide_entry = [nhklplanes_max](std::size_t i) { return ( i>30 && i+10 < nhklplanes_max ) ? HideResult::Hide : HideResult::Show; };
          } else {
            hide_entry = [](std::size_t) { return HideResult::Show; };
          }
        }
      } else {
        //Full list is not available and we only need a few planes. For speedup,
        //try to create smaller partial list.
        double d1 = c.hklDLower();
        double d2 = c.hklDUpper();
        if ( std::isinf(d2) ) {
          //upper scale not available, try to estimate it:
          if ( !c.hasStructureInfo() ) {
            d2 = d1 + 1.0;//pure guess
          } else {
            auto& si = c.getStructureInfo();
            d2 = ncmax( 1.2 * d1, 2.0 * ncmax( si.lattice_a, si.lattice_b, si.lattice_c ) );
          }
        }
        nc_assert(!std::isinf(d1) && !std::isinf(d2) && d2>d1 && d1 >= 0.0);
        double dmin_quick = d2 < 15 ? 0.5 : 0.8;
        if ( dmin_quick <= 1.1 * c.hklDLower() )
          dmin_quick = 2.0*c.hklDLower();
        if ( dmin_quick >= 0.9 * d2 )
          dmin_quick = 0.5*d2;
        if ( dmin_quick > c.hklDLower()*1.1 && dmin_quick < d2*0.9 )
          localHKLList = c.hklListPartialCalc( dmin_quick, c.hklDUpper() );
        if ( localHKLList.has_value() && localHKLList.value().size() > nhklplanes_max ) {//> not >= so we know if some were hidden
          thelist = &localHKLList.value();
          some_hidden = true;
        } else {
          //Didn't work for whatever reason, fall back to asking for the full list:
          localHKLList.reset();
          thelist = &c.hklList();
        }
      }

      //NB: Be careful what methods we call here, so we don't trigger full list build by mistake!
      const HKLList& hklList = *thelist;
      os<<hr;
      os << "HKL info type: " << c.hklInfoType() << '\n';
      ExpandHKLHelper hklExpander( c );
      const bool doExpandHKL = verbose && hklExpander.canExpand( c.hklInfoType() );
      const bool doExpandHKLSnipSome = ( verbose == 1 );
      static bool allow_utf8 = ncgetenv_bool("UTF8");

      enum class Encoding { ASCII, UTF8 };
      Encoding encoding = ( allow_utf8 ? Encoding::UTF8 : Encoding::ASCII );

      auto prettyPrintHKL = []( HKL v, Encoding enc )
      {
        auto encDigit = [enc] ( int d ) {
          constexpr auto unicode_combined_overline = "\xcc\x85";
          if ( enc == Encoding::ASCII ) {
            return std::to_string( d );
          } else {
            auto sa = std::to_string( std::abs(d) );
            if ( d>=0 )
              return sa;
            std::string r;
            for ( auto& ch : sa ) {
              r += ch;
              r += unicode_combined_overline;
            }
            return r;
          }
        };
        auto isSingleDigit = [](int d) { return d >= -9 && d <= 9; };
        std::string res;
        if ( enc == Encoding::UTF8 && isSingleDigit(v.h) && isSingleDigit(v.k) && isSingleDigit(v.l) ) {
          res += encDigit(v.h);
          res += encDigit(v.k);
          res += encDigit(v.l);
        } else {
          res += '(';
          res += encDigit(v.h);
          res += ',';
          res += encDigit(v.k);
          res += ',';
          res += encDigit(v.l);
          res += ')';
        }
        return res;
      };
      auto prHKL = [&prettyPrintHKL,encoding](HKL v) { return prettyPrintHKL(v,encoding); };

      os<<hr;
      os<<"HKL planes (d_lower = "<<fmtg(c.hklDLower())
        <<" Aa, d_upper = "<<fmtg(c.hklDUpper())<<" Aa):\n";
      os<<"  H   K   L  d_hkl[Aa] Mult. FSquared[barn]"
        <<(doExpandHKL?" Expanded-HKL-list (+sign flips)":"")<<'\n';
      std::size_t idx(0);
      for ( auto& hkl : hklList ) {
        auto hideres = hide_entry( idx++ );
        if ( hideres != HideResult::Show ) {
          if ( (verbose==1) && hideres == HideResult::Hide && !some_hidden && hklList.size() != nhklplanes_max )
            os<<"...\n";
          some_hidden = true;
          if ( hideres == HideResult::HideRest )
            break;
          continue;
        }
        if (nhklplanes_max==0)
          break;
        --nhklplanes_max;
        os << str_rjust(std::to_string(hkl.hkl.h),3) << ' '
           << str_rjust(std::to_string(hkl.hkl.k),3) << ' '
           << str_rjust(std::to_string(hkl.hkl.l),3) << ' '
           << fmt(hkl.dspacing,"%10g") << ' '
           << str_rjust(std::to_string(hkl.multiplicity),4) << ' '
           << fmt(hkl.fsquared,"%11g")
           << (doExpandHKL?"":"\n");

        if ( doExpandHKL ) {

          auto v = hklExpander.expand(hkl);

          unsigned nEqvMaxIfSnip;
          if ( encoding == Encoding::UTF8 ) {
            nEqvMaxIfSnip = 12;
            for ( auto& e : v )
              if ( std::max(std::max(std::abs(e.h),std::abs(e.k)),std::abs(e.l)) > 9 )
                nEqvMaxIfSnip = 6;
          } else {
            nEqvMaxIfSnip = 8;
          }
          nc_assert(nEqvMaxIfSnip%2 == 0);

          auto nEqv = v.size();
          os<<"     ";
          for ( auto i : ncrange(nEqv) ) {
            if ( doExpandHKLSnipSome && nEqv > nEqvMaxIfSnip && ( i >= nEqvMaxIfSnip/2 && i+nEqvMaxIfSnip/2 < nEqv ) ) {
              continue;
            }
            auto& e = v[i];
            nc_assert(e.h||e.k||e.l);
            os<<prHKL(e);
            if ( i+1 != nEqv )
              os<<' ';
          }
          if ( doExpandHKLSnipSome && nEqv > nEqvMaxIfSnip )
            os<<" ...";
          os<<'\n';
        }
      }
      if (some_hidden)
        os<<"  (some planes left out for brevity, increase verbosity to show all)\n";
    }
  }

  os<<hr;
}

std::string NCrystal::dump_str( const Info&c, DumpVerbosity verbosityVal )
{
  std::ostringstream os;
  dump( os, c, verbosityVal );
  return std::move(os).str();
}

void NCrystal::dump( const Info&c, DumpVerbosity verbosityVal )
{
  Msg::outputMsg( dump_str(c,verbosityVal).c_str(), MsgType::RawOutput );
}
