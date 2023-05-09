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

#include "NCrystal/NCDump.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/NCMatCfg.hh"
#include "NCrystal/internal/NCCfgManip.hh"
#include "NCrystal/internal/NCPlaneProvider.hh"

#include <cstdio>
#include <sstream>
#include <iostream>
#include <iomanip>

void NCrystal::dump( const Info&c, DumpVerbosity verbosityVal )
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
  std::cout << hr <<"------------------------------------   NCrystal Material Info   ------------------------------------\n";

  if ( !c.getDataSourceName().str().empty() && c.getDataSourceName().str() != "<unknown>" ) {
    printf("%s", hr);
    std::cout<<"Data source: "<<c.getDataSourceName()<<std::endl;
  }

  if (!Cfg::CfgManip::empty(c.getCfgData())) {
    printf("%s", hr);
    std::cout<<"Process cfg hints: ";
    Cfg::CfgManip::stream( c.getCfgData(), std::cout );
    std::cout<<std::endl;
  }

  //Densities and composition-related information:
  printf("%s", hr);
  const auto averageMass = c.getAverageAtomMass();
  std::cout<<"Density : "<<fmtg(c.getDensity().dbl())<<" g/cm3, "<<fmtg(c.getNumberDensity().dbl())<<" atoms/Aa^3"<<std::endl;

  {//Composition
    printf("%s", hr);
    printf("Composition (by mole):");
    for (auto& e : c.getComposition()) {
      auto lbl = safeDisplayLabel(e.atom);
      printf(" %g%% %s",e.fraction*100.0,lbl.c_str());
    }
    printf("\n%s", hr);
    printf("Composition (by mass):");
    for (auto& e : c.getComposition()) {
      auto lbl = safeDisplayLabel(e.atom);
      printf(" %g%% %s",e.fraction*100.0 * e.atom.data().averageMassAMU().dbl()/averageMass.dbl(),lbl.c_str());
    }
    printf("\n%s", hr);
    printf("Atom data:\n");
    for (auto& e : c.getComposition()) {
      auto lbl = safeDisplayLabel(e.atom);
      nc_assert(longestDisplayLabel>0);
      std::cout<<"   "<<std::left << std::setw(longestDisplayLabel) << lbl<<" = "<<*e.atom.atomDataSP<<std::endl;
    }
  }
  printf("%s", hr);
  printf("Averaged quantities:\n");
  printf("   Atomic mass               : %gu\n",averageMass.dbl());
  printf("   Absorption XS at 2200m/s  : %g barn\n", c.getXSectAbsorption().dbl());
  printf("   Free scattering XS        : %g barn\n", c.getXSectFree().dbl());
  auto sld = c.getSLD();
  printf("   Scattering length density : %g %s\n", c.getSLD().dbl(), decltype(sld)::unit());

  if (c.hasTemperature()) {
    printf("%s", hr);
    printf("Temperature : %g kelvin\n",c.getTemperature().dbl());
  }

  if ( c.stateOfMatter() != Info::StateOfMatter::Unknown ) {
    printf("%s", hr);
    const char * subtype="";
    if ( c.stateOfMatter() == Info::StateOfMatter::Solid )
      subtype = isSinglePhase ? ( c.isCrystalline() ? " (crystalline)" : " (amorphous)" ) : "";
    printf("State of matter: %s%s\n",Info::toString(c.stateOfMatter()).c_str(),subtype);
  }

  if ( !isSinglePhase ) {
    ///////////////////
    //  Multi phase  //
    ///////////////////

    auto& phases = c.getPhases();
    std::cout<<hr<<"Multi-phase material with "<<phases.size()<<" phases:\n";
    for ( auto& phidx : ncrange( phases.size() ) ) {
      const auto& ph = phases.at(phidx);
      const auto& volfrac = ph.first;
      const auto& phinfo = *ph.second;
      const double mole_frac = volfrac*phinfo.getNumberDensity().dbl()/c.getNumberDensity().dbl();
      std::cout<<"  "<<phidx<<") "<< fmt(100.*volfrac,"%7g") <<"% (volume) "
               <<fmt(100.0*mole_frac,"%7g")<<"% (mole) "
               <<fmt(100.*mole_frac*phinfo.getAverageAtomMass().dbl()/averageMass.dbl(),"%7g")<<"% (mass): ";
      bool firstph(true);
      for ( auto& phc : phinfo.getComposition() ) {
        if (firstph)
          firstph=false;
        else
          std::cout<<" / ";
        std::cout<<safeDisplayLabelImpl(phinfo,phc.atom);
      }
      if ( phinfo.isMultiPhase() ) {
        std::cout<<"  {"<<phinfo.getPhases().size()<<" sub-phases}";
      }

      std::cout<<"\n";
    }
  } else {
    ////////////////////
    //  Single phase  //
    ////////////////////

    if ( c.hasStructureInfo() ) {
      const StructureInfo& si = c.getStructureInfo();
      printf("%s", hr);
      if (si.spacegroup!=0)
        printf("Space group number      : %i\n",si.spacegroup);
      printf("Lattice spacings   [Aa] : %g %g %g\n",si.lattice_a,si.lattice_b,si.lattice_c);
      printf("Lattice angles    [deg] : %g %g %g\n",si.alpha,si.beta,si.gamma);
      printf("Unit cell volume [Aa^3] : %g\n",si.volume);
      printf("Atoms / unit cell       : %i\n",si.n_atoms);
    };

    if ( c.hasAtomInfo() ) {
      printf("%s", hr);
      unsigned ntot = 0;
      for ( auto& ai : c.getAtomInfos() )
        ntot += ai.numberPerUnitCell();
      printf("Atoms in unit cell (total %i):\n",ntot);
      for ( auto& ai : c.getAtomInfos() ) {
        auto lbl = safeDisplayLabel( ai.atom() );
        std::ostringstream s;
        nc_assert(longestDisplayLabel>0);
        s << "     "<< ai.numberPerUnitCell() <<" "
          << std::left << std::setw(longestDisplayLabel)
          <<lbl
          <<" atom"<<(ai.numberPerUnitCell()==1?"":"s");
        if ( ai.debyeTemp().has_value() || ai.msd().has_value() ) {
          s <<" [";
          if ( ai.debyeTemp().has_value() ) {
            s <<"T_Debye="<<ai.debyeTemp().value();
            if ( ai.msd().has_value() )
              s << ", ";
          }
          if ( ai.msd().has_value() )
            s <<"MSD="<<ai.msd().value()<<"Aa^2";
          s<<"]";
        }
        printf("%s\n",s.str().c_str());
      }
      {
        printf("%s", hr);
        printf("Atomic coordinates:\n");
        auto& atomlist = c.getAtomInfos();
        if ( !verbose && ntot>30) {
          printf("  (suppressed due to their large number, increase verbosity to show)\n");
        } else {
          auto prettyPrintValue2Str = [](double value)->std::string
          {
            std::ostringstream ss;
            ss << fmtg_frac(value);
            return ss.str();
          };
          for ( auto& ai : atomlist ) {
            auto lbl = safeDisplayLabel(ai.atom());
            for ( auto& pos : ai.unitCellPositions() ) {
              std::ostringstream ss;
              nc_assert(longestDisplayLabel>0);
              ss << std::left << std::setw(longestDisplayLabel) << lbl;
              printf("     %s   %10s   %10s   %10s\n",
                     ss.str().c_str(),
                     prettyPrintValue2Str(pos[0]).c_str(),
                     prettyPrintValue2Str(pos[1]).c_str(),
                     prettyPrintValue2Str(pos[2]).c_str());
            }
          }
        }
      }
    }

    if ( c.hasDynamicInfo()) {
      printf("%s", hr);
      for (auto& di: c.getDynamicInfoList()) {
        auto lbl = safeDisplayLabel(di->atom());
        printf("Dynamic info for %s (%g%%):\n",lbl.c_str(),di->fraction()*100.0);
        auto di_knl = dynamic_cast<const DI_ScatKnl*>(di.get());
        if (di_knl) {
          auto di_skd = dynamic_cast<const DI_ScatKnlDirect*>(di_knl);
          auto di_vdos = dynamic_cast<const DI_VDOS*>(di_knl);
          auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di_knl);
          printf("   type: S(alpha,beta)%s\n",(di_vdos?" [from VDOS]":(di_vdosdebye?" [from VDOSDebye]":"")));
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
              printf("   egrid: ");
              if ( eg_emin == 0.0 )
                printf("auto");
              else
                printf("%g eV",eg_emin);
              printf(" -> ");
              if ( eg_emax == 0.0 )
                printf("auto");
              else
                printf("%g eV",eg_emax);
              if ( eg_npts > 0)
                printf(" (%g points)",eg_npts);
              printf("\n");
            } else {
              nc_assert(eg.size()>=4);
              printf("   egrid: %g eV -> %g eV (%llu points)\n",eg.front(),eg.back(), (unsigned long long)eg.size());
            }
          }

          if (di_skd) {
            const auto& sabdata = *(di_skd->ensureBuildThenReturnSAB());
            const auto& ag = sabdata.alphaGrid();
            const auto& bg = sabdata.betaGrid();
            const auto& sab = sabdata.sab();
            printf("   alpha-grid   : %g -> %g (%llu points)\n",ag.front(),ag.back(), (unsigned long long)ag.size());
            printf("   beta-grid    : %g -> %g (%llu points)\n",bg.front(),bg.back(), (unsigned long long)bg.size());
            printf("   S(alpha,beta): %llu points, S_max = %g\n",(unsigned long long)sab.size(), *std::max_element(sab.begin(),sab.end()));
          }
          if (di_vdos) {
            printf("   VDOS Source: %llu points\n",(unsigned long long)di_vdos->vdosData().vdos_density().size());
            printf("   VDOS E_max: %g meV\n",di_vdos->vdosData().vdos_egrid().second*1000.0);
          } else if (di_vdosdebye) {
            printf("   VDOS E_max: %g meV\n",di_vdosdebye->debyeTemperature().kT()*1000.0);
          }
        } else if (dynamic_cast<const DI_Sterile*>(di.get())) {
          printf("   type: sterile\n");
        } else if (dynamic_cast<const DI_FreeGas*>(di.get())) {
          printf("   type: freegas\n");
        } else {
          nc_assert_always(false);
        }
      }
    }

    if ( c.providesNonBraggXSects() ) {
      printf("%s", hr);
      printf("Provides non-bragg cross section calculations. A few values are:\n");
      printf("   lambda[Aa]  sigma_scat[barn]\n");
      double ll[] = {0.5, 1.0, 1.798, 2.5, 5, 10, 20 };
      for (unsigned i = 0; i < sizeof(ll)/sizeof(ll[0]); ++i)
        printf("%13g %17g\n",ll[i],c.xsectScatNonBragg(NeutronWavelength{ll[i]}).dbl());
    }

    {
      const auto& customsections = c.getAllCustomSections();
      if ( !customsections.empty() ) {
        printf("%s", hr);
        printf("Custom data sections:\n");
        for (const auto& e: customsections) {
          printf("  %s:\n", e.first.c_str());
          for (const auto& line : e.second) {
            printf("    ");
            for (std::size_t i = 0; i < line.size(); ++i ) {
              printf("%s",line.at(i).c_str());
              if ( i+1 != line.size() )
                printf(" ");
            }
            printf("\n");
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
      printf("%s", hr);
      std::cout << "HKL info type: " << c.hklInfoType() << std::endl;
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

      printf("%s", hr);
      printf("HKL planes (d_lower = %g Aa, d_upper = %g Aa):\n",c.hklDLower(),c.hklDUpper());
      printf("  H   K   L  d_hkl[Aa] Mult. FSquared[barn]%s\n",
             (doExpandHKL?" Expanded-HKL-list (+sign flips)":""));
      std::size_t idx(0);
      for ( auto& hkl : hklList ) {
        auto hideres = hide_entry( idx++ );
        if ( hideres != HideResult::Show ) {
          if ( (verbose==1) && hideres == HideResult::Hide && !some_hidden && hklList.size() != nhklplanes_max )
            printf("...\n");
          some_hidden = true;
          if ( hideres == HideResult::HideRest )
            break;
          continue;
        }
        if (nhklplanes_max==0)
          break;
        --nhklplanes_max;
        printf("%3i %3i %3i %10g %4i %11g%s",hkl.hkl.h,hkl.hkl.k,hkl.hkl.l,hkl.dspacing,
               hkl.multiplicity,hkl.fsquared,(doExpandHKL?"":"\n"));

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
          printf("     ");
          for ( auto i : ncrange(nEqv) ) {
            if ( doExpandHKLSnipSome && nEqv > nEqvMaxIfSnip && ( i >= nEqvMaxIfSnip/2 && i+nEqvMaxIfSnip/2 < nEqv ) ) {
              continue;
            }
            auto& e = v[i];
            nc_assert(e.h||e.k||e.l);
            printf( "%s%s",
                    prHKL(e).c_str(),
                    (i+1==nEqv?"":" ") );
          }
          if ( doExpandHKLSnipSome && nEqv > nEqvMaxIfSnip )
            printf(" ...");
          printf("\n");
        }
      }
      if (some_hidden)
        printf("  (some planes left out for brevity, increase verbosity to show all)\n");
    }
  }

  std::cout<<hr<<std::flush;
}
