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

#include "NCrystal/NCDump.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCString.hh"
#include <cstdio>
#include <sstream>
#include <iostream>
#include <iomanip>

void NCrystal::dump(const Info&c)
{
  //Figure out max display label width for column alignment:
  unsigned longestDisplayLabel(0);
  if (c.hasComposition()) {
    nc_assert_always(c.getComposition().size()<std::numeric_limits<unsigned>::max());
    auto ncomps = static_cast<unsigned>(c.getComposition().size());
    for ( unsigned i = 0; i < ncomps; ++i )
      longestDisplayLabel = std::max<unsigned>(longestDisplayLabel,
                                               static_cast<unsigned>(c.displayLabel(AtomIndex{i}).size()));
    nc_assert_always(longestDisplayLabel<1000);
  }

  const char * hr = "---------------------------------------------------------\n";

  if ( c.stateOfMatter() != Info::StateOfMatter::Unknown ) {
    printf("%s", hr);
    const char * subtype="";
    if ( c.stateOfMatter() == Info::StateOfMatter::Solid )
      subtype = c.isCrystalline() ? " (crystalline)" : " (amorphous)";
    printf("State of matter: %s%s\n",Info::toString(c.stateOfMatter()).c_str(),subtype);
  }

  if (c.hasStructureInfo()) {
    const StructureInfo& si = c.getStructureInfo();
    printf("%s", hr);
    if (si.spacegroup!=0)
      printf("Space group number      : %i\n",si.spacegroup);
    printf("Lattice spacings   [Aa] : %g %g %g\n",si.lattice_a,si.lattice_b,si.lattice_c);
    printf("Lattice angles    [deg] : %g %g %g\n",si.alpha,si.beta,si.gamma);
    printf("Unit cell volume [Aa^3] : %g\n",si.volume);
    printf("Atoms / unit cell       : %i\n",si.n_atoms);
  };

  if (c.hasAtomInfo()) {
    printf("%s", hr);
    unsigned ntot = 0;
    for ( auto& ai : c.getAtomInfos() )
      ntot += ai.numberPerUnitCell();
    printf("Atoms in unit cell (total %i):\n",ntot);
    for ( auto& ai : c.getAtomInfos() ) {
      auto lbl = c.displayLabel( ai.atom().index );
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
      for ( auto& ai : c.getAtomInfos() ) {
        auto lbl = c.displayLabel(ai.atom().index);
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

  if (c.hasDensity()) {
    printf("%s", hr);
    printf("Density : %g g/cm3\n",c.getDensity().dbl());
  }

  if (c.hasNumberDensity()) {
    printf("%s", hr);
    printf("NumberDensity : %g atoms/Aa3\n",c.getNumberDensity().dbl());
  }

  if (c.hasComposition()) {
    printf("%s", hr);
    printf("Composition:\n");
    for (auto& e : c.getComposition()) {
      auto lbl = c.displayLabel(e.atom.index);
      printf(" %20g%% %s\n",e.fraction*100.0,lbl.c_str());
    }
    printf("%s", hr);
    printf("Atom data:\n");
    for (auto& e : c.getComposition()) {
      auto lbl = c.displayLabel(e.atom.index);
      nc_assert(longestDisplayLabel>0);
      std::cout<<"   "<<std::left << std::setw(longestDisplayLabel) << lbl<<" = "<<*e.atom.atomDataSP<<std::endl;
    }
  }

  if (c.hasTemperature()) {
    printf("%s", hr);
    printf("Temperature : %g kelvin\n",c.getTemperature().dbl());
  }

  if (c.hasDynamicInfo()) {
    printf("%s", hr);
    for (auto& di: c.getDynamicInfoList()) {
      auto lbl = c.displayLabel(di->atom().index);
      printf("Dynamic info for %s (%g%%):\n",lbl.c_str(),di->fraction()*100.0);
      auto di_knl = dynamic_cast<const DI_ScatKnl*>(di.get());
      if (di_knl) {
        auto di_skd = dynamic_cast<const DI_ScatKnlDirect*>(di_knl);
        auto di_vdos = dynamic_cast<const DI_VDOS*>(di_knl);
        auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di_knl);
        printf("   type: S(alpha,beta)%s\n",(di_vdos?" [from VDOS]":(di_vdosdebye?" [from VDOSDebye]":"")));
        auto sp_egrid = di_knl->energyGrid();
        if (!!sp_egrid)
          printf("   egrid: %g -> %g (%llu points)\n",sp_egrid->front(),sp_egrid->back(), (unsigned long long)sp_egrid->size());
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

  if (c.hasXSectAbsorption()||c.hasXSectFree()) {
    printf("%s", hr);
    printf("Neutron cross-sections:\n");
    if (c.hasXSectAbsorption())
      printf("   Absorption at 2200m/s : %g barn\n", c.getXSectAbsorption().dbl());
    if (c.hasXSectFree())
      printf("   Free scattering       : %g barn\n", c.getXSectFree().dbl());
  }

  if (c.providesNonBraggXSects()) {
    printf("%s", hr);
    printf("Provides non-bragg cross-section calculations. A few values are:\n");
    printf("   lambda[Aa]  sigma_scat[barn]\n");
    double ll[] = {0.5, 1.0, 1.798, 2.5, 5, 10, 20 };
    for (unsigned i = 0; i < sizeof(ll)/sizeof(ll[0]); ++i)
      printf("%13g %17g\n",ll[i],c.xsectScatNonBragg(NeutronWavelength{ll[i]}).dbl());
  }

  const auto& customsections = c.getAllCustomSections();
  if (!customsections.empty()) {
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

  if (c.hasHKLInfo()) {
    printf("%s", hr);
    printf("HKL planes (d_lower = %g Aa, d_upper = %g Aa):\n",c.hklDLower(),c.hklDUpper());
    printf("  H   K   L  d_hkl[Aa] Multiplicity FSquared[barn]%s\n",
           (c.hasExpandedHKLInfo()?" Expanded-HKL-list":""));
    HKLList::const_iterator itE = c.hklEnd();
    for (HKLList::const_iterator it = c.hklBegin();it!=itE;++it) {
      printf("%3i %3i %3i %10g %12i %14g%s",it->h,it->k,it->l,it->dspacing,
             it->multiplicity,it->fsquared,(it->eqv_hkl?"":"\n"));
      if (it->eqv_hkl!=nullptr) {
        const short * eqv_hkl = &(it->eqv_hkl[0]);
        const size_t nEqv = it->demi_normals.size();
        nc_assert_always( nEqv );
        printf(" ");
        for (size_t i = 0 ; i < nEqv; ++i) {
          const short h(eqv_hkl[i*3]), k(eqv_hkl[i*3+1]), l(eqv_hkl[i*3+2]);
          nc_assert(h||k||l);
          printf("%i,%i,%i | %i,%i,%i%s",h,k,l,-h,-k,-l,(i+1==nEqv?"":" | "));
        }
        printf("\n");
      }
    }
  }
  printf("%s", hr);
}
