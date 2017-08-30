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

#include "NCrystal/NCDump.hh"
#include "NCrystal/NCInfo.hh"
#include <cstdio>

void NCrystal::dump(const NCrystal::Info*c)
{
  const char * hr = "---------------------------------------------------------\n";
  if (c->hasStructureInfo()) {
    const StructureInfo& si = c->getStructureInfo();
    printf("%s", hr);
    if (si.spacegroup!=0)
      printf("Space group number      : %i\n",si.spacegroup);
    printf("Lattice spacings   [Aa] : %g %g %g\n",si.lattice_a,si.lattice_b,si.lattice_c);
    printf("Lattice angles    [deg] : %g %g %g\n",si.alpha,si.beta,si.gamma);
    printf("Unit cell volume [Aa^3] : %g\n",si.volume);
    printf("Atoms / unit cell       : %i\n",si.n_atoms);
  };

  if (c->hasAtomInfo()) {
    printf("%s", hr);
    AtomList::const_iterator itE = c->atomInfoEnd();
    unsigned ntot = 0;
    for (AtomList::const_iterator it = c->atomInfoBegin();it!=itE;++it)
      ntot += it->number_per_unit_cell;
    printf("Atoms per unit cell (total %i):\n",ntot);
    for (AtomList::const_iterator it = c->atomInfoBegin();it!=itE;++it) {
      if (c->hasPerElementDebyeTemperature()) {
        printf("     %i Z=%i atoms [T_Debye=%gK]\n",
               it->number_per_unit_cell,it->atomic_number,it->debye_temp);
      } else {
        printf("     %i Z=%i atoms\n",it->number_per_unit_cell,it->atomic_number);
      }
    }
  }

  if (c->hasDensity()) {
    printf("%s", hr);
    printf("   Density : %g g/cm3\n",c->getDensity());
  }

  if (c->hasTemperature()) {
    printf("%s", hr);
    printf("   Temperature : %g kelvin\n",c->getTemperature());
  }

  if (c->hasDebyeTemperature()) {
    printf("%s", hr);
    printf("   Debye temperature (global) : %g kelvin\n",c->getDebyeTemperature());
  }

  if (c->hasXSectAbsorption()||c->hasXSectFree()) {
    printf("%s", hr);
    printf("Neutron cross-sections:\n");
    if (c->hasXSectAbsorption())
      printf("   Absorption at 2200m/s : %g barn\n", c->getXSectAbsorption());
    if (c->hasXSectFree())
      printf("   Free scattering       : %g barn\n", c->getXSectFree());
  }

  if (c->providesNonBraggXSects()) {
    printf("%s", hr);
    printf("Provides non-bragg cross-section calculations. A few values are:\n");
    printf("   lambda[Aa]  sigma_scat[barn]\n");
    double ll[] = {0.5, 1.0, 1.798, 2.5, 5, 10, 20 };
    for (unsigned i = 0; i < sizeof(ll)/sizeof(ll[0]); ++i)
      printf("%13g %17g\n",ll[i],c->xsectScatNonBragg(ll[i]));
  }
  if (c->hasHKLInfo()) {
    printf("%s", hr);
    printf("HKL planes (d_lower = %g Aa, d_upper = %g Aa):\n",c->hklDLower(),c->hklDUpper());
    printf("  H   K   L  d_hkl[Aa] Multiplicity FSquared[barn]%s\n",
           (c->hasExpandedHKLInfo()?" Expanded-HKL-list":""));
    HKLList::const_iterator itE = c->hklEnd();
    for (HKLList::const_iterator it = c->hklBegin();it!=itE;++it) {
      printf("%3i %3i %3i %10g %12i %14g%s",it->h,it->k,it->l,it->dspacing,
             it->multiplicity,it->fsquared,(it->eqv_hkl?"":"\n"));
      if (it->eqv_hkl) {
        const short * eqv_hkl = it->eqv_hkl;
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
