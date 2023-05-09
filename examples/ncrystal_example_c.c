
/******************************************************************************/
/*                                                                            */
/*  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   */
/*                                                                            */
/*  Copyright 2015-2023 NCrystal developers                                   */
/*                                                                            */
/*  Licensed under the Apache License, Version 2.0 (the "License");           */
/*  you may not use this file except in compliance with the License.          */
/*  You may obtain a copy of the License at                                   */
/*                                                                            */
/*      http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                            */
/*  Unless required by applicable law or agreed to in writing, software       */
/*  distributed under the License is distributed on an "AS IS" BASIS,         */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/*  See the License for the specific language governing permissions and       */
/*  limitations under the License.                                            */
/*                                                                            */
/******************************************************************************/

/* Include NCrystal C-interface: */
#include "NCrystal/ncrystal.h"

#include <stdio.h>
#include <stdlib.h>

double my_randgen() {
  /* Using rand from stdlib.h for this small example to generate numbers     */
  /* uniformly in (0,1] (note, this is NOT recommended for scientific work). */
  return 1.0-rand() / (RAND_MAX + 1.);
}

int main() {
  /* declarations first (to support ancient pedantic C compilers) */
  ncrystal_scatter_t pc, sc;
  ncrystal_process_t pc_proc, sc_proc;
  double wl, ekin, xsect;
  double ekin_final, cos_scat_angle;
  unsigned i;
  const double dir1[3] = { 0., 1., 1. };
  const double dir2[3] = { 1., 1., 0. };

  /* Setup random generator */
  ncrystal_setrandgen(my_randgen);

  /* Create and use polycrystalline aluminium */
  pc =  ncrystal_create_scatter( "Al_sg225.ncmat;temp=25C;dcutoff=0.5" );

  /* Cast to ncrystal_process_t for usage in cross-section functions: */
  pc_proc = ncrystal_cast_scat2proc(pc);

  wl = 2.5;
  ekin = ncrystal_wl2ekin(wl);
  ncrystal_crosssection_nonoriented(pc_proc,ekin,&xsect);
  printf("polycrystal Al x-sect at %g Aa is %g barn\n",wl,xsect);

  for (i=0;i<20;++i) {
    ncrystal_samplescatterisotropic( pc, ekin, &ekin_final, &cos_scat_angle);
    printf( "polycrystal random cos(scatangle) and delta_energy at %g Aa is %g and %g meV\n",
            wl,
            cos_scat_angle,
            (ekin_final-ekin)*1000.0 );
  }

  /* Create and use single-crystal germanium */
  wl = 1.540;/*angstrom*/
  ekin = ncrystal_wl2ekin(wl);

  sc = ncrystal_create_scatter("Ge_sg227.ncmat;dcutoff=0.5;mos= 40.0 arcsec"
                               ";dir1=@crys_hkl:5,1,1@lab:0,0,1"
                               ";dir2=@crys_hkl:0,-1,1@lab:0,1,0");
  sc_proc = ncrystal_cast_scat2proc(sc);

  ncrystal_crosssection(sc_proc,ekin,&dir1,&xsect);
  printf("singlecrystal Ge x-sect at %g Aa is %g barn (orientation 1)\n",wl,xsect);

  ncrystal_crosssection(sc_proc,ekin,&dir2,&xsect);
  printf("singlecrystal Ge x-sect at %g Aa is %g barn (orientation 2)\n",wl,xsect);

  /* Unref to release memory (just to be tidy): */
  ncrystal_unref(&pc);
  ncrystal_unref(&sc);

  return 0;
}
