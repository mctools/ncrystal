
/******************************************************************************/
/*                                                                            */
/*  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   */
/*                                                                            */
/*  Copyright 2015-2026 NCrystal developers                                   */
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

/* Test of the C API for filter tables (ncrystal_filtertable), for a Be     */
/* filter at 80 K.                                                           */

#include "NCrystal/ncrystal.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static void require( int cond, const char * what )
{
  if ( !cond ) {
    printf( "FAILED: %s\n", what );
    exit( 1 );
  }
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  ncrystal_sethaltonerror( 0 );
  ncrystal_setquietonerror( 1 );

  unsigned n = 0, i, iedge;
  double * wl = NULL;
  double * macroxs = NULL;
  ncrystal_filtertable( "stdlib::Be_sg194.ncmat;temp=80K", &n, &wl,
                        &macroxs, NULL );
  require( !ncrystal_error(), "no error" );
  require( n > 10 && wl && macroxs, "table returned" );
  require( wl[0] == 0.0 && wl[n-1] > 20.0, "range from 0 to long wavelengths" );
  for ( i = 1; i < n; ++i )
    require( wl[i] >= wl[i-1] && macroxs[i] >= 0.0, "sorted and non-negative" );
  printf( "Be 80K: %u points, macroxs(0) = %.4g/cm\n", n, macroxs[0] );

  /* The Be filter cut-off: the last Bragg edge (a pair of points with the  */
  /* same wavelength) is at 3.96 Aa, with a large drop of the cross section: */
  iedge = n;
  for ( i = 0; i + 1 < n; ++i )
    if ( wl[i] == wl[i+1] )
      iedge = i;
  require( iedge < n, "Bragg edges in table" );
  require( fabs( wl[iedge] - 3.96 ) < 0.01, "last Bragg edge at 3.96 Aa" );
  require( macroxs[iedge] > 20.0 * macroxs[iedge+1], "Bragg cut-off" );
  printf( "  last Bragg edge at %.4f Aa: %.4g/cm -> %.4g/cm\n",
          wl[iedge], macroxs[iedge], macroxs[iedge+1] );
  ncrystal_dealloc_doubleptr( wl );
  ncrystal_dealloc_doubleptr( macroxs );

  /* Errors: unsupported options, and an oriented material. No arrays are    */
  /* returned. */
  const char * cfgs[] = { "stdlib::Al_sg225.ncmat",
                          "stdlib::Al_sg225.ncmat;mos=0.3deg"
                          ";dir1=@crys_hkl:0,0,1@lab:0,0,1"
                          ";dir2=@crys_hkl:0,1,0@lab:0,1,0" };
  const char * opts[] = { "foo=1", NULL };
  for ( i = 0; i < 2; ++i ) {
    ncrystal_filtertable( cfgs[i], &n, &wl, &macroxs, opts[i] );
    require( ncrystal_error() && n == 0 && !wl && !macroxs, "error expected" );
    printf( "Expected error: %s\n", ncrystal_lasterror() );
    ncrystal_clearerror();
  }
  return 0;
}
