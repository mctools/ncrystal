
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

#include "NCrystal/ncrystal.h"
#include <stdint.h>
#include <stdio.h>

//A small dummy stateful RNG in a form suitable for passing to
//ncrystal_samplescatter_rs.

typedef struct
{
  uint32_t data;
} my_custom_rng_state_t;

my_custom_rng_state_t my_custom_rng_createstream( uint32_t seed )
{
  my_custom_rng_state_t state;
  state.data = 1789569706 + seed;
  return state;
}

double my_custom_rngfct( void * state_raw )
{
  my_custom_rng_state_t * state = (my_custom_rng_state_t*)state_raw;
  state->data = (1103515245 * (state->data) + 12345) % 2147483648;
  const double f = 1.0 / 2147483648;
  const double res = state->data * f;
  printf("...custom_rng produces %g\n",res);
  return res;
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  ncrystal_scatter_t sc
    = ncrystal_create_scatter("stdlib::Al_sg225.ncmat;comp=bragg;dcutoff=2.2");

  const double indir[3] = { 0.0, 0.0, 1.0 };
  double ekin_final;
  double outdir[3];

  printf("reset stream1 state\n");
  my_custom_rng_state_t stream1 = my_custom_rng_createstream( 12345 );

  printf("stream1 state: %lu\n",(unsigned long)(stream1.data));
  printf("Requesting scatter (stream1)\n");
  ncrystal_samplescatter_rs( my_custom_rngfct, &stream1, sc,
                             0.025,
                             (const double (*)[3]) & indir,
                             &ekin_final,
                             (double (*)[3]) & outdir );
  printf("stream1 state: %lu\n",(unsigned long)(stream1.data));
  printf("Requesting scatter (stream1)\n");
  ncrystal_samplescatter_rs( my_custom_rngfct, &stream1, sc,
                             0.025,
                             (const double (*)[3]) & indir,
                             &ekin_final,
                             (double (*)[3]) & outdir );
  printf("stream1 state: %lu\n",(unsigned long)(stream1.data));

  printf("reset stream1 state\n");
  stream1 = my_custom_rng_createstream( 12345 );
  printf("create stream2 as copy of stream1 state\n");
  my_custom_rng_state_t stream2 = stream1;

  printf("stream1 state: %lu\n",(unsigned long)(stream1.data));
  printf("stream2 state: %lu\n",(unsigned long)(stream2.data));

  printf("Requesting scatter (stream1)\n");
  ncrystal_samplescatter_rs( my_custom_rngfct, &stream1, sc,
                             0.025,
                             (const double (*)[3]) & indir,
                             &ekin_final,
                             (double (*)[3]) & outdir );
  printf("stream1 state: %lu\n",(unsigned long)(stream1.data));
  printf("stream2 state: %lu\n",(unsigned long)(stream2.data));

  printf("Requesting scatter (stream2)\n");
  ncrystal_samplescatter_rs( my_custom_rngfct, &stream2, sc,
                             0.025,
                             (const double (*)[3]) & indir,
                             &ekin_final,
                             (double (*)[3]) & outdir );
  printf("stream1 state: %lu\n",(unsigned long)(stream1.data));
  printf("stream2 state: %lu\n",(unsigned long)(stream2.data));

  printf("Requesting scatter (stream1)\n");
  ncrystal_samplescatter_rs( my_custom_rngfct, &stream1, sc,
                             0.025,
                             (const double (*)[3]) & indir,
                             &ekin_final,
                             (double (*)[3]) & outdir );
  printf("stream1 state: %lu\n",(unsigned long)(stream1.data));
  printf("stream2 state: %lu\n",(unsigned long)(stream2.data));

  ncrystal_unref(&sc);

  return 0;
}
