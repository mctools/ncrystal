
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

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  /* Trigger various obsoletion msgs */
  ncrystal_decodecfg_packfact("stdlib::Al_sg225.ncmat");
  ncrystal_clear_info_caches();
  /* Only trigger warning the first time: */
  ncrystal_clear_info_caches();
  ncrystal_clear_info_caches();
  ncrystal_info_t info = ncrystal_create_info("stdlib::Al_sg225.ncmat");
  ncrystal_scatter_t scatter
    = ncrystal_create_scatter("stdlib::Al_sg225.ncmat");
  ncrystal_info_hasatompos(info);

  ncrystal_info_hasanydebyetemp( info );
  ncrystal_info_getdebyetempbyelement( info, 0 );
  ncrystal_info_getglobaldebyetemp( info );
  char* foo = ncrystal_get_file_contents( "stdlib::Al_sg225.ncmat" );
  ncrystal_dealloc_string(foo);

  double angle, dekin;
  ncrystal_genscatter_nonoriented( scatter,0.025,&angle,&dekin );

  ncrystal_unref(&info);
  ncrystal_unref(&scatter);

  return 0;
}
