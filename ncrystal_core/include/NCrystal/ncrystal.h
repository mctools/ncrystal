#ifndef ncrystal_h
#define ncrystal_h

/******************************************************************************/
/*                                                                            */
/*  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   */
/*                                                                            */
/*  Copyright 2015-2019 NCrystal developers                                   */
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

/****************************/
/* C-interface for NCrystal */
/****************************/

#include "NCrystal/ncapi.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Object handle types. All are pointer-sized, thus small enough to pass    == */
  /*== around by value. The internal data-structures are reference-counted, so  == */
  /*== users should call ncrystal_ref/ncrystal_unref as appropriate if keeping  == */
  /*== such objects around. All ncrystal_create_xxx functions returns handles   == */
  /*== which have already have their reference counts increased, so users only  == */
  /*== need to call ncrystal_unref or invalidate on them after usage            == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  typedef struct { void * internal; } ncrystal_info_t;
  typedef struct { void * internal; } ncrystal_process_t;
  typedef struct { void * internal; } ncrystal_scatter_t;
  typedef struct { void * internal; } ncrystal_absorption_t;

  NCRYSTAL_API int  ncrystal_refcount( void* object );
  NCRYSTAL_API void ncrystal_ref( void* object );
  NCRYSTAL_API void ncrystal_unref( void* object );/*unrefs and deletes if count reaches 0*/
  NCRYSTAL_API void ncrystal_unrefnodelete( void* object ); /*unrefs but never deletes*/
  NCRYSTAL_API int  ncrystal_valid( void* object );
  NCRYSTAL_API void ncrystal_invalidate( void* object );/*invalidates handle (does not unref!)*/

  /*Casts might be needed to use shared interfaces (doesn't increase the           */
  /*ref-count of the underlying object):                                           */
  NCRYSTAL_API ncrystal_process_t ncrystal_cast_scat2proc(ncrystal_scatter_t);
  NCRYSTAL_API ncrystal_process_t ncrystal_cast_abs2proc(ncrystal_absorption_t);

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Factory functions. The cfgstr arguments will be passed directly to the   == */
  /*== MatCfg constructor. The file NCMatCfg.hh contains more info about the    == */
  /*== format, which always starts with the name of a datafile.                 == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  NCRYSTAL_API ncrystal_info_t ncrystal_create_info( const char * cfgstr );
  NCRYSTAL_API ncrystal_scatter_t ncrystal_create_scatter( const char * cfgstr );
  NCRYSTAL_API ncrystal_absorption_t ncrystal_create_absorption( const char * cfgstr );

  /* Fine tuning factory availability and caching                                  */
  NCRYSTAL_API void ncrystal_clear_info_caches();
  NCRYSTAL_API void ncrystal_disable_caching();
  NCRYSTAL_API void ncrystal_enable_caching();
  NCRYSTAL_API void ncrystal_clear_factory_registry();
  NCRYSTAL_API int ncrystal_has_factory( const char * name );

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Query ncrystal process handles (see previous section for casting         == */
  /*== scatter and absorption handles to ncrystal_process_t                     == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /*Name of object:                                                                */
  NCRYSTAL_API const char * ncrystal_name(ncrystal_process_t);

  /*Determine if process is non-oriented (normally) or not (single-crystal):       */
  NCRYSTAL_API int ncrystal_isnonoriented(ncrystal_process_t);

  /*Access cross-sections [barn] by neutron kinetic energy [eV]:                   */
  NCRYSTAL_API void ncrystal_crosssection_nonoriented( ncrystal_process_t,
                                                       double ekin,
                                                       double* result);
  NCRYSTAL_API void ncrystal_crosssection( ncrystal_process_t,
                                           double ekin,
                                           const double (*direction)[3],
                                           double* result );
  NCRYSTAL_API void ncrystal_domain( ncrystal_process_t,
                                     double* ekin_low, double* ekin_high);

  /*Generate random scatterings (radians, eV) by neutron kinetic energy [eV].      */
  NCRYSTAL_API void ncrystal_genscatter_nonoriented( ncrystal_scatter_t,
                                                     double ekin,
                                                     double* result_angle,
                                                     double* result_dekin );

  NCRYSTAL_API void ncrystal_genscatter( ncrystal_scatter_t,
                                         double ekin,
                                         const double (*direction)[3],
                                         double (*result_direction)[3],
                                         double* result_deltaekin );

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Query ncrystal info handles (see NCInfo.hh for more details)             == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /*Access crystal structure. Returns 0 if structure information is unavailable,   */
  /*otherwise the passed parameters are set and 1 is returned:                     */
  NCRYSTAL_API int ncrystal_info_getstructure( ncrystal_info_t,
                                               unsigned* spacegroup,
                                               double* lattice_a, double* lattice_b, double* lattice_c,
                                               double* alpha, double* beta, double* gamma,
                                               double* volume, unsigned* n_atoms );

  /*Access various scalar information (return values of -1 means unavailable:      */
  NCRYSTAL_API double ncrystal_info_gettemperature( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getxsectabsorption( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getxsectfree( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getdebyetemp( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getdensity( ncrystal_info_t );

  /*Access HKL info:                                                               */
  NCRYSTAL_API int ncrystal_info_nhkl( ncrystal_info_t ); /* -1 when not available */
  NCRYSTAL_API double ncrystal_info_hkl_dlower( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_hkl_dupper( ncrystal_info_t );
  NCRYSTAL_API void ncrystal_info_gethkl( ncrystal_info_t, int idx,
                                          int* h, int* k, int* l, int* multiplicity,
                                          double * dspacing, double* fsquared );

  /* TODO for NC2: more NCInfo data available here (including AtomInfo)            */

  /* Convenience: */
  NCRYSTAL_API double ncrystal_info_dspacing_from_hkl( ncrystal_info_t, int h, int k, int l );

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Error handling                                                           == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /*By default, any error encountered will result in a message printed to stdout   */
  /*and termination of the programme. This behaviour can be changed by calling     */
  /*ncrystal_sethaltonerror(0) and ncrystal_setquietonerror(1) respectively.       */
  NCRYSTAL_API int ncrystal_sethaltonerror(int);/* returns old value */
  NCRYSTAL_API int ncrystal_setquietonerror(int);/* returns old value */

  /*If not halting on error, these functions can be used to access information     */
  /*about errors encountered:                                                      */
  NCRYSTAL_API int ncrystal_error();/* returns 1 if an error condition occurred. */
  NCRYSTAL_API const char * ncrystal_lasterror();/* returns description of last error (NULL if none) */
  NCRYSTAL_API const char * ncrystal_lasterrortype();/* returns description of last error (NULL if none) */
  /* TODO for NC2: error file/line-no as well? */

  NCRYSTAL_API void ncrystal_clearerror();/* clears previous error if any */

  /*Another option is to provide a custom error handler which will be called on    */
  /*each error:                                                                    */
  NCRYSTAL_API void ncrystal_seterrhandler(void (*handler)(char*,char*));

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Miscellaneous                                                            == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /* Dump info to stdout:                                                          */
  NCRYSTAL_API void ncrystal_dump(ncrystal_info_t);

  /* Utility converting between neutron wavelength [Aa] to kinetic energy [eV]:    */
  NCRYSTAL_API double ncrystal_wl2ekin( double wl );
  NCRYSTAL_API double ncrystal_ekin2wl( double ekin );

  /* Extract NCMatCfg variables which can not be inferred from an ncrystal_info_t  */
  /* object and which might be needed in plugins (to be expanded as needed):       */
  NCRYSTAL_API double ncrystal_decodecfg_packfact( const char * cfgstr );

  /* For serious scientific usage, users should register their own random          */
  /* generator function before using the genscatter functions. It must return      */
  /* numbers uniformly in [0,1):                                                   */
  NCRYSTAL_API void ncrystal_setrandgen( double (*rg)() );

  /* For special uses it is possible to trigger save/restore of the rng            */
  NCRYSTAL_API void ncrystal_save_randgen();    /* save & ref current randgen                   */
  NCRYSTAL_API void ncrystal_restore_randgen(); /* restore from and clear+unref last saved      */
  NCRYSTAL_API void ncrystal_setbuiltinrandgen(); /* for reproducibility use NCrystals own rng   */

  /* NCrystal version info:                                                        */
#define NCRYSTAL_VERSION_MAJOR 1
#define NCRYSTAL_VERSION_MINOR 0
#define NCRYSTAL_VERSION_PATCH 0
#define NCRYSTAL_VERSION   1000000 /* (1000000*MAJOR+1000*MINOR+PATCH)                */
#define NCRYSTAL_VERSION_STR "1.0.0"
  NCRYSTAL_API int ncrystal_version(); /* returns NCRYSTAL_VERSION                              */
  NCRYSTAL_API const char * ncrystal_version_str(); /* returns NCRYSTAL_VERSION_STR             */

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Specialised array versions of some functions                             == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  NCRYSTAL_API void ncrystal_genscatter_nonoriented_many( ncrystal_scatter_t,
                                                          const double * ekin,
                                                          unsigned long n_ekin,
                                                          unsigned long repeat,
                                                          double* results_angle,
                                                          double* results_dekin );

  NCRYSTAL_API void ncrystal_crosssection_nonoriented_many( ncrystal_process_t,
                                                            const double * ekin,
                                                            unsigned long n_ekin,
                                                            unsigned long repeat,
                                                            double* results );

#ifdef __cplusplus
}
#endif

#endif
