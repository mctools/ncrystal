#ifndef ncrystal_h
#define ncrystal_h

/******************************************************************************/
/*                                                                            */
/*  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   */
/*                                                                            */
/*  Copyright 2015-2020 NCrystal developers                                   */
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
  typedef struct { void * internal; } ncrystal_atomdata_t;

  NCRYSTAL_API int  ncrystal_refcount( void* object );
  NCRYSTAL_API void ncrystal_ref( void* object );
  NCRYSTAL_API void ncrystal_unref( void* object );/*unrefs and deletes if count reaches 0*/
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
  NCRYSTAL_API void ncrystal_clear_info_caches(); /*NB: ncrystal_clear_caches below clears more! */
  NCRYSTAL_API void ncrystal_disable_caching(); /*NB: this concerns Info object caching only! */
  NCRYSTAL_API void ncrystal_enable_caching();  /*NB: this concerns Info object caching only! */
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
  NCRYSTAL_API double ncrystal_info_getdensity( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getnumberdensity( ncrystal_info_t );

  /*Access HKL info:                                                               */
  NCRYSTAL_API int ncrystal_info_nhkl( ncrystal_info_t ); /* -1 when not available */
  NCRYSTAL_API double ncrystal_info_hkl_dlower( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_hkl_dupper( ncrystal_info_t );
  NCRYSTAL_API void ncrystal_info_gethkl( ncrystal_info_t, int idx,
                                          int* h, int* k, int* l, int* multiplicity,
                                          double * dspacing, double* fsquared );

  /*Access AtomInfo:                                                               */
  NCRYSTAL_API unsigned ncrystal_info_natominfo( ncrystal_info_t );/* 0=unavail    */
  NCRYSTAL_API int ncrystal_info_hasatompos( ncrystal_info_t );
  NCRYSTAL_API int ncrystal_info_hasatommsd( ncrystal_info_t );
  NCRYSTAL_API void ncrystal_info_getatominfo( ncrystal_info_t, unsigned iatom,
                                               unsigned* atomdataindex,
                                               unsigned* number_per_unit_cell,
                                               double* debye_temp, double* msd );
  NCRYSTAL_API void ncrystal_info_getatompos( ncrystal_info_t,
                                              unsigned iatom, unsigned ipos,
                                              double* x, double* y, double* z );

  /*Debye temperatures (per-element also via atominfo):                            */
  NCRYSTAL_API int ncrystal_info_hasanydebyetemp( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getdebyetempbyelement( ncrystal_info_t,
                                                           unsigned atomdataindex );
  NCRYSTAL_API double ncrystal_info_getglobaldebyetemp( ncrystal_info_t );/* -1=unavail */

  /*Access dynamic info:                                                           */
  /*ditypeid: 0->nonscat, 1:freegas, 2:scatknl 3:vdos, 4:vdosdebye, 99:unknown     */
  NCRYSTAL_API unsigned ncrystal_info_ndyninfo( ncrystal_info_t );
  NCRYSTAL_API void ncrystal_dyninfo_base( ncrystal_info_t,
                                           unsigned idyninfo,
                                           double* fraction,
                                           unsigned* atomdataindex,
                                           double* temperature,
                                           unsigned* ditypeid );

  /* Extract scattering kernel for ditype 2,3,4 (vdoslux ignored for type 2).      */
  NCRYSTAL_API void ncrystal_dyninfo_extract_scatknl( ncrystal_info_t,
                                                      unsigned idyninfo,
                                                      unsigned vdoslux,
                                                      double* suggestedEmax,
                                                      unsigned* negrid,
                                                      unsigned* nalpha,
                                                      unsigned* nbeta,
                                                      const double** egrid,
                                                      const double** alphagrid,
                                                      const double** betagrid,
                                                      const double** sab );

  /* Access vdos data for ditype 3.                                                */
  NCRYSTAL_API void ncrystal_dyninfo_extract_vdos( ncrystal_info_t,
                                                   unsigned idyninfo,
                                                   double * egridMin,
                                                   double * egridMax,
                                                   unsigned * vdos_ndensity,
                                                   const double ** vdos_density );

  /* Access vdos debye temperature for ditype 4.                                   */
  NCRYSTAL_API void ncrystal_dyninfo_extract_vdosdebye( ncrystal_info_t,
                                                        unsigned idyninfo,
                                                        double * debye_temp );

  /* Access input curve ditype 3 (returns vdos_negrid=0 if not available).         */
  NCRYSTAL_API void ncrystal_dyninfo_extract_vdos_input( ncrystal_info_t,
                                                         unsigned idyninfo,
                                                         unsigned* vdos_negrid,
                                                         const double ** vdos_egrid,
                                                         unsigned* vdos_ndensity,
                                                         const double ** vdos_density );

  /* Convenience:                                                                  */
  NCRYSTAL_API double ncrystal_info_dspacing_from_hkl( ncrystal_info_t, int h, int k, int l );


  /* Composition (ncomponents=0 means composition unavailable):                    */
  NCRYSTAL_API unsigned ncrystal_info_ncomponents( ncrystal_info_t );
  NCRYSTAL_API void ncrystal_info_getcomponent( ncrystal_info_t, unsigned icomponent,
                                                unsigned* atomdataindex,
                                                double* fraction );


  /* Turn returned atomdata_idx's from calls above into actual ncrystal_atomdata_t */
  /* objects. The returned objects are ref-counted and the calling code should     */
  /* eventually unref them (with a call to ncrystal_unref) in order to prevent     */
  /* resource leaks.                                                               */
  NCRYSTAL_API ncrystal_atomdata_t ncrystal_create_atomdata( ncrystal_info_t,
                                                             unsigned atomdataindex );

  /* Get atom data fields. Each object falls in one of three categories:           */
  /* 1) Natural elements (ncomponents=A=0,Z>0)                                     */
  /* 2) Single isotope (ncomponents=0, Z>0, A>=Z)                                  */
  /* 3) Composite (A=0,ncomponents>1,Z>0 if all components share Z, otherwise Z=0) */
  /* Note that displaylabel=0 if atomdata object was sub-component.                */
  NCRYSTAL_API void ncrystal_atomdata_getfields( ncrystal_atomdata_t,
                                                 const char** displaylabel,
                                                 const char** description,
                                                 double* mass, double *incxs,
                                                 double* cohsl_fm, double* absxs,
                                                 unsigned* ncomponents,
                                                 unsigned* zval, unsigned* aval );

  /* Get atomdata and fraction of component:                                       */
  /* NB: Returned object should eventually be unref'ed by calling code.            */
  NCRYSTAL_API ncrystal_atomdata_t ncrystal_create_atomdata_subcomp( ncrystal_atomdata_t,
                                                                     unsigned icomponent,
                                                                     double* fraction );

  /* Custom data section:                                                         */
  NCRYSTAL_API unsigned ncrystal_info_ncustomsections( ncrystal_info_t );
  NCRYSTAL_API const char* ncrystal_info_customsec_name( ncrystal_info_t, unsigned isection );
  NCRYSTAL_API unsigned ncrystal_info_customsec_nlines( ncrystal_info_t, unsigned isection );
  NCRYSTAL_API unsigned ncrystal_info_customline_nparts( ncrystal_info_t, unsigned isection, unsigned iline );
  NCRYSTAL_API const char* ncrystal_info_customline_getpart( ncrystal_info_t, unsigned isection, unsigned iline, unsigned ipart );

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
  /* TODO: error file/line-no as well? */

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

  /* Access internal DB for isotopes and natural elements.                         */
  /* NB: Will return invalid handle in case lookup failed. Otherwise, the          */
  /* returned object should eventually be unref'ed by calling code:                */
  NCRYSTAL_API ncrystal_atomdata_t ncrystal_create_atomdata_fromdb( unsigned z,
                                                                    unsigned a );
  /* Version which accepts strings like "Al", "H2", "D" ...:                       */
  NCRYSTAL_API ncrystal_atomdata_t ncrystal_create_atomdata_fromdbstr( const char* );

  /* Get all (Z,A) values in internal DB (A=0 means natural element). The second   */
  /* fct accepts two preallocated arrays with length given by the first fct:       */
  NCRYSTAL_API unsigned ncrystal_atomdatadb_getnentries();
  NCRYSTAL_API void ncrystal_atomdatadb_getallentries( unsigned* zvals,
                                                       unsigned* avals );



  /* Extract NCMatCfg variables which can not be inferred from an ncrystal_info_t  */
  /* object and which might be needed in plugins (to be expanded as needed):       */
  NCRYSTAL_API double ncrystal_decodecfg_packfact( const char * cfgstr );
  NCRYSTAL_API unsigned ncrystal_decodecfg_vdoslux( const char * cfgstr );

  /* For serious scientific usage, users should register their own random          */
  /* generator function before using the genscatter functions. It must return      */
  /* numbers uniformly in [0,1):                                                   */
  NCRYSTAL_API void ncrystal_setrandgen( double (*rg)() );

  /* For special uses it is possible to trigger save/restore of the rng            */
  NCRYSTAL_API void ncrystal_save_randgen();    /* save & ref current randgen                   */
  NCRYSTAL_API void ncrystal_restore_randgen(); /* restore from and clear+unref last saved      */
  NCRYSTAL_API void ncrystal_setbuiltinrandgen(); /* for reproducibility use NCrystal's own rng */

  /* Clear various caches employed inside NCrystal:                                */
  NCRYSTAL_API void ncrystal_clear_caches();

  /* Register in-memory file data:                                                 */
  NCRYSTAL_API void ncrystal_register_in_mem_file_data(const char* virtual_filename,
                                                       const char* data);

  /* Browse (some) available files. Resulting string list must be deallocated by a */
  /* call to ncrystal_dealloc_stringlist by, and contains entries in the format    */
  /* name0,src0,hid0,name1,src1,hid1...:                                           */
  NCRYSTAL_API void ncrystal_get_file_list( const char* extension,
                                            unsigned* nstrs, char*** strs );

  /* Get list of plugins. Resulting string list must be deallocated by a call to   */
  /* ncrystal_dealloc_stringlist by, and contains entries in the format            */
  /* pluginname0,filename0,plugintype0,pluginname1,filename1,plugintype1,...:      */
  NCRYSTAL_API void ncrystal_get_plugin_list( unsigned* nstrs, char*** strs );

  /* Access file contents (must free afterwards with ncrystal_dealloc_string)      */
  NCRYSTAL_API char* ncrystal_get_file_contents( const char * name );

  /* Deallocate strings:                                                           */
  NCRYSTAL_API void ncrystal_dealloc_stringlist( unsigned len, char** );
  NCRYSTAL_API void ncrystal_dealloc_string( char* );

  /* NCrystal version info:                                                        */
#define NCRYSTAL_VERSION_MAJOR 2
#define NCRYSTAL_VERSION_MINOR 3
#define NCRYSTAL_VERSION_PATCH 1
#define NCRYSTAL_VERSION   2003001 /* (1000000*MAJOR+1000*MINOR+PATCH)             */
#define NCRYSTAL_VERSION_STR "2.3.1"
  NCRYSTAL_API int ncrystal_version(); /* returns NCRYSTAL_VERSION                  */
  NCRYSTAL_API const char * ncrystal_version_str(); /* returns NCRYSTAL_VERSION_STR */

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

  NCRYSTAL_API void ncrystal_genscatter_many( ncrystal_scatter_t,
                                              double ekin,
                                              const double (*direction)[3],
                                              unsigned long repeat,
                                              double * results_dirx,
                                              double * results_diry,
                                              double * results_dirz,
                                              double * results_dekin );

  NCRYSTAL_API void ncrystal_crosssection_nonoriented_many( ncrystal_process_t,
                                                            const double * ekin,
                                                            unsigned long n_ekin,
                                                            unsigned long repeat,
                                                            double* results );

#ifdef __cplusplus
}
#endif

#endif
