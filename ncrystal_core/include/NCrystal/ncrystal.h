#ifndef ncrystal_h
#define ncrystal_h

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

  /*Can down-cast as well (returns invalid object with .internal null ptr in case  */
  /*object is not of that type:                                                    */
  NCRYSTAL_API ncrystal_scatter_t ncrystal_cast_proc2scat(ncrystal_process_t);
  NCRYSTAL_API ncrystal_absorption_t ncrystal_cast_proc2abs(ncrystal_process_t);

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

  /* Notice: ncrystal_scatter_t objects contain RNG streams, which a lot of other  */
  /* functions in this file are dedicated to handling.                             */

  /* Alternative creation method with new RNG stream (WARNING: Using this is       */
  /* intended for unit-tests only, as it is hard to guarantee two RNG streams are  */
  /* truly independent solely based on the seed value).                            */
  NCRYSTAL_API ncrystal_scatter_t ncrystal_create_scatter_builtinrng( const char * cfgstr,
                                                                      unsigned long seed );

  /* Cheaply clone scatter and absorption instances. The cloned objects will be    */
  /* using the same physics models and sharing any read-only data, but will be     */
  /* using their own copy of caches. For the case of scatter handles they will     */
  /* also get their own independent RNG stream. All in all, this means that        */
  /* the objects are safe to use concurrently in multi-threaded programming, as    */
  /* long as each thread gets its own clone. Cloned objects must still be cleaned  */
  /* up by calling ncrystal_unref. */
  NCRYSTAL_API ncrystal_scatter_t ncrystal_clone_scatter( ncrystal_scatter_t );
  NCRYSTAL_API ncrystal_absorption_t ncrystal_clone_absorption( ncrystal_absorption_t );

  /* Clone function where resulting object will use a specific rngstream index.    */
  /* All objects with the same indeed will share the same RNG state, so a sensible */
  /* strategy is to use the same index for all scatter objects which are to be     */
  /* used in the same thread:                                                      */
  NCRYSTAL_API ncrystal_scatter_t ncrystal_clone_scatter_rngbyidx( ncrystal_scatter_t,
                                                                   unsigned long rngstreamidx );

  /* Clone function where resulting object will use specific rngstream which has   */
  /* been set aside for the current thread. Thus, this function can be called      */
  /* from a given work-thread, in order to get a thread-safe scatter handle.       */
  NCRYSTAL_API ncrystal_scatter_t ncrystal_clone_scatter_rngforcurrentthread( ncrystal_scatter_t );

  /* Convenience function which creates objects directly from a data string        */
  /* rather than an on-disk or in-memory file. Such usage obviously precludes      */
  /* proper caching behind the scenes, and is intended for scenarios where the     */
  /* same data should not be used repeatedly. The ncrystal_xxx_t* arguments will   */
  /* be overriden with new handles (nullptrs results in no such object created).   */
  NCRYSTAL_API void ncrystal_multicreate_direct( const char* data,
                                                 const char* dataType,/*NULL => determine from data */
                                                 const char* cfg_params,/*e.g. "temp=300K;dcutoff=1"*/
                                                 ncrystal_info_t*,
                                                 ncrystal_scatter_t*,
                                                 ncrystal_absorption_t* );

  /* Factory availablity:                                                          */
  NCRYSTAL_API int ncrystal_has_factory( const char * name );

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Query ncrystal process handles (see previous section for casting         == */
  /*== scatter and absorption handles to ncrystal_process_t                     == */
  /*== NB: Also notice the "_many" versions of functions further below.         == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /*Name and UID of underlying ProcImpl::Process object:                           */
  NCRYSTAL_API const char * ncrystal_name(ncrystal_process_t);

  /*Determine if process is non-oriented (normally) or not (single-crystal):       */
  NCRYSTAL_API int ncrystal_isnonoriented(ncrystal_process_t);

  /*Access cross sections [barn] by neutron kinetic energy [eV]:                   */
  NCRYSTAL_API void ncrystal_crosssection_nonoriented( ncrystal_process_t,
                                                       double ekin,
                                                       double* result);
  NCRYSTAL_API void ncrystal_crosssection( ncrystal_process_t,
                                           double ekin,
                                           const double (*direction)[3],
                                           double* result );
  NCRYSTAL_API void ncrystal_domain( ncrystal_process_t,
                                     double* ekin_low, double* ekin_high);

  /*Generate random scatterings (neutron kinetic energy is in eV). The isotropic   */
  /*functions can only be called whe ncrystal_isnonoriented returns true (1).      */
  NCRYSTAL_API void ncrystal_samplescatterisotropic( ncrystal_scatter_t,
                                                     double ekin,
                                                     double* ekin_final,
                                                     double* cos_scat_angle );

  NCRYSTAL_API void ncrystal_samplescatter( ncrystal_scatter_t,
                                            double ekin,
                                            const double (*direction)[3],
                                            double* ekin_final,
                                            double (*direction_final)[3] );

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

  /*Access various scalar information:                                             */
  NCRYSTAL_API double ncrystal_info_getxsectabsorption( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getxsectfree( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getdensity( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getnumberdensity( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getsld( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_gettemperature( ncrystal_info_t );/*-1 if N/A. */

  /*State of matter (Unknown = 0, Solid = 1, Gas = 2, Liquid = 3)                  */
  NCRYSTAL_API int ncrystal_info_getstateofmatter( ncrystal_info_t );

  /* Access phase information (nphases=0 means single phase)                       */
  NCRYSTAL_API int ncrystal_info_nphases( ncrystal_info_t );
  NCRYSTAL_API ncrystal_info_t ncrystal_info_getphase( ncrystal_info_t,
                                                       int iphase,
                                                       double* fraction );

  /*Access HKL info:                                                               */
  NCRYSTAL_API int ncrystal_info_nhkl( ncrystal_info_t ); /* -1 when not available */
  NCRYSTAL_API double ncrystal_info_hkl_dlower( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_hkl_dupper( ncrystal_info_t );
  NCRYSTAL_API void ncrystal_info_gethkl( ncrystal_info_t, int idx,
                                          int* h, int* k, int* l, int* multiplicity,
                                          double * dspacing, double* fsquared );
  /*All HKL indices in a given group (returns first value h[0]==k[0]==l[0]==0 if not possible). */
  NCRYSTAL_API void ncrystal_info_gethkl_allindices( ncrystal_info_t, int idx,
                                                     int* h, int* k, int* l );/* arrays of length multiplicity/2 */

  NCRYSTAL_API double ncrystal_info_braggthreshold( ncrystal_info_t ); /* [Aa], -1 when not available */
  NCRYSTAL_API int ncrystal_info_hklinfotype( ncrystal_info_t ); /* integer casted value of HKLInfoType */

  /*Access AtomInfo:                                                               */
  NCRYSTAL_API unsigned ncrystal_info_natominfo( ncrystal_info_t );/* 0=unavail    */
  NCRYSTAL_API int ncrystal_info_hasatommsd( ncrystal_info_t );
  NCRYSTAL_API int ncrystal_info_hasatomdebyetemp( ncrystal_info_t );
  NCRYSTAL_API int ncrystal_info_hasdebyetemp( ncrystal_info_t );/* alias of hasatomdebyetemp */
  NCRYSTAL_API void ncrystal_info_getatominfo( ncrystal_info_t, unsigned iatom,
                                               unsigned* atomdataindex,
                                               unsigned* number_per_unit_cell,
                                               double* debye_temp, double* msd );
  NCRYSTAL_API void ncrystal_info_getatompos( ncrystal_info_t,
                                              unsigned iatom, unsigned ipos,
                                              double* x, double* y, double* z );

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

  /* Extract Sjolander Gn functions directly from VDOS curves. The res_gn_vals     */
  /* array should be freed with ncrystal_dealloc_doubleptr after usage.            */
  NCRYSTAL_API void ncrystal_raw_vdos2gn( const double* vdos_egrid,
                                          const double* vdos_density,
                                          unsigned vdos_egrid_npts,
                                          unsigned vdos_density_npts,
                                          double scattering_xs,
                                          double mass_amu,
                                          double temperature,
                                          unsigned nvalue,
                                          double* res_gn_xmin,
                                          double* res_gn_xmax,
                                          unsigned* res_gn_npts,
                                          double** res_gn_vals );

  /* Expand VDOS to scattering kernel with Sjolander's method, functions directly */
  /* from VDOS curves. The alpha, beta, and sab arrays should be freed with       */
  /* ncrystal_dealloc_doubleptr after usage. The order_weight_fct can optionally  */
  /* be used to change the contribution of a given Gn order to the kernel.        */
  NCRYSTAL_API void ncrystal_raw_vdos2knl( const double* vdos_egrid,
                                           const double* vdos_density,
                                           unsigned vdos_egrid_npts,
                                           unsigned vdos_density_npts,
                                           double scattering_xs,
                                           double mass_amu,
                                           double temperature,
                                           unsigned vdoslux,
                                           double (*order_weight_fct)( unsigned order ),
                                           unsigned* nalpha,
                                           unsigned* nbeta,
                                           double** alpha,
                                           double** beta,
                                           double** sab );

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


  /* Composition (always >=1 component). Note that for multiphase objects, the     */
  /* provided atomdataidx will be invalid, so it is important in general to use    */
  /* ncrystal_create_component_atomdata(..) to access the atomdata object          */
  /* associated with the i'th component, and NOT ncrystal_create_atomdata(..):     */
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

  /* Same, but via index in composition vector:                                    */
  NCRYSTAL_API ncrystal_atomdata_t ncrystal_create_component_atomdata( ncrystal_info_t,
                                                                       unsigned icomponent );

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

  /* Get flattened breakdown of composition i.e. to set up a structureless base    */
  /* material with same atoms in another application). Set prefernatelem to 1 or   */
  /* 0 to indicate whether the breakdowns should prefer natural elements (A=0)     */
  /* where possible, or always expand such to individial isotopes. The             */
  /* natelemprovider argument is a callback to a function which provides such      */
  /* natural abundances, and this is needed when prefernatelem==0 or when the      */
  /* composition contains both natural elements and isotopes of a particular       */
  /* element. If prefernatelem=1 the natelemprovider can be left out by passing    */
  /* a NULL ptr, which will then still work for all materials except for this with */
  /* the mentioned mixed compositions. Results are returned as a JSON string (list */
  /* of (Z,A,frac) entries where A=0 indicates natural elements), which must be    */
  /* cleaned up with ncrystal_dealloc_string).                                     */
  /*                                                                               */
  /* To keep the interface simple (if ugly), the natelemprovider function is       */
  /* passed the Z value along with a preallocated int and double array of length   */
  /* 128, which it must fill with corresponding A and abundance values. The return */
  /* value is used to indicate the number of isotopes present. Thus, at most 128   */
  /* isotopes can be returned (which should be plenty for all known elements).     */
  /* It should simply return 0 in case of an element with no nat. abundance data.  */
  char * ncrystal_get_flatcompos( ncrystal_info_t, int prefernatelem,
                                  unsigned (*natelemprovider)(unsigned,unsigned*,double*) );

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
  /*== Random number stream handling                                            == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /* Register custom RNG stream (it must return numbers uniformly in [0,1)). This  */
  /* RNG will be used for any subsequent calls to ncrystal_create_scatter, and it  */
  /* will NOT be changed when cloning ncrystal_scatter_t objects. Thus, for multi- */
  /* threaded usage, the caller should ensure that the function is thread-safe and */
  /* returns numbers from different streams in each thread (through suitable usage */
  /* of thread_local objects perhaps).                                             */
  NCRYSTAL_API void ncrystal_setrandgen( double (*rg)() );

  /* It is also possible to (re) set the RNG to the builtin generator (optionally  */
  /* by state or integer seed) */
  NCRYSTAL_API void ncrystal_setbuiltinrandgen();
  NCRYSTAL_API void ncrystal_setbuiltinrandgen_withseed(unsigned long seed);
  NCRYSTAL_API void ncrystal_setbuiltinrandgen_withstate(const char*);

  /* If supported (which it will NOT be if the RNG was set using the C API and the */
  /* ncrystal_setrandgen function), the state of the RNG stream can be accessed    */
  /* and manipulated via the following functions (returned strings must be free'd  */
  /* by calling ncrystal_dealloc_string). Note that setting the rng state will     */
  /* affect all objects sharing the RNG stream with the given scatter object (and  */
  /* those subsequently cloned from any of those). Note that if the provided state */
  /* originates in (the current version of) NCrystal's builtin RNG algorithm, it   */
  /* can always be used in the ncrystal_setrngstate_ofscatter function, even if    */
  /* the current RNG uses a different algorithm (it will simply be replaced).      */
  /* Finally note that getrngstate returns NULL if state manip. is not supported.  */
  NCRYSTAL_API int ncrystal_rngsupportsstatemanip_ofscatter( ncrystal_scatter_t );
  NCRYSTAL_API char* ncrystal_getrngstate_ofscatter( ncrystal_scatter_t );
  NCRYSTAL_API void ncrystal_setrngstate_ofscatter( ncrystal_scatter_t, const char* );

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Data sources                                                             == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /* Access TextData. Returns a string list of length 5:                           */
  /* [contents, uid(as string), sourcename, datatype, resolvedphyspath].           */
  /* The last entry is optional and will be an empty str if absent.                */
  /* Must free list with call to ncrystal_dealloc_stringlist.                      */
  NCRYSTAL_API char** ncrystal_get_text_data( const char * name );

  /* Register in-memory file data (as a special case data can be "ondisk://<path>" */
  /* to instead create a virtual alias for an on-disk file).:                      */
  NCRYSTAL_API void ncrystal_register_in_mem_file_data(const char* virtual_filename,
                                                       const char* data);

  /* Browse (some) available files. Resulting string list must be deallocated by a */
  /* call to ncrystal_dealloc_stringlist by, and contains entries in the format    */
  /* name0,src0,fact0,priority0,name1,src1,fact1,priority1,..:                     */
  NCRYSTAL_API void ncrystal_get_file_list( unsigned* nstrs, char*** strs );

  /* Add directory to custom search path or clear custom search path:              */
  NCRYSTAL_API void ncrystal_add_custom_search_dir( const char * dir );
  NCRYSTAL_API void ncrystal_remove_custom_search_dirs();

  /* Enable/disable standard search mechanisms (all enabled by default). For the   */
  /* case of ncrystal_enable_stddatalib, provide dir=NULL if you wish to use the   */
  /* usual stddatalib path (determined by NCRYSTAL_DATA_DIR environment and        */
  /* compile-time variables: */
  NCRYSTAL_API void ncrystal_enable_abspaths( int doEnable );
  NCRYSTAL_API void ncrystal_enable_relpaths( int doEnable );
  NCRYSTAL_API void ncrystal_enable_stddatalib( int doEnable, const char * dir );
  NCRYSTAL_API void ncrystal_enable_stdsearchpath( int doEnable );

  /* Remove all current data sources (supposedly in order to subsequently enable   */
  /* sources selectively). This also removes virtual files and caches.             */
  NCRYSTAL_API void ncrystal_remove_all_data_sources();


  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Miscellaneous                                                            == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /* Dump info to stdout (1st fct is same as calling 2nd one with verbosity_lvl=0):*/
  NCRYSTAL_API void ncrystal_dump(ncrystal_info_t);
  NCRYSTAL_API void ncrystal_dump_verbose(ncrystal_info_t, unsigned verbosity_lvl );

  /* Utility converting between neutron wavelength [Aa] to kinetic energy [eV]:    */
  NCRYSTAL_API double ncrystal_wl2ekin( double wl );
  NCRYSTAL_API double ncrystal_ekin2wl( double ekin );

  /* Extract extra debug information about objects (as JSON string which must be   */
  /* cleaned up with ncrystal_dealloc_string).                                     */
  NCRYSTAL_API char * ncrystal_dbg_process( ncrystal_process_t );

  /*UID of underlying Info or ProcImpl::Process object as string (must free with   */
  /*call to ncrystal_dealloc_string:                                               */
  NCRYSTAL_API char * ncrystal_process_uid(ncrystal_process_t);
  NCRYSTAL_API char * ncrystal_info_uid(ncrystal_info_t);

  /*Generate cfg-str variable documentation as string (must free with call to      */
  /*ncrystal_dealloc_string). Mode 0 (full), 1 (short), 2 (json):                  */
  NCRYSTAL_API char * ncrystal_gencfgstr_doc(int mode);

  /*Underlying UID (in case density value or cfgdata  was overridden):             */
  NCRYSTAL_API char * ncrystal_info_underlyinguid(ncrystal_info_t);


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

  /* Convert between atomic mean squared displacements and Debye temperatures.     */
  /* Units are kelvin, AMU, and Angstrom^2:                                        */
  double ncrystal_debyetemp2msd( double debyetemp, double temperature, double mass );
  double ncrystal_msd2debyetemp( double msd, double temperature, double mass );

  /* Extract information from VDOS curve (see NCVDOSEval.hh for details):          */
  void ncrystal_vdoseval( double vdos_emin, double vdos_emax,
                          unsigned vdos_ndensity, const double* vdos_density,
                          double temperature, double atom_mass_amu,
                          double* msd, double* debye_temp, double* gamma0,
                          double* temp_eff, double* origIntegral );

  /* Extract NCMatCfg variables which can not be inferred from an ncrystal_info_t  */
  /* object and which might be needed in plugins (to be expanded as needed).       */
  /* Returned strings must be cleaned up with ncrystal_dealloc_string.             */
  NCRYSTAL_API unsigned ncrystal_decodecfg_vdoslux( const char * cfgstr );
  NCRYSTAL_API char* ncrystal_decodecfg_json( const char * cfgstr );
  /* Parse and reencode cfg (as NCrystal::MatCfg(cfgstr).toStrCfg()):              */
  NCRYSTAL_API char* ncrystal_normalisecfg( const char * cfgstr );

  /* Clear various caches employed inside NCrystal:                                */
  NCRYSTAL_API void ncrystal_clear_caches();

  /* Get list of plugins. Resulting string list must be deallocated by a call to   */
  /* ncrystal_dealloc_stringlist by, and contains entries in the format            */
  /* pluginname0,filename0,plugintype0,pluginname1,filename1,plugintype1,...:      */
  NCRYSTAL_API void ncrystal_get_plugin_list( unsigned* nstrs, char*** strs );

  /* Deallocate strings / double arrays:                                           */
  NCRYSTAL_API void ncrystal_dealloc_stringlist( unsigned len, char** );
  NCRYSTAL_API void ncrystal_dealloc_string( char* );
  NCRYSTAL_API void ncrystal_dealloc_doubleptr( double* );

  /* NCrystal version info:                                                        */
#ifdef NCRYSTAL_VERSION_MAJOR
#  undef NCRYSTAL_VERSION_MAJOR
#endif
#ifdef NCRYSTAL_VERSION_MINOR
#  undef NCRYSTAL_VERSION_MINOR
#endif
#ifdef NCRYSTAL_VERSION_PATCH
#  undef NCRYSTAL_VERSION_PATCH
#endif
#ifdef NCRYSTAL_VERSION
#  undef NCRYSTAL_VERSION
#endif
#ifdef NCRYSTAL_VERSION_STR
#  undef NCRYSTAL_VERSION_STR
#endif
#define NCRYSTAL_VERSION_MAJOR 3
#define NCRYSTAL_VERSION_MINOR 6
#define NCRYSTAL_VERSION_PATCH 80
#define NCRYSTAL_VERSION   3006080 /* (1000000*MAJOR+1000*MINOR+PATCH)             */
#define NCRYSTAL_VERSION_STR "3.6.80"
  NCRYSTAL_API int ncrystal_version(); /* returns NCRYSTAL_VERSION                  */
  NCRYSTAL_API const char * ncrystal_version_str(); /* returns NCRYSTAL_VERSION_STR */

  /* Load raw NCMAT data into JSON structures. Must deallocate with call to        */
  /* ncrystal_dealloc_string as usual. (WARNING: JSON is incomplete for now!!!!!)  */
  NCRYSTAL_API char * ncrystal_ncmat2json( const char * textdataname );


  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Specialised array versions of some functions                             == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  NCRYSTAL_API void ncrystal_samplescatterisotropic_many( ncrystal_scatter_t,
                                                          const double * ekin,
                                                          unsigned long n_ekin,
                                                          unsigned long repeat,
                                                          double* results_ekin,
                                                          double* results_cos_scat_angle );

  NCRYSTAL_API void ncrystal_samplescatter_many( ncrystal_scatter_t,
                                                 double ekin,
                                                 const double (*direction)[3],
                                                 unsigned long repeat,
                                                 double* results_ekin,
                                                 double * results_dirx,
                                                 double * results_diry,
                                                 double * results_dirz );

  NCRYSTAL_API void ncrystal_crosssection_nonoriented_many( ncrystal_process_t,
                                                            const double * ekin,
                                                            unsigned long n_ekin,
                                                            unsigned long repeat,
                                                            double* results );

  /*============================================================================== */
  /*============================================================================== */
  /*==                                                                          == */
  /*== Various obsolete functions which are bound to be removed in future       == */
  /*== releases of NCrystal.                                                    == */
  /*==                                                                          == */
  /*============================================================================== */
  /*============================================================================== */

  /*Obsolete function which now always returns 1.0. Packing factors are now        */
  /*instead absorbed into the material densities:                                  */
  NCRYSTAL_API double ncrystal_decodecfg_packfact( const char * cfgstr );

  /*Obsolete function which now just is an alias for ncrystal_clear_caches above:  */
  NCRYSTAL_API void ncrystal_clear_info_caches();

  /*Obsolete function. Atom positions are now always available when                */
  /*ncrystal_info_natominfo returns a value greater than 0:                        */
  NCRYSTAL_API int ncrystal_info_hasatompos( ncrystal_info_t );

  /*Obsolete functions. Debye temperatures are now always available via the        */
  /*AtomInfo objects, and there is no longer a concept of a "global" Debye temp.:  */
  NCRYSTAL_API int ncrystal_info_hasanydebyetemp( ncrystal_info_t );
  NCRYSTAL_API double ncrystal_info_getdebyetempbyelement( ncrystal_info_t,
                                                           unsigned atomdataindex );
  NCRYSTAL_API double ncrystal_info_getglobaldebyetemp( ncrystal_info_t );/* -1=unavail */

  /* Obsolete function, the ncrystal_get_text_data function should be used instead */
  /* The returned content should be deallocated with ncrystal_dealloc_string.      */
  NCRYSTAL_API char* ncrystal_get_file_contents( const char * name );

  /*Obsolete genscatter functions. Users should use the ncrystal_samplescatter     */
  /*functions above instead:                                                       */
  NCRYSTAL_API void ncrystal_genscatter_nonoriented( ncrystal_scatter_t,
                                                     double ekin,
                                                     double* result_angle,
                                                     double* result_dekin );
  NCRYSTAL_API void ncrystal_genscatter( ncrystal_scatter_t,
                                         double ekin,
                                         const double (*direction)[3],
                                         double (*result_direction)[3],
                                         double* result_deltaekin );
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

#ifdef __cplusplus
}
#endif

#endif
