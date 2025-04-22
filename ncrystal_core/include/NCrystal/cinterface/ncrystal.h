#ifndef ncrystal_h
#define ncrystal_h

/******************************************************************************/
/*                                                                            */
/*  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   */
/*                                                                            */
/*  Copyright 2015-2025 NCrystal developers                                   */
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
  /*== First some tedious definitions to effectively implement a namespace in C == */
  /*== symbols, if requested through the NCRYSTAL_NAMESPACE_PROTECTION variable == */
  /*== (cf. also ncapi.h). Normally, NCRYSTAL_NAMESPACE_PROTECTION is not set   == */
  /*== and a define like ncrystal_create_info simply expands to the same name,  == */
  /*== resulting in the same symbol in the compiled binaries. However, if for   == */
  /*== instance =DNCRYSTAL_NAMESPACE_PROTECTION=foobar, then defines such as    == */
  /*== ncrystal_create_info will instead expand to ncrystalfoobar_create_info,  == */
  /*== resulting in namespaced symbols in the libraries. Due to the following   == */
  /*== list of defines, user code using the C API here should almost always     == */
  /*== write ncrystal_create_info(...) in their code.                           == */
  /*============================================================================== */
  /*============================================================================== */

  /* Note: when there is no namespace we get self-referential macro's below, but   */
  /* those are in fact allowed by the standard and won't be evaluated recursively. */

#ifdef ncrystal_absorption_t
#  undef ncrystal_absorption_t
#endif
#define ncrystal_absorption_t NCRYSTAL_APPLY_C_NAMESPACE(absorption_t)
#ifdef ncrystal_add_custom_search_dir
#  undef ncrystal_add_custom_search_dir
#endif
#define ncrystal_add_custom_search_dir NCRYSTAL_APPLY_C_NAMESPACE(add_custom_search_dir)
#ifdef ncrystal_atomdata_getfields
#  undef ncrystal_atomdata_getfields
#endif
#define ncrystal_atomdata_getfields NCRYSTAL_APPLY_C_NAMESPACE(atomdata_getfields)
#ifdef ncrystal_atomdata_t
#  undef ncrystal_atomdata_t
#endif
#define ncrystal_atomdata_t NCRYSTAL_APPLY_C_NAMESPACE(atomdata_t)
#ifdef ncrystal_atomdatadb_getallentries
#  undef ncrystal_atomdatadb_getallentries
#endif
#define ncrystal_atomdatadb_getallentries NCRYSTAL_APPLY_C_NAMESPACE(atomdatadb_getallentries)
#ifdef ncrystal_atomdatadb_getnentries
#  undef ncrystal_atomdatadb_getnentries
#endif
#define ncrystal_atomdatadb_getnentries NCRYSTAL_APPLY_C_NAMESPACE(atomdatadb_getnentries)
#ifdef ncrystal_cast_abs2proc
#  undef ncrystal_cast_abs2proc
#endif
#define ncrystal_cast_abs2proc NCRYSTAL_APPLY_C_NAMESPACE(cast_abs2proc)
#ifdef ncrystal_cast_proc2abs
#  undef ncrystal_cast_proc2abs
#endif
#define ncrystal_cast_proc2abs NCRYSTAL_APPLY_C_NAMESPACE(cast_proc2abs)
#ifdef ncrystal_cast_proc2scat
#  undef ncrystal_cast_proc2scat
#endif
#define ncrystal_cast_proc2scat NCRYSTAL_APPLY_C_NAMESPACE(cast_proc2scat)
#ifdef ncrystal_cast_scat2proc
#  undef ncrystal_cast_scat2proc
#endif
#define ncrystal_cast_scat2proc NCRYSTAL_APPLY_C_NAMESPACE(cast_scat2proc)
#ifdef ncrystal_clear_caches
#  undef ncrystal_clear_caches
#endif
#define ncrystal_clear_caches NCRYSTAL_APPLY_C_NAMESPACE(clear_caches)
#ifdef ncrystal_clear_info_caches
#  undef ncrystal_clear_info_caches
#endif
#define ncrystal_clear_info_caches NCRYSTAL_APPLY_C_NAMESPACE(clear_info_caches)
#ifdef ncrystal_clearerror
#  undef ncrystal_clearerror
#endif
#define ncrystal_clearerror NCRYSTAL_APPLY_C_NAMESPACE(clearerror)
#ifdef ncrystal_clone_absorption
#  undef ncrystal_clone_absorption
#endif
#define ncrystal_clone_absorption NCRYSTAL_APPLY_C_NAMESPACE(clone_absorption)
#ifdef ncrystal_clone_scatter
#  undef ncrystal_clone_scatter
#endif
#define ncrystal_clone_scatter NCRYSTAL_APPLY_C_NAMESPACE(clone_scatter)
#ifdef ncrystal_clone_scatter_rngbyidx
#  undef ncrystal_clone_scatter_rngbyidx
#endif
#define ncrystal_clone_scatter_rngbyidx NCRYSTAL_APPLY_C_NAMESPACE(clone_scatter_rngbyidx)
#ifdef ncrystal_clone_scatter_rngforcurrentthread
#  undef ncrystal_clone_scatter_rngforcurrentthread
#endif
#define ncrystal_clone_scatter_rngforcurrentthread NCRYSTAL_APPLY_C_NAMESPACE(clone_scatter_rngforcurrentthread)
#ifdef ncrystal_create_absorption
#  undef ncrystal_create_absorption
#endif
#define ncrystal_create_absorption NCRYSTAL_APPLY_C_NAMESPACE(create_absorption)
#ifdef ncrystal_create_atomdata
#  undef ncrystal_create_atomdata
#endif
#define ncrystal_create_atomdata NCRYSTAL_APPLY_C_NAMESPACE(create_atomdata)
#ifdef ncrystal_create_atomdata_fromdb
#  undef ncrystal_create_atomdata_fromdb
#endif
#define ncrystal_create_atomdata_fromdb NCRYSTAL_APPLY_C_NAMESPACE(create_atomdata_fromdb)
#ifdef ncrystal_create_atomdata_fromdbstr
#  undef ncrystal_create_atomdata_fromdbstr
#endif
#define ncrystal_create_atomdata_fromdbstr NCRYSTAL_APPLY_C_NAMESPACE(create_atomdata_fromdbstr)
#ifdef ncrystal_create_atomdata_subcomp
#  undef ncrystal_create_atomdata_subcomp
#endif
#define ncrystal_create_atomdata_subcomp NCRYSTAL_APPLY_C_NAMESPACE(create_atomdata_subcomp)
#ifdef ncrystal_create_component_atomdata
#  undef ncrystal_create_component_atomdata
#endif
#define ncrystal_create_component_atomdata NCRYSTAL_APPLY_C_NAMESPACE(create_component_atomdata)
#ifdef ncrystal_create_info
#  undef ncrystal_create_info
#endif
#define ncrystal_create_info NCRYSTAL_APPLY_C_NAMESPACE(create_info)
#ifdef ncrystal_create_scatter
#  undef ncrystal_create_scatter
#endif
#define ncrystal_create_scatter NCRYSTAL_APPLY_C_NAMESPACE(create_scatter)
#ifdef ncrystal_create_scatter_builtinrng
#  undef ncrystal_create_scatter_builtinrng
#endif
#define ncrystal_create_scatter_builtinrng NCRYSTAL_APPLY_C_NAMESPACE(create_scatter_builtinrng)
#ifdef ncrystal_crosssection
#  undef ncrystal_crosssection
#endif
#define ncrystal_crosssection NCRYSTAL_APPLY_C_NAMESPACE(crosssection)
#ifdef ncrystal_crosssection_nonoriented
#  undef ncrystal_crosssection_nonoriented
#endif
#define ncrystal_crosssection_nonoriented NCRYSTAL_APPLY_C_NAMESPACE(crosssection_nonoriented)
#ifdef ncrystal_crosssection_nonoriented_many
#  undef ncrystal_crosssection_nonoriented_many
#endif
#define ncrystal_crosssection_nonoriented_many NCRYSTAL_APPLY_C_NAMESPACE(crosssection_nonoriented_many)
#ifdef ncrystal_dbg_process
#  undef ncrystal_dbg_process
#endif
#define ncrystal_dbg_process NCRYSTAL_APPLY_C_NAMESPACE(dbg_process)
#ifdef ncrystal_enable_factory_threadpool
#  undef ncrystal_enable_factory_threadpool
#endif
#define ncrystal_enable_factory_threadpool NCRYSTAL_APPLY_C_NAMESPACE(enable_factory_threadpool)
#ifdef ncrystal_dealloc_doubleptr
#  undef ncrystal_dealloc_doubleptr
#endif
#define ncrystal_dealloc_doubleptr NCRYSTAL_APPLY_C_NAMESPACE(dealloc_doubleptr)
#ifdef ncrystal_dealloc_string
#  undef ncrystal_dealloc_string
#endif
#define ncrystal_dealloc_string NCRYSTAL_APPLY_C_NAMESPACE(dealloc_string)
#ifdef ncrystal_dealloc_stringlist
#  undef ncrystal_dealloc_stringlist
#endif
#define ncrystal_dealloc_stringlist NCRYSTAL_APPLY_C_NAMESPACE(dealloc_stringlist)
#ifdef ncrystal_debyetemp2msd
#  undef ncrystal_debyetemp2msd
#endif
#define ncrystal_debyetemp2msd NCRYSTAL_APPLY_C_NAMESPACE(debyetemp2msd)
#ifdef ncrystal_decodecfg_json
#  undef ncrystal_decodecfg_json
#endif
#define ncrystal_decodecfg_json NCRYSTAL_APPLY_C_NAMESPACE(decodecfg_json)
#ifdef ncrystal_decodecfg_packfact
#  undef ncrystal_decodecfg_packfact
#endif
#define ncrystal_decodecfg_packfact NCRYSTAL_APPLY_C_NAMESPACE(decodecfg_packfact)
#ifdef ncrystal_decodecfg_vdoslux
#  undef ncrystal_decodecfg_vdoslux
#endif
#define ncrystal_decodecfg_vdoslux NCRYSTAL_APPLY_C_NAMESPACE(decodecfg_vdoslux)
#ifdef ncrystal_domain
#  undef ncrystal_domain
#endif
#define ncrystal_domain NCRYSTAL_APPLY_C_NAMESPACE(domain)
#ifdef ncrystal_dump
#  undef ncrystal_dump
#endif
#define ncrystal_dump NCRYSTAL_APPLY_C_NAMESPACE(dump)
#ifdef ncrystal_dump_verbose
#  undef ncrystal_dump_verbose
#endif
#define ncrystal_dump_verbose NCRYSTAL_APPLY_C_NAMESPACE(dump_verbose)
#ifdef ncrystal_dump_tostr
#  undef ncrystal_dump_tostr
#endif
#define ncrystal_dump_tostr NCRYSTAL_APPLY_C_NAMESPACE(dump_tostr)
#ifdef ncrystal_dyninfo_base
#  undef ncrystal_dyninfo_base
#endif
#define ncrystal_dyninfo_base NCRYSTAL_APPLY_C_NAMESPACE(dyninfo_base)
#ifdef ncrystal_dyninfo_extract_scatknl
#  undef ncrystal_dyninfo_extract_scatknl
#endif
#define ncrystal_dyninfo_extract_scatknl NCRYSTAL_APPLY_C_NAMESPACE(dyninfo_extract_scatknl)
#ifdef ncrystal_dyninfo_extract_vdos
#  undef ncrystal_dyninfo_extract_vdos
#endif
#define ncrystal_dyninfo_extract_vdos NCRYSTAL_APPLY_C_NAMESPACE(dyninfo_extract_vdos)
#ifdef ncrystal_dyninfo_extract_vdos_input
#  undef ncrystal_dyninfo_extract_vdos_input
#endif
#define ncrystal_dyninfo_extract_vdos_input NCRYSTAL_APPLY_C_NAMESPACE(dyninfo_extract_vdos_input)
#ifdef ncrystal_dyninfo_extract_vdosdebye
#  undef ncrystal_dyninfo_extract_vdosdebye
#endif
#define ncrystal_dyninfo_extract_vdosdebye NCRYSTAL_APPLY_C_NAMESPACE(dyninfo_extract_vdosdebye)
#ifdef ncrystal_ekin2wl
#  undef ncrystal_ekin2wl
#endif
#define ncrystal_ekin2wl NCRYSTAL_APPLY_C_NAMESPACE(ekin2wl)
#ifdef ncrystal_enable_abspaths
#  undef ncrystal_enable_abspaths
#endif
#define ncrystal_enable_abspaths NCRYSTAL_APPLY_C_NAMESPACE(enable_abspaths)
#ifdef ncrystal_enable_relpaths
#  undef ncrystal_enable_relpaths
#endif
#define ncrystal_enable_relpaths NCRYSTAL_APPLY_C_NAMESPACE(enable_relpaths)
#ifdef ncrystal_enable_stddatalib
#  undef ncrystal_enable_stddatalib
#endif
#define ncrystal_enable_stddatalib NCRYSTAL_APPLY_C_NAMESPACE(enable_stddatalib)
#ifdef ncrystal_enable_stdsearchpath
#  undef ncrystal_enable_stdsearchpath
#endif
#define ncrystal_enable_stdsearchpath NCRYSTAL_APPLY_C_NAMESPACE(enable_stdsearchpath)
#ifdef ncrystal_error
#  undef ncrystal_error
#endif
#define ncrystal_error NCRYSTAL_APPLY_C_NAMESPACE(error)
#ifdef ncrystal_gencfgstr_doc
#  undef ncrystal_gencfgstr_doc
#endif
#define ncrystal_gencfgstr_doc NCRYSTAL_APPLY_C_NAMESPACE(gencfgstr_doc)
#ifdef ncrystal_genscatter
#  undef ncrystal_genscatter
#endif
#define ncrystal_genscatter NCRYSTAL_APPLY_C_NAMESPACE(genscatter)
#ifdef ncrystal_genscatter_many
#  undef ncrystal_genscatter_many
#endif
#define ncrystal_genscatter_many NCRYSTAL_APPLY_C_NAMESPACE(genscatter_many)
#ifdef ncrystal_genscatter_nonoriented
#  undef ncrystal_genscatter_nonoriented
#endif
#define ncrystal_genscatter_nonoriented NCRYSTAL_APPLY_C_NAMESPACE(genscatter_nonoriented)
#ifdef ncrystal_genscatter_nonoriented_many
#  undef ncrystal_genscatter_nonoriented_many
#endif
#define ncrystal_genscatter_nonoriented_many NCRYSTAL_APPLY_C_NAMESPACE(genscatter_nonoriented_many)
#ifdef ncrystal_get_file_contents
#  undef ncrystal_get_file_contents
#endif
#define ncrystal_get_file_contents NCRYSTAL_APPLY_C_NAMESPACE(get_file_contents)
#ifdef ncrystal_get_file_list
#  undef ncrystal_get_file_list
#endif
#define ncrystal_get_file_list NCRYSTAL_APPLY_C_NAMESPACE(get_file_list)
#ifdef ncrystal_get_flatcompos
#  undef ncrystal_get_flatcompos
#endif
#define ncrystal_get_flatcompos NCRYSTAL_APPLY_C_NAMESPACE(get_flatcompos)
#ifdef ncrystal_get_plugin_list
#  undef ncrystal_get_plugin_list
#endif
#define ncrystal_get_plugin_list NCRYSTAL_APPLY_C_NAMESPACE(get_plugin_list)
#ifdef ncrystal_get_text_data
#  undef ncrystal_get_text_data
#endif
#define ncrystal_get_text_data NCRYSTAL_APPLY_C_NAMESPACE(get_text_data)
#ifdef ncrystal_getrngstate_ofscatter
#  undef ncrystal_getrngstate_ofscatter
#endif
#define ncrystal_getrngstate_ofscatter NCRYSTAL_APPLY_C_NAMESPACE(getrngstate_ofscatter)
#ifdef ncrystal_has_factory
#  undef ncrystal_has_factory
#endif
#define ncrystal_has_factory NCRYSTAL_APPLY_C_NAMESPACE(has_factory)
#ifdef ncrystal_info_braggthreshold
#  undef ncrystal_info_braggthreshold
#endif
#define ncrystal_info_braggthreshold NCRYSTAL_APPLY_C_NAMESPACE(info_braggthreshold)
#ifdef ncrystal_info_customline_getpart
#  undef ncrystal_info_customline_getpart
#endif
#define ncrystal_info_customline_getpart NCRYSTAL_APPLY_C_NAMESPACE(info_customline_getpart)
#ifdef ncrystal_info_customline_nparts
#  undef ncrystal_info_customline_nparts
#endif
#define ncrystal_info_customline_nparts NCRYSTAL_APPLY_C_NAMESPACE(info_customline_nparts)
#ifdef ncrystal_info_customsec_name
#  undef ncrystal_info_customsec_name
#endif
#define ncrystal_info_customsec_name NCRYSTAL_APPLY_C_NAMESPACE(info_customsec_name)
#ifdef ncrystal_info_customsec_nlines
#  undef ncrystal_info_customsec_nlines
#endif
#define ncrystal_info_customsec_nlines NCRYSTAL_APPLY_C_NAMESPACE(info_customsec_nlines)
#ifdef ncrystal_info_dspacing_from_hkl
#  undef ncrystal_info_dspacing_from_hkl
#endif
#define ncrystal_info_dspacing_from_hkl NCRYSTAL_APPLY_C_NAMESPACE(info_dspacing_from_hkl)
#ifdef ncrystal_info_getatominfo
#  undef ncrystal_info_getatominfo
#endif
#define ncrystal_info_getatominfo NCRYSTAL_APPLY_C_NAMESPACE(info_getatominfo)
#ifdef ncrystal_info_getatompos
#  undef ncrystal_info_getatompos
#endif
#define ncrystal_info_getatompos NCRYSTAL_APPLY_C_NAMESPACE(info_getatompos)
#ifdef ncrystal_info_getcomponent
#  undef ncrystal_info_getcomponent
#endif
#define ncrystal_info_getcomponent NCRYSTAL_APPLY_C_NAMESPACE(info_getcomponent)
#ifdef ncrystal_info_getdebyetempbyelement
#  undef ncrystal_info_getdebyetempbyelement
#endif
#define ncrystal_info_getdebyetempbyelement NCRYSTAL_APPLY_C_NAMESPACE(info_getdebyetempbyelement)
#ifdef ncrystal_info_getdensity
#  undef ncrystal_info_getdensity
#endif
#define ncrystal_info_getdensity NCRYSTAL_APPLY_C_NAMESPACE(info_getdensity)
#ifdef ncrystal_info_getglobaldebyetemp
#  undef ncrystal_info_getglobaldebyetemp
#endif
#define ncrystal_info_getglobaldebyetemp NCRYSTAL_APPLY_C_NAMESPACE(info_getglobaldebyetemp)
#ifdef ncrystal_info_gethkl
#  undef ncrystal_info_gethkl
#endif
#define ncrystal_info_gethkl NCRYSTAL_APPLY_C_NAMESPACE(info_gethkl)
#ifdef ncrystal_info_gethkl_allindices
#  undef ncrystal_info_gethkl_allindices
#endif
#define ncrystal_info_gethkl_allindices NCRYSTAL_APPLY_C_NAMESPACE(info_gethkl_allindices)
#ifdef ncrystal_info_getnumberdensity
#  undef ncrystal_info_getnumberdensity
#endif
#define ncrystal_info_getnumberdensity NCRYSTAL_APPLY_C_NAMESPACE(info_getnumberdensity)
#ifdef ncrystal_info_getphase
#  undef ncrystal_info_getphase
#endif
#define ncrystal_info_getphase NCRYSTAL_APPLY_C_NAMESPACE(info_getphase)
#ifdef ncrystal_info_getsld
#  undef ncrystal_info_getsld
#endif
#define ncrystal_info_getsld NCRYSTAL_APPLY_C_NAMESPACE(info_getsld)
#ifdef ncrystal_info_getstateofmatter
#  undef ncrystal_info_getstateofmatter
#endif
#define ncrystal_info_getstateofmatter NCRYSTAL_APPLY_C_NAMESPACE(info_getstateofmatter)
#ifdef ncrystal_info_getstructure
#  undef ncrystal_info_getstructure
#endif
#define ncrystal_info_getstructure NCRYSTAL_APPLY_C_NAMESPACE(info_getstructure)
#ifdef ncrystal_info_gettemperature
#  undef ncrystal_info_gettemperature
#endif
#define ncrystal_info_gettemperature NCRYSTAL_APPLY_C_NAMESPACE(info_gettemperature)
#ifdef ncrystal_info_getxsectabsorption
#  undef ncrystal_info_getxsectabsorption
#endif
#define ncrystal_info_getxsectabsorption NCRYSTAL_APPLY_C_NAMESPACE(info_getxsectabsorption)
#ifdef ncrystal_info_getxsectfree
#  undef ncrystal_info_getxsectfree
#endif
#define ncrystal_info_getxsectfree NCRYSTAL_APPLY_C_NAMESPACE(info_getxsectfree)
#ifdef ncrystal_info_hasanydebyetemp
#  undef ncrystal_info_hasanydebyetemp
#endif
#define ncrystal_info_hasanydebyetemp NCRYSTAL_APPLY_C_NAMESPACE(info_hasanydebyetemp)
#ifdef ncrystal_info_hasatomdebyetemp
#  undef ncrystal_info_hasatomdebyetemp
#endif
#define ncrystal_info_hasatomdebyetemp NCRYSTAL_APPLY_C_NAMESPACE(info_hasatomdebyetemp)
#ifdef ncrystal_info_hasatommsd
#  undef ncrystal_info_hasatommsd
#endif
#define ncrystal_info_hasatommsd NCRYSTAL_APPLY_C_NAMESPACE(info_hasatommsd)
#ifdef ncrystal_info_hasatompos
#  undef ncrystal_info_hasatompos
#endif
#define ncrystal_info_hasatompos NCRYSTAL_APPLY_C_NAMESPACE(info_hasatompos)
#ifdef ncrystal_info_hasdebyetemp
#  undef ncrystal_info_hasdebyetemp
#endif
#define ncrystal_info_hasdebyetemp NCRYSTAL_APPLY_C_NAMESPACE(info_hasdebyetemp)
#ifdef ncrystal_info_hkl_dlower
#  undef ncrystal_info_hkl_dlower
#endif
#define ncrystal_info_hkl_dlower NCRYSTAL_APPLY_C_NAMESPACE(info_hkl_dlower)
#ifdef ncrystal_info_hkl_dupper
#  undef ncrystal_info_hkl_dupper
#endif
#define ncrystal_info_hkl_dupper NCRYSTAL_APPLY_C_NAMESPACE(info_hkl_dupper)
#ifdef ncrystal_info_hklinfotype
#  undef ncrystal_info_hklinfotype
#endif
#define ncrystal_info_hklinfotype NCRYSTAL_APPLY_C_NAMESPACE(info_hklinfotype)
#ifdef ncrystal_info_natominfo
#  undef ncrystal_info_natominfo
#endif
#define ncrystal_info_natominfo NCRYSTAL_APPLY_C_NAMESPACE(info_natominfo)
#ifdef ncrystal_info_ncomponents
#  undef ncrystal_info_ncomponents
#endif
#define ncrystal_info_ncomponents NCRYSTAL_APPLY_C_NAMESPACE(info_ncomponents)
#ifdef ncrystal_info_ncustomsections
#  undef ncrystal_info_ncustomsections
#endif
#define ncrystal_info_ncustomsections NCRYSTAL_APPLY_C_NAMESPACE(info_ncustomsections)
#ifdef ncrystal_info_ndyninfo
#  undef ncrystal_info_ndyninfo
#endif
#define ncrystal_info_ndyninfo NCRYSTAL_APPLY_C_NAMESPACE(info_ndyninfo)
#ifdef ncrystal_info_nhkl
#  undef ncrystal_info_nhkl
#endif
#define ncrystal_info_nhkl NCRYSTAL_APPLY_C_NAMESPACE(info_nhkl)
#ifdef ncrystal_info_nphases
#  undef ncrystal_info_nphases
#endif
#define ncrystal_info_nphases NCRYSTAL_APPLY_C_NAMESPACE(info_nphases)
#ifdef ncrystal_info_t
#  undef ncrystal_info_t
#endif
#define ncrystal_info_t NCRYSTAL_APPLY_C_NAMESPACE(info_t)
#ifdef ncrystal_info_uid
#  undef ncrystal_info_uid
#endif
#define ncrystal_info_uid NCRYSTAL_APPLY_C_NAMESPACE(info_uid)
#ifdef ncrystal_info_underlyinguid
#  undef ncrystal_info_underlyinguid
#endif
#define ncrystal_info_underlyinguid NCRYSTAL_APPLY_C_NAMESPACE(info_underlyinguid)
#ifdef ncrystal_invalidate
#  undef ncrystal_invalidate
#endif
#define ncrystal_invalidate NCRYSTAL_APPLY_C_NAMESPACE(invalidate)
#ifdef ncrystal_isnonoriented
#  undef ncrystal_isnonoriented
#endif
#define ncrystal_isnonoriented NCRYSTAL_APPLY_C_NAMESPACE(isnonoriented)
#ifdef ncrystal_lasterror
#  undef ncrystal_lasterror
#endif
#define ncrystal_lasterror NCRYSTAL_APPLY_C_NAMESPACE(lasterror)
#ifdef ncrystal_lasterrortype
#  undef ncrystal_lasterrortype
#endif
#define ncrystal_lasterrortype NCRYSTAL_APPLY_C_NAMESPACE(lasterrortype)
#ifdef ncrystal_msd2debyetemp
#  undef ncrystal_msd2debyetemp
#endif
#define ncrystal_msd2debyetemp NCRYSTAL_APPLY_C_NAMESPACE(msd2debyetemp)
#ifdef ncrystal_multicreate_direct
#  undef ncrystal_multicreate_direct
#endif
#define ncrystal_multicreate_direct NCRYSTAL_APPLY_C_NAMESPACE(multicreate_direct)
#ifdef ncrystal_name
#  undef ncrystal_name
#endif
#define ncrystal_name NCRYSTAL_APPLY_C_NAMESPACE(name)
#ifdef ncrystal_ncmat2json
#  undef ncrystal_ncmat2json
#endif
#ifdef ncrystal_namespace
#  undef ncrystal_namespace
#endif
#define ncrystal_namespace NCRYSTAL_APPLY_C_NAMESPACE(namespace)
#define ncrystal_ncmat2json NCRYSTAL_APPLY_C_NAMESPACE(ncmat2json)
#ifdef ncrystal_normalisecfg
#  undef ncrystal_normalisecfg
#endif
#define ncrystal_normalisecfg NCRYSTAL_APPLY_C_NAMESPACE(normalisecfg)
#ifdef ncrystal_process_t
#  undef ncrystal_process_t
#endif
#define ncrystal_process_t NCRYSTAL_APPLY_C_NAMESPACE(process_t)
#ifdef ncrystal_process_uid
#  undef ncrystal_process_uid
#endif
#define ncrystal_process_uid NCRYSTAL_APPLY_C_NAMESPACE(process_uid)
#ifdef ncrystal_raw_vdos2gn
#  undef ncrystal_raw_vdos2gn
#endif
#define ncrystal_raw_vdos2gn NCRYSTAL_APPLY_C_NAMESPACE(raw_vdos2gn)
#ifdef ncrystal_raw_vdos2knl
#  undef ncrystal_raw_vdos2knl
#endif
#define ncrystal_raw_vdos2knl NCRYSTAL_APPLY_C_NAMESPACE(raw_vdos2knl)
#ifdef ncrystal_raw_vdos2kernel
#  undef ncrystal_raw_vdos2kernel
#endif
#define ncrystal_raw_vdos2kernel NCRYSTAL_APPLY_C_NAMESPACE(raw_vdos2kernel)
#ifdef ncrystal_ref
#  undef ncrystal_ref
#endif
#define ncrystal_ref NCRYSTAL_APPLY_C_NAMESPACE(ref)
#ifdef ncrystal_refcount
#  undef ncrystal_refcount
#endif
#define ncrystal_refcount NCRYSTAL_APPLY_C_NAMESPACE(refcount)
#ifdef ncrystal_register_in_mem_file_data
#  undef ncrystal_register_in_mem_file_data
#endif
#define ncrystal_register_in_mem_file_data NCRYSTAL_APPLY_C_NAMESPACE(register_in_mem_file_data)
#ifdef ncrystal_remove_all_data_sources
#  undef ncrystal_remove_all_data_sources
#endif
#define ncrystal_remove_all_data_sources NCRYSTAL_APPLY_C_NAMESPACE(remove_all_data_sources)
#ifdef ncrystal_remove_custom_search_dirs
#  undef ncrystal_remove_custom_search_dirs
#endif
#define ncrystal_remove_custom_search_dirs NCRYSTAL_APPLY_C_NAMESPACE(remove_custom_search_dirs)
#ifdef ncrystal_rngsupportsstatemanip_ofscatter
#  undef ncrystal_rngsupportsstatemanip_ofscatter
#endif
#define ncrystal_rngsupportsstatemanip_ofscatter NCRYSTAL_APPLY_C_NAMESPACE(rngsupportsstatemanip_ofscatter)
#ifdef ncrystal_samplescatter
#  undef ncrystal_samplescatter
#endif
#define ncrystal_samplescatter NCRYSTAL_APPLY_C_NAMESPACE(samplescatter)
#ifdef ncrystal_samplescatter_many
#  undef ncrystal_samplescatter_many
#endif
#define ncrystal_samplescatter_many NCRYSTAL_APPLY_C_NAMESPACE(samplescatter_many)
#ifdef ncrystal_samplescatterisotropic
#  undef ncrystal_samplescatterisotropic
#endif
#define ncrystal_samplescatterisotropic NCRYSTAL_APPLY_C_NAMESPACE(samplescatterisotropic)
#ifdef ncrystal_samplescatterisotropic_many
#  undef ncrystal_samplescatterisotropic_many
#endif
#define ncrystal_samplescatterisotropic_many NCRYSTAL_APPLY_C_NAMESPACE(samplescatterisotropic_many)
#ifdef ncrystal_scatter_t
#  undef ncrystal_scatter_t
#endif
#define ncrystal_scatter_t NCRYSTAL_APPLY_C_NAMESPACE(scatter_t)
#ifdef ncrystal_setbuiltinrandgen
#  undef ncrystal_setbuiltinrandgen
#endif
#define ncrystal_setbuiltinrandgen NCRYSTAL_APPLY_C_NAMESPACE(setbuiltinrandgen)
#ifdef ncrystal_setbuiltinrandgen_withseed
#  undef ncrystal_setbuiltinrandgen_withseed
#endif
#define ncrystal_setbuiltinrandgen_withseed NCRYSTAL_APPLY_C_NAMESPACE(setbuiltinrandgen_withseed)
#ifdef ncrystal_setbuiltinrandgen_withstate
#  undef ncrystal_setbuiltinrandgen_withstate
#endif
#define ncrystal_setbuiltinrandgen_withstate NCRYSTAL_APPLY_C_NAMESPACE(setbuiltinrandgen_withstate)
#ifdef ncrystal_seterrhandler
#  undef ncrystal_seterrhandler
#endif
#define ncrystal_seterrhandler NCRYSTAL_APPLY_C_NAMESPACE(seterrhandler)
#ifdef ncrystal_sethaltonerror
#  undef ncrystal_sethaltonerror
#endif
#define ncrystal_sethaltonerror NCRYSTAL_APPLY_C_NAMESPACE(sethaltonerror)
#ifdef ncrystal_setquietonerror
#  undef ncrystal_setquietonerror
#endif
#define ncrystal_setquietonerror NCRYSTAL_APPLY_C_NAMESPACE(setquietonerror)
#ifdef ncrystal_setrandgen
#  undef ncrystal_setrandgen
#endif
#define ncrystal_setrandgen NCRYSTAL_APPLY_C_NAMESPACE(setrandgen)
#ifdef ncrystal_setrngstate_ofscatter
#  undef ncrystal_setrngstate_ofscatter
#endif
#define ncrystal_setrngstate_ofscatter NCRYSTAL_APPLY_C_NAMESPACE(setrngstate_ofscatter)
#ifdef ncrystal_unref
#  undef ncrystal_unref
#endif
#define ncrystal_unref NCRYSTAL_APPLY_C_NAMESPACE(unref)
#ifdef ncrystal_valid
#  undef ncrystal_valid
#endif
#define ncrystal_valid NCRYSTAL_APPLY_C_NAMESPACE(valid)
#ifdef ncrystal_vdoseval
#  undef ncrystal_vdoseval
#endif
#define ncrystal_vdoseval NCRYSTAL_APPLY_C_NAMESPACE(vdoseval)
#ifdef ncrystal_version
#  undef ncrystal_version
#endif
#define ncrystal_version NCRYSTAL_APPLY_C_NAMESPACE(version)
#ifdef ncrystal_version_str
#  undef ncrystal_version_str
#endif
#define ncrystal_version_str NCRYSTAL_APPLY_C_NAMESPACE(version_str)
#ifdef ncrystal_wl2ekin
#  undef ncrystal_wl2ekin
#endif
#define ncrystal_wl2ekin NCRYSTAL_APPLY_C_NAMESPACE(wl2ekin)
#ifdef ncrystal_benchloadcfg
#  undef ncrystal_benchloadcfg
#endif
#define ncrystal_benchloadcfg NCRYSTAL_APPLY_C_NAMESPACE(benchloadcfg)
#ifdef ncrystal_setmsghandler
#  undef ncrystal_setmsghandler
#endif
#define ncrystal_setmsghandler NCRYSTAL_APPLY_C_NAMESPACE(setmsghandler)
#ifdef ncrystal_runmmcsim_stdengine
#  undef ncrystal_runmmcsim_stdengine
#endif
#define ncrystal_runmmcsim_stdengine NCRYSTAL_APPLY_C_NAMESPACE(runmmcsim_stdengine)


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

  /* Expand VDOS to scattering kernel with Sjolander's method, directly from  */
  /* VDOS curves. The alpha, beta, and sab arrays should be freed with        */
  /* ncrystal_dealloc_doubleptr after usage. The order_weight_fct can         */
  /* optionally be used to change the contribution of a given Gn order to the */
  /* kernel (otherwise pass in a NULL ptr). The target_emax parameter will be */
  /* used to specifically request a target_emax value for the expansion (set  */
  /* to 0 to use a value based on the vdoslux parameter). The suggested_emax  */
  /* parameter will return the actual emax value supported by the returned    */
  /* kernel, unless an order_weight_fct is provided in which case 0 is always */
  /* returned. It might be lower than the target value, which can normally    */
  /* happen if the temperature is very low or the target emax was very high.  */
  NCRYSTAL_API void ncrystal_raw_vdos2kernel( const double* vdos_egrid,
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
                                              double** sab,
                                              double target_emax,
                                              double* suggested_emax );

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
  NCRYSTAL_API char * ncrystal_get_flatcompos( ncrystal_info_t, int prefernatelem,
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
  /*== Error and message handling                                               == */
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
  NCRYSTAL_API int ncrystal_error(void);/* returns 1 if an error condition occurred. */
  NCRYSTAL_API const char * ncrystal_lasterror(void);/* returns description of last error (NULL if none) */
  NCRYSTAL_API const char * ncrystal_lasterrortype(void);/* returns description of last error (NULL if none) */
  /* TODO: error file/line-no as well? */

  NCRYSTAL_API void ncrystal_clearerror(void);/* clears previous error if any */

  /*Another option is to provide a custom error handler which will be called on    */
  /*each error:                                                                    */
  NCRYSTAL_API void ncrystal_seterrhandler(void (*handler)(char*,char*));

  /* By default, NCrystal output is emitted on stdout. The following function can  */
  /* be used to instead provide a custom handling of these via a call-back         */
  /* function. The first argument is the message, and the second is the message    */
  /* class (0: info, 1: warning, 2: raw output):                                   */
  NCRYSTAL_API void ncrystal_setmsghandler(void (*handler)(const char*,unsigned));

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
  NCRYSTAL_API void ncrystal_setrandgen( double (*rg)(void) );

  /* It is also possible to (re) set the RNG to the builtin generator (optionally  */
  /* by state or integer seed) */
  NCRYSTAL_API void ncrystal_setbuiltinrandgen(void);
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
  /* call to ncrystal_dealloc_stringlist, and contains entries in the format       */
  /* name0,src0,fact0,priority0,name1,src1,fact1,priority1,..:                     */
  NCRYSTAL_API void ncrystal_get_file_list( unsigned* nstrs, char*** strs );

  /* Add directory to custom search path or clear custom search path:              */
  NCRYSTAL_API void ncrystal_add_custom_search_dir( const char * dir );
  NCRYSTAL_API void ncrystal_remove_custom_search_dirs(void);

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
  NCRYSTAL_API void ncrystal_remove_all_data_sources(void);


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

  /* Dump info to a string (must be cleaned up with ncrystal_dealloc_string):      */
  NCRYSTAL_API char* ncrystal_dump_tostr(ncrystal_info_t, unsigned verbosity_lvl );

  /* Utility converting between neutron wavelength [Aa] to kinetic energy [eV]:    */
  NCRYSTAL_API double ncrystal_wl2ekin( double wl );
  NCRYSTAL_API double ncrystal_ekin2wl( double ekin );

  /* Enable factory thread pool. Assuming that user code runs in single thread     */
  /* (at least while initialising materials), this requested value is the TOTAL    */
  /* number of threads utilised INCLUDING that user thread. Thus, a value of 0     */
  /* or 1 will disable this thread pool, while a value of 8 will result in 7       */
  /* secondary worker threads being allocated.  Supply a value >= 9999 to simply   */
  /* use a number of threads appropriate for the system:                           */
  NCRYSTAL_API void ncrystal_enable_factory_threadpool( unsigned );

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
  NCRYSTAL_API unsigned ncrystal_atomdatadb_getnentries(void);
  NCRYSTAL_API void ncrystal_atomdatadb_getallentries( unsigned* zvals,
                                                       unsigned* avals );

  /* Convert between atomic mean squared displacements and Debye temperatures.     */
  /* Units are kelvin, AMU, and Angstrom^2:                                        */
  NCRYSTAL_API double ncrystal_debyetemp2msd( double debyetemp, double temperature, double mass );
  NCRYSTAL_API double ncrystal_msd2debyetemp( double msd, double temperature, double mass );

  /* Extract information from VDOS curve (see NCVDOSEval.hh for details):          */
  NCRYSTAL_API void ncrystal_vdoseval( double vdos_emin, double vdos_emax,
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
  NCRYSTAL_API void ncrystal_clear_caches(void);

  /* Get list of plugins. Resulting string list must be deallocated by a call to   */
  /* ncrystal_dealloc_stringlist by, and contains entries in the format            */
  /* pluginname0,filename0,plugintype0,pluginname1,filename1,plugintype1,...:      */
  NCRYSTAL_API void ncrystal_get_plugin_list( unsigned* nstrs, char*** strs );

  /* Deallocate strings / double arrays:                                           */
  NCRYSTAL_API void ncrystal_dealloc_stringlist( unsigned len, char** );
  NCRYSTAL_API void ncrystal_dealloc_string( char* );
  NCRYSTAL_API void ncrystal_dealloc_doubleptr( double* );

  /* NCrystal version info:                                                             */
  /* Note that through the ncapi.h include above you also have access to macros:        */
  /* NCRYSTAL_VERSION_MAJOR (integer)                                                   */
  /* NCRYSTAL_VERSION_MINOR (integer)                                                   */
  /* NCRYSTAL_VERSION_PATCH (integer)                                                   */
  /* NCRYSTAL_VERSION (integer == 1000000*MAJOR+1000*MINOR+PATCH)                       */
  /* NCRYSTAL_VERSION_STRING (string, "MAJOR.MINOR.PATCH" )                             */
  NCRYSTAL_API int ncrystal_version(void); /* returns NCRYSTAL_VERSION                  */
  NCRYSTAL_API const char * ncrystal_version_str(void); /* returns NCRYSTAL_VERSION_STR */

  /* If compiled with NCRYSTAL_NAMESPACE_PROTECTION, return the namespace here:    */
  /* (will be an empty string in default installations):                           */
  NCRYSTAL_API const char * ncrystal_namespace(void);

  /* Load raw NCMAT data into JSON structures. Must deallocate with call to        */
  /* ncrystal_dealloc_string as usual. (WARNING: JSON is incomplete for now!!!!!)  */
  NCRYSTAL_API char * ncrystal_ncmat2json( const char * textdataname );

  /* Get time in seconds to load the cfg in question (if do_scatter=0 it will only */
  /* create Info objects). Caches are cleared as a side effect: */
  NCRYSTAL_API double ncrystal_benchloadcfg( const char * cfgstr, int do_scat, int repeat );


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
  NCRYSTAL_API void ncrystal_clear_info_caches(void);

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


  /* Same as ncrystal_raw_vdos2knl but without target_emax/suggested_emax: */
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

  /* Run the experimental embedded simulation engine for diffraction patterns. */
  /* Depending on the tally_detail_lvl, various results are returned in the    */
  /* output variables:                                                         */
  /*     0 : only exit angle hist (contents + errors, no tally_json)           */
  /*     1 : also running stats of exit angle dist in tally_json.              */
  /*     2 : many other hists as well (completely contained in tally_json).    */
  /* tally_json must be deallocated with ncrystal_dealloc_string if not NULL.  */
  /* tally_exitangle_contents and _errsq must be passed to                     */
  /* ncrystal_dealloc_doubleptr after usage. Set nthreads>=9999 for            */
  /* auto-detection of a suitable number of threads:                           */
  NCRYSTAL_API void ncrystal_runmmcsim_stdengine( unsigned nthreads,
                                                  unsigned tally_detail_lvl,
                                                  const char * mat_cfgstr,
                                                  const char * mmc_geomcfg,
                                                  const char * mmc_srccfg,
                                                  char ** tally_json,
                                                  unsigned * tally_exitangle_nbins,
                                                  double ** tally_exitangle_contents,
                                                  double ** tally_exitangle_errsq );

#ifdef __cplusplus
}
#endif

#endif
