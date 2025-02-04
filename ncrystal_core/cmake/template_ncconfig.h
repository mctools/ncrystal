
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

//Note about this template: At CMake configuration time, the AT-SIGN enclosed
//CMake variables are expanded. Then at CMake build time, the Generator
//expressions are expanded.

const char * nccfg_const_bin2libdir(void) { return "@NCrystal_relpath_BINDIR2LIBDIR@"; }
const char * nccfg_const_bin2shlibdir(void) { return "@NCrystal_relpath_BINDIR2SHLIBDIR@"; }
const char * nccfg_const_libname(void) { return "@NCrystal_libname_genexp@"; }//NB: generator expression
const char * nccfg_const_shlibname(void) { return "@NCrystal_shlibname_genexp@"; }//NB: generator expression
const char * nccfg_const_bin2libpath(void) { return "@NCrystal_relpath_BINDIR2LIBDIR@/@NCrystal_libname_genexp@"; }
const char * nccfg_const_bin2shlibpath(void) { return "@NCrystal_relpath_BINDIR2SHLIBDIR@/@NCrystal_shlibname_genexp@"; }
const char * nccfg_const_bin2datadir(void) { return "@NCrystal_relpath_BINDIR2DATADIR@"; }
const char * nccfg_const_bin2incdir(void) { return "@NCrystal_relpath_BINDIR2INCDIR@"; }
const char * nccfg_const_bin2cmakedir(void) { return "@NCrystal_relpath_BINDIR2CMAKEDIR@"; }
const char * nccfg_const_version(void) { return "@NCrystal_VERSION@"; }
const char * nccfg_const_intversion(void) { return "@nccfgapp_intversion@"; }
const char * nccfg_const_namespace(void) { return "@NCRYSTAL_NAMESPACE@"; }
const char * nccfg_const_cmakebuildtype(void) { return "$<CONFIG>"; }//NB: generator expression
int nccfg_boolopt_data(void) { return @nccfgapp_has_data_01@; }
int nccfg_boolopt_dynamic_plugins(void) { return @nccfgapp_has_dynload_01@; }
int nccfg_boolopt_embed_data(void) { return @nccfgapp_has_dataembed_01@; }
int nccfg_boolopt_examples(void) { return @nccfgapp_has_examples_01@; }
int nccfg_boolopt_modify_rpath(void) { return @nccfgapp_has_modrpath_01@; }
int nccfg_boolopt_threads(void) { return @nccfgapp_has_threads_01@; }
int nccfg_boolopt_expects_shlibdir_override(void) { return @nccfgapp_expect_shlibdir_override_01@; }
