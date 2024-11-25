
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

//Note about this template: At CMake configuration time, the AT-SIGN enclosed
//CMake variables are expanded. Then at CMake build time, the Generator
//expressions are expanded.

const char * nccfg_const_bin2libdir() { return "@NCrystal_relpath_BINDIR2LIBDIR@"; }
const char * nccfg_const_bin2shlibdir() { return "@NCrystal_relpath_BINDIR2SHLIBDIR@"; }
const char * nccfg_const_libname() { return "@NCrystal_libname_genexp@"; }//NB: generator expression
const char * nccfg_const_shlibname() { return "@NCrystal_shlibname_genexp@"; }//NB: generator expression
const char * nccfg_const_bin2libpath() { return "@NCrystal_relpath_BINDIR2LIBDIR@/@NCrystal_libname_genexp@"; }
const char * nccfg_const_bin2shlibpath() { return "@NCrystal_relpath_BINDIR2SHLIBDIR@/@NCrystal_shlibname_genexp@"; }
const char * nccfg_const_bin2datadir() { return "@NCrystal_relpath_BINDIR2DATADIR@"; }
const char * nccfg_const_bin2incdir() { return "@NCrystal_relpath_BINDIR2INCDIR@"; }
const char * nccfg_const_bin2cmakedir() { return "@NCrystal_relpath_BINDIR2CMAKEDIR@"; }
const char * nccfg_const_version() { return "@NCrystal_VERSION@"; }
const char * nccfg_const_intversion() { return "@nccfgapp_intversion@"; }
const char * nccfg_const_builtinplugins() { return "@NCrystal_builtin_plugin_names@"; }
const char * nccfg_const_namespace() { return "@NCRYSTAL_NAMESPACE@"; }
const char * nccfg_const_cmakebuildtype() { return "$<CONFIG>"; }//NB: generator expression
int nccfg_boolopt_data() { return @nccfgapp_has_data_01@; }
int nccfg_boolopt_dynamic_plugins() { return @nccfgapp_has_dynload_01@; }
int nccfg_boolopt_embed_data() { return @nccfgapp_has_dataembed_01@; }
int nccfg_boolopt_examples() { return @nccfgapp_has_examples_01@; }
int nccfg_boolopt_modify_rpath() { return @nccfgapp_has_modrpath_01@; }
int nccfg_boolopt_threads() { return @nccfgapp_has_threads_01@; }
