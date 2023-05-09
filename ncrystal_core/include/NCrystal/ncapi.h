#ifndef ncrystal_api_h
#define ncrystal_api_h

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

/********************************************************************************/
/* API macros for NCrystal, used to mark classes and functions exported to API  */
/* clients. This is needed to support Windows DLL builds as well as Unix builds */
/* with -fvisibility-hidden. Macros inspired by suggestions found at            */
/* https://gcc.gnu.org/wiki/Visibility                                          */
/*                                                                              */
/* On windows, the NCrystal_EXPORTS macro must be defined during compilation of */
/* the NCrystal library (dll), but not when client applications and libraries   */
/* are being compiled. If using CMake and the CMakeLists.txt distributed with   */
/* NCrystal, this should be taken care of automatically.                        */
/********************************************************************************/

#ifdef NCRYSTAL_API
#  undef NCRYSTAL_API
#endif
#ifdef NCRYSTAL_LOCAL
#  undef NCRYSTAL_LOCAL
#endif
#if defined (_WIN32) || defined (__CYGWIN__) || defined (WIN32)
#  ifdef NCrystal_EXPORTS
#    define NCRYSTAL_API __declspec(dllexport)
#  else
#    define NCRYSTAL_API __declspec(dllimport)
#  endif
#  define NCRYSTAL_LOCAL
#else
#  if ( defined(__GNUC__) && __GNUC__ >= 4) || defined(__clang__)
#    define NCRYSTAL_API    __attribute__ ((visibility ("default")))
#    define NCRYSTAL_LOCAL  __attribute__ ((visibility ("hidden")))
#  else
#    define NCRYSTAL_API
#    define NCRYSTAL_LOCAL
#  endif
#endif

#ifdef __cplusplus
/* For decorating with constexpr only in C++17 and later, or noexcept only in */
/* non-dbg builds:                                                            */
#  ifdef ncconstexpr17
#    undef ncconstexpr17
#  endif
#  if __cplusplus >= 201703L
#    define ncconstexpr17 constexpr
#  else
#    define ncconstexpr17
#  endif
#  ifdef ncnoexceptndebug
#    undef ncnoexceptndebug
#  endif
#  ifndef NDEBUG
#    define ncnoexceptndebug
#  else
#    define ncnoexceptndebug noexcept
#  endif
#  ifdef ncnodiscard17
#    undef ncnodiscard17
#  endif
#  if __cplusplus >= 201703L
#    define ncnodiscard17 [[nodiscard]]
#  else
#    define ncnodiscard17
#  endif
#  ifdef nclikely
#    undef nclikely
#  endif
#  if __cplusplus >= 202002L
#    define nclikely [[likely]]
#  else
#    define nclikely
#  endif
#  ifdef ncunlikely
#    undef ncunlikely
#  endif
#  if __cplusplus >= 202002L
#    define ncunlikely [[unlikely]]
#  else
#    define ncunlikely
#  endif
#endif

#endif
