////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

//No-op unless built as the special NCRYSTAL_WIN_FMA-flagged /arch:AVX2 twin
//of NCFastConvolve.cc (MSVC/x86 only, per the NC*_WINFMA.cc convention in
//ncrystal_core/CMakeLists.txt) -- see
//NCrystal/internal/utils/NCFMADispatch.hh for why this file exists and how
//it fits together with NCFastConvolve.cc/_FMA.hh.
#ifdef NCRYSTAL_WIN_FMA
#include "NCFastConvolve_FMA.hh"
#endif
