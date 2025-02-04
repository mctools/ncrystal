#ifndef NCrystal_hh
#define NCrystal_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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

////////////////////////////////////////////////////
// Convenience header for all NCrystal interfaces //
////////////////////////////////////////////////////

#ifndef NCrystal_Dump_hh
#  include "NCrystal/dump/NCDump.hh"
#endif
#ifndef NCrystal_Defs_hh
#  include "NCrystal/core/NCDefs.hh"
#endif
#ifndef NCrystal_Exception_hh
#  include "NCrystal/core/NCException.hh"
#endif
#ifndef NCrystal_Mem_hh
#  include "NCrystal/core/NCMem.hh"
#endif
#ifndef NCrystal_SmallVector_hh
#  include "NCrystal/core/NCSmallVector.hh"
#endif
#ifndef NCrystal_Variant_hh
#  include "NCrystal/core/NCVariant.hh"
#endif
#ifndef NCrystal_Fmt_hh
#  include "NCrystal/core/NCFmt.hh"
#endif
#ifndef NCrystal_PluginMgmt_hh
#  include "NCrystal/plugins/NCPluginMgmt.hh"
#endif
#ifndef NCrystal_AtomData_hh
#  include "NCrystal/interfaces/NCAtomData.hh"
#endif
#ifndef NCrystal_TextData_hh
#  include "NCrystal/text/NCTextData.hh"
#endif
#ifndef NCrystal_Info_hh
#  include "NCrystal/interfaces/NCInfo.hh"
#endif
#ifndef NCrystal_InfoTypes_hh
#  include "NCrystal/interfaces/NCInfoTypes.hh"
#endif
#ifndef NCrystal_SABData_hh
#  include "NCrystal/interfaces/NCSABData.hh"
#endif
#ifndef NCrystal_MatCfg_hh
#  include "NCrystal/factories/NCMatCfg.hh"
#endif
#ifndef NCrystal_RNG_hh
#  include "NCrystal/interfaces/NCRNG.hh"
#endif
#ifndef NCrystal_SCOrientation_hh
#  include "NCrystal/interfaces/NCSCOrientation.hh"
#endif
#ifndef NCrystal_Version_hh
#  include "NCrystal/interfaces/NCVersion.hh"
#endif
#ifndef NCrystal_CompositionUtils_hh
#  include "NCrystal/misc/NCCompositionUtils.hh"
#endif
#ifndef NCrystal_ProcImpl_hh
#  include "NCrystal/interfaces/NCProcImpl.hh"
#endif
#ifndef NCrystal_Proc_hh
#  include "NCrystal/interfaces/NCProc.hh"
#endif
#ifndef NCrystal_MsgCtrl_hh
#  include "NCrystal/misc/NCMsgCtrl.hh"
#endif
#ifndef NCrystal_FactTypes_hh
#  include "NCrystal/factories/NCFactTypes.hh"
#endif
#ifndef NCrystal_FactImpl_hh
#  include "NCrystal/factories/NCFactImpl.hh"
#endif
#ifndef NCrystal_FactRequests_hh
#  include "NCrystal/factories/NCFactRequests.hh"
#endif
#ifndef NCrystal_FactRequestsImpl_hh
#  include "NCrystal/factories/NCFactRequestsImpl.hh"
#endif
#ifndef NCrystal_DataSources_hh
#  include "NCrystal/factories/NCDataSources.hh"
#endif
#ifndef NCrystal_Fact_hh
#  include "NCrystal/factories/NCFact.hh"
#endif
#ifndef NCrystal_Types_hh
#  include "NCrystal/core/NCTypes.hh"
#endif
#ifndef NCrystal_ImmutBuf_hh
#  include "NCrystal/core/NCImmutBuf.hh"
#endif
#ifndef NCrystal_FactThreads_hh
#  include "NCrystal/threads/NCFactThreads.hh"
#endif

#endif
