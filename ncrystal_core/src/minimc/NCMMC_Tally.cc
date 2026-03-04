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

#include "NCrystal/internal/minimc/NCMMC_Tally.hh"
namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

NCMMC::TallyMgr::TallyMgr( TallyPtr tally_template )
  : m_template(std::move(tally_template))
{
}

NCMMC::TallyPtr NCMMC::TallyMgr::getIndependentTallyPtr() const
{
  return m_template->cloneSetup();
}

void NCMMC::TallyMgr::addResult( TallyPtr res_SO )
{
  NCRYSTAL_DEBUGMMCMSG("TallyMgr::addResult");
  std::shared_ptr<TallyBase> res = std::move(res_SO);
  std::shared_ptr<TallyBase> to_merge;
  {
    NCRYSTAL_LOCK_GUARD(m_final_mutex);//hold it only briefly!
    if ( m_final == nullptr ) {
      //first time, or someone else is currently merging.
      m_final = std::move(res);
      return;
    } else {
      to_merge.swap(m_final);
    }
  }
  if ( to_merge != nullptr ) {
    //perform the merging in the current thread without holding a lock,
    //and then put the merged result back:
    to_merge->merge( std::move( *res.get() ) );
    //Put result back:
    res = nullptr;//to be safe (not strictly needed)
    this->addResult( std::move(to_merge) );
  }
}

NCMMC::TallyPtr NCMMC::TallyMgr::getFinalResult()
{
  NCRYSTAL_LOCK_GUARD(m_final_mutex);//should not really be needed if used correctly
  nc_assert_always(m_final!=nullptr);
  return std::move(m_final);
}
