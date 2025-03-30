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

#include "NCrystal/internal/extd_utils/NCProcCompBldr.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include <list>

namespace NC = NCrystal;

struct NC::Utils::ProcCompBldr::Impl final : NoCopyMove {
  //We use data structures here that do not invalidate existing entries when new
  //ones are appended. That way, each thread can get a pointer to its entry, and
  //write to it freely without any locking mechanism required. For efficiency,
  //we use a combination of an array for usual cases, and a std::list for
  //overflow (std::list's do not invalidate existing iterators when new entries
  //are appended).
  FactoryJobs m_jobMgr;
  //Fixed output buffer:
  std::mutex m_mutex;
  static constexpr unsigned comps_buf_maxsize = 32;
  unsigned m_comps_buf_size = 0;
  Optional<ComponentList> m_comps[comps_buf_maxsize];
  //Overflow buffer:
  std::list<Optional<ComponentList>> m_extra_comps;
};


NC::Utils::ProcCompBldr::ProcCompBldr() = default;
NC::Utils::ProcCompBldr::~ProcCompBldr() = default;

void NC::Utils::ProcCompBldr::addfct( std::function<ProcPtr()> fct )
{
  this->addfct_cl( [fct]()
  {
    ProcImpl::ProcComposition::ComponentList complist;
    complist.emplace_back(fct());
    return complist;
  });
}

void NC::Utils::ProcCompBldr::add( ProcPtr procptr, double scale )
{
  ComponentList cl;
  cl.emplace_back( scale, std::move(procptr) );
  this->add_cl( std::move(cl) );
}

void NC::Utils::ProcCompBldr::addfct_cl( std::function<ComponentList()> fct )
{
  NCRYSTAL_LOCK_GUARD(m_impl->m_mutex);
  unsigned buf_idx = m_impl->m_comps_buf_size++;
  Optional<ComponentList>* bufptr;
  if ( buf_idx < Impl::comps_buf_maxsize ) {
    bufptr = &m_impl->m_comps[buf_idx];
  } else {
    //Overflow:
    bufptr = &(*m_impl->m_extra_comps.emplace(m_impl->m_extra_comps.end()));
  }
  m_impl->m_jobMgr.queue( [fct,bufptr](){
    ComponentList res = fct();
    *bufptr = std::move(res);
  } );
}

void NC::Utils::ProcCompBldr::add_cl( ComponentList cl )
{
  NCRYSTAL_LOCK_GUARD(m_impl->m_mutex);
  unsigned buf_idx = m_impl->m_comps_buf_size++;
  if ( buf_idx < Impl::comps_buf_maxsize ) {
    m_impl->m_comps[buf_idx] = std::move(cl);
  } else {
    //Overflow:
    m_impl->m_extra_comps.emplace_back( std::move(cl) );
  }
}

NC::Utils::ProcCompBldr::ComponentList NC::Utils::ProcCompBldr::finalise()
{
  NCRYSTAL_LOCK_GUARD(m_impl->m_mutex);
  m_impl->m_jobMgr.waitAll();
  ComponentList res;
  auto transferEntries = [&res]( ComponentList& cl )
  {
    for ( auto& e : cl )
      if ( ( e.scale>0.0 && !e.process->isNull() ) )
        res.emplace_back(std::move(e));
  };

  const char * errmsg = ( "ProcCompBldr did not receive expected"
                          " component list from job" );
  unsigned stdsize = ( m_impl->m_comps_buf_size > Impl::comps_buf_maxsize
                       ? Impl::comps_buf_maxsize : m_impl->m_comps_buf_size );

  for ( unsigned i = 0; i < stdsize; ++i ) {
    if ( !m_impl->m_comps[i].has_value() )
      NCRYSTAL_THROW(LogicError,errmsg);
    nc_assert_always( m_impl->m_comps[i].has_value() );
    transferEntries( m_impl->m_comps[i].value() );
  }
  //Overflow:
  for ( auto& cl : m_impl->m_extra_comps ) {
    if (!cl.has_value())
      NCRYSTAL_THROW(LogicError,errmsg);
    transferEntries( cl.value() );
  }
  m_impl->m_comps_buf_size = 0;
  m_impl->m_extra_comps.clear();
  return res;
}

NC::ProcImpl::ProcPtr NC::Utils::ProcCompBldr::finalise_scatter()
{
  return ProcImpl::ProcComposition::consumeAndCombine( this->finalise(),
                                                       ProcessType::Scatter );
}

NC::ProcImpl::ProcPtr NC::Utils::ProcCompBldr::finalise_absorption()
{
  return ProcImpl::ProcComposition::consumeAndCombine( this->finalise(),
                                                       ProcessType::Absorption );
}
