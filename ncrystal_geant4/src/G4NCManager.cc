////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "G4NCrystal/G4NCManager.hh"

#include "G4Material.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <sstream>
#include <iomanip>
#include <cassert>

namespace NC = NCrystal;
namespace NCG4 = G4NCrystal;

NCG4::Manager::Manager()
  : m_key("NCScat")
{}

NCG4::Manager::~Manager() = default;

NCG4::Manager * NCG4::Manager::s_mgr = 0;

NCG4::Manager * NCG4::Manager::getInstance()
{
  //TODO: Clearly this is not MT safe.
  if (!s_mgr)
    s_mgr = new Manager;
  return s_mgr;
}

NC::CachePtr& NCG4::Manager::getCachePtrForCurrentThreadAndProcess( unsigned scatter_idx ) const {
  static
#ifdef G4MULTITHREADED //protect thread_local keyword to avoid potential headaches in ST builds
 thread_local
#endif
    std::unique_ptr<std::vector<NCrystal::CachePtr>> cacheptrs;
  if ( cacheptrs == nullptr ) {
    cacheptrs = std::make_unique<decltype(cacheptrs)::element_type>();
    cacheptrs->resize( m_scatters.size() );
  }
  assert( scatter_idx < cacheptrs->size() );
  return (*cacheptrs)[scatter_idx];
}

NCrystal::ProcImpl::OptionalProcPtr NCG4::Manager::getScatterPropertyPtr(G4Material*mat) const
{
  //Returns numeric_limits<unsigned>::max() if not available:
  unsigned scatidx = lookupScatterPropertyIndex(mat);
  if ( scatidx == std::numeric_limits<unsigned>::max() )
    return nullptr;
  return m_scatters.at(scatidx);
}

void NCG4::Manager::addScatterProperty(G4Material* mat,NCrystal::ProcImpl::ProcPtr&&scat)
{
  if ( mat == nullptr || scat == nullptr )
    G4Exception ("NCG4::Manager::addScatterProperty", "NCAddingNull",
                 JustWarning, "Got nullptr argument.");

  if ( scat->processType() != NC::ProcessType::Scatter )
    G4Exception ("NCG4::Manager::addScatterProperty", "NCAddNonScatter",
                 JustWarning, "Can only add scattering process (processType() is not NCrystal::Process::Scatter).");

  auto matprop = mat->GetMaterialPropertiesTable();
  if (matprop) {
    if (matprop->ConstPropertyExists(m_key.c_str())) {
      std::stringstream ss;
      ss<<"Setting property on material "<<mat->GetName()<<" more than once overrides previous content!";
      G4Exception ("NCG4::Manager::addScatterProperty", "NCAddTwice",
                   JustWarning, ss.str().c_str());
    }
  } else {
    matprop = new G4MaterialPropertiesTable();
    mat->SetMaterialPropertiesTable(matprop);
  }

  unsigned idx = std::numeric_limits<unsigned>::max();
  uint64_t scatuid = scat->getUniqueID().value;
  auto it = m_scat2idx.find(scatuid);
  if ( it == m_scat2idx.end() ) {
    idx = m_scatters.size();
    m_scat2idx[scatuid] = idx;
    m_scatters.push_back(std::move(scat));
  } else {
    //already known:
    idx = it->second;
  }
  assert( unsigned(double(idx)) == idx );//make sure we can get the idx back out
  matprop->AddConstProperty(m_key.c_str(), idx
#if G4VERSION_NUMBER >= 1100
                            , true//createNewKey=true is required explicitly since G4 11.0
#endif
                            );
}

void NCG4::Manager::cleanup()
{
  if (s_mgr) {
    delete s_mgr;
    s_mgr=0;
  }
  NC::clearCaches();
}

void NCG4::Manager::clearCaches()
{
  NC::clearCaches();
}

void NCG4::Manager::handleError(const char* origin, unsigned id, NC::Error::Exception& e)
{
  std::ostringstream s_code;
  s_code << "G4NCrystal::"<<e.getTypeName()<<std::setfill('0') << std::setw(3)<<id;
  G4Exception(origin, s_code.str().c_str(),FatalException, e.what());
}
