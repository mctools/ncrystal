#ifndef G4NCrystal_Manager_hh
#define G4NCrystal_Manager_hh

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

#include <map>
#include <vector>
#include "G4Material.hh"
#include "G4Version.hh"
#include "G4MaterialPropertiesTable.hh"
#include "NCrystal/NCProcImpl.hh"
//Manager class tracking indices of NCrystal::Scatter instances associated to
//G4Materials, via entries in the G4MaterialPropertiesTable's on the materials.

class G4Material;
namespace NCrystal {
  namespace Error {
    class Exception;
  }
}

namespace G4NCrystal {

  class Manager {
  public:

    //Methods needed to add NCrystal::Scatter* properties to G4Materials (via
    //the property tables):
    static Manager * getInstance();//Get the singleton
    void addScatterProperty(G4Material*, NCrystal::ProcImpl::ProcPtr&&);

    //Methods for framework implementers:
    const NCrystal::ProcImpl::Process* getScatterProperty(G4Material*) const;//returns nullptr when absent.
    NCrystal::ProcImpl::OptionalProcPtr getScatterPropertyPtr(G4Material*) const;//returns nullptr when absent.

    //Same but with per-thread per-scatter CachePtr:
    using ProcAndCache = std::pair<const NCrystal::ProcImpl::Process*,NCrystal::CachePtr*>;
    ProcAndCache getScatterPropertyWithThreadSafeCache(G4Material*) const;

    //Thoroughly clear caches, manager singleton, and possibly NCrystal
    //factories. It is NOT safe to use the Scatter properties of already created
    //G4Materials after this.
    static void cleanup();

    //Safer cleanup which should only remove objects with no current usage:
    static void clearCaches();

    unsigned nMaterialsWithProperties() const { return m_scatters.size(); }

    //Translate thrown NCrystal exceptions to G4Exception(..) calls (id should
    //be unique and fixed for each call location):
    static void handleError( const char*origin, unsigned id,
                             NCrystal::Error::Exception& );

  private:
    Manager( const Manager & );
    Manager & operator= ( const Manager & );
    Manager();
    ~Manager();
    static Manager * s_mgr;
    std::vector<NCrystal::ProcImpl::ProcPtr> m_scatters;
    std::map<uint64_t,unsigned> m_scat2idx;
    G4String m_key;
    NCrystal::CachePtr& getCachePtrForCurrentThreadAndProcess( unsigned scatter_idx ) const;
    //Returns numeric_limits<unsigned>::max() if not available:
    unsigned lookupScatterPropertyIndex(G4Material*) const;

  };

  ///////////////////////////////////////////////
  // Inline for fast access during event loop: //
  ///////////////////////////////////////////////

  inline std::pair<const NCrystal::ProcImpl::Process*, NCrystal::CachePtr*>
  Manager::getScatterPropertyWithThreadSafeCache(G4Material* mat) const
  {
    //Returns numeric_limits<unsigned>::max() if not available:
    unsigned scatidx = lookupScatterPropertyIndex(mat);
    if ( scatidx == std::numeric_limits<unsigned>::max() )
      return {nullptr,nullptr};
    assert(scatidx<m_scatters.size());
    const NCrystal::ProcImpl::Process* sp = m_scatters[scatidx].get();
    return { sp, &getCachePtrForCurrentThreadAndProcess( scatidx ) };
  }

  inline unsigned Manager::lookupScatterPropertyIndex(G4Material*mat) const
  {
    G4MaterialPropertiesTable* matprop = mat->GetMaterialPropertiesTable();
    constexpr unsigned not_found = std::numeric_limits<unsigned>::max();
    if (!matprop)
      return not_found;
#if G4VERSION_NUMBER < 1040
    //Property maps pre-Geant4 10.4
    using PropMap_T = std::map< G4String, G4double, std::less<G4String> > ;
    const PropMap_T* propcmap = matprop->GetPropertiesCMap();
    PropMap_T::size_type size = propcmap->size();
    //Access property like this instead of via
    //ConstPropertyExists+GetConstProperty, to avoid needless string allocations
    //and double map lookup:
    if (size<=1) {
      //Optimise to quickly check in small maps, since we are likely the only
      //user of constant material properties (this check can be removed if the
      //situation changes drastically):
      PropMap_T::const_iterator it;
      return (size && (it=propcmap->begin())->first == m_key) ? static_cast<unsigned>(it->second) : not_found;
    } else {
      //Full map search:
      PropMap_T::const_iterator it = propcmap->find(m_key);
      return it==propcmap->end() ? not_found : static_cast<unsigned>(it->second);
    }
#else
    //Property maps post Geant4 10.4. Here the internal keys are integers, but
    //we can't be sure it is always the same in different materials. Unfortunately
    //GetConstProperty triggers a G4 exception in case the key is not present,
    //so for safety we have to trigger one extra map search via ConstPropertyExists:

#  if G4VERSION_NUMBER >= 1100
    if ( matprop->GetMaterialConstPropertyNames().empty() || !matprop->ConstPropertyExists(m_key) )
      return not_found;
#  else
    if ( matprop->GetConstPropertyMap()->empty() || !matprop->ConstPropertyExists(m_key) )
      return not_found;
#endif
    return static_cast<unsigned>(matprop->GetConstProperty(m_key));
#endif
  }

  inline const NCrystal::ProcImpl::Process* Manager::getScatterProperty(G4Material*mat) const
  {
    //Returns numeric_limits<unsigned>::max() if not available:
    unsigned scatidx = lookupScatterPropertyIndex(mat);
    if ( scatidx == std::numeric_limits<unsigned>::max() )
      return nullptr;
    return m_scatters.at(scatidx).get();
  }

}

#endif
