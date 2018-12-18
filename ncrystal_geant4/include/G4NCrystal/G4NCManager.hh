#ifndef G4NCrystal_Manager_hh
#define G4NCrystal_Manager_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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

//Manager class tracking indices of NCrystal::Scatter instances associated to
//G4Materials, via entries in the G4MaterialPropertiesTable's on the materials.

class G4Material;
namespace NCrystal {
  namespace Error {
    class Exception;
  }
  class Scatter;
}

namespace G4NCrystal {

  class Manager {
  public:

    //Methods needed to add NCrystal::Scatter* properties to G4Materials (via
    //the property tables):
    static Manager * getInstance();//Get the singleton
    void addScatterProperty(G4Material*,const NCrystal::Scatter*);

    //Methods for framework implementers:
    const NCrystal::Scatter* getScatterProperty(G4Material*);//returns 0 when absent.
    static Manager * getInstanceNoInit();//get singleton if created, else 0
    static void cleanup();//delete singleton, unref kept scatter instances (for valgrind).
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
    std::vector<const NCrystal::Scatter*> m_scatters;
    std::map<const NCrystal::Scatter*,unsigned> m_scat2idx;
    G4String m_key;
#if G4VERSION_NUMBER < 1040
    typedef std::map< G4String, G4double, std::less<G4String> > PropMap_T;
#endif
  };

  ///////////////////////////////////////////////
  // Inline for fast access during event loop: //
  ///////////////////////////////////////////////

  inline const NCrystal::Scatter* Manager::getScatterProperty(G4Material*mat)
  {
    G4MaterialPropertiesTable* matprop = mat->GetMaterialPropertiesTable();
    if (!matprop)
      return 0;
    //Access property like this instead of via
    //ConstPropertyExists+GetConstProperty, to avoid needless string allocations
    //and double map lookup:
#if G4VERSION_NUMBER < 1040
    //Property maps pre-Geant4 10.4
    const PropMap_T* propcmap = matprop->GetPropertiesCMap();
    PropMap_T::size_type size = propcmap->size();
    if (size<=1) {
      //Optimise to quickly check in small maps, since we are likely the only
      //user of constant material properties (this check can be removed if the
      //situation changes drastically):
      PropMap_T::const_iterator it;
      return (size && (it=propcmap->begin())->first == m_key) ?
        m_scatters.at(unsigned(it->second)) : 0;
    } else {
      //Full map search:
      PropMap_T::const_iterator it = propcmap->find(m_key);
      return it==propcmap->end() ? 0 : m_scatters.at(unsigned(it->second));
    }
#else
    //Property maps post Geant4 10.4. Here the internal keys are integers, but
    //we can't be sure it is always the same in different materials. Unfortunately
    //GetConstProperty triggers a G4 exception in case the key is not present,
    //so for safety we have to trigger one extra map search via ConstPropertyExists:
    if ( matprop->GetConstPropertyMap()->empty() || !matprop->ConstPropertyExists(m_key) )
      return 0;
    return m_scatters.at((unsigned)matprop->GetConstProperty(m_key));
#endif

  }

}

#endif
