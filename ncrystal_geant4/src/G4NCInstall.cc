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

#include "G4NCrystal/G4NCInstall.hh"
#include "G4NCrystal/G4NCManager.hh"
#include "G4NCProcWrapper.hh"

#include "G4Material.hh"
#include "G4ProcessManager.hh"
#include "G4Neutron.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronElasticProcess.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "globals.hh"

#ifdef G4MULTITHREADED
#  include "G4MTRunManager.hh"
#endif

namespace G4NCrystal {
  class ProcWrapper;

  static ProcWrapper * s_proc = 0;
  void doInstall(bool onDemand) {
    if (s_proc)
      return;

    Manager * mgr = Manager::getInstance();
    if (onDemand) {
      bool foundany (false);
      if (mgr->nMaterialsWithProperties()) {
        //NB: This includes base materials, i.e. might turn up materials only
        //indirectly in the geometry.
        const G4RegionStore* region_store = G4RegionStore::GetInstance();
        for (std::vector<G4Region*>::const_iterator itRegion = region_store->begin();itRegion!=region_store->end();++itRegion) {
          (*itRegion)->UpdateMaterialList();
          std::vector<G4Material*>::const_iterator itMat = (*itRegion)->GetMaterialIterator();
          const unsigned nmat = (*itRegion)->GetNumberOfMaterials();
          for (unsigned imat=0;imat<nmat;++itMat,++imat) {
            assert(*itMat);
            if (mgr->getScatterProperty(*itMat)) {
              foundany = true;
              break;
            }
          }
          if (foundany)
            break;
        }
      }
      if (!foundany) {
        G4cout<<"G4NCInstall :: No materials with \"NCrystal\" property found in active geometry."<<G4endl;
        G4cout<<"G4NCInstall :: Will not touch existing processes for neutrons"<<G4endl;
        return;
      }
    }


#ifdef G4MULTITHREADED
    if (dynamic_cast<G4MTRunManager*>(G4RunManager::GetRunManager())) {
      G4cout<<G4endl;
      G4cout<<G4endl;
      G4cout<<G4endl;
      G4cout<<"G4NCrystal WARNING :: Detected usage of G4MTRunManager - but NCrystal code is not currently supporting multi-threaded mode!!"<<G4endl;
      G4cout<<G4endl;
      G4cout<<G4endl;
      G4cout<<G4endl;
    }
#endif

    G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
    if (!pmanager) {
      G4Exception("G4NCrystal::doInstall","Error",FatalException,
                  "Could not get process manager for neutron");
      return;//for static analysis code that do not consider G4Exception non-returning
    }

    const G4ProcessVector* pl = pmanager->GetProcessList();
    G4HadronElasticProcess* pHadElastic(0);

    auto nprocs = pl->size();
    for ( decltype(nprocs) i = 0; i < nprocs; ++i ) {
      if (!pmanager->GetProcessActivation(i))
        continue;
      G4HadronElasticProcess* pHadElasticTest = dynamic_cast<G4HadronElasticProcess*>((*pl)[i]);
      if (pHadElasticTest) {
        if (pHadElastic) {
          G4Exception("G4NCrystal::doInstall","Error",FatalException,
                      "More than one process derived from G4HadronElasticProcess"
                      " found in active process list for Neutrons");
        } else {
          pHadElastic = pHadElasticTest;
        }
      }
    }
    if (!pHadElastic) {
      G4Exception("G4NCrystal::doInstall","Error",FatalException,
                  "No process derived from G4HadronElasticProcess found in active process list for Neutrons");
    } else {
      G4cout<<"G4NCrystal :: Wrapping and replacing existing "<<
        pHadElastic->GetProcessName()<<" process for neutrons"<<G4endl;
      pmanager->AddDiscreteProcess(s_proc = new G4NCrystal::ProcWrapper(pHadElastic));
      if(!pmanager->SetProcessActivation(pHadElastic,false))
        G4Exception("G4NCrystal::doInstall","Error",FatalException,
                    "Encountered error when deactivating the neutron elastic hadronic process");
    }
  }
}

void G4NCrystal::install()
{
  doInstall(false);
}

void G4NCrystal::installOnDemand()
{
  doInstall(true);
}
