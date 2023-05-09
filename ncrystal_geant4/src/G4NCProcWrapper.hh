#ifndef G4NCrystal_ProcWrapper_hh
#define G4NCrystal_ProcWrapper_hh

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

#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"

class G4HadronElasticProcess;
class G4Material;

namespace G4NCrystal {

  class Manager;

  class ProcWrapper : public G4VDiscreteProcess
  {
    // Wrapper process used by G4NCInstall to dynamically support NCrystal
    // physics with any physics model.

  public:
    ProcWrapper(G4HadronElasticProcess * procToWrap,
                const G4String& processName = "");
    virtual ~ProcWrapper();

    //Intercepted methods in which NCrystal scatter physics is supplied when
    //appropriate, passing the call through to the wrapped process when not:
    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& ) final;
    G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) final;

    void BuildPhysicsTable(const G4ParticleDefinition&) final;
    G4bool IsApplicable(const G4ParticleDefinition& pd) final;

  private:
    G4ParticleChange m_particleChange;
    G4HadronElasticProcess * m_wrappedProc;
    Manager * m_mgr;
  };

}

#endif
