////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "G4NCProcWrapper.hh"
#include "G4NCrystal/G4NCManager.hh"
#include "NCrystal/NCScatter.hh"
#include "NCrystal/NCException.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronElasticProcess.hh"
#include "G4ParticleChange.hh"
#include "G4VTouchable.hh"
#include "G4NavigationHistory.hh"

G4NCrystal::ProcWrapper::ProcWrapper(G4HadronElasticProcess * procToWrap,
    const G4String& processName)
: G4VDiscreteProcess(processName.empty()?procToWrap->GetProcessName():processName),
  m_wrappedProc(procToWrap), m_mgr(G4NCrystal::Manager::getInstance())
{
  if (verboseLevel>1)
    G4cout << "G4NCrystal ProcWrapper named "<<GetProcessName() << " is created"<< G4endl;
}

G4NCrystal::ProcWrapper::~ProcWrapper()
{
}

G4VParticleChange* G4NCrystal::ProcWrapper::PostStepDoIt(const G4Track& trk, const G4Step& step)
{
  nc_assert(trk.GetParticleDefinition()->GetPDGEncoding()==2112);

  const NCrystal::Scatter* scat;
  const double ekin = trk.GetKineticEnergy();

  //Important to always clear the interaction lengths here, even when we pass on
  //the call to the wrapped process below (since G4 does not directly interact
  //with the now disabled wrapped process):
  ClearNumberOfInteractionLengthLeft();


  //TODO for NC2: Revisit the 5eV global threshold. Issue is perhaps tungsten
  //which has resonances below 5eV. If possible, auto-detect the threshold by
  //comparing G4 and NCrystal xsects around the threshold. Otherwise allow
  //per-material specifications.
  if ( ekin > 5*CLHEP::eV || !(scat=m_mgr->getScatterProperty(trk.GetMaterial()) ) || !ekin )
    return m_wrappedProc->PostStepDoIt(trk,step);


  //Let the NCrystal::Scatter instance do its work!

  double ekin_eV = ekin * (1.0/CLHEP::eV);//NCrystal unit is eV
  const G4ThreeVector& indir = trk.GetMomentumDirection();
  G4ThreeVector outdir;
  double delta_ekin_eV;
  try {
    if(!scat->isOriented()) {
      scat->generateScattering( ekin_eV, NC_CVECTOR_CAST(indir), NC_VECTOR_CAST(outdir), delta_ekin_eV );
    } else {
      const G4AffineTransform &trf = step.GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform();
      G4ThreeVector inlocalDir = trf.TransformAxis(indir);
      G4ThreeVector outlocaldir;
      scat->generateScattering( ekin_eV, NC_CVECTOR_CAST(inlocalDir), NC_VECTOR_CAST(outlocaldir), delta_ekin_eV );
      outdir = trf.Inverse().TransformAxis(outlocaldir);
    }
  } catch ( NCrystal::Error::Exception& e ) {
    Manager::handleError("G4NCrystal::ProcWrapper::PostStepDoIt",101,e);
  }
  m_particleChange.Clear();
  m_particleChange.Initialize(trk);
  m_particleChange.ProposeWeight(trk.GetWeight());//a bit weird that we have to do this.
  m_particleChange.ProposeMomentumDirection(outdir);
  m_particleChange.ProposeEnergy( ekin + delta_ekin_eV * CLHEP::eV );
  return &m_particleChange;
}


G4double G4NCrystal::ProcWrapper::GetMeanFreePath(const G4Track& trk, G4double p, G4ForceCondition* f)
{
  nc_assert(trk.GetParticleDefinition()->GetPDGEncoding()==2112);
  const NCrystal::Scatter* scat;
  const double ekin = trk.GetKineticEnergy();
  if ( ekin > 5*CLHEP::eV || !ekin || !(scat=m_mgr->getScatterProperty(trk.GetMaterial())) )
    return m_wrappedProc->GetMeanFreePath(trk,p,f);

  double ekin_eV = ekin * (1.0/CLHEP::eV);//NCrystal unit is eV
  const G4ThreeVector& indir = trk.GetMomentumDirection();

  double xs(0.0);
  try {
    if(!scat->isOriented()) {
      xs = scat->crossSection(ekin_eV, NC_CVECTOR_CAST(indir)) * CLHEP::barn;
    } else {
      G4ThreeVector inlocalDir = trk.GetStep()->GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform().TransformAxis(indir);
      xs = scat->crossSection(ekin_eV, NC_CVECTOR_CAST(inlocalDir)) * CLHEP::barn;
    }
  } catch ( NCrystal::Error::Exception& e ) {
    Manager::handleError("G4NCrystal::ProcWrapper::GetMeanFreePath",102,e);
  }

  return xs
      ? 1.0 / ( trk.GetMaterial()->GetTotNbOfAtomsPerVolume() * xs )
          : std::numeric_limits<double>::infinity() ;
}

void G4NCrystal::ProcWrapper::BuildPhysicsTable(const G4ParticleDefinition&)
{
}

G4bool G4NCrystal::ProcWrapper::IsApplicable(const G4ParticleDefinition& pd)
{
  return pd.GetPDGEncoding()==2112;
}
