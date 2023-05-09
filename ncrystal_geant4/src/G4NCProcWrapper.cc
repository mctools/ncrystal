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

#include "G4NCProcWrapper.hh"
#include "G4NCrystal/G4NCManager.hh"
#include "NCrystal/NCRNG.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VParticleChange.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronElasticProcess.hh"
#include "G4ParticleChange.hh"
#include "G4VTouchable.hh"
#include "G4NavigationHistory.hh"

namespace NC = NCrystal;
namespace NCG4 = G4NCrystal;

namespace G4NCrystal {
  namespace {
    class RNG_G4Wrapper : public NC::RNGStream {
      CLHEP::HepRandomEngine * m_engine;
    public:
      //Can be cheaply created on the stack just before being used in calls to
      //ProcImpl::Scatter objects, like:
      //
      //  RNG_G4Wrapper rng(G4Random::getTheEngine());
      //
      //Since G4Random::getTheEngine() always returns a thread-local engine, this is then MT-safe!
      constexpr RNG_G4Wrapper(CLHEP::HepRandomEngine * e) noexcept : m_engine(e) {}
    protected:
      double actualGenerate() override { return m_engine->flat(); }
    };

  }
}

NCG4::ProcWrapper::ProcWrapper(G4HadronElasticProcess * procToWrap,
    const G4String& processName)
: G4VDiscreteProcess(processName.empty()?procToWrap->GetProcessName():processName),
  m_wrappedProc(procToWrap), m_mgr(Manager::getInstance())
{
  if (verboseLevel>1)
    G4cout << "G4NCrystal ProcWrapper named "<<GetProcessName() << " is created"<< G4endl;
}

NCG4::ProcWrapper::~ProcWrapper()
{
}

G4VParticleChange* NCG4::ProcWrapper::PostStepDoIt(const G4Track& trk, const G4Step& step)
{
  assert(trk.GetParticleDefinition()->GetPDGEncoding()==2112);

  const double ekin = trk.GetKineticEnergy();

  //Important to always clear the interaction lengths here, even when we pass on
  //the call to the wrapped process below (since G4 does not directly interact
  //with the now disabled wrapped process):
  ClearNumberOfInteractionLengthLeft();

  //TODO: Revisit the 5eV global threshold. Issue is perhaps tungsten
  //which has resonances below 5eV. If possible, auto-detect the threshold by
  //comparing G4 and NCrystal xsects around the threshold. Otherwise allow
  //per-material specifications.

  Manager::ProcAndCache procandcache;
  if ( !(ekin>0.0)
       || ekin > 5*CLHEP::eV
       || ( procandcache = m_mgr->getScatterPropertyWithThreadSafeCache( trk.GetMaterial() ) ).first == nullptr )
    return m_wrappedProc->PostStepDoIt(trk,step);

 auto& ncscat = *procandcache.first;
 assert(procandcache.second!=nullptr);
 auto& cacheptr = *procandcache.second;

  G4ThreeVector g4outcome_dir;
  double g4outcome_ekin;
  try {
    auto random_engine = G4Random::getTheEngine();
    nc_assert(random_engine!=nullptr);
    RNG_G4Wrapper rng(random_engine);
    constexpr double inv_eV = 1.0/CLHEP::eV;
    NC::NeutronEnergy nc_ekin_in{ekin * inv_eV};//NCrystal unit is eV
    const G4ThreeVector& indir = trk.GetMomentumDirection();
    if( ! ncscat.isOriented() ) {
      //Orientation of material does not matter:
      auto outcome = ncscat.sampleScatter(cacheptr,rng,nc_ekin_in, NC::NeutronDirection{indir.x(),indir.y(),indir.z()} );
      g4outcome_ekin = outcome.ekin.get() * CLHEP::eV;
      g4outcome_dir.set(outcome.direction[0],outcome.direction[1],outcome.direction[2]);
    } else {
      //Orientation of material matters, need to transform to-and-from the frame of the volume (touchable):
      const G4AffineTransform &trf = step.GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform();
      G4ThreeVector indir_local = trf.TransformAxis(indir);
      auto outcome = ncscat.sampleScatter(cacheptr,rng,nc_ekin_in, NC::NeutronDirection{indir_local.x(),indir_local.y(),indir_local.z()} );
      g4outcome_ekin = outcome.ekin.get() * CLHEP::eV;
      g4outcome_dir = trf.Inverse().TransformAxis(G4ThreeVector{outcome.direction[0],outcome.direction[1],outcome.direction[2]});
    }
  } catch ( NC::Error::Exception& e ) {
    Manager::handleError("G4NCrystal::ProcWrapper::PostStepDoIt",101,e);
  }
  m_particleChange.Clear();
  m_particleChange.Initialize(trk);
  m_particleChange.ProposeWeight(trk.GetWeight());//a bit weird that we have to do this.
  m_particleChange.ProposeMomentumDirection( g4outcome_dir );
  m_particleChange.ProposeEnergy( g4outcome_ekin );
  return &m_particleChange;
}


G4double NCG4::ProcWrapper::GetMeanFreePath(const G4Track& trk, G4double p, G4ForceCondition* f)
{
  assert(trk.GetParticleDefinition()->GetPDGEncoding()==2112);
  const double ekin = trk.GetKineticEnergy();

  Manager::ProcAndCache procandcache;
  if ( !(ekin>0.0)
       || ekin > 5*CLHEP::eV
       || ( procandcache = m_mgr->getScatterPropertyWithThreadSafeCache( trk.GetMaterial() ) ).first == nullptr )
    return m_wrappedProc->GetMeanFreePath(trk,p,f);


  auto& ncscat = *procandcache.first;
  assert(procandcache.second!=nullptr);
  auto& cacheptr = *procandcache.second;

  double xs(0.0);
  try {
    constexpr double inv_eV = 1.0/CLHEP::eV;
    NC::NeutronEnergy nc_ekin_in{ekin * inv_eV};//NCrystal unit is eV
    const G4ThreeVector& indir = trk.GetMomentumDirection();
    if( ! ncscat.isOriented() ) {
      xs = ncscat.crossSection( cacheptr, nc_ekin_in, NC::NeutronDirection{indir.x(),indir.y(),indir.z()}).get() * CLHEP::barn;
    } else {
      G4ThreeVector indir_local = trk.GetStep()->GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform().TransformAxis(indir);
      xs = ncscat.crossSection( cacheptr, nc_ekin_in, NC::NeutronDirection{indir_local.x(),indir_local.y(),indir_local.z()} ).get() * CLHEP::barn;
    }
  } catch ( NC::Error::Exception& e ) {
    Manager::handleError("G4NCrystal::ProcWrapper::GetMeanFreePath",102,e);
  }

  return xs
      ? 1.0 / ( trk.GetMaterial()->GetTotNbOfAtomsPerVolume() * xs )
      : NC::kInfinity ;
}

void NCG4::ProcWrapper::BuildPhysicsTable(const G4ParticleDefinition&)
{
}

G4bool NCG4::ProcWrapper::IsApplicable(const G4ParticleDefinition& pd)
{
  return pd.GetPDGEncoding()==2112;
}

