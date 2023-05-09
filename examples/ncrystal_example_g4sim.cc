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


// Small self-contained example application of G4NCrystal. We simulate an
// idealistic neutron scattering experiment, where a very small spherical
// Aluminium sample is surrounded by a much larger spherical detector where hits
// are recorded and, from their position, the scattering angle is inferred and
// printed. Due to the usage of G4NCrystal, the scattering in the sample will
// correctly be dominated by diffraction in the polycrystalline Aluminium, and
// the resulting spectrum will be dominated by scatterings at certain angles
// corresponding to the Bragg edges of Aluminium.

#include "G4NCrystal/G4NCrystal.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysListFactory.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4PVPlacement.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include "G4ios.hh"

class MySD : public G4VSensitiveDetector {
  //////////////////////////////////////////////////////////////////////////////
  // Sensitive detector for monitoring neutron hits in the spherical detector
  // volume and printing out the "detected" scattering angle.
  //////////////////////////////////////////////////////////////////////////////
public:
  MySD() : G4VSensitiveDetector("MySD") {}
  virtual ~MySD(){}

  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) final
  {
    if (step->GetPreStepPoint()->GetStepStatus() != fGeomBoundary)
      return true;//must have just entered the volume
    if (step->GetTrack()->GetDynamicParticle()->GetPDGcode()!=2112)
      return true;//must be neutron
    G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
    double r = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
    if (pos.z()>0&&r<0.001*mm)
      return true;//No scattering took place
    G4cout << "Hit detected at theta = "<<atan2(r,pos.z())*NCrystal::kToDeg<<" deg"<<G4endl;
    return true;
  }
};

class MyGeo : public G4VUserDetectorConstruction {
  //////////////////////////////////////////////////////////////////////////////
  // Constructs an r=1*mm spherical sample of polycrystalline Aluminium inside an
  // r=100*cm spherical vacuum inside a 1*mm thick spherical counting volume,
  // inside a 110*cm world box. The sample is small enough for multiple neutron
  // scattering events to be negligible and the detector is far enough from the
  // sample to make sample size effects on the angular measurement equally
  // negligible.
  //////////////////////////////////////////////////////////////////////////////
public:
  MyGeo() = default;
  virtual ~MyGeo() = default;
  G4VPhysicalVolume* Construct() final
  {
    G4Material * mat_vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic",true);
    G4Material * mat_aluminium = G4NCrystal::createMaterial("Al_sg225.ncmat");

    G4LogicalVolume * world_log = new G4LogicalVolume(new G4Box("world",110*cm,110*cm,110*cm),
                                                      mat_vacuum,"world",0,0,0);
    G4PVPlacement * world_phys = new G4PVPlacement(0,G4ThreeVector(),world_log,"world",0,false,0);
    G4LogicalVolume * det_log = new G4LogicalVolume(new G4Sphere("detector",0,100.1*cm,0,CLHEP::twopi,0,CLHEP::pi),
                                                      mat_vacuum,"detector",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(),det_log,"detector",world_log,false,0);
    G4LogicalVolume * vacuum_log = new G4LogicalVolume(new G4Sphere("vacuum",0,100.0*cm,0,CLHEP::twopi,0,CLHEP::pi),
                                                       mat_vacuum,"vacuum",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(),vacuum_log,"vacuum",det_log,false,0);
    G4LogicalVolume * sample_log = new G4LogicalVolume(new G4Sphere("sample",0,0.1*cm,0,CLHEP::twopi,0,CLHEP::pi),
                                                       mat_aluminium,"sample",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(),sample_log,"sample",vacuum_log,false,0);
    MySD * sd = new MySD();
    G4SDManager::GetSDMpointer()->AddNewDetector( sd );
    det_log->SetSensitiveDetector(sd);
    return world_phys;
  }
};

class MyGun : public G4VUserPrimaryGeneratorAction {
  //////////////////////////////////////////////////////////////////////////////
  // Setup monochromatic source of neutrons, hitting the sample with initial direction (0,0,1)
  //////////////////////////////////////////////////////////////////////////////
public:

  MyGun(double neutron_wavelength) : m_particleGun(new G4ParticleGun(1))
  {
    m_particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("neutron"));
    m_particleGun->SetParticleEnergy(neutronWavelengthToEKin(neutron_wavelength));
    m_particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -1.0*cm));
    m_particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  }

  virtual ~MyGun() = default;

  void GeneratePrimaries(G4Event* evt) final
  {
    m_particleGun->GeneratePrimaryVertex(evt);
  }

  double neutronWavelengthToEKin(double l) {
    return 0.5 * CLHEP::h_Planck * CLHEP::h_Planck * CLHEP::c_squared / (l*l*CLHEP::neutron_mass_c2);
  }

private:
  std::unique_ptr<G4ParticleGun> m_particleGun;
};

int main(int,char**) {

  NCrystal::libClashDetect();//Detect broken installation

  //Set seed:
  CLHEP::HepRandom::setTheSeed(123);

  //G4 Run manager:
  std::unique_ptr<G4RunManager> runManager(new G4RunManager);

  //Setup geometry and physics-list:
  runManager->SetUserInitialization(new MyGeo);
  runManager->SetUserInitialization(G4PhysListFactory().GetReferencePhysList("QGSP_BIC_HP"));

  //Setup monochromatic source of 4.0Aa neutrons. Note that at 4.0 Aangstrom,
  //more than 90% of scattering events in polycrystalline aluminium are coherent
  //with peaks only at theta = 118 deg and theta = 162 deg:
  const double neutron_wavelength = 4.0*angstrom;
  runManager->SetUserAction(new MyGun(neutron_wavelength));

  //Initialize g4 run manager:
  runManager->Initialize();

  //Install G4NCrystal:
  G4NCrystal::installOnDemand();

  //Run 1000 events (should give us ~20 scattering events on average for 4.0 aangstrom):
  runManager->BeamOn(1000);

  //Cleanup:
  G4NCrystal::Manager::cleanup();//delete manager singleton, unref cached ncrystal objects (for valgrind).

  return 0;
}
