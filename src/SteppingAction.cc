//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
// #include "DetectorAnalysis.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Initialize autolock for multiple threads writing into a single file
namespace { G4Mutex myParticleLog = G4MUTEX_INITIALIZER; }


SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  backgroundFileName("../data/hits.csv"),
  signalFileName("../data/signal.csv")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // Allocate variable for particle logging checks
  G4bool isEnteringDetector, isInDetector;
  G4String volName, nextVolName;
  G4int flag;
  
  G4Track* track = aStep->GetTrack();

  // Particles are background by default
  G4bool isBackground = true;

  // Get pre and post step logical volume names
  if (track->GetVolume()) {volName = track->GetVolume()->GetName();}
  if (track->GetNextVolume()) {nextVolName = track->GetNextVolume()->GetName();}

  // Logical check between pre- and post- steps to determine if 
  // the particle is either in or entering the detector volume 
  isInDetector = (volName == "Detector" && nextVolName == "Detector");
  isEnteringDetector = (volName != "Detector" && nextVolName == "Detector");

  // Get particle name string, either "e-" or "gamma" 
  G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();

  // Creator process corresponds to particle generation via ParticleGun or 
  // a physics process, in this case electron bremsstrahlung (eBrem) 
  // or electron ionization (eIoni)
  // **N.B.: need to check CreatorProcess is not null since primary tracks have 
  // no process name and will segfault if dereferenced**
  if(particleName == "gamma")
  {
    
    isBackground = false;
    
    if(track->GetCreatorProcess() != NULL){
      const G4String& processName = track->GetCreatorProcess()->GetProcessName();
  
      // If the particle is not a primary track, then checked if it's created
      // via electron bremsstrahlung (eBrem); if so, reset to isBackgroung == true
      isBackground = (processName == "eBrem");
  
      } 
  }

  if(isInDetector){flag = 0;}
  else if(isEnteringDetector){flag = 1;}

  if(isInDetector || isEnteringDetector)
  {
    const G4StepPoint* postPoint = aStep->GetPostStepPoint();
    G4double ene = postPoint->GetKineticEnergy();
    G4ThreeVector pos = postPoint->GetPosition();
  
    if(ene > 50.*keV)
    {
      // write to background hits file	    
      if(isBackground) LogParticle(pos, ene, backgroundFileName, flag, particleName);
      
      // write to signal hits file
      else LogParticle(pos, ene, signalFileName, flag, particleName);
    }
  }    
    
}


void SteppingAction::LogParticle(G4ThreeVector pos, G4double ene, G4String detectorFileName, G4int flag, G4String PID)
{
    G4AutoLock lock(&myParticleLog);

    std::ofstream hitFile_detector;
    hitFile_detector.open(detectorFileName, std::ios_base::app);

    hitFile_detector << flag << "," << pos.x()/cm << "," << pos.y()/cm 
	    << "," << pos.z()/cm << "," << ene/keV << "," << PID << "\n";

    hitFile_detector.close();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
