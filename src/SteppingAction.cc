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

#include "SteppingMessenger.hh"

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
  fSteppingMessenger(0)
{

  fSteppingMessenger = new SteppingMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
  delete fSteppingMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // Allocate variable for particle logging checks
  G4bool check1, check2;
  check1 = check2 = false;
  G4String volName, nextVolName;
  
  G4Track* track = aStep->GetTrack();

  // Get pre and post step logical volume names
  if (track->GetVolume()) {volName = track->GetVolume()->GetName();}
  if (track->GetNextVolume()) {nextVolName = track->GetNextVolume()->GetName();}


  // Below is the logical operation that determines if a particle is 
  // either in or entering the detector
  

  // Searching for match to "av_1_impr_X_Detector_pv_Y", starting from
  // character 12, where X is the impression copy 
  // (detector assembly) [1,2,3] and Y is the detector number 
  // within each assembly [2,4,6,8]
  check1 = (std::string::npos == volName.find("P", 12));
  check2 = (std::string::npos != nextVolName.find("P", 12));
  

  // Creator process corresponds to particle generation via ParticleGun or
  // a physics process, in this case electron bremsstrahlung (eBrem) 
  // or electron ionization (eIoni)
  // **N.B.: need to check CreatorProcess is not null since primary tracks
  // have no process name and will segfault if dereferenced**

  if(check1 && check2) // particle is in detector
  {
    // Position of hit
    //const G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition();
    
    // Particle energy immediately before it enters detector material 
    const G4double ene = aStep->GetPreStepPoint()->GetKineticEnergy(); 

    // Get generated particle position, energy, and momentum direction
    const G4ThreeVector vtx = track->GetVertexPosition();

    // Redlen lower energy detection threshold
    if(ene > 20.*keV)
    {
  
      // Get particle name string, either "e-" or "gamma" 
      G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
     
      if(particleName == "gamma" && track->GetCreatorProcess() == NULL)
      // Write to signal hits file (no part. name since all are gammas)
        {LogParticle(vtx, ene, signalFileName, nextVolName);} 
      else
     // Write to background hits file	   
        {LogParticle(vtx, ene, backgroundFileName, nextVolName);}

    }
  }    
  

}

void SteppingAction::LogParticle(G4ThreeVector vtx, G4double ene, G4String detectorFileName, G4String volName)
{
    // locks program so that multiple threads cannot write to file
    // at once, unlocks when current scope (i.e. this method) is left
    G4AutoLock lock(&myParticleLog);
    
    G4int loc1  = volName.find('P', 10);
    G4int loc2  = volName.find('_', loc1);
    G4int loc2p = volName.find('i', loc1);
    G4int loc3  = volName.find('j', loc2);
    G4int loc4  = volName.find('_', loc3);
 
    
    // Find detector and pixel number from volume name string
    G4String avPlacementNum = volName.substr(3, volName.find('_',3)-3);
    G4String detNumber   = volName.substr(10, volName.find('_', 10)-10);
    G4String iNum   	 = volName.substr(loc2p+1, loc3-(loc2p+2));
    G4String jNum 	 = volName.substr(loc3+1, loc4-(loc3+1));

    //std::cout << volName << std::endl; 
    //std::cout << avPlacementNum << std::endl; 

    // Correction needed after dark detector was implemented
    if(std::stoi(avPlacementNum) == 2)
    {
      switch(std::stoi(detNumber))
      {
	case 2:
	  detNumber = "12";
	  break;
	case 3:
	  detNumber = "13";
	  break;
	case 4:
	  detNumber = "14";
	  break;
	case 5:
	  detNumber = "15";
          break;
	default:
	  throw std::invalid_argument("Something wrong with detector switch method!");
      }
    }


    std::ofstream hitFile_detector;
    hitFile_detector.open(detectorFileName, std::ios_base::app);

    hitFile_detector << ene/keV << "," 
		     << detNumber << "," 
	    	     << iNum << "," 
		     << jNum << ","
		     << vtx.x()/cm << ","
		     << vtx.y()/cm << ","
		     << vtx.z()/cm << "\n";

    hitFile_detector.close();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

