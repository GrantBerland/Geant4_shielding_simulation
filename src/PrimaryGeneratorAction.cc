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
// $Id: PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{

  G4double PI = 3.14159265358979323846;
  G4int n_particle = 10;
  fParticleGun  = new G4ParticleGun();
  G4ParticleTable* particleType = G4ParticleTable::GetParticleTable();

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(particleType->FindParticle("e-")); 


  G4double xPos,yPos,zPos,xDir,yDir,zDir;
  G4double sphereR = 30.*cm;
  G4double theta_exclusion = 64.*PI/180.*rad;

  for(G4int i = 0; i<n_particle; i++){

    // Reset position and direction of particle
    xPos = 0; xDir = 0;
    yPos = 0; yDir = 0;
    zPos = 0; zDir = 0;

    // Rand on u = [-1, 1)
    G4double u = G4UniformRand()*2.-1.;

    // Rand on theta = [0, 2*pi)
    G4double theta = G4UniformRand()*2.*PI;
    
    // Calculate random particle position on sphere, excluding spherical cap    
    do{
      xPos = sphereR * sqrt(1 - u * u) * cos(theta);
      yPos = sphereR * sqrt(1 - u * u) * sin(theta);
      zPos = sphereR * u;
      }
    while(yPos < sphereR * cos(theta_exclusion));


    xDir = G4UniformRand();
    yDir = G4UniformRand();
    zDir = G4UniformRand();

    // Enforces inward directionality to particles
    if(xPos > 0) xDir = -xDir;
    if(yPos > 0) yDir = -yDir;
    if(zPos > 0) zDir = -zDir;

    G4cout << "Generating particle " << i << G4endl;

    // Selects random energy according to exponential distribution
    G4double randEnergy = G4UniformRand()*100.*keV;


    fParticleGun->SetParticlePosition(G4ThreeVector(xPos, yPos, zPos)); 
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xDir, yDir, zDir));
    fParticleGun->SetParticleEnergy(randEnergy);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

