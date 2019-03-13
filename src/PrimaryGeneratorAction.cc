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
#include <random>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  
  fParticleGun  = new G4ParticleGun();
  G4ParticleDefinition* electronParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(electronParticle); 

  // Default parameters
  //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.5*m,0.));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
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

  G4double PI = 3.14159265358979323846;

  // N particles generated per simulation run
  G4int n_particle = 1000;

  // Allocate variables for random position, direction
  G4double xPos,yPos,zPos,xDir,yDir,zDir;

  // Radius of sphere surface where particles are generated
  G4double sphereR = 30.*cm;

  // Loss cone angle (half angle w.r.t. sphere angle) at 500 km in radians
  G4double theta_exclusion = 64.*PI/180.;

  // E-folding (E0 energy)
  double E0 = 300.;

  std::default_random_engine generator;
  std::exponential_distribution<double> exp_distribution(1/E0);
  
  // f0 flux
  //G4double f0 = 2e8;


  for(G4int i = 0; i<n_particle; i++){


    // Reset position and direction of particle
    xPos = 0; xDir = 0;
    yPos = 0; yDir = 0;
    zPos = 0; zDir = 0;

    
    // Calculate random particle position on sphere, excluding spherical cap    
    do{
      // Rand on [-1, 1)
      G4double u = G4UniformRand()*2.-1.;

      // Rand on [0, 2*pi)
      G4double theta = G4UniformRand()*2.*PI;

      xPos = sphereR * std::sqrt(1 - u * u) * std::cos(theta);
      yPos = sphereR * std::sqrt(1 - u * u) * std::sin(theta);
      zPos = sphereR * u;
      }
    while(yPos > sphereR * std::cos(theta_exclusion)); 
// exits when y position falls below spherical cap

    // Uniform random numbers on [0, 1)
    xDir = G4UniformRand();
    yDir = G4UniformRand();
    zDir = G4UniformRand();

    // Enforces inward directionality to particles
    if(xPos > 0) xDir = -xDir;
    if(yPos > 0) yDir = -yDir;
    if(zPos > 0) zDir = -zDir;

    // Selects random energy according to exponential distribution
    G4double randEnergy = exp_distribution(generator)*300.*keV;

    fParticleGun->SetParticlePosition(G4ThreeVector(xPos, yPos, zPos));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xDir, yDir, zDir));
    fParticleGun->SetParticleEnergy(randEnergy);

    fParticleGun->GeneratePrimaryVertex(anEvent);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

