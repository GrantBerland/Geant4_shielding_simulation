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

  fParticleGun  = new G4ParticleGun();
  G4ParticleDefinition* electronParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(electronParticle);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateLossConeSample(LossConeSample* r)
{
  G4int dataSize = 32;
  G4double lossConeData[32][2] = {
  { 1.0000 , 0.0000 },
  { 3.0000 , 0.0019 },
  { 5.0000 , 0.0056 },
  { 7.0000 , 0.0110 },
  { 9.0000 , 0.0182 },
  { 11.0000 , 0.0272 },
  { 13.0000 , 0.0379 },
  { 15.0000 , 0.0503 },
  { 17.0000 , 0.0647 },
  { 19.0000 , 0.0817 },
  { 21.0000 , 0.1007 },
  { 23.0000 , 0.1210 },
  { 25.0000 , 0.1430 },
  { 27.0000 , 0.1671 },
  { 29.0000 , 0.1938 },
  { 31.0000 , 0.2222 },
  { 33.0000 , 0.2527 },
  { 35.0000 , 0.2852 },
  { 37.0000 , 0.3189 },
  { 39.0000 , 0.3548 },
  { 41.0000 , 0.3925 },
  { 43.0000 , 0.4320 },
  { 45.0000 , 0.4743 },
  { 47.0000 , 0.5189 },
  { 49.0000 , 0.5653 },
  { 51.0000 , 0.6147 },
  { 53.0000 , 0.6676 },
  { 55.0000 , 0.7238 },
  { 57.0000 , 0.7834 },
  { 59.0000 , 0.8479 },
  { 61.0000 , 0.9186 },
  { 63.0000 , 1.0000 }};

  //G4double fluxSum = 18.155822371325151;
  G4double randomNumber = G4UniformRand();
  G4int angleIndex = -1;

  // Find the random angle that corresponds to loss cone flux
  for(G4int angle = 0; angle<dataSize; angle++){
    if(randomNumber < lossConeData[angle][1]){
      angleIndex = angle;
      break; }
  }


  // Selects exponential folding energy E0 based on backscattered pitch angle range
  G4double E0 = -1;
  if(     lossConeData[angleIndex][0] < 20) {E0 = 159.;}
  else if(lossConeData[angleIndex][0] < 30) {E0 = 141.;}
  else if(lossConeData[angleIndex][0] < 40) {E0 = 177.;}
  else if(lossConeData[angleIndex][0] < 50) {E0 = 194.;}
  else if(lossConeData[angleIndex][0] < 64) {E0 = 230.;}
  
  // Mathematics spherical coordinates definition!!!
  G4double sphereR = 15.*cm;
  G4double PI = 3.14159265358979323846;

  // Angle about the field line theta on [0 , 2pi)
  G4double theta = G4UniformRand()*2.*PI; 
  
  // Zenith angle (pitch angle)
  G4double phi = lossConeData[angleIndex][0] * PI / 180.;
  
  // We want our Y direction to be "up"
  r->x = sphereR * std::cos(theta) * std::sin(phi);
  r->y = sphereR * std::cos(phi);
  r->z = sphereR * std::sin(theta) * std::sin(phi);
  

  // Uniform random numbers on [0, 1)
  r->xDir = G4UniformRand();
  r->yDir = G4UniformRand();
  r->zDir = G4UniformRand();

  //G4double norm = std::sqrt(r->xDir * r->xDir + r->yDir * r->yDir + r->zDir * r->zDir);

  // Enforces inward directionality to particles
  if(r->x > 0) {r->xDir = -r->xDir;}// else{r->xDir = r->xDir/norm;}
  if(r->y > 0) {r->yDir = -r->yDir;}// else{r->yDir = r->yDir/norm;}
  if(r->z > 0) {r->zDir = -r->zDir;}// else{r->zDir = r->zDir/norm;}


  randomNumber = G4UniformRand();

  // Inverse CDF sampling for exponential RV
  r->energy = (std::log(1 - randomNumber)*-E0)*keV;
}
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event

  G4double PI = 3.14159265358979323846;

  // N particles generated per simulation run
  G4int nParticles = 100;	// trapped particles
  G4int nLCparticles = std::floor(0.1851*nParticles);	
  // loss cone particles (backscattered)

  // Allocate variables for random position, direction
  G4double xPos,yPos,zPos,xDir,yDir,zDir;

  // Constant sphere offsets
  //G4double xShift = 5.*cm;
  //G4double yShift = 5.*cm;
  //G4double zShift = 5.*cm;
  //TODO: why do these lines fuck everything up?-> try no cm multiplication
  // it worked! but why?

  G4double xShift = 0.;
  G4double yShift = 0.;
  G4double zShift = 0.;
  
  // Radius of sphere surface where particles are generated
  G4double sphereR = 15.*cm;

  // Loss cone angle (same as polar angle, phi) at 500 km, in radians
  G4double theta_exclusion = 64.*PI/180.;

  // E-folding (E0 energy) in keV (Wei)
  G4double E0 = 150.;
  for(G4int i = 0; i<nParticles; i++){


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
      yPos = sphereR * u;
      zPos = sphereR * std::sqrt(1 - u * u) * std::sin(theta); 
      }
    while(yPos > sphereR * std::cos(theta_exclusion));
    // exits when y position falls below spherical cap

    // Uniform random numbers on [0, 1)
    xDir = G4UniformRand();
    yDir = G4UniformRand();
    zDir = G4UniformRand();

    //G4double norm = std::sqrt(xDir * xDir + yDir * yDir + zDir * zDir);
    
    // Enforces inward directionality to particles
    if(xPos > 0) {xDir = -xDir;}// else {xDir = xDir/norm;}
    if(yPos > 0) {yDir = -yDir;}// else {yDir = yDir/norm;}
    if(zPos > 0) {zDir = -zDir;}// else {zDir = zDir/norm;}


    // Selects random energy according to exponential distribution
    G4double randomNumber = G4UniformRand();
    G4double randEnergy = std::log(1 - randomNumber)*-E0*keV;


    fParticleGun->SetParticlePosition(G4ThreeVector(xPos+xShift, yPos+yShift, zPos+zShift));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xDir, yDir, zDir));
    fParticleGun->SetParticleEnergy(randEnergy);

    fParticleGun->GeneratePrimaryVertex(anEvent);

    // for debuggin
    /*    
    std::ofstream particleSource;

    particleSource.open("./detector_trappedSourceParts.txt", std::ios_base::app);

    particleSource << (xPos+xShift)/cm << "," << (yPos+yShift)/cm << "," << (zPos+zShift)/cm
    << "," << xDir << "," << yDir << "," << zDir << "," << randEnergy/keV << "," << sphereR <<"\n";
    particleSource.close();
    */
  }


    LossConeSample* r = new LossConeSample();
  

    //std::ofstream particleLCSource;
    //particleLCSource.open("./detector_LCSourceParts.txt", std::ios_base::app);

  for(G4int i = 0; i<nLCparticles; i++){

    GenerateLossConeSample(r);
    
    // For debugging
    /*
    fParticleGun->SetParticlePosition(G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);

    fParticleGun->GeneratePrimaryVertex(anEvent);

    particleLCSource << (r->x+xShift)/cm << "," << (r->y+yShift)/cm << "," << (r->z+zShift)/cm
    << "," << r->xDir << "," << r->yDir << "," << r->zDir << "," << r->energy << "\n";
    */
}

    //particleLCSource.close();
    delete r;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
