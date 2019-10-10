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

#include "PrimaryGeneratorMessenger.hh"

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

#include <fstream>
#include <stdexcept>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  E_folding(100.),
  fPI(3.14159265358979323846),
  fDeg2Rad(3.14159265358979323846 / 180.),
  sphereR(25.*cm),
  lossConeAngleDeg(64.),
  photonPhiLimitDeg(40.), 
  fWhichParticle(0),
  electronParticle(0),
  photonParticle(0),
  fPrimaryGeneratorMessenger(0)
{

  fParticleGun  = new G4ParticleGun();
 
  fPrimaryGeneratorMessenger = new PrimaryGeneratorMessenger(this);

  electronParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  
  photonParticle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryGeneratorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateLossConeElectrons(ParticleSample* r)
{
  // This method generates electrons that were in the loss cone, but have
  // backscattered off of the atmosphere and are directed anti-Earthward
  // and towards the satellite.


  // Discrete CDF of particles within the loss cone and backscattered
  // from the atmosphere with 0.6527 degree resolution from [0, 63] 
  // degrees (i.e. the loss cone at approximately 500 km altitude)
  const G4int dataSize = 96;
  static const G4double lossConeData[96][2] = {
{ 1.0000 , 0.0000 },{ 1.6526 , 0.0006 },{ 2.3053 , 0.0012 },
{ 2.9579 , 0.0019 },{ 3.6105 , 0.0030 },{ 4.2632 , 0.0042 },
{ 4.9158 , 0.0054 },{ 5.5684 , 0.0071 },{ 6.2211 , 0.0089 },
{ 6.8737 , 0.0107 },{ 7.5263 , 0.0129 },{ 8.1789 , 0.0153 },
{ 8.8316 , 0.0176 },{ 9.4842 , 0.0204 },{ 10.1368 , 0.0233 },
{ 10.7895 , 0.0263 },{ 11.4421 , 0.0296 },{ 12.0947 , 0.0331 },
{ 12.7474 , 0.0366 },{ 13.4000 , 0.0404 },{ 14.0526 , 0.0444 },
{ 14.7053 , 0.0485 },{ 15.3579 , 0.0529 },{ 16.0105 , 0.0575 },
{ 16.6632 , 0.0622 },{ 17.3158 , 0.0673 },{ 17.9684 , 0.0729 },
{ 18.6211 , 0.0785 },{ 19.2737 , 0.0843 },{ 19.9263 , 0.0905 },
{ 20.5789 , 0.0967 },{ 21.2316 , 0.1030 },{ 21.8842 , 0.1097 },
{ 22.5368 , 0.1163 },{ 23.1895 , 0.1231 },{ 23.8421 , 0.1303 },
{ 24.4947 , 0.1375 },{ 25.1474 , 0.1448 },{ 25.8000 , 0.1526 },
{ 26.4526 , 0.1605 },{ 27.1053 , 0.1685 },{ 27.7579 , 0.1772 },
{ 28.4105 , 0.1859 },{ 29.0632 , 0.1947 },{ 29.7158 , 0.2040 },
{ 30.3684 , 0.2133 },{ 31.0211 , 0.2226 },{ 31.6737 , 0.2325 },
{ 32.3263 , 0.2425 },{ 32.9789 , 0.2524 },{ 33.6316 , 0.2630 },
{ 34.2842 , 0.2736 },{ 34.9368 , 0.2842 },{ 35.5895 , 0.2951 },
{ 36.2421 , 0.3062 },{ 36.8947 , 0.3172 },{ 37.5474 , 0.3288 },
{ 38.2000 , 0.3405 },{ 38.8526 , 0.3522 },{ 39.5053 , 0.3643 },
{ 40.1579 , 0.3766 },{ 40.8105 , 0.3889 },{ 41.4632 , 0.4017 },
{ 42.1158 , 0.4145 },{ 42.7684 , 0.4274 },{ 43.4211 , 0.4409 },
{ 44.0737 , 0.4547 },{ 44.7263 , 0.4685 },{ 45.3789 , 0.4827 },
{ 46.0316 , 0.4973 },{ 46.6842 , 0.5119 },{ 47.3368 , 0.5267 },
{ 47.9895 , 0.5419 },{ 48.6421 , 0.5570 },{ 49.2947 , 0.5726 },
{ 49.9474 , 0.5887 },{ 50.6000 , 0.6048 },{ 51.2526 , 0.6213 },
{ 51.9053 , 0.6386 },{ 52.5579 , 0.6559 },{ 53.2105 , 0.6735 },
{ 53.8632 , 0.6919 },{ 54.5158 , 0.7102 },{ 55.1684 , 0.7288 },
{ 55.8211 , 0.7483 },{ 56.4737 , 0.7677 },{ 57.1263 , 0.7875 },
{ 57.7789 , 0.8085 },{ 58.4316 , 0.8296 },{ 59.0842 , 0.8509 },
{ 59.7368 , 0.8740 },{ 60.3895 , 0.8970 },{ 61.0421 , 0.9203 },
{ 61.6947 , 0.9469 },{ 62.3474 , 0.9734 },{ 63.0000 , 1.0000 }};
  
  G4double randomNumber = G4UniformRand();
  G4int angleIndex = -1;

  // Discrete inverse CDF sampling to determine the spatial distribution
  // of backscattered particles from the atmosphere
  for(G4int angle = 0; angle<dataSize; angle++)
  {    
    if(randomNumber < lossConeData[angle][1])
    {
      angleIndex = angle;
      break; 
    }
  }


  // Selects exponential folding energy E0 based on backscattered 
  // pitch angle range
  G4double E0 = -1;
  if(     lossConeData[angleIndex][0] < 20) {E0 = 159.;}
  else if(lossConeData[angleIndex][0] < 30) {E0 = 141.;}
  else if(lossConeData[angleIndex][0] < 40) {E0 = 177.;}
  else if(lossConeData[angleIndex][0] < 50) {E0 = 194.;}
  else if(lossConeData[angleIndex][0] < 64) {E0 = 230.;}
  

  // NB: Mathematics spherical coordinates definition used below
  
  // Uniform sampling of angle about the field line theta on [0 , 2pi)
  G4double theta = G4UniformRand() * 2. * fPI; 
  
  // Zenith angle (pitch angle), converted to radians
  G4double phi = lossConeData[angleIndex][0] * fDeg2Rad;
  
  // We want our Y direction to be "up," otherwise standard
  // spherical to cartesian coordinate transform
  
  r->x = sphereR * std::cos(theta) * std::sin(phi);
  r->y = sphereR * std::cos(phi);
  r->z = sphereR * std::sin(theta) * std::sin(phi);
  


  // Uniform random numbers on [0, 1) for particle direction
  r->xDir = G4UniformRand();
  r->yDir = G4UniformRand();
  r->zDir = G4UniformRand();


  // Enforces inward directionality to particles
  if(r->x > 0) {r->xDir = -r->xDir;}
  if(r->y > 0) {r->yDir = -r->yDir;}
  if(r->z > 0) {r->zDir = -r->zDir;}


  // Continous inverse CDF sampling for exponential energy distribution
  randomNumber = G4UniformRand();
  r->energy = ((std::log(1 - randomNumber)*-E0))*keV;

}

void PrimaryGeneratorAction::GenerateSignalSource(ParticleSample* r)
{
  // This method generates a photon signal from a energetic particle
  // precipitation event that occurs at approximately 50 - 100 km 
  // altitude. Variables are photon energy and the 
  // zenith angle distribution of the photons.

  G4double theta, randomNumber;
  
  //Sample from exponential, energy dist. fitted to Wei's results
  // (valid for E0,source = 100 keV) 
  G4double shiftThreshold = 50.;
  G4double meanEnergy     = 241.4;

  // Rejection sampling on particle energy 
  do
  {
    randomNumber = G4UniformRand();
    
    r->energy = -(meanEnergy - shiftThreshold) *
	    std::log(1 - randomNumber) * keV;
  
  } while(r->energy < shiftThreshold*keV);
  

  // Uniformly distributed around azimuthal direction ab. fieldline 
  // theta ~ U[0, 2 pi]
  theta = G4UniformRand() * 2. * fPI;  
  
  // Phi (half angle) takes values in a cone determined by 
  // spacecraft altitude, precipitation event altitude and size
  G4double photonPhiLimitRad = photonPhiLimitDeg * fDeg2Rad;       
  
  // u ~ U[0, 1-1/2*cos(phi_limit)] s.t. phi is distributed uniformly
  G4double u = G4UniformRand()*(1.-std::cos(photonPhiLimitRad))/2.;
  
  // phi ~ U[0, phi_limit] = cos^-1(1 - 2 u)
  G4double phi = std::acos(1 - 2 * u);


  // We want our Y direction to be "up"
  r->x = (sphereR + 15.*cm) * std::sin(phi) * std::sin(theta);
  r->y = (sphereR + 15.*cm) * std::cos(phi);
  r->z = (sphereR + 15.*cm) * std::sin(phi) * std::cos(theta);


  // Inward unit normal
  r->xDir = -std::sin(theta) * std::sin(phi);
  r->zDir = -std::cos(theta) * std::sin(phi); 
  r->yDir = -std::cos(phi);


  // Uniform nudges s.t. each particle isn't aimed at the origin
  G4double nudgeFactor = 0.3;

  r->xDir += (G4UniformRand()*2. - 1.) * nudgeFactor; 
  r->yDir += (G4UniformRand()*2. - 1.) * nudgeFactor; 
  r->zDir += (G4UniformRand()*2. - 1.) * nudgeFactor; 
}


void PrimaryGeneratorAction::GenerateTrappedElectrons(ParticleSample* r)
{
  
    // Loss cone angle (same as polar angle, phi) at 500 km, in radians
    G4double theta_exclusion = 64.*fPI/180.;

    // Calculate random particle position on sphere via rejection 
    // sampling, excluding spherical cap that makes up the loss cone
    do{
      // Rand on [-1, 1)
      G4double u = G4UniformRand()*2.-1.;

      // Rand on [0, 2*pi)
      G4double theta = G4UniformRand()*2.*fPI;
      r->x = sphereR * std::sqrt(1 - u * u) * std::cos(theta); 
      r->y = sphereR * u; 
      r->z = sphereR * std::sqrt(1 - u * u) * std::sin(theta); 
      }
    while(r->y > sphereR * std::cos(theta_exclusion));
    // exits when y position falls below spherical cap

    // Uniform random numbers on [0, 1)
    r->xDir = G4UniformRand();
    r->yDir = G4UniformRand();
    r->zDir = G4UniformRand();

    
    // Enforces inward directionality to particles
    if(r->x > 0) {r->xDir = -r->xDir;}
    if(r->y > 0) {r->yDir = -r->yDir;}
    if(r->z > 0) {r->zDir = -r->zDir;}


    // Selects random energy according to exponential distribution
    G4double randomNumber = G4UniformRand();
    r->energy = ((std::log(1 - randomNumber)*-E_folding))*keV;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Method called to populate member variables with number of 
  // particles to generate
  //CalculateParticlesToGenerate();

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(electronParticle);

  // Constant sphere offsets
  G4double xShift = 0.;
  G4double yShift = 1.*cm;
  G4double zShift = 2.5*cm;

  
  // Struct that holds position, momentum direction, and energy
  ParticleSample* r = new ParticleSample();
 
  switch(fWhichParticle){
    case(0): // Background electron, outside loss cone

    GenerateTrappedElectrons(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    break;

    case(1): // Background electron, inside loss cone

    GenerateLossConeElectrons(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    break;

    case(2): // Signal photon

    // Selects photon for particle type
    fParticleGun->SetParticleDefinition(photonParticle);
  
    GenerateSignalSource(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    break;

    default: 
       throw std::invalid_argument("Need to chose particle type!");
  }

  // Free particle utility struct
  delete r;

}

