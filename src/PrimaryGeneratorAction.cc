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

#include <fstream>
#include <stdexcept>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  E_folding(100.),
  E_shift(0.),
  fPI(3.14159265358979323846),
  sphereR(24.*cm),
  lossConeAngleDeg(64.),
  photonPhiLimitDeg(40.),
  multModifier(0.),
  electronParticle(0),
  photonParticle(0),
  nBackgroundElectrons(0),
  nLossConeElectrons(0),
  nSignalPhotons(0),
  photonEnergyProb100keV{
0.03656131,0.25322138,0.41933641,0.54648931,
0.64401902,0.71941948,0.77817282,0.82384164,
0.85938710,0.88718833,0.90880856,0.92570210,
0.93902379,0.94978780,0.95852530,0.96564760,
0.97141985,0.97609667,0.97988323,0.98306097,
0.98572875,0.98799985,0.98983837,0.99134501,
0.99263692,0.99369583,0.99459389,0.99533782,
0.99598452,0.99654148,0.99702781,0.99745401,
0.99781161,0.99810876,0.99835896,0.99858048,
0.99876756,0.99891555,0.99904548,0.99915724,
0.99925246,0.99933265,0.99939845,0.99945309,
0.99950227,0.99955484,0.99960124,0.99964836,
0.99969242,0.99972854,0.99975926,0.99978796,
0.99981340,0.99983692,0.99986030,0.99988102,
0.99989542,0.99990494,0.99991505,0.99992416,
0.99993166,0.99993927,0.99994694,0.99995436,
},
  photonEnergyProb200keV{
0.03066625,0.21750923,0.36495406,0.48218012,
0.57618225,0.65148647,0.71196796,0.76104236,
0.80088337,0.83345108,0.85997623,0.88170227,
0.89957654,0.91442320,0.92685580,0.93737503,
0.94619513,0.95361777,0.95986459,0.96517422,
0.96969009,0.97357224,0.97689726,0.97971784,
0.98216937,0.98430482,0.98616823,0.98779024,
0.98920731,0.99047318,0.99159633,0.99257494,
0.99342281,0.99416446,0.99482191,0.99540139,
0.99591614,0.99636694,0.99676320,0.99712006,
0.99743749,0.99771491,0.99795808,0.99817147,
0.99836162,0.99853478,0.99868415,0.99882250,
0.99894752,0.99905448,0.99914777,0.99923378,
0.99931322,0.99938551,0.99945173,0.99951425,
0.99956924,0.99961657,0.99966119,0.99970238,
0.99973941,0.99977286,0.99980319,0.99983036,
},
  photonEnergyProb300keV{
0.02818334,0.20200644,0.34211434,0.45510159,
0.54695364,0.62176885,0.68278631,0.73297487,
0.77433868,0.80858439,0.83694915,0.86060972,
0.88040922,0.89705892,0.91117350,0.92316284,
0.93331228,0.94196522,0.94931855,0.95558576,
0.96096569,0.96561659,0.96965985,0.97314387,
0.97618687,0.97888122,0.98126791,0.98335913,
0.98519104,0.98680357,0.98824378,0.98952201,
0.99063815,0.99160880,0.99246592,0.99323222,
0.99392045,0.99453989,0.99508919,0.99558152,
0.99601587,0.99640281,0.99675101,0.99706151,
0.99734280,0.99759477,0.99781584,0.99802117,
0.99821001,0.99837820,0.99852953,0.99866966,
0.99879825,0.99891640,0.99902618,0.99912726,
0.99921862,0.99930154,0.99937869,0.99944781,
0.99950865,0.99956353,0.99961266,0.99965665,
}
{

  fParticleGun  = new G4ParticleGun();
  
  electronParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  
  photonParticle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::CalculateParticlesToGenerate()
{

  // Convert sphere radius from Geant internal units to cm
  G4double sphereRcm = sphereR / cm;  
  
  // Surface area of sphere in cm^2
  G4double sphereAreacm2 = 4. * fPI * sphereRcm * sphereRcm; 
  
  // Approximate flux at L = 3.6, 500 km altitude at high latitudes
  G4double trappedElectronFlux = 1e5;      // el/cm^2/s/str

  // Converts loss cone angle from degrees to radians
  G4double lossConeAngleRad = lossConeAngleDeg * fPI / 180.;
  
  // Calculates pi * r^2 for photon flux calculations, using
  // r = R * sin(phi), where R is the radius of the sphere
  G4double sphereCrossSectionalArea = 
	  (sphereRcm  * std::sin(photonPhiLimitDeg * fPI / 180.)) *
	  (sphereRcm  * std::sin(photonPhiLimitDeg * fPI / 180.)) * fPI;

  // Calculation of solid angle based on cone half angle (loss cone angle)
  G4double trappedElectronSolidAngle =      // str  
	  2 * fPI * (1 + std::cos(lossConeAngleRad));    

  // Trapped electrons to generate from flux calculation 
  nBackgroundElectrons =                                 //[el/s]
	  (unsigned long long int)(trappedElectronFlux * //[el/cm^2/s/str]
			           sphereAreacm2 *       //[cm^2]
			           trappedElectronSolidAngle); //[str]


  // Loss cone particles (backscattered), fractional flux derived 
  // from Marshall, Bortnik work
  // N_{loss cone particles}/N_{total particles}
  G4double lossConeParticleFraction = 0.316;
 
  // gamma = ^ above fraction
  // N_{loss cone parts} = gamma * N_{trapped parts}/(1 - gamma)
  nLossConeElectrons = 
	  (unsigned long long int)(lossConeParticleFraction*nBackgroundElectrons/(1-lossConeParticleFraction));	


  // Photon flux [ph/cm^2/s] per each E_0 folding energy range (in keV)
  if      (E_folding <= 100.) nSignalPhotons = 67;
  else if (E_folding <= 200.) nSignalPhotons = 302;
  else if (E_folding <= 300.) nSignalPhotons = 690;
  else throw std::invalid_argument("Non-realizable E_0");

  // Converts [ph/cm^2/s] to [ph/s] through a circle the size of the 
  // generation area
  nSignalPhotons *= sphereCrossSectionalArea;
  
  // Optional particle generation numbers read in from file
  unsigned long long int source1, source2, source3;
  int debugFlag;
  // Reads in multiplicative factor to evenly reduce the number of 
  // particles generated, as well as optional source counts
  std::fstream particleNumberMultiplierFile;
  particleNumberMultiplierFile.open("numberOfParticles.txt", std::ios_base::in);
  particleNumberMultiplierFile >> multModifier >> source1
	  >> source2 >> source3 >> debugFlag;
  particleNumberMultiplierFile.close();

  // Runs imposed number of particles if nonzero integer in file 
  // else, runs the calculated number of particles
  if(source1 != 0){
    nBackgroundElectrons = source1;
  }
  if(source2 != 0){
    nLossConeElectrons   = source2;
  }  
  if(source3 != 0){
    nSignalPhotons       = source3;
  }


  // Reduces number of particles by modifier in order to be simulatable 
  nBackgroundElectrons *= multModifier;
  nLossConeElectrons   *= multModifier;
  nSignalPhotons       *= multModifier;

  G4cout << "Background Electrons: " << nBackgroundElectrons 
	  << "\nLoss Cone Electrons: " << nLossConeElectrons 
	  << "\nSignal Photons: " << nSignalPhotons << G4endl;

  if(debugFlag == 1){
    // For debugging
    throw std::invalid_argument("Stop here!");
  }

}


void PrimaryGeneratorAction::GenerateLossConeElectrons(ParticleSample* r)
{
  // This method generates electrons that were in the loss cone, but have
  // backscattered off of the atmosphere and are direction anti-Earthward
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
  for(G4int angle = 0; angle<dataSize; angle++){
    if(randomNumber < lossConeData[angle][1]){
      angleIndex = angle;
      break; }
  }


  // Selects exponential folding energy E0 based on backscattered 
  // pitch angle range
  G4double E0 = -1;
  if(     lossConeData[angleIndex][0] < 20) {E0 = 159.;}
  else if(lossConeData[angleIndex][0] < 30) {E0 = 141.;}
  else if(lossConeData[angleIndex][0] < 40) {E0 = 177.;}
  else if(lossConeData[angleIndex][0] < 50) {E0 = 194.;}
  else if(lossConeData[angleIndex][0] < 64) {E0 = 230.;}
  

  // N.B.: Mathematics spherical coordinates definition used below
  
  // Uniform sampling of angle about the field line theta on [0 , 2pi)
  G4double theta = G4UniformRand()*2.*fPI; 
  
  // Zenith angle (pitch angle), converted to radians
  G4double phi = lossConeData[angleIndex][0] * fPI / 180.;
  
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
  r->energy = ((std::log(1 - randomNumber)*-E0 + E_shift))*keV;
}

void PrimaryGeneratorAction::GenerateSignalSource(ParticleSample* r)
{
  // This method generates a photon signal from a energetic particle
  // precipitation event that occurs at approximately 50 - 100 km 
  // altitude. Variables are photon energy and the 
  // pitch, or zenith angle distribution of the photons.

  G4double theta, phi, u, randomNumber;
  static const G4double photonEnergyArray[64] = 
{50.00000,57.14286,64.28571,71.42857,
78.57143,85.71429,92.85714,100.00000,
107.14286,114.28571,121.42857,128.57143,
135.71429,142.85714,150.00000,157.14286,
164.28571,171.42857,178.57143,185.71429,
192.85714,200.00000,207.14286,214.28571,
221.42857,228.57143,235.71429,242.85714,
250.00000,257.14286,264.28571,271.42857,
278.57143,285.71429,292.85714,300.00000,
307.14286,314.28571,321.42857,328.57143,
335.71429,342.85714,350.00000,357.14286,
364.28571,371.42857,378.57143,385.71429,
392.85714,400.00000,407.14286,414.28571,
421.42857,428.57143,435.71429,442.85714,
450.00000,457.14286,464.28571,471.42857,
478.57143,485.71429,492.85714,500.00000,
};
  const G4int dataSize = 64;
  const G4double* photonEnergyTablePointer;
  randomNumber = G4UniformRand();
 
  // If statements to determine which photon energy probability table to 
  // use, based on the E_0 folding energy of the background electrons
  if(E_folding <= 100.){
    photonEnergyTablePointer = photonEnergyProb100keV;}
  else if(E_folding <= 200.){
    photonEnergyTablePointer = photonEnergyProb200keV;}
  else if(E_folding <= 300.){
    photonEnergyTablePointer = photonEnergyProb300keV;}
  else {throw std::invalid_argument("Source folding energy not in {100,200,300} keV");}

    // Discrete inverse CDF lookup of probabilites, which are then
    // linked to an energy
    for(G4int energyIndex=0; energyIndex<dataSize; energyIndex++)
    {
      if(randomNumber < photonEnergyTablePointer[energyIndex])
      {
        r->energy = photonEnergyArray[energyIndex]*keV;
        break;
      }
    }

  
  // Uniformly distributed around field line 
  theta = G4UniformRand()*2.*fPI;  
  
  // Phi (half angle) takes values in a cone determined by 
  // spacecraft altitude, precipitation event altitude and size
  G4double photonPhiLimitRad   = photonPhiLimitDeg * fPI / 180.;       
  u = G4UniformRand()*(1.-std::cos(photonPhiLimitRad))/2.;
  phi = std::acos(1 - 2 * u);

  // We want our Y direction to be "up"
  /*
  r->x = sphereR * std::sin(phi) * std::cos(theta);
  r->y = sphereR * std::cos(phi);
  r->z = sphereR * std::sin(phi) * std::sin(theta);
  */
  G4double radius = sphereR * std::sin(photonPhiLimitRad) 
	             * std::sqrt(G4UniformRand());
  r->x = radius * std::sin(theta);
  r->y = sphereR;
  r->z = radius * std::cos(theta);


  // Uniform distributed downwards, towards detector
  r->xDir = 0;
  r->yDir = -1;
  r->zDir = 0;

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
    r->energy = ((std::log(1 - randomNumber)*-E_folding + E_shift))*keV;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Method called to populate member varaibles with number of 
  // particles to generate
  CalculateParticlesToGenerate();

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(electronParticle);

  // Constant sphere offsets
  G4double xShift = 0.;
  G4double yShift = 0.;
  G4double zShift = 0.;

  
  // Struct that holds position, momentum direction, and energy
  ParticleSample* r = new ParticleSample();
  
  for(unsigned long long int i = 0; i<nBackgroundElectrons; i++){

    GenerateTrappedElectrons(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
  
  for(unsigned long long int i = 0; i<nLossConeElectrons; i++){

    GenerateLossConeElectrons(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

  // Selects photon for particle type
  fParticleGun->SetParticleDefinition(photonParticle);
  
  for(unsigned long long int i = 0; i<nSignalPhotons; i++){

    GenerateSignalSource(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}


  delete r;

}

