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
// $Id: DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
//#include "G4Tubs.hh"
//#include "G4Cons.hh"
//#include "G4Orb.hh"
//#include "G4Sphere.hh"
//#include "G4GenericPolycone.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
//#include "G4IntersectionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
#include "B2TrackerSD.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 50.*cm, env_sizeZ = 50.*cm;

  G4bool checkOverlaps    = true;
  G4bool validateSchema   = false;

  fParser.Read("./top.gdml", validateSchema);
  G4SolidStore* solids = G4SolidStore::GetInstance();


    // Material: Vacuum
    //TODO: check pressures, environment for Van Allen belt altitudes
  G4Material* vacuum_material = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );
  

  G4LogicalVolume* logicWorld = new G4LogicalVolume((*solids)[0],
							     vacuum_material,
							     "World");



  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        vacuum_material,     //its material
                        "Envelope");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
checkOverlaps); //overlaps checking




  // Option to switch on/off checking of volumes overlaps
  //



  //G4cout << *(G4Material::GetMaterialTable() ) << G4endl;


  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size


  // CZT for detector
  G4Element* Cd = new G4Element("Cadmium","Cd",48., 112.41*g/mole);
  G4Element* Zn = new G4Element("Zinc","Zn", 30., 65.38*g/mole);
  G4Element* Te = new G4Element("Tellurium","Te", 52., 127.60*g/mole);
  G4Material* CZT = new G4Material("CZT", 5.8*g/cm3, 3);
  CZT->AddElement(Cd, 48*perCent);
  CZT->AddElement(Zn, 2*perCent);
  CZT->AddElement(Te, 50*perCent);


  G4ThreeVector ele_pos  = G4ThreeVector(40.*mm,40.*mm,28.*mm);
  G4VSolid* electronics_box = new G4Box("Electronics",
                                    1.5*cm, 1.5*cm, 0.1*cm);




  G4LogicalVolume* electronics_LV =
  new G4LogicalVolume(electronics_box,                     //its solid
                      CZT,  //its material
                      "Detector1");                        //its name
  G4RotationMatrix* rotmatrix = new G4RotationMatrix();
  rotmatrix->rotateX(90.*deg);


  new G4PVPlacement(rotmatrix,                     //no rotation
                  ele_pos,                 //at position
                  electronics_LV,          //its logical volume
		  "Detector1",           //its name
		  logicEnv,                //its mother  volume
		  false,                   //no boolean operation:u
	 	  0,                       //copy number
		  checkOverlaps);          //overlaps checking




	  for(unsigned int i = 0; i < solids->size(); i++){
	    G4VSolid* psol = (*solids)[i];
	    std::cout << "Solid ID: " << i << " Name: " << psol->GetName() << std::endl;
	  };




  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  G4LogicalVolume* logicShielding = new G4LogicalVolume((*solids)[1],
					nist->FindOrBuildMaterial("G4_Al"),
							     "Shielding");

  G4LogicalVolume* logicWindows = new G4LogicalVolume((*solids)[2],
					nist->FindOrBuildMaterial("G4_Al"),
							     "Windows");

  G4LogicalVolume* logicDetector = new G4LogicalVolume((*solids)[3],
								CZT,
							      "Detector");

  G4LogicalVolume* logicBaffles = new G4LogicalVolume((*solids)[4],
					nist->FindOrBuildMaterial("G4_Al"),
							      "Detector");



  G4LogicalVolume* shielding_LV = new G4LogicalVolume((*solids)[0],
							     CZT,
							     "Shielding");


  new G4PVPlacement(0,
                    G4ThreeVector(-14.*cm,-10.*cm,-11.*cm),
                    shielding_LV,
                    "Shielding",
                    logicEnv,
                    false,
                    0,
                    checkOverlaps);

  fScoringVolume = electronics_LV;

  // Register the detector as a sensitive detector
  B2TrackerSD* aTrackerSD = new B2TrackerSD("shielding_sim/Detector1","TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  SetSensitiveDetector("Detector1", aTrackerSD, true);


  // always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
