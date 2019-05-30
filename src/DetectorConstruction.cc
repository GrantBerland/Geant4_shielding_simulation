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

  fParser.Read("./top_v1.gdml", validateSchema);
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




  // CZT for detector
  G4Element* Cd = new G4Element("Cadmium","Cd",48., 112.41*g/mole);
  G4Element* Zn = new G4Element("Zinc","Zn", 30., 65.38*g/mole);
  G4Element* Te = new G4Element("Tellurium","Te", 52., 127.60*g/mole);
  G4Material* CZT = new G4Material("CZT", 5.8*g/cm3, 3);
  CZT->AddElement(Cd, 48*perCent);
  CZT->AddElement(Zn, 2*perCent);
  CZT->AddElement(Te, 50*perCent);

  G4int natoms;
  G4double z, a;

  G4Element* O  = new G4Element("Oxygen"  ,"O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silicon","Si" , z= 14., a= 28.09*g/mole);
  G4Element* H  = new G4Element("Hydrogen","H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,"C" , z= 6., a= 12.01*g/mole);

  // Si02 for silicon
  G4Material* SiO2 = new G4Material("quartz", 2.200*g/cm3, 2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);


  //Epoxy (for FR4 )
  G4Material* Epoxy = new G4Material("Epoxy" , 1.2*g/cm3,2);
  Epoxy->AddElement(H, natoms=2);
  Epoxy->AddElement(C, natoms=2);
  
  //FR4 (Glass + Epoxy)
  G4Material* FR4 = new G4Material("FR4",1.86*g/cm3, 2);
  FR4->AddMaterial(SiO2, 0.528);
  FR4->AddMaterial(Epoxy, 0.472);




  for(unsigned int i = 0; i < solids->size(); i++){
    G4VSolid* psol = (*solids)[i];
    G4cout << "Solid ID: " << i << " Name: " << psol->GetName() << G4endl;
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


  new G4PVPlacement(0,                     //no rotation
                  G4ThreeVector(0.,0.,-10.*cm),                 //at position
                  logicShielding,          //its logical volume
		  "Shielding",           //its name
		  logicEnv,                //its mother  volume
		  false,                   //no boolean operation:u
	 	  0,                       //copy number
		  checkOverlaps);          //overlaps checking


  G4LogicalVolume* logicWindows = new G4LogicalVolume((*solids)[2],
					nist->FindOrBuildMaterial("G4_Be"),
							     "Windows");
  new G4PVPlacement(0,                       //no rotation
	            G4ThreeVector(0.,0.,0.), //its logical volume
	            logicWindows,            //at position
	   	    "Windows",               //its name
	 	    logicShielding,                //its mother  volume
		    false,                   //no boolean operation:u
		    0,                       //copy number
		    checkOverlaps);          //overlaps checking


  G4LogicalVolume* logicDetectorTop = new G4LogicalVolume((*solids)[3],
								CZT,
							      "DetectorCZT");
  new G4PVPlacement(0,                     //no rotation
	 	   G4ThreeVector(0.,0.,0.),   //at position
	 	   logicDetectorTop,          //its logical volume
		   "DetectorCZT",              //its name
		   logicShielding,         //its mother  volume
		   false,                  //no boolean operation:u
		    0,                     //copy number
	  	   checkOverlaps);         //overlaps checking
	 
  G4LogicalVolume* logicDetectorBottom = new G4LogicalVolume((*solids)[4],
								FR4,
							      "DetectorFR4");
  new G4PVPlacement(0,                     //no rotation
	 	   G4ThreeVector(0.,0.,0.),   //at position
	 	   logicDetectorBottom,          //its logical volume
		   "DetectorFR4",              //its name
		   logicShielding,
		   false,                  //no boolean operation:u
		    0,                     //copy number
	  	   checkOverlaps);         //overlaps checking
 
  G4LogicalVolume* logicBaffles = new G4LogicalVolume((*solids)[5],
						nist->FindOrBuildMaterial("G4_W"),
								      "Baffles");


  new G4PVPlacement(0,                     //no rotation
		  G4ThreeVector(0.,0.,0.),                 //at position
		  logicBaffles,          //its logical volume
		  "Baffles",           //its name
		  logicShielding,                //its mother  volume
		  false,                   //no boolean operation:u
		  0,                       //copy number
		  checkOverlaps);          //overlaps checking



  fScoringVolume = logicDetectorTop;

  // Register the detector as a sensitive detector
  B2TrackerSD* aTrackerSD = new B2TrackerSD("DetectorCZT","TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  SetSensitiveDetector("DetectorCZT", aTrackerSD, true);


  // always return the physical World
  return physWorld;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
