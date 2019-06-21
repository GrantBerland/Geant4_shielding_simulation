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
#include "G4SubtractionSolid.hh"
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
  G4double env_sizeXY = 40.*cm, env_sizeZ = 40.*cm;

  G4bool checkOverlaps    = true;
  
  // Material: Vacuum
  G4Material* vacuum_material = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );
  
  G4Box* worldBox = new G4Box("World", 0.6*env_sizeXY, 0.6*env_sizeXY, 0.6*env_sizeZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(worldBox,
						     vacuum_material,
						     "World");

  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,
		  			     G4ThreeVector(),
					     logicWorld,
					     "World",
					     0,
					     false,
					     0);

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

  //////////////////////////////////////////
  /////////// Detector Construction ////////
  //////////////////////////////////////////

  // 0.5 factor due to full height definition
  G4double outerShieldingThickness = 1.5*mm * 0.5;  
  G4double innerShieldingThickness = 3.0*mm * 0.5;
  G4double detectorXY      = 40.*mm;
  G4double detectorZ       = 5.*mm;
  G4double detectorElectronicsZ = 10.185*mm;
  G4double boxInnerSizeXY  = 90.*mm * 0.5;
  G4double boxInnerSizeZ   = 50.*mm * 0.5;
  G4double windowThickness = 1.*mm;
  G4double baffleHeight    = 20.*mm;
  G4double baffleThickness = 0.5*mm;
  G4double frontEndBoardThickness = 2.86*mm;


  /////////////////////////////////////////
  //////////////// Geometry ///////////////
  /////////////////////////////////////////

  G4VSolid* detectorBox = new G4Box("Detector",
		  		    0.5*detectorXY,
				    0.5*detectorZ,
				    0.5*detectorXY);

  G4VSolid* detectorElectronics = new G4Box("DetectorFR4",
		  		    0.5*detectorXY,
				    0.5*detectorElectronicsZ,
				    0.5*detectorXY);

  G4VSolid* frontEndBoard = new G4Box("Electronics",
		  		    boxInnerSizeXY,
				    0.5*frontEndBoardThickness,
				    boxInnerSizeXY);


  G4VSolid* berylliumWindow = new G4Box("Be_Window",
		  		    0.5*15.*mm,
				    0.5*windowThickness,
				    boxInnerSizeXY);
  
  G4VSolid* topWindow = new G4Box("top_Be_Window",
		  		    boxInnerSizeXY,
				    0.5*windowThickness,
				    boxInnerSizeXY);

  G4VSolid* baffles = new G4Box("Baffles",
		  		    0.5*baffleThickness,
				    0.5*baffleHeight,
				    boxInnerSizeXY);

  G4VSolid* collimeterOuterBox = new G4Box("Outer-collimeter",
boxInnerSizeXY+2*innerShieldingThickness+2*outerShieldingThickness,
25.*mm * 0.5,
boxInnerSizeXY+2*innerShieldingThickness+2*outerShieldingThickness);

  G4VSolid* outerShieldingBox = new G4Box("Outer-shielding",
boxInnerSizeXY+2*innerShieldingThickness+2*outerShieldingThickness,
boxInnerSizeZ+2*innerShieldingThickness+2*outerShieldingThickness,
boxInnerSizeXY+2*innerShieldingThickness+2*outerShieldingThickness);


  G4VSolid* innerShieldingBox = new G4Box("Inner-shielding",
	  boxInnerSizeXY+2*innerShieldingThickness,
	  boxInnerSizeZ+2*innerShieldingThickness,
	  boxInnerSizeXY+2*innerShieldingThickness);


  G4VSolid* slit1 = new G4Box("Slit",
		  	1.1*mm,
			outerShieldingThickness+1.*cm,
			boxInnerSizeXY-0.5*mm);
  G4VSolid* slit2 = new G4Box("Slit",
		  	1.1*mm,
			innerShieldingThickness+1.*cm,
			boxInnerSizeXY-0.5*mm);

  G4VSolid* subtractionBox = new G4Box("Subtraction-box",
    		        boxInnerSizeXY,
			boxInnerSizeZ,
			boxInnerSizeXY);
  
  
  //////////////////////////////////////////
  ///////////// Subtractions ///////////////
  //////////////////////////////////////////


  // empty rotation matrix for SubtractionSolid constructor
  G4RotationMatrix* rotm = new G4RotationMatrix();   
  
  G4double windowPlacement = boxInnerSizeZ + innerShieldingThickness + outerShieldingThickness + windowThickness*1.5;

  G4double windowSlitShiftX = 2.5*cm;

  // Hollows out collimeter
  G4SubtractionSolid* collimeterOuter = 
	  new G4SubtractionSolid("Outer-collimeter",
	  collimeterOuterBox,
	  subtractionBox);

  // Hollows out outer shielding (tungsten)
  G4SubtractionSolid* shieldingBoxOuter = 
	  new G4SubtractionSolid("W_Shielding",
	  outerShieldingBox,
	  innerShieldingBox);


  // Outer shielding window depressions 1 & 2 and slits 1 & 2
  G4SubtractionSolid* shieldingBoxOuter_slit1 = 
	  new G4SubtractionSolid("W_Shielding",
	  shieldingBoxOuter,
	  slit1,
	  rotm,
	  G4ThreeVector(windowSlitShiftX, boxInnerSizeZ + outerShieldingThickness, 0.));
  
  G4SubtractionSolid* shieldingBoxOuter_slit2 = 
	  new G4SubtractionSolid("W_Shielding",
	  shieldingBoxOuter_slit1,
	  slit1,
	  rotm,
	  G4ThreeVector(-windowSlitShiftX, boxInnerSizeZ + outerShieldingThickness, 0.));


  G4SubtractionSolid* shieldingBoxOuter_slit_window1 = 
	  new G4SubtractionSolid("W_Shielding",
	  shieldingBoxOuter_slit2,
	  berylliumWindow,
	  rotm,
	  G4ThreeVector(windowSlitShiftX, windowPlacement, 0.));

  G4SubtractionSolid* shieldingBoxOuter_slit_window2 = 
	  new G4SubtractionSolid("W_Shielding",
	  shieldingBoxOuter_slit_window1,
	  berylliumWindow,
	  rotm,
	  G4ThreeVector(-windowSlitShiftX, windowPlacement, 0.));

  // Hollows out inner shielding (aluminum)
  G4SubtractionSolid* shieldingBoxInner = 
	  new G4SubtractionSolid("Al_Shielding",
	  innerShieldingBox,
	  subtractionBox);

  // Inner shielding slits 1 & 2
  G4SubtractionSolid* shieldingBoxInner_slit1 = 
	  new G4SubtractionSolid("Al_Shielding",
	  shieldingBoxInner,
	  slit2,
	  rotm,
	  G4ThreeVector(windowSlitShiftX, boxInnerSizeZ + innerShieldingThickness, 0.));
  
  G4SubtractionSolid* shieldingBoxInner_slit2 = 
	  new G4SubtractionSolid("Al_Shielding",
	  shieldingBoxInner_slit1,
	  slit2,
	  rotm,
	  G4ThreeVector(-windowSlitShiftX, boxInnerSizeZ + innerShieldingThickness, 0.));
  

  ////////////////////////////////////////////
  ////////////// Logical Volumes /////////////
  ////////////////////////////////////////////

  G4LogicalVolume* logicCollimeter = new G4LogicalVolume(collimeterOuter,
				   nist->FindOrBuildMaterial("G4_Al"),
							"Collimeter");
  

  G4LogicalVolume* logicDetector = new G4LogicalVolume(detectorBox,
							CZT,
							"Detector");
  
  G4LogicalVolume* logicDetectorElectronics = new G4LogicalVolume(detectorElectronics,
							FR4,
							"DetectorFR4");
  
  G4LogicalVolume* logicFrontEndBoard = new G4LogicalVolume(frontEndBoard,
							FR4,
							"Electronics");
  
  G4LogicalVolume* logicalOuterShielding = new G4LogicalVolume(shieldingBoxOuter_slit_window2,
		  nist->FindOrBuildMaterial("G4_W"),
		  "W_Shielding");
  
  G4LogicalVolume* logicalInnerShielding = new G4LogicalVolume(shieldingBoxInner_slit2,
		  nist->FindOrBuildMaterial("G4_Al"),
		  "Al_Shielding");
  
  G4LogicalVolume* logicalWindow = new G4LogicalVolume(berylliumWindow,
		  nist->FindOrBuildMaterial("G4_Be"),
		  "Be_Window");
  
  G4LogicalVolume* logicalTopWindow = new G4LogicalVolume(topWindow,
		  nist->FindOrBuildMaterial("G4_Be"),
		  "top_Be_Window");
  
  G4LogicalVolume* logicalBaffles = new G4LogicalVolume(baffles,
		  nist->FindOrBuildMaterial("G4_W"),
		  "Baffles");


  ////////////////////////////////////////////////
  ///////////////// Placements ///////////////////
  ////////////////////////////////////////////////
  
  new G4PVPlacement(0,
		    G4ThreeVector(0.,boxInnerSizeZ+16.15*mm+855*um,0.),
		    logicCollimeter,
		    "Collimeter",
		    logicEnv,
		    false,
		    checkOverlaps);
 
  G4double detectorPosX = -22.5*mm;
  G4double detectorPosZ = -21.*mm;
  G4double detectorPosY = frontEndBoardThickness + 0.5*cm;

  // position multiplier
  G4int pm1[4] = {1, -1, 1, -1};
  G4int pm2[4] = {1, 1, -1, -1};


  for(G4int nDet=0; nDet<4;nDet++)
  {

    new G4PVPlacement(0,
		    G4ThreeVector(pm1[nDet]*detectorPosX,
			    -1.25*cm+0.1*mm+detectorPosY,
			    pm2[nDet]*detectorPosZ),
		    logicDetector,
		    "Detector",
		    logicEnv,
		    false,
		    nDet,
		    checkOverlaps);

    new G4PVPlacement(0,
		    G4ThreeVector(pm1[nDet]*detectorPosX,
			    -2.*cm+detectorPosY,
			    pm2[nDet]*detectorPosZ),
		    logicDetectorElectronics,
		    "DetectorFR4",
		    logicEnv,
		    false,
		    nDet,
		    checkOverlaps);
  } 

  new G4PVPlacement(0,
		    G4ThreeVector(0.,-2.0*cm,0.),
		    logicFrontEndBoard,
		    "Electronics",
		    logicEnv,
		    false,
		    checkOverlaps);
  
  new G4PVPlacement(0,
		  G4ThreeVector(),
		  logicalOuterShielding,
		  "W_Shielding",
		  logicEnv,
		  true,
		  checkOverlaps);
  
  new G4PVPlacement(0,
		  G4ThreeVector(),
		  logicalInnerShielding,
		  "Al_Shielding",
		  logicEnv,
		  true,
		  checkOverlaps);
  
  new G4PVPlacement(0,
		  G4ThreeVector(windowSlitShiftX, windowPlacement,0.),
		  logicalWindow,
		  "Be_Window",
		  logicEnv,
		  false,
		  0,
		  checkOverlaps);

  new G4PVPlacement(0,
		  G4ThreeVector(-windowSlitShiftX, windowPlacement,0.),
		  logicalWindow,
		  "Be_Window",
		  logicEnv,
		  false,
		  1,
		  checkOverlaps);


  new G4PVPlacement(0,
		  G4ThreeVector(0., windowPlacement+20.*mm+2.*mm,0.),
		  logicalTopWindow,
		  "top_Be_Window",
		  logicEnv,
		  false,
		  checkOverlaps);
  
  // Baffle parameterisation
  G4int numberBaffles = 28;

  G4double bafflePlacement = windowPlacement+0.5*baffleHeight+outerShieldingThickness+0.3*mm;
  
  G4double axialDistance;
  rotm->rotateY(90.*deg); 

  for(G4int i = 0; i<numberBaffles; i++)
  {
    // baffles are spaced 2.53 mm apart
    axialDistance = (i-14) * (2.53 + 0.5) * mm;

    new G4PVPlacement(rotm,
		     G4ThreeVector(0., bafflePlacement, axialDistance),
		     logicalBaffles,
		     "Baffles",
		     logicEnv,
		     false,
		     i,
		     false);
  }


  // always return the physical World
  return physWorld;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
