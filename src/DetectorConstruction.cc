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
#include "G4Sphere.hh"
//#include "G4GenericPolycone.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
//#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4AssemblyVolume.hh"


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
  G4double env_sizeXY = 85.*cm; 

  G4bool checkOverlaps    = true;
  
  // Material: Vacuum
  G4Material* vacuum_material = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );
  
  G4Sphere* worldBox = new G4Sphere("World", 
		  0.*cm, 0.6*env_sizeXY, 
		  0.*deg , 360.*deg,
		  0.*deg , 180.*deg);

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

  G4Sphere* solidEnv =
    new G4Sphere("Envelope",                    //its name
		  0.*cm, 0.5*env_sizeXY, 
		  0.*deg , 360.*deg,
		  0.*deg , 180.*deg);

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
  G4double outerShieldingThickness = 15.*mm * 0.5; // PE 
  G4double shieldingThickness2     = 3.5*mm * 0.5; // W
  G4double shieldingThickness3     = 2.5*mm * 0.5; // Sn
  G4double detectorXY      = 40.*mm;
  G4double detectorZ       = 5.*mm;
  G4double detectorElectronicsZ = 10.185*mm;
  G4double boxInnerSizeXY  = 90.*mm * 0.5;
  G4double boxInnerSizeZ   = 50.*mm * 0.5;
  G4double windowThickness = 1.*mm;
  G4double frontEndBoardThickness = 2.86*mm;
  G4double detectorApertureSpacing = 20.*mm;

  /////////////////////////////////////////
  //////////////// Geometry ///////////////
  /////////////////////////////////////////

  G4VSolid* detectorBox = new G4Box("Detector",
		  		    0.5*detectorXY,
				    0.5*detectorZ,
				    0.5*detectorXY);

  G4VSolid* detectorElectronics = new G4Box("BottomFR4",
		  		    0.5*detectorXY,
				    0.5*detectorElectronicsZ,
				    0.5*detectorXY);

  G4VSolid* frontEndBoard = new G4Box("Electronics",
		  		    boxInnerSizeXY,
				    0.5*frontEndBoardThickness,
				    boxInnerSizeXY);

  G4VSolid* topWindow = new G4Box("top_Be_Window",
		  		    boxInnerSizeXY,
				    0.5*windowThickness,
				    boxInnerSizeXY);


  // Shielding boxes
  G4VSolid* outerShieldingBox = new G4Box("PE_Shielding",
boxInnerSizeXY+2*outerShieldingThickness+2*shieldingThickness2+2*shieldingThickness3,
boxInnerSizeZ+2*outerShieldingThickness+2*shieldingThickness2+2*shieldingThickness3,
boxInnerSizeXY+2*outerShieldingThickness+2*shieldingThickness2+2*shieldingThickness3);


  G4VSolid* shieldingBox2 = new G4Box("W_Shielding",
boxInnerSizeXY+2*shieldingThickness2+2*shieldingThickness3,
boxInnerSizeZ+2*shieldingThickness2+2*shieldingThickness3,
boxInnerSizeXY+2*shieldingThickness2+2*shieldingThickness3);

  G4VSolid* shieldingBox3 = new G4Box("Sn_Shielding",
boxInnerSizeXY+2*shieldingThickness3,
boxInnerSizeZ+2*shieldingThickness3,
boxInnerSizeXY+2*shieldingThickness3);



  G4VSolid* subtractionBox = new G4Box("Subtraction-box",
    		        boxInnerSizeXY,
			boxInnerSizeZ,
			boxInnerSizeXY);
  
  // Bus structure boxes
  G4VSolid* busBackPlate = new G4Box("Back-plate",
		        14.5*cm,
			8.5*mm/2.,
			14.5*cm);

  G4VSolid* busFrontPlate = new G4Box("Front-plate",
		        7.*cm,
			8.5*mm/2.,
			10.65*cm+0.25*cm/2+4.*cm);

  G4VSolid* busSidePlate = new G4Box("Side-plate",
		        14.*cm,
			6.35*mm/2.,
			7.*cm);
  
  G4VSolid* busThickPlate = new G4Box("Thick-plate",
		        14.*cm,
			59.*mm/2.,
			7.*cm);

  //////////////////////////////////////////
  ///////////// Subtractions ///////////////
  //////////////////////////////////////////


  // empty rotation matrix for SubtractionSolid constructor
  G4RotationMatrix* rotm = new G4RotationMatrix();   
  
  G4double windowPlacement = boxInnerSizeZ + outerShieldingThickness + windowThickness*1.5;


  // Hollows out outer shielding 
  G4SubtractionSolid* shieldingBox1sub = 
	  new G4SubtractionSolid("PE_Shielding",
	  outerShieldingBox,
	  shieldingBox2);
  
  G4SubtractionSolid* shieldingBox1sub_t = 
	  new G4SubtractionSolid("PE_Shielding",
	  shieldingBox1sub,
	  subtractionBox,
	  rotm,
	  G4ThreeVector(0.,50.*mm,0.));

  // Hollows out shielding 2
  G4SubtractionSolid* shieldingBox2sub = 
	  new G4SubtractionSolid("W_Shielding",
	  shieldingBox2,
	  shieldingBox3);
  
  G4SubtractionSolid* shieldingBox2sub_t = 
	  new G4SubtractionSolid("W_Shielding",
	  shieldingBox2sub,
	  subtractionBox,
	  rotm,
	  G4ThreeVector(0.,50.*mm,0.));
  
  // Hollows out shielding 3
  G4SubtractionSolid* shieldingBox3sub = 
	  new G4SubtractionSolid("Sn_Shielding",
	  shieldingBox3,
	  subtractionBox);
  
  G4SubtractionSolid* shieldingBox3sub_t = 
	  new G4SubtractionSolid("Sn_Shielding",
	  shieldingBox3sub,
	  subtractionBox,
	  rotm,
	  G4ThreeVector(0.,50.*mm,0.));
  
  

  ////////////////////////////////////////////
  ////////////// Logical Volumes /////////////
  ////////////////////////////////////////////

  
  G4LogicalVolume* logicDetector = new G4LogicalVolume(detectorBox,
							CZT,
							"Detector");
  
  G4LogicalVolume* logicDetectorElectronics = new G4LogicalVolume(detectorElectronics,
							FR4,
							"BottomFR4");
  
  G4LogicalVolume* logicFrontEndBoard = new G4LogicalVolume(frontEndBoard,
							FR4,
							"Electronics");
  
  G4LogicalVolume* logicalOuterShielding = new G4LogicalVolume(shieldingBox1sub_t,
		  nist->FindOrBuildMaterial("G4_POLYETHYLENE"),
		  "PE_Shielding");

  G4LogicalVolume* logicalShielding2 = new G4LogicalVolume(shieldingBox2sub_t,
		  nist->FindOrBuildMaterial("G4_W"),
		  "W_Shielding");
  
  G4LogicalVolume* logicalShielding3 = new G4LogicalVolume(shieldingBox3sub_t,
		  nist->FindOrBuildMaterial("G4_Sn"),
		  "Sn_Shielding");
  
  G4LogicalVolume* logicalTopWindow = new G4LogicalVolume(topWindow,
		  nist->FindOrBuildMaterial("G4_Be"),
		  "top_Be_Window");
  

  // Bus structures

  G4LogicalVolume* logicalBusBackPlate = new G4LogicalVolume(busBackPlate,
		  nist->FindOrBuildMaterial("G4_Al"),
		  "Back-plate");

  G4LogicalVolume* logicalBusFrontPlate = new G4LogicalVolume(busFrontPlate,
		  nist->FindOrBuildMaterial("G4_Al"),
		  "Front-plate");

  G4LogicalVolume* logicalBusSidePlate = new G4LogicalVolume(busSidePlate,
		  nist->FindOrBuildMaterial("G4_Al"),
		  "Side-plate");

  G4LogicalVolume* logicalBusThickPlate = new G4LogicalVolume(busThickPlate,
		  nist->FindOrBuildMaterial("G4_Al"),
		  "Thick-plate");

  ////////////////////////////////////////////////
  ///////////////// Placements ///////////////////
  ////////////////////////////////////////////////
  
  G4AssemblyVolume* detectorAssembly = new G4AssemblyVolume();

  G4RotationMatrix Rm;
  G4ThreeVector    Tm;
  G4Transform3D    Tr;

  // Front end electronics board
  Tm.setX(0.); Tm.setY(-2.*cm); Tm.setZ(0.);
  Tr = G4Transform3D(Rm, Tm); 
  detectorAssembly->AddPlacedVolume(logicFrontEndBoard, Tr);
  
  // Call method to return coded aperture sub solid, logical volume below
  G4SubtractionSolid* logicAp1 = CreateCodedAperture();
  
  G4LogicalVolume* logic_aperature_base =
    new G4LogicalVolume(logicAp1,            //its solid
                        nist->FindOrBuildMaterial("G4_W"), // material
                        "Aperature-base");         //its name
  
  G4double detectorPosX = -21.*mm;
  G4double detectorPosZ = -21.*mm;
  G4double detectorPosY = frontEndBoardThickness + 0.5*cm;

  // position multiplier arrays
  G4int pm1[4] = {1, -1, 1, -1};
  G4int pm2[4] = {1, 1, -1, -1};

  // Create all Redlen detectors and apertures per assembly
  for(G4int nDet=0; nDet<4;nDet++)
  {
 
    // CZT Detector
    Tm.setX(pm1[nDet]*detectorPosX); 
    Tm.setY(-1.25*cm+0.1*mm+detectorPosY); 
    Tm.setZ(pm2[nDet]*detectorPosZ);
    
    Tr = G4Transform3D(Rm, Tm); 
  
    detectorAssembly->AddPlacedVolume(logicDetector, Tr);
    

    // FR4 Detector
    Tm.setX(pm1[nDet]*detectorPosX); 
    Tm.setY(-2.*cm+detectorPosY); 
    Tm.setZ(pm2[nDet]*detectorPosZ);
    
    Tr = G4Transform3D(Rm, Tm); 

    detectorAssembly->AddPlacedVolume(logicDetectorElectronics, Tr);
  
    // Coded Aperture
    Tm.setX(pm1[nDet]*detectorPosX);
    Tm.setY(detectorPosY + detectorApertureSpacing + detectorZ/2.);
    Tm.setZ(pm2[nDet]*detectorPosZ);
    
    Rm.rotateX(90.*deg);
    Tr = G4Transform3D(Rm, Tm); 

    detectorAssembly->AddPlacedVolume(logic_aperature_base, Tr);
    Rm.rotateX(-90.*deg);
    
   } 

  // Polyethylene shielding
  Tm.setX(0.); Tm.setY(0.); Tm.setZ(0.);
  Tr = G4Transform3D(Rm, Tm); 

  detectorAssembly->AddPlacedVolume(logicalOuterShielding, Tr);


  // Tungsten shielding
  Tm.setX(0.); Tm.setY(0.); Tm.setZ(0.);
  Tr = G4Transform3D(Rm, Tm); 

  detectorAssembly->AddPlacedVolume(logicalShielding2, Tr);

  // Tin shielding
  Tm.setX(0.); Tm.setY(0.); Tm.setZ(0.);
  Tr = G4Transform3D(Rm, Tm); 

  detectorAssembly->AddPlacedVolume(logicalShielding3, Tr);


  // Top beryllium window
  Tm.setX(0.); Tm.setY(windowPlacement-5.*mm); Tm.setZ(0.);
  Tr = G4Transform3D(Rm, Tm); 

  detectorAssembly->AddPlacedVolume(logicalTopWindow, Tr);
  
  
 
  // Bus structure placements
  
  G4double busHeight = 1.*cm;

  G4RotationMatrix* wallRotm = new G4RotationMatrix();
  wallRotm->rotateZ(90.*deg);

  new G4PVPlacement(0,
		  G4ThreeVector(0., -7.0*cm+busHeight, 0.),
		  logicalBusBackPlate,
		  "Back-plate",
		  logicEnv,
		  false,
		  checkOverlaps);

  
  new G4PVPlacement(wallRotm,
		  G4ThreeVector(-14.*cm-8.5*mm/2, busHeight+0.65*cm, 0.25*cm/2),
		  logicalBusFrontPlate,
		  "Front-plate",
		  logicEnv,
		  false,
		  0,
		  checkOverlaps);

  
  new G4PVPlacement(wallRotm,
		  G4ThreeVector(14.*cm+8.5*mm/2, busHeight+0.65*cm, 0.25*cm/2),
		  logicalBusFrontPlate,
		  "Front-plate",
		  logicEnv,
		  false,
		  1,
		  checkOverlaps);

  G4RotationMatrix* sideWallRotm = new G4RotationMatrix();
  sideWallRotm->rotateX(90.*deg);

  new G4PVPlacement(sideWallRotm,
		  G4ThreeVector(0., busHeight+0.5*cm, -14.5*cm),
		  logicalBusSidePlate,
		  "Side-plate",
		  logicEnv,
		  false,
		  checkOverlaps);

  new G4PVPlacement(sideWallRotm,
		  G4ThreeVector(0., busHeight+5.*mm, 14.5*cm+59.*mm/2),
		  logicalBusThickPlate,
		  "Thick-plate",
		  logicEnv,
		  false,
		  checkOverlaps);

  // Place the 3 copies of the detector assemblies using the position 
  // multiplier arrays from above
  unsigned int numDetectorAssemblies = 3;
  G4double dimX = -7.0*cm;
  G4double dimZ = -7.0*cm;
  Rm.rotateY(0.*deg);
 
  G4double detectorHeight = 50.*mm;
  for(unsigned int i=0; i<numDetectorAssemblies; i++){	 
    Tm.setX(pm1[i]*dimX); Tm.setY(detectorHeight); Tm.setZ(pm2[i]*dimZ);
    Tr = G4Transform3D(Rm, Tm);

    detectorAssembly->MakeImprint(logicEnv, Tr);

  }


  // always return the physical World
  return physWorld;
}

G4SubtractionSolid* DetectorConstruction::CreateCodedAperture()
{
  G4RotationMatrix Rm;
  G4ThreeVector    Tm;
  G4Transform3D    Tr;
  
  G4double boxXY 	   = 4.*cm;
  G4double boxZ  	   = 1.5*mm;
  G4double aperatureSquare = 0.2*cm;

  // added dimension to "fill the gap" between detectors
  G4Box* aperature_base = new G4Box("Aperature-base",
		   		    (boxXY+2.*mm)/2.,
				    (boxXY+2.*mm)/2.,
				    boxZ/2.);

  G4RotationMatrix* rotm = new G4RotationMatrix();   

  G4Box* coded_box = new G4Box("Coded-box",
		  		aperatureSquare/2.,
				aperatureSquare/2.,
				boxZ+5.*mm);
  
  
  G4UnionSolid* swapSolid;
  G4String placementXY_str; 
  G4double placementX, placementY; 
  G4String token;
  std::ifstream placementFile("coded_aperture_array.txt", std::ios_base::in);
  
  // Get number of lines in file
  int numberOfBoxes = 0;
  while(getline(placementFile, placementXY_str, '\n'))
    { numberOfBoxes++; }
  
  placementFile.close();


  // Reopen file to start from first line
  placementFile.open("coded_aperture_array.txt", std::ios_base::in);
  getline(placementFile, placementXY_str, '\n');
  
  token = placementXY_str.substr(
		  0, 
  		  placementXY_str.find(',')); 
  
  placementX = std::stod(token);
  
  token = placementXY_str.substr(
		  placementXY_str.find(',')+1, 
		  placementXY_str.find('\n'));
  
  placementY = std::stod(token);
  
  
  G4UnionSolid* coded_boxes = new G4UnionSolid("Combined-boxes",
		  				coded_box,
						coded_box,
						rotm,
						G4ThreeVector(
							placementX*cm,
							placementY*cm,
							0.)); 
  
  // starts at 1 since logicAp1 uses first line of file 
  for(int i=1; i<numberOfBoxes; i++)
  {
    getline(placementFile, placementXY_str, '\n');

    token = placementXY_str.substr(
		  0, 
  		  placementXY_str.find(',')); 
    
    placementX = std::stod(token); 
 
    token = placementXY_str.substr(
		  placementXY_str.find(',')+1, 
		  placementXY_str.find('\n'));
    
    placementY = std::stod(token); 

    swapSolid = new G4UnionSolid("Aperature-base",
	  			   coded_boxes,
	  			   coded_box,
	  			   rotm,
	  			   G4ThreeVector(placementX*cm,
					         placementY*cm,
						 0.));
 
    coded_boxes = swapSolid;
  }
  
  placementFile.close();

  G4SubtractionSolid* logicAp1 = 
	    new G4SubtractionSolid("Aperature-base",
	  			   aperature_base,
	  			   coded_boxes,
	  			   rotm,
	  			   G4ThreeVector(0.,0.,0.));
  
  return logicAp1; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
