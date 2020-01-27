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
#include "G4PVReplica.hh"

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

  G4double outerShieldingThickness = 15.*mm; // PE 
  G4double shieldingThickness2     = 3.5*mm; // W
  G4double shieldingThickness3     = 2.5*mm; // Sn
  G4double detectorXY      	   = 40.*mm;
  G4double detectorZ       	   = 6.*mm;
  G4double detectorElectronicsZ    = 10.185*mm;
  G4double boxInnerSizeXY  	   = 90.*mm * 0.5;
  G4double windowThickness 	   = 1.*mm;// 2 windows, each 1 mm
  G4double frontEndBoardThickness  = 2.86*mm;
  G4double detectorApertureSpacing = 12.*mm;
  G4double detectorHeight 	   = 50.*mm;
  G4double collimatorHeight 	   = 13.*mm;

  G4double pixelSize      = 2.5*mm;
  /////////////////////////////////////////
  //////////////// Geometry ///////////////
  /////////////////////////////////////////
  
  G4VSolid* pixelBlock = new G4Box("Pixel",
		  		0.5*pixelSize,
				0.5*detectorZ,
		  		0.5*pixelSize);


  G4VSolid* detectorElectronics = new G4Box("BottomFR4",
		  		    0.5*detectorXY,
				    0.5*detectorElectronicsZ,
				    0.5*detectorXY);

  G4VSolid* frontEndBoard = new G4Box("Electronics",
		  		    boxInnerSizeXY-5.*mm,
				    0.5*frontEndBoardThickness,
				    boxInnerSizeXY-5.*mm);

  G4VSolid* collimatorBlock = new G4Box("Collimator",
		  		   0.5*collimatorHeight,
				   0.5*1.*mm,
				   0.5*2*detectorXY); 
  
  G4RotationMatrix* rotm = new G4RotationMatrix();   
  rotm->rotateX(90.*deg);

  G4VSolid* collimatorUnion = new G4UnionSolid("Collimator",
		  		   collimatorBlock,
				   collimatorBlock,
				   rotm,
				   G4ThreeVector()); 
  delete rotm;

  G4double shieldingHeight = 5.*cm;
  G4double shieldingXZ     = -4.30*cm; 
  // Polyethylene shielding
  G4String name1 = "pe_Shielding";  

  G4double box1OuterDim = 86.5*mm+3.5*mm;
  G4double boxDepth1    = 70.*mm;
  G4double aBit 	= 4.*mm;
  G4double aLittleBit 	= 0.6*mm;
  G4double depthCorrection = 8.*mm;

  G4LogicalVolume* logic_shielding1 = CreateLshielding(box1OuterDim,
		  boxDepth1-27.*mm-depthCorrection/2.,
		  outerShieldingThickness,
		  -10.5*mm,
		  nist->FindOrBuildMaterial("G4_POLYETHYLENE"),
		  name1); 
  
  new G4PVPlacement(0,
		G4ThreeVector(shieldingXZ+2.*mm-0.6*mm, 
			      shieldingHeight-7.25*mm-depthCorrection/2., 
			      shieldingXZ+2.*mm-0.6*mm),
		  logic_shielding1,
		  name1,
		  logicEnv,
		  false,
		  checkOverlaps);
  
  // Tungsten shielding
  G4String name2 = "W_Shielding"; 
  
  G4double box2OuterDim = 86.5*mm;
  G4double boxDepth2    = boxDepth1 - outerShieldingThickness;
  
  G4LogicalVolume* logic_shielding2 = CreateLshielding(box2OuterDim,
			boxDepth2+aBit-12.*mm-depthCorrection+2.*mm,
		  	shieldingThickness2,
			-3.*mm,
		  	nist->FindOrBuildMaterial("G4_W"),
		  	name2); 
  
  
  new G4PVPlacement(0,
		  G4ThreeVector(shieldingXZ-aLittleBit+2.*mm,
			  shieldingHeight-3.5*mm-depthCorrection/2.,
			  shieldingXZ-aLittleBit+2.*mm),
		  logic_shielding2,
		  name2,
		  logicEnv,
		  false,
		  checkOverlaps);
 
  // Tin shielding
  G4String name3 = "Sn_Shielding"; 
  
  G4double box3OuterDim = 84.*mm;
  G4double boxDepth3    = boxDepth2 - shieldingThickness2; 
  G4LogicalVolume* logic_shielding3 = CreateLshielding(box3OuterDim,
			boxDepth3 - 5.*mm-depthCorrection+2.*mm,
		  	shieldingThickness3,
			-1.25*mm,
		  	nist->FindOrBuildMaterial("G4_Sn"),
		  	name3); 
  
  new G4PVPlacement(0,
		  G4ThreeVector(shieldingXZ-aLittleBit+2.*mm,
			  	shieldingHeight+0.25*cm-5.*mm-depthCorrection/2.,
			  	shieldingXZ-aLittleBit+2.*mm),
		  logic_shielding3,
		  name3,
		  logicEnv,
		  false,
		  checkOverlaps);
  
  ////////////////////////////////////////////
  ////////////// Logical Volumes /////////////
  ////////////////////////////////////////////

  G4LogicalVolume* logicPixel; 
 

  G4LogicalVolume* logicCollimator = new G4LogicalVolume(collimatorUnion,
		  		      nist->FindOrBuildMaterial("G4_W"),
				      "Collimator");

  G4LogicalVolume* logicDetectorElectronics = new G4LogicalVolume(detectorElectronics,
							FR4,
							"BottomFR4");
  
  G4LogicalVolume* logicFrontEndBoard = new G4LogicalVolume(frontEndBoard,
							FR4,
							"Electronics");
  
  G4LogicalVolume* logicalTopWindow = CreateBerylliumWindow(
		84.*mm, 
		windowThickness,
		nist->FindOrBuildMaterial("G4_Be"),
		"top_Be_Window");
  

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
 

  // Collimator placements
  Tm.setX(0.*cm); Tm.setY(0.60*cm); Tm.setZ(0.);
  Rm.rotateX(90.*deg);
  Rm.rotateZ(90.*deg);
  Tr = G4Transform3D(Rm, Tm); 
  detectorAssembly->AddPlacedVolume(logicCollimator, Tr);
  
  Rm.rotateZ(-90.*deg);
  Rm.rotateX(-90.*deg);


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

  
  G4AssemblyVolume* pixelAssembly = new G4AssemblyVolume();
  
  G4RotationMatrix RmT;
  G4ThreeVector    TmT;
  G4Transform3D    TrT;
  
  G4String pixelName;

  // Create pixelated detector
  G4int nameCounter = 0;
  G4int numPixels = 16;
  for(G4int i=0; i<numPixels; i++){
    for(G4int j=0; j<numPixels; j++){
  
	pixelName = "P";
	pixelName += std::to_string(nameCounter);
	pixelName += "_i";
	pixelName += std::to_string(i);
	pixelName += "_j";
	pixelName += std::to_string(j);
	nameCounter++;

	logicPixel = new G4LogicalVolume(pixelBlock,
		  		   CZT,
				   pixelName);

	TmT.setX(pixelSize*(i-numPixels));
	TmT.setZ(pixelSize*(j-numPixels));
	TmT.setY(-1.25*cm+0.1*mm+detectorPosY);
	
	TrT = G4Transform3D(RmT, TmT);	
        
	pixelAssembly->AddPlacedVolume(logicPixel, TrT);
    }
  }
  
  // Create all Redlen detectors and apertures per assembly
  for(G4int nDet=0; nDet<4;nDet++)
  {
 
    // CZT Detector
    Tm.setX(pm1[nDet]*detectorPosX + 2.1*cm); 
    Tm.setY(-1.25*cm+detectorPosY+detectorZ); 
    Tm.setZ(pm2[nDet]*detectorPosZ + 2.1*cm);
    
    Tr = G4Transform3D(Rm, Tm); 
  
    detectorAssembly->AddPlacedAssembly(pixelAssembly, Tr);

    // FR4 Detector
    Tm.setX(pm1[nDet]*detectorPosX); 
    Tm.setY(-2.*cm+detectorPosY); 
    Tm.setZ(pm2[nDet]*detectorPosZ);
    
    Tr = G4Transform3D(Rm, Tm); 

    detectorAssembly->AddPlacedVolume(logicDetectorElectronics, Tr);
  
    // Coded Aperture
    Tm.setX(pm1[nDet]*detectorPosX + 0.5*mm-100.*um);
    Tm.setY(detectorPosY + detectorApertureSpacing + detectorZ-1.25*cm);
    Tm.setZ(pm2[nDet]*detectorPosZ + 0.5*mm-100.*um);
    
    Rm.rotateX(90.*deg);
    Tr = G4Transform3D(Rm, Tm); 
    detectorAssembly->AddPlacedVolume(logic_aperature_base, Tr);
    Rm.rotateX(-90.*deg);
    
   } 



  // Top beryllium window
  new G4PVPlacement(0,
		  G4ThreeVector(shieldingXZ-aLittleBit+2.*mm,
    detectorHeight+detectorPosY+detectorApertureSpacing
    +detectorZ-1.25*cm+(windowThickness+1.5*mm)/2.,
			  	shieldingXZ-aLittleBit+2.*mm),
		  logicalTopWindow,
		  "top_Be_Window",
		  logicEnv,
		  false,
		  checkOverlaps);
  
  // Bottom beryllium window 
  new G4PVPlacement(0,
		  G4ThreeVector(shieldingXZ-aLittleBit+2.*mm,
    detectorHeight+detectorPosY+detectorApertureSpacing
    +detectorZ-1.25*cm-(windowThickness+1.5*mm)/2.,
			  	shieldingXZ-aLittleBit+2.*mm),
		  logicalTopWindow,
		  "bottom_Be_Window",
		  logicEnv,
		  false,
		  checkOverlaps);

  // Place the 3 copies of the detector assemblies using the position 
  // multiplier arrays from above
  unsigned int numDetectorAssemblies = 3;
  G4double dimX = -4.2*cm;
  G4double dimZ = -4.2*cm;
  Rm.rotateY(0.*deg);
 
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
  // FIXME
  G4double aperatureSquare = 0.22*cm;

  // added dimension to "fill the gap" between detectors
  G4double fillTheGap = 2.*mm;
  G4Box* aperature_base = new G4Box("Aperature-base",
		   		    (boxXY+fillTheGap)/2.,
				    (boxXY+fillTheGap)/2.,
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

G4LogicalVolume* DetectorConstruction::CreateLshielding(G4double outerDim,
		G4double    boxDepth,				
		G4double    shieldingThickness,
		G4double    finiteThicknessOffset, 
		G4Material* shieldingMaterial,
		G4String    shieldingName)
{
  // Units assigned at function call!!
  
  G4Box* mainBlock = new G4Box("B1",
		   	    (outerDim+shieldingThickness)/2.,
			    (boxDepth+shieldingThickness)/2.,
			    (outerDim+shieldingThickness)/2.);
  
  G4Box* sideBlock = new G4Box("B1_s",
		(outerDim+shieldingThickness)/2.+finiteThicknessOffset,
			    (boxDepth+shieldingThickness)/2.,
			    (outerDim+shieldingThickness)/2.);


  // Subtraction box, made arbitrarily tall to remove 
  // top wall of shielding
  G4Box* subtraction_block = new G4Box("B1",
		   	    outerDim/2.,
			    boxDepth/2. + shieldingThickness,
			    outerDim/2.);

  G4Box* sideSubtraction_block = new G4Box("B1_s",
		   	    outerDim/2. + finiteThicknessOffset,
			    boxDepth/2. + shieldingThickness,
			    outerDim/2.);
  
  
  G4Box* middleWallRemover = new G4Box("MWR",
		 shieldingThickness*2.5, 			
		 boxDepth/2. + shieldingThickness,
		 outerDim/2.);
  
  
  // Empty rotation matrix for union and subtraction solids
  G4RotationMatrix* rotm = new G4RotationMatrix();
  
  // Main L shape of shielding, union of two blocks
  G4UnionSolid* solid_L = new G4UnionSolid("solid-L",
		  		mainBlock,
				sideBlock,
				rotm,
				G4ThreeVector(
			outerDim+shieldingThickness+finiteThicknessOffset,
		  			0,
					0));
  // Union of 3rd block
  // union on Z-side
  rotm->rotateY(90.*deg);
  solid_L = new G4UnionSolid("solid-L",
		  		solid_L,
				sideBlock,
				rotm,
				G4ThreeVector(0,
		  			      0,
		 outerDim+shieldingThickness+finiteThicknessOffset));

  // (unrotate)
  rotm->rotateY(-90.*deg);
  
  // L shape to subtraction from main L shielding, union of 2 blocks
  // union occurs on X-side
  G4UnionSolid* sub_L = new G4UnionSolid("sub-L",
		  			subtraction_block,
					sideSubtraction_block,
					rotm,
	G4ThreeVector(outerDim+shieldingThickness+finiteThicknessOffset,
		  			0,
					0));
  
  // Union of 3rd block fo subtracting off from the L-shape
  rotm->rotateY(90.*deg);
  sub_L = new G4UnionSolid("sub-L",
	  		   sub_L,
			   sideSubtraction_block,
			   rotm,
			   G4ThreeVector(0,
		  		         0,
		 		         outerDim+shieldingThickness+finiteThicknessOffset));
  
  rotm->rotateY(-90.*deg);

  // Upward shift to ensure top of shielding is removed,
  // and bottom has thickness shieldingThickness
  G4ThreeVector subtraction_shift = G4ThreeVector(0.,
		  				  shieldingThickness,
		  				  0.);

  // Main subtraction
  G4SubtractionSolid* hollow_L = new G4SubtractionSolid(shieldingName,
		  					solid_L,
							sub_L,
							rotm,
						       subtraction_shift);


  // Removes middle walls created by finite subtraction box thickness
  hollow_L = new G4SubtractionSolid(shieldingName,
 					hollow_L,
					middleWallRemover,
					rotm,
 				        G4ThreeVector(outerDim/2.,
						shieldingThickness,
						0.)); 
  
  // Removes second inner wall
  rotm->rotateY(90.*deg);
  hollow_L = new G4SubtractionSolid(shieldingName,
 					hollow_L,
					middleWallRemover,
					rotm,
 				        G4ThreeVector(0.,
						shieldingThickness,
						outerDim/2.)); 
  
 
  // Final logical volume to be returned
  G4LogicalVolume* logic_L = new G4LogicalVolume(hollow_L,
						shieldingMaterial,
						shieldingName);  


  return logic_L;
}

G4LogicalVolume* DetectorConstruction::CreateBerylliumWindow(
		G4double windowDimension, 
		G4double windowThickness,
		G4Material* windowMaterial,
		G4String windowName)
{

  G4Box* windowBox = new G4Box("Window-box",
		  	       windowDimension/2.,
			       windowThickness/2.,
			       windowDimension/2.);

  G4Box* sideWindowBox = new G4Box("Window-box",
		  	       windowDimension/2.,
			       windowThickness/2.,
			       windowDimension/2.);
  
  G4RotationMatrix* rotm = new G4RotationMatrix();

  
  G4UnionSolid* L_window = new G4UnionSolid("Window-concat",
		  	      windowBox,
			      sideWindowBox,
			      rotm,
			      G4ThreeVector(windowDimension,
				      	    0.,
					    0.));
  
 
  rotm->rotateY(90.*deg);
  L_window = new G4UnionSolid("Window-concat",
		  	      L_window,
			      sideWindowBox,
			      rotm,
			      G4ThreeVector(0.,
				      	    0.,
					    windowDimension));
  


  G4LogicalVolume* window = new G4LogicalVolume(L_window,
		  				windowMaterial,
						windowName);

  return window;
}


