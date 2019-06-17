

#include "PhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4DecayPhysics.hh"
///#include "G4EmStandardPhysics72.hh"
///#include "G4EmStandardPhysics71.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics.hh"
//#include "G4QStoppingPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"
#include "G4SystemOfUnits.hh"


PhysicsList::PhysicsList(): G4VModularPhysicsList()
{
	 // Particles
  particleList = new G4DecayPhysics("decays");

  // EM physics
  ///emPhysicsList = new G4EmStandardPhysics71("standard");
  emPhysicsList = new G4EmStandardPhysics_option1();
	
	
	}

PhysicsList::~PhysicsList()
{;}

void PhysicsList::ConstructParticle()
{
	// In this method, static member functions should be called
	// for all particles which you want to use.
	// This ensures that objects of these particle types will be
	// created in the program. 

	//G4Geantino::GeantinoDefinition();
	//G4Electron::ElectronDefinition();
	//G4Proton::ProtonDefinition();
	particleList->ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
	// Define transportation process

	AddTransportation();
    emPhysicsList->ConstructProcess();
    particleList->ConstructProcess();
	
}

void PhysicsList::SetCuts()
{
	// uppress error messages even in case e/gamma/proton do not exist            
	G4int temp = GetVerboseLevel();                                                
	SetVerboseLevel(0);                                                           
	// " G4VUserPhysicsList::SetCutsWithDefault" method sets 
	// the default cut value for all particle types 
	//SetCutsWithDefault();   
	SetCutValue(0.1*mm, "electron");

	// Retrieve verbose level
	SetVerboseLevel(temp);  
}

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0)
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    
//     if (name == "LHEP") {
//
//    hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
//    hadronPhys.push_back( new G4HadronElasticPhysics("LElastic",
//						     verboseLevel,false));
//    hadronPhys.push_back( new HadronPhysicsLHEP());
//    hadronPhys.push_back( new G4IonPhysics("ion"));

 // } else
 
   if (name == "G4standard_exp") {

    delete emPhysicsList;
    ///emPhysicsList = new G4EmStandardPhysics72("standard");
	emPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "G4standard_fast") {

    delete emPhysicsList;
    ///emPhysicsList = new G4EmStandardPhysics71("standard");
    emPhysicsList = new G4EmStandardPhysics_option1();

  }  else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

void PhysicsList::SetStandardList(G4bool flagHP, G4bool glauber)
{
//  hadronPhys.push_back( new G4EmExtraPhysics("gamma_nuc"));
//  hadronPhys.push_back( new G4HadronHElasticPhysics("elastic",verboseLevel,
//						    flagHP,glauber));
//  hadronPhys.push_back( new G4QStoppingPhysics("stopping"));
//  hadronPhys.push_back( new G4IonBinaryCascadePhysics("binary_ion"));
//  hadronPhys.push_back( new G4NeutronTrackingCut());
}
