

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIdirectory.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* prim):G4UImessenger(),fPrimaryGenerator(prim) 
{
  fPrimDir = new G4UIdirectory("/particleSource/");
  fPrimDir->SetGuidance("Select trapped background, loss cone background, or photon signal source.");

  fcmd = new G4UIcmdWithAnInteger("/particleSource/setBackgroundType",this);
  fcmd->SetParameterName("Background Type {0,1,2}",true);
  fcmd->SetDefaultValue(0);
  fcmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDcmd = new G4UIcmdWithADouble("/particleSource/setFoldingEnergy",this);
  fDcmd->SetParameterName("Folding Energy {100, 200, 300} keV",true);
  fDcmd->SetDefaultValue(100.);
  fDcmd->AvailableForStates(G4State_PreInit, G4State_Idle);


}



PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPrimDir;
  delete fcmd;
  delete fDcmd;
}


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{

  if(command == fcmd){
    fPrimaryGenerator->SetWhichParticle(std::stoi(newValue));
  }    	  

  if(command == fDcmd){
    fPrimaryGenerator->SetFoldingEnergy(std::stod(newValue));
  }

}
