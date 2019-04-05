/*
written by R. Ian Crocker on 3 Mar 2008
*/

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;

class PhysicsList:  public G4VModularPhysicsList
{
  public:
    PhysicsList();
    ~PhysicsList();

  protected:
    // Construct particle and physics process
    void ConstructParticle();

    void SetCuts();

    //the following is from drews PhysicsList.hh
    //These methods construct particles
	//void ConstructLeptons();
	//void ConstructBaryons();

    //These methods consruct physics processes and register them
    //void ConstructGeneral();
    //void ConstructEM();
    
   void AddPhysicsList(const G4String& name);
   void ConstructProcess();
   void List();
   
 private:

  void SetStandardList(G4bool flagHP = false, G4bool glauber = false);

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  G4VPhysicsConstructor*  emPhysicsList;
  G4VPhysicsConstructor*  particleList;
  std::vector<G4VPhysicsConstructor*>  hadronPhys;
    
  G4bool dump;  
   
};

#endif







