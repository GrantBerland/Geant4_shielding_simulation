

#include "RedlenSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"


RedlenSensitiveDetector::RedlenSensitiveDetector(const G4String& name,
		                 const G4String& hitsCollectionName):
	G4VSensitiveDetector(name),
	fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

RedlenSensitiveDetector::~RedlenSensitiveDetector()
{}


void RedlenSensitiveDetector::Initialize(G5HCofThisEvent* hce)
{
  fHitsCollection = 
  new RedlenHitsCollection(SensitiveDetectorName, collectionName[0]);


  G4int hcID = 
	G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  hce->AddHitsCollection(hcID, fHitsCollection);

}

G4bool RedlenSensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  RedlenHit* newHit = new RedlenHit();
  
  
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());

  fHitsCollection->insert(newHit);

  return true;
}

void RedlenSensitiveDetector::EndOfEvent()
{
  // Write results to file

  G4int nOfHits = fHitsCollection->entires();

  for(G4int i=0; i<nOfHits; i++)
  {
    std::cout << "Hit!" << std::endl;
    // (*fHitsCollection)[i]-> (access methods)

  }
}

