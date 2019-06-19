
#ifndef RedlenSensitiveDetector_h
#define RedlenSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"

#include "RedlenHit.hh"

class G4Step;
class G4HCofThisEvent;

class RedlenSensitiveDetector : public G4VSensitiveDetector
{
  public:
    RedlenSensitiveDetector(const G4String& SDname);
    virtual ~RedlenSensitiveDetector();

  public:
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void Initialize(G4HCofThisEvent* HCE);
    virtual void EndOfEvent(G4HCofThisEvent* HCE);

  private:
    RedlenHitsCollection* fHitsCollection;

};

#endif
