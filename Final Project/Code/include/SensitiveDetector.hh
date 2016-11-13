#ifndef sensitivedetector_hh
#define sensitivedetector_hh

#include "G4VSensitiveDetector.hh"
#include "Hit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class SensitiveDetector : public G4VSensitiveDetector
{
public:
    SensitiveDetector(const G4String& name, const G4String& hitsCollectionName);
    virtual ~SensitiveDetector();

    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

private:
    HitsCollection* fHitsCollection;
};

#endif
