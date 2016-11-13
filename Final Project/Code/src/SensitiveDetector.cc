#include "SensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

using namespace CLHEP;

//==================================================================================

SensitiveDetector::SensitiveDetector(const G4String& name,const G4String& hitsCollectionName)
: G4VSensitiveDetector(name),
  fHitsCollection(NULL)
{
    collectionName.insert(hitsCollectionName);
}

//==================================================================================

SensitiveDetector::~SensitiveDetector(){}

//==================================================================================

void SensitiveDetector::Initialize(G4HCofThisEvent* HCE){

    static int HCID = -1;

    // Create hits collection
    fHitsCollection= new HitsCollection(SensitiveDetectorName, collectionName[0]);

    // Add this collection in hce
    if(HCID < 0){HCID = GetCollectionID(0);}         //G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    HCE->AddHitsCollection(HCID, fHitsCollection);

}

//==================================================================================

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*){

    class Hit* newHit = new class Hit();

    // Track ID
    newHit->SetTrackID(aStep->GetTrack()->GetTrackID());

    // Incident energy
    newHit->SetEnergy(G4float(aStep->GetPreStepPoint()->GetTotalEnergy()));

    // Detector ID
    newHit->SetVol(atoi(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()));

    // Interaction Process
    G4String proc = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    if (proc == "photo"){
        newHit->SetProcess(1);
    }
    else if (proc == "compt"){
        newHit->SetProcess(2);
    }
    else if (proc == "rayl"){
        newHit->SetProcess(3);
    }

    // HEALPix index
    newHit->SetHPindex(PrimaryGeneratorAction::Instance()->GetHP_index());


    // Depth of interaction
    G4ThreeVector detcent = DetectorConstruction::Instance()->GetDetCenters()[atoi(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()) - 1];
    G4ThreeVector intpos = aStep->GetPostStepPoint()->GetPosition();
    G4double doi = (intpos-detcent).dot(detcent)/detcent.mag() + 0.500*cm;
    if      (doi < 0.) {newHit->SetDOI(0. *mm);}
    else if (doi > 10.){newHit->SetDOI(10.*mm);}
    else {newHit->SetDOI(doi);}

    // Add the newHit to the HitCollection
    fHitsCollection->insert( newHit );

    return true;
}

//==================================================================================

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*){}

//==================================================================================
