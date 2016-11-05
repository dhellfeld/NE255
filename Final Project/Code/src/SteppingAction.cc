#include "SteppingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "Analysis.hh"


#include "G4Step.hh"
#include "G4RunManager.hh"

using namespace CLHEP;

// ====================================================================================================

SteppingAction::SteppingAction(): G4UserSteppingAction(){}

// ====================================================================================================

SteppingAction::~SteppingAction(){}

// ====================================================================================================

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    
    //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // Kill track if a compton scatter
    if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "compt"){
        
        step->GetTrack()->SetTrackStatus(fStopAndKill);
    }
    
    // Fill histograms if photoabsoprtion
    else if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "phot"){
        
        det = atoi(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName());
        hpi = PrimaryGeneratorAction::Instance()->GetHP_index();
        
        if (PrimaryGeneratorAction::Instance()->GetDOIsetting()){ doi = CalculateDOI(step); }
        else { doi = 1; }
        
        //analysisManager->FillH2(1, hpi,det);
        
        RunAction::Instance()->FillResponse(hpi-1, (det-1 + 192*(doi-1)));
    }
    
    
    
}

// ====================================================================================================

int SteppingAction::CalculateDOI(const G4Step * step){
    
    G4ThreeVector detcent = DetectorConstruction::Instance()->GetDetCenters()[atoi(step->GetPreStepPoint()->GetPhysicalVolume()->GetName()) - 1];
    G4ThreeVector intpos = step->GetPostStepPoint()->GetPosition();
    G4double _doi = (intpos-detcent).dot(detcent)/detcent.mag() + 0.50000*cm;

    int doibin = 0;
    if      (_doi < 0 ){doibin = 1;}
    else if (_doi > 10){doibin = 10;}
    else               {doibin = floor(_doi) + 1;}

    return doibin;
    
    
/*  
    // A little slower ....
 
    int doibin = 0;
    if      (_doi <  1.0*mm)                 {doibin = 1; }
    else if (_doi >= 1.0*mm && _doi < 2.0*mm){doibin = 2; }
    else if (_doi >= 2.0*mm && _doi < 3.0*mm){doibin = 3; }
    else if (_doi >= 3.0*mm && _doi < 4.0*mm){doibin = 4; }
    else if (_doi >= 4.0*mm && _doi < 5.0*mm){doibin = 5; }
    else if (_doi >= 5.0*mm && _doi < 6.0*mm){doibin = 6; }
    else if (_doi >= 6.0*mm && _doi < 7.0*mm){doibin = 7; }
    else if (_doi >= 7.0*mm && _doi < 8.0*mm){doibin = 8; }
    else if (_doi >= 8.0*mm && _doi < 9.0*mm){doibin = 9; }
    else if (_doi >  9.0*mm)                 {doibin = 10;}

    return doibin;
*/
    
}

// ====================================================================================================
