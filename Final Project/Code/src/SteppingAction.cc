#include "SteppingAction.hh"

#include "G4RunManager.hh"
#include "RunAction.hh"
#include "G4Step.hh"


//==================================================================================================

SteppingAction* SteppingAction::fgInstance = 0;

//==================================================================================================

SteppingAction* SteppingAction::Instance() {
    
	return fgInstance;
}

//==================================================================================================

SteppingAction::SteppingAction() {
    
	fgInstance = this;
}

//==================================================================================================

SteppingAction::~SteppingAction() {
    
	fgInstance = 0;
}

//==================================================================================================

void SteppingAction::UserSteppingAction(const G4Step* /*step*/) {}

//==================================================================================================


