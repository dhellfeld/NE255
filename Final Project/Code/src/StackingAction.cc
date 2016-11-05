#include "StackingAction.hh"
#include "G4Track.hh"
#include "G4Electron.hh"

//==================================================================================

StackingAction* StackingAction::fgInstance = 0;

//==================================================================================

StackingAction* StackingAction::Instance() {
    
    return fgInstance;
    
}

//==================================================================================

StackingAction::StackingAction(){
    
    fgInstance = this;
    
}

//==================================================================================

StackingAction::~StackingAction(){

    fgInstance = 0;
}

//==================================================================================

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* track){
  
    // Keep primary particle
    if (track->GetParentID() == 0){
        return fUrgent;
    }

    // Kill secondary electrons
    if (track->GetDefinition() == G4Electron::Electron()){
        
        return fKill;
    }
    
    else { return fUrgent; }

}

//==================================================================================
