#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>


//==================================================================================

EventAction* EventAction::fgInstance = 0;

//==================================================================================

EventAction* EventAction::Instance() {
    
    return fgInstance;
    
}

// ====================================================================================================

EventAction::EventAction(): G4UserEventAction(){

    fgInstance = this;

}

// ====================================================================================================

EventAction::~EventAction(){
 
    fgInstance = 0;
}

// ====================================================================================================

void EventAction::BeginOfEventAction(const G4Event* event)
{
    
    RunAction::Instance()->SetEventID(event->GetEventID());
    
}

// ====================================================================================================

void EventAction::EndOfEventAction(const G4Event* /*event*/)
{
    
  //G4int eventID = event->GetEventID();

}  

// ====================================================================================================
