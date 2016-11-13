#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"

//============================================================================//

ActionInitialization::ActionInitialization(){}

//============================================================================//

ActionInitialization::~ActionInitialization(){}

//============================================================================//

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction);
}

//============================================================================//

void ActionInitialization::Build() const {

    // for thread-local class instances
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
    SetUserAction(new EventAction);
    SetUserAction(new SteppingAction);
    SetUserAction(new StackingAction);

}

//============================================================================//
