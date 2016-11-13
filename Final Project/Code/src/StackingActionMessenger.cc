#include "StackingActionMessenger.hh"
#include "StackingAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

#include <istream>
#include <sstream>


//==================================================================================================

StackingActionMessenger::StackingActionMessenger(StackingAction * stacking)
: G4UImessenger(), fStackingAction(stacking)
{

    fTrackDirectory = new G4UIdirectory("/PRISM/tracking/");
    fTrackDirectory->SetGuidance("Track settings");

    fTrackElectronCmd = new G4UIcmdWithAString("/PRISM/tracking/trackElectrons",this);
    fTrackElectronCmd->SetGuidance("Turn electron tracking on/off");
    fTrackElectronCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//==================================================================================================

StackingActionMessenger::~StackingActionMessenger()
{
    delete fTrackElectronCmd;
    delete fTrackDirectory;

}

//==================================================================================================

void StackingActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

    if (command == fTrackElectronCmd) {

        G4String choice = newValue;

        if (choice == "on" || choice == "ON" || choice == "On"){

            // if its turned off, turn it on
            if (!fStackingAction->GetElectronTracking()){
                fStackingAction->SetElectronTracking(true);
            }

        }

        else if (choice == "off" || choice == "OFF" || choice == "Off"){

            // if its turn on, turn it off
            if (fStackingAction->GetElectronTracking()){
                fStackingAction->SetElectronTracking(false);
            }

        }
    }
}

//==================================================================================================


// Write GetCurrentValue functions...
