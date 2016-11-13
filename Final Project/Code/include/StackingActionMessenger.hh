#ifndef stackingactionmessenger_hh
#define stackingactionmessenger_hh 1

#include "globals.hh"
#include "G4UIcommand.hh"
#include "G4UImessenger.hh"

class StackingAction;
class G4UIdirectory;
class G4UIcmdWithAString;

class StackingActionMessenger : public G4UImessenger {
    
public:
    StackingActionMessenger(StackingAction* stacking);
    ~StackingActionMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValue);
    
    
private:
    StackingAction* fStackingAction;
    
    G4UIdirectory* fTrackDirectory;
    G4UIcmdWithAString* fTrackElectronCmd;

};

#endif
