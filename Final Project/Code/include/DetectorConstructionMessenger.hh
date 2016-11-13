
#ifndef detectorconstructionnmessenger_hh
#define detectorconstructionnmessenger_hh 1

#include "globals.hh"
#include "G4UIcommand.hh"
#include "G4UImessenger.hh"

// Initialize classes
class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWith3Vector;
class G4UIcmdWithoutParameter;

class DetectorConstructionMessenger : public G4UImessenger {
    
public:
    DetectorConstructionMessenger(DetectorConstruction* detector);
    ~DetectorConstructionMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValue);
    
    
private:
    DetectorConstruction* fDetectorConstruction;
    
    G4UIdirectory* fDirectory;
    
    G4UIcmdWithAString*      fSetMaskCmd;
    G4UIcmdWithoutParameter* fRandMaskCmd;
    G4UIcmdWithoutParameter* fFullMaskCmd;
    G4UIcmdWith3Vector*      fDetDimCmd;
    G4UIcmdWithAString*      fSetDetIndexCmd;
    
};

#endif
