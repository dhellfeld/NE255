
#include "DetectorConstructionMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithoutParameter.hh"

#include <istream>
#include <sstream>

using namespace std;

//==================================================================================================

DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction * detector)
: G4UImessenger(), fDetectorConstruction(detector)
{

    
    fDirectory = new G4UIdirectory("/PRISM_SIM/");
    fDirectory->SetGuidance("Top directory for sim");
    
    fDirectory = new G4UIdirectory("/PRISM_SIM/geom/");
    fDirectory->SetGuidance("Manipulate geometry");
    
    fSetMaskCmd = new G4UIcmdWithAString("/PRISM_SIM/geom/setMask",this);
    fSetMaskCmd->SetGuidance("Set the mask with hexidecimal code");
    fSetMaskCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fRandMaskCmd = new G4UIcmdWithoutParameter("/PRISM_SIM/geom/randomMask",this);
    fRandMaskCmd->SetGuidance("Generate random mask and update geometry");
    fRandMaskCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fFullMaskCmd = new G4UIcmdWithoutParameter("/PRISM_SIM/geom/fullMask",this);
    fFullMaskCmd->SetGuidance("Fully populate the mask and update geometry");
    fFullMaskCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    
}

//==================================================================================================

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{

    delete fDirectory;
    delete fSetMaskCmd;
    delete fRandMaskCmd;

}

//==================================================================================================

void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

    if (command == fSetMaskCmd) {
        
        vector<int> maskarray = fDetectorConstruction->HexToBin(newValue);
        fDetectorConstruction->SetMask(maskarray);
        fDetectorConstruction->UpdateGeometry();
    }
    
    else if (command == fRandMaskCmd) {
        
        vector<int> maskarray = fDetectorConstruction->GetRandomMask();
        fDetectorConstruction->SetMask(maskarray);
        fDetectorConstruction->UpdateGeometry();
    }
    
    else if (command == fFullMaskCmd) {
        
        vector<int> maskarray(192, 1);
        fDetectorConstruction->SetMask(maskarray);
        fDetectorConstruction->UpdateGeometry();
    }

    
}

//==================================================================================================
