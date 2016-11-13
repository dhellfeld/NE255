
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

    fDirectory = new G4UIdirectory("/PRISM/");
    fDirectory->SetGuidance("Top directory for sim");

    fDirectory = new G4UIdirectory("/PRISM/geo/");
    fDirectory->SetGuidance("Manipulate geometry");

    fSetMaskCmd = new G4UIcmdWithAString("/PRISM/geo/setMask",this);
    fSetMaskCmd->SetGuidance("Set the mask with hexidecimal code");
    fSetMaskCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fRandMaskCmd = new G4UIcmdWithoutParameter("/PRISM/geo/randomMask",this);
    fRandMaskCmd->SetGuidance("Generate random mask and update geometry");
    fRandMaskCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fFullMaskCmd = new G4UIcmdWithoutParameter("/PRISM/geo/fullMask",this);
    fFullMaskCmd->SetGuidance("Fully populate the mask and update geometry");
    fFullMaskCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fDetDimCmd = new G4UIcmdWith3Vector("/PRISM/geo/cubedim",this);
    fDetDimCmd->SetGuidance("Set the x,y,z lengths of the detectors (in cm) (full lengths)");
    fDetDimCmd->SetParameterName("xdim","ydim","zdim", false);
    fDetDimCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSetDetIndexCmd = new G4UIcmdWithAString("/PRISM/geo/detIndexing",this);
    fSetDetIndexCmd->SetGuidance("Set the way the detectors are indexed (ring or nested)");
    fSetDetIndexCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//==================================================================================================

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{
    delete fDirectory;
    delete fSetMaskCmd;
    delete fRandMaskCmd;
    delete fDetDimCmd;
    delete fSetDetIndexCmd;
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

    else if (command == fDetDimCmd) {

        fDetectorConstruction->SetDetDim(fDetDimCmd->GetNew3VectorValue(newValue));
        fDetectorConstruction->CheckOverlapsOn();
        fDetectorConstruction->UpdateGeometry();
    }

    else if (command == fSetDetIndexCmd) {

        G4String indexing = newValue;

        if (indexing == "Ring" || indexing == "ring" || indexing == "RING"){

            if (fDetectorConstruction->GetDetIndexing() == "Nested"){
                fDetectorConstruction->SetDetIndexing("Ring");
                fDetectorConstruction->UpdateGeometry();
            }
        }

        else if (indexing == "Nested" || indexing == "nested" || indexing == "NESTED"){

            if (fDetectorConstruction->GetDetIndexing() == "Ring"){
                fDetectorConstruction->SetDetIndexing("Nested");
                fDetectorConstruction->UpdateGeometry();
            }

        }

        else {G4cerr << "\n\nI don't recognize the command /PRISM_SIM/geom/detIndexing " << indexing << ". Setting to default (ring)\n\n";}

    }


}

//==================================================================================================
