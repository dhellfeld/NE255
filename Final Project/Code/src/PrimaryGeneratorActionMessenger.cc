
#include "PrimaryGeneratorActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"

#include <istream>
#include <sstream>

using namespace CLHEP;
using namespace std;


//==================================================================================================

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction * generator)
: G4UImessenger(), fPrimaryGeneratorAction(generator)
{

    fDirDirectory = new G4UIdirectory("/PRISM_SIM/direction/");
    fDirDirectory->SetGuidance("Direction of incoming rays");
    
    fRunDirectory = new G4UIdirectory("/PRISM_SIM/run/");
    fRunDirectory->SetGuidance("Run controls");
    
    fHPCmd = new G4UIcmdWithAnInteger("/PRISM_SIM/direction/setHPindex",this);
    fHPCmd->SetGuidance("Set the HEALPix index far field souce");
    fHPCmd->SetParameterName("HPindex",false);
    fHPCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fHPNsideCmd = new G4UIcmdWithAnInteger("/PRISM_SIM/direction/setHPNside",this);
    fHPNsideCmd->SetGuidance("Set nside for HEALPix (8 or 16)");
    fHPNsideCmd->SetParameterName("HPNside",false);
    fHPNsideCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fHPindexingCmd = new G4UIcmdWithAString("/PRISM_SIM/direction/setHPindexscheme",this);
    fHPindexingCmd->SetGuidance("Set indexing scheme for HEALPix (nested or ring)");
    fHPindexingCmd->SetParameterName("HPNside",false);
    fHPindexingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetUpHEALPixCmd = new G4UIcmdWithoutParameter("/PRISM_SIM/direction/SetUpHEALPix",this);
    fSetUpHEALPixCmd->SetGuidance("Set up HEALPix (command to follow the Nside and IndexingScheme commands)");
    fSetUpHEALPixCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    
    fNumPartCmd = new G4UIcmdWithAnInteger("/PRISM_SIM/run/particles",this);
    fNumPartCmd->SetGuidance("Tell me how many particles");
    fNumPartCmd->SetParameterName("NumPart",false);
    fNumPartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fStartHPCmd = new G4UIcmdWithAnInteger("/PRISM_SIM/run/HPstart",this);
    fStartHPCmd->SetGuidance("Start HP index");
    fStartHPCmd->SetParameterName("NumPart",false);
    fStartHPCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    
    fFinishHPCmd = new G4UIcmdWithAnInteger("/PRISM_SIM/run/HPfinish",this);
    fFinishHPCmd->SetGuidance("Finish HP index");
    fFinishHPCmd->SetParameterName("NumPart",false);
    fFinishHPCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    
    fDOIsettingCmd = new G4UIcmdWithABool("/PRISM_SIM/run/DOI",this);
    fDOIsettingCmd->SetGuidance("Turn DOI on/off");
    fDOIsettingCmd->SetParameterName("DOIsetting",false);
    fDOIsettingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//==================================================================================================

PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{

    delete fHPCmd;
    delete fHPNsideCmd;
    delete fHPindexingCmd;
    delete fSetUpHEALPixCmd;
    delete fDirDirectory;
    delete fRunDirectory;
    delete fNumPartCmd;
    delete fStartHPCmd;
    delete fFinishHPCmd;
    delete fDOIsettingCmd;

}

//==================================================================================================

void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

    
    if (command == fHPCmd){
                
        G4int HPindex = fHPCmd->GetNewIntValue(newValue);
        
        G4int nside = fPrimaryGeneratorAction->GetHPNside();
        G4int maxindex;
        if (nside == 8){maxindex = 768;}
        else if (nside == 16){maxindex = 3072;}
        else {maxindex = 0;}
        
        if (HPindex <= maxindex){
            
            vector<struct Angles> HPangles = fPrimaryGeneratorAction->GetHPangles();

            fPrimaryGeneratorAction->SetTheta(HPangles[HPindex-1].theta*(180./pi));
            fPrimaryGeneratorAction->SetPhi(HPangles[HPindex-1].phi*(180./pi));
        
            fPrimaryGeneratorAction->SetHP_index(HPindex);

        }
        
        else {
            G4cerr << "\n\nThat HEALPix index is out of bounds for the current HPNside value... HPindex will be set to 0\n\n";
        }
                              
    }
    
    else if (command == fHPNsideCmd){
        
        G4int HPNsideindex = fHPCmd->GetNewIntValue(newValue);
        
        if (HPNsideindex == 8 || HPNsideindex == 16 || HPNsideindex == 32){
            fPrimaryGeneratorAction->SetHPNside(HPNsideindex);
        }
        
        else {G4cerr << "\n\nI do not recognize your command: /setHPNside " << HPNsideindex << "\nStopping...\n\n";}

    }
    
    else if (command == fHPindexingCmd){
        
        
        G4String HPindexing = newValue;
        
        if  (HPindexing == "nested" || HPindexing == "Nested" || HPindexing == "NESTED" ||
             HPindexing == "ring"   || HPindexing == "Ring"   || HPindexing == "RING"){
            
            fPrimaryGeneratorAction->SetHPindexing(HPindexing);
        }
        
        else {G4cerr << "\n\nI do not recognize your command: '/setHPindexscheme " << HPindexing << "'\nStopping...\n\n";}

    }
    
    else if (command == fSetUpHEALPixCmd){
        
        G4String filename;
        
        G4String indexing = fPrimaryGeneratorAction->GetHPindexing();
        G4int nside = fPrimaryGeneratorAction->GetHPNside();
        
        if (nside == 8 && (indexing == "nested" || indexing == "Nested" || indexing == "NESTED")){filename = "HEALpix/n=8/HEALPix_PhiTheta_Nested.txt";}
        else if (nside == 8 && (indexing == "ring" || indexing == "Ring" || indexing == "RING")){filename = "HEALpix/n=8/HEALPix_PhiTheta_Ring.txt";}
        else if (nside == 16 && (indexing == "nested" || indexing == "Nested" || indexing == "NESTED")){filename = "HEALpix/n=16/HEALPix_PhiTheta_Nested.txt";}
        else if (nside == 16 && (indexing == "ring" || indexing == "Ring" || indexing == "RING")){filename = "HEALpix/n=16/HEALPix_PhiTheta_Ring.txt";}
        else if (nside == 32 && (indexing == "nested" || indexing == "Nested" || indexing == "NESTED")){filename = "HEALpix/n=32/HEALPix_PhiTheta_Nested.txt";}
        else if (nside == 32 && (indexing == "ring" || indexing == "Ring" || indexing == "RING")){filename = "HEALpix/n=32/HEALPix_PhiTheta_Ring.txt";}

        fPrimaryGeneratorAction->ReadInHEALPixAngles(filename);
        
    }
    
    else if (command == fNumPartCmd){
        
        fPrimaryGeneratorAction->SetNumPart(fNumPartCmd->GetNewIntValue(newValue));
        
    }

    else if (command == fStartHPCmd){
        
        fPrimaryGeneratorAction->SetHPStart(fStartHPCmd->GetNewIntValue(newValue));
        
    }
    
    else if (command == fFinishHPCmd){
        
        fPrimaryGeneratorAction->SetHPFinish(fFinishHPCmd->GetNewIntValue(newValue));
        
    }
    
    else if (command == fDOIsettingCmd){
        
        fPrimaryGeneratorAction->SetDOIsetting(fDOIsettingCmd->GetNewBoolValue(newValue));
        
    }

}

//==================================================================================================


// Write GetCurrentValue functions...



