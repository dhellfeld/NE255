
#include "PrimaryGeneratorActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4UImanager.hh"

#include <istream>
#include <sstream>

using namespace CLHEP;
using namespace std;


//==================================================================================================

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction * generator)
: G4UImessenger(), fPrimaryGeneratorAction(generator)
{

    UI = G4UImanager::GetUIpointer();

    fHPDirectory = new G4UIdirectory("/PRISM/healpix/");
    fHPDirectory->SetGuidance("Setting up Healpix");

    fSrcDirectory = new G4UIdirectory("/PRISM/source/");
    fSrcDirectory->SetGuidance("Distribution/direction/position of incoming rays");

    fOutDirectory = new G4UIdirectory("/PRISM/output/");
    fOutDirectory->SetGuidance("Output file");

    fHPCmd = new G4UIcmdWithAnInteger("/PRISM/source/setHPindex",this);
    fHPCmd->SetGuidance("Set the HEALPix index far field souce");
    fHPCmd->SetParameterName("HPindex",false);
    fHPCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fHPNsideCmd = new G4UIcmdWithAnInteger("/PRISM/healpix/setHPNside",this);
    fHPNsideCmd->SetGuidance("Set nside for HEALPix (8 or 16)");
    fHPNsideCmd->SetParameterName("HPNside",false);
    fHPNsideCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fHPindexingCmd = new G4UIcmdWithAString("/PRISM/healpix/setHPindexscheme",this);
    fHPindexingCmd->SetGuidance("Set indexing scheme for HEALPix (nested or ring)");
    fHPindexingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSetUpHEALPixCmd = new G4UIcmdWithoutParameter("/PRISM/healpix/SetUpHEALPix",this);
    fSetUpHEALPixCmd->SetGuidance("Set up HEALPix (command to follow the Nside and IndexingScheme commands)");
    fSetUpHEALPixCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fOutputFileNameCmd = new G4UIcmdWithAString("/PRISM/output/filename",this);
    fOutputFileNameCmd->SetGuidance("Name the output file");
    fOutputFileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPrintTextCmd = new G4UIcmdWithAString("/PRISM/output/printText",this);
    fPrintTextCmd->SetGuidance("Print to text (on/off)");
    fPrintTextCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPrintBinaryCmd = new G4UIcmdWithAString("/PRISM/output/printBinary",this);
    fPrintBinaryCmd->SetGuidance("Print to binary (on/off)");
    fPrintBinaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fFarFieldSourceCmd = new G4UIcmdWithoutParameter("/PRISM/source/farfieldsource", this);
    fFarFieldSourceCmd->SetGuidance("Source is at infinity, parallel rays");
    fFarFieldSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fFarFieldRingSourceCmd = new G4UIcmdWithoutParameter("/PRISM/source/farfieldringsource", this);
    fFarFieldRingSourceCmd->SetGuidance("Ring source is at infinity, parallel rays");
    fFarFieldRingSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fNearFieldSourceCmd = new G4UIcmdWithADoubleAndUnit("/PRISM/source/nearfieldsource",this);
    fNearFieldSourceCmd->SetGuidance("Near field point source, provide distance in cm");
    fNearFieldSourceCmd->SetParameterName("NearFieldDist",false);
    fNearFieldSourceCmd->SetDefaultUnit("cm");
    fNearFieldSourceCmd->SetUnitCandidates("cm m");
    fNearFieldSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    //fSourceStrengthCmd = new G4UIcmdWithADoubleAndUnit("/PRISM/source/sourcestrength",this);
    fSourceStrengthCmd = new G4UIcmdWithADouble("/PRISM/source/sourcestrength",this);
    fSourceStrengthCmd->SetGuidance("Source strength (in mCi)");
    fSourceStrengthCmd->SetParameterName("SrcStrength",false);
    //fSourceStrengthCmd->SetDefaultUnit("Ci");
    //fSourceStrengthCmd->SetUnitCandidates("Ci");
    fSourceStrengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAcqTimeCmd = new G4UIcmdWithAnInteger("/PRISM/source/acqtime",this);
    fAcqTimeCmd->SetGuidance("Acquisition time (in seconds)");
    fAcqTimeCmd->SetParameterName("Acqtime",false);
    fAcqTimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fRunCmd = new G4UIcmdWithoutParameter("/PRISM/source/run",this);
    fRunCmd->SetGuidance("Calculate beamOn for specified scenario");
    fRunCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//==================================================================================================

PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{
    delete fHPCmd;
    delete fHPNsideCmd;
    delete fHPindexingCmd;
    delete fSetUpHEALPixCmd;
    delete fHPDirectory;
    delete fSrcDirectory;
    delete fOutDirectory;
    delete fOutputFileNameCmd;
    delete fPrintTextCmd;
    delete fPrintBinaryCmd;
    delete fFarFieldSourceCmd;
    delete fFarFieldRingSourceCmd;
    delete fNearFieldSourceCmd;
    delete fSourceStrengthCmd;
    delete fAcqTimeCmd;
    delete fRunCmd;

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

        if (HPNsideindex == 8 || HPNsideindex == 16){
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

        if (nside == 16 && (indexing == "nested" || indexing == "Nested" || indexing == "NESTED")){filename = "HEALpix/n=16/HEALPix_PhiTheta_Nested.txt";}
        if (nside == 16 && (indexing == "ring" || indexing == "Ring" || indexing == "RING")){filename = "HEALpix/n=16/HEALPix_PhiTheta_Ring.txt";}
        if (nside == 32 && (indexing == "nested" || indexing == "Nested" || indexing == "NESTED")){filename = "HEALpix/n=32/HEALPix_PhiTheta_Nested.txt";}
        if (nside == 32 && (indexing == "ring" || indexing == "Ring" || indexing == "RING")){filename = "HEALpix/n=32/HEALPix_PhiTheta_Ring.txt";}

        fPrimaryGeneratorAction->ReadInHEALPixAngles(filename);

    }

    else if (command == fOutputFileNameCmd){

        RunAction::Instance()->SetOutputFilename(newValue);

    }

    else if (command == fPrintTextCmd){

        if (newValue == "On" || newValue == "ON" || newValue == "on"){
            RunAction::Instance()->SetPrintText(true);
        }
        else if (newValue == "Off" || newValue == "OFF" || newValue == "off"){
            RunAction::Instance()->SetPrintText(false);
        }


    }

    else if (command == fPrintBinaryCmd){

        if (newValue == "On" || newValue == "ON" || newValue == "on"){
            RunAction::Instance()->SetPrintBinary(true);
        }
        else if (newValue == "Off" || newValue == "OFF" || newValue == "off"){
            RunAction::Instance()->SetPrintBinary(false);
        }

    }

    else if (command == fFarFieldSourceCmd){

        fPrimaryGeneratorAction->SetFarFieldSource(true);
        fPrimaryGeneratorAction->SetNearFieldSource(false);
        fPrimaryGeneratorAction->SetFarFieldRingSource(false);

    }

    else if (command == fFarFieldRingSourceCmd){

        fPrimaryGeneratorAction->SetFarFieldSource(false);
        fPrimaryGeneratorAction->SetNearFieldSource(false);
        fPrimaryGeneratorAction->SetFarFieldRingSource(true);

    }

    else if (command == fNearFieldSourceCmd){

        fPrimaryGeneratorAction->SetFarFieldSource(false);
        fPrimaryGeneratorAction->SetNearFieldSource(true);
        fPrimaryGeneratorAction->SetFarFieldRingSource(false);
        fPrimaryGeneratorAction->SetNearFieldSourceDist(fNearFieldSourceCmd->GetNewDoubleValue(newValue));

    }

    else if (command == fSourceStrengthCmd){

        fPrimaryGeneratorAction->SetSourceStrength(fSourceStrengthCmd->GetNewDoubleValue(newValue));

    }

    else if (command == fAcqTimeCmd){

        fPrimaryGeneratorAction->SetAcquisitionTime(fAcqTimeCmd->GetNewIntValue(newValue));

    }

    else if (command == fRunCmd){

        G4double strength = fPrimaryGeneratorAction->GetSourceStrength();
        G4int acqtime     = fPrimaryGeneratorAction->GetAcquisitionTime();

        G4double relsolidang = 1.;
        if (fPrimaryGeneratorAction->GetNearFieldSource()){
            G4double srcdist  = fPrimaryGeneratorAction->GetNearFieldSourceDist();
            relsolidang = pow((10.*cm / (2. * srcdist)), 2);
        }
        else if (fPrimaryGeneratorAction->GetFarFieldSource()){
            relsolidang = 1. / (4.* CLHEP::pi);
        }

        G4double numparticles = floor((strength * 3.7*pow(10,7)) * acqtime * relsolidang);

        // beam on command
        cout << "Beam On! (" << numparticles << " particles!)\n";
        stringstream cmd_;
        cmd_ << "/run/beamOn " << numparticles;
        UI->ApplyCommand( cmd_.str() );

    }
}

//==================================================================================================


// Write GetCurrentValue functions...
