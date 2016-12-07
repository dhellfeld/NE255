
#ifndef primarygeneratoractionmessenger_hh
#define primarygeneratoractionmessenger_hh

#include "globals.hh"
#include "G4UIcommand.hh"
#include "G4UImessenger.hh"

#include "G4UImanager.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class PrimaryGeneratorActionMessenger : public G4UImessenger {

public:
    PrimaryGeneratorActionMessenger(PrimaryGeneratorAction* generator);
    ~PrimaryGeneratorActionMessenger();

    void SetNewValue(G4UIcommand * command,G4String newValue);

    G4UImanager* UI;

private:
    PrimaryGeneratorAction* fPrimaryGeneratorAction;

    G4UIdirectory* fHPDirectory;
    G4UIdirectory* fSrcDirectory;
    G4UIdirectory* fOutDirectory;

    G4UIcmdWithAnInteger*      fHPCmd;
    G4UIcmdWithAnInteger*      fHPNsideCmd;
    G4UIcmdWithAString*        fHPindexingCmd;
    G4UIcmdWithoutParameter*   fSetUpHEALPixCmd;
    G4UIcmdWithAString*        fOutputFileNameCmd;
    G4UIcmdWithAString*        fPrintTextCmd;
    G4UIcmdWithAString*        fPrintBinaryCmd;
    G4UIcmdWithoutParameter*   fFarFieldSourceCmd;
    G4UIcmdWithoutParameter*   fFarFieldRingSourceCmd;
    G4UIcmdWithADoubleAndUnit* fNearFieldSourceCmd;
    G4UIcmdWithADouble*        fSourceStrengthCmd;
    G4UIcmdWithAnInteger*      fAcqTimeCmd;
    G4UIcmdWithoutParameter*   fRunCmd;

};

#endif
