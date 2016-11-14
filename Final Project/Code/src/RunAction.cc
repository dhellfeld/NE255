#include "RunAction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "SensitiveDetector.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"

using namespace CLHEP;
using namespace std;


//==================================================================================================

RunAction* RunAction::fgInstance = 0;

//==================================================================================================

RunAction* RunAction::Instance() {

	return fgInstance;
}

//==================================================================================================

RunAction::RunAction(): printtext(false), printbin(true), outputfilename("output") {

	fgInstance = this;

}

//==================================================================================================

RunAction::~RunAction() {

	fgInstance = 0;
}

//==================================================================================================

void RunAction::BeginOfRunAction(const G4Run* /*aRun*/) {}

//==================================================================================================

void RunAction::EndOfRunAction(const G4Run* /*aRun*/) {

    // Dump data to file
    if (printbin){PrintToBinaryFile();}
    if (printtext){PrintToTextFile();}

    // Clear data tuples for next run
    ClearTuples();

}

//==================================================================================================
// Fill functions

void RunAction::FillDOIbintuple     (G4int doi)         {DOIbintuple.push_back(doi); }
void RunAction::FillDetIDtuple      (G4int detid)       {DetIDtuple.push_back(detid); }
void RunAction::FillEvtNtuple       (G4int evtN)        {EvtNtuple.push_back(evtN); }
void RunAction::FillHitNumtuple     (G4int hitnum)      {HitNumtuple.push_back(hitnum); }
void RunAction::FillTrackIDtuple    (G4int trackid)     {TrackIDtuple.push_back(trackid); }
void RunAction::FillEnergytuple     (G4float energy)    {Energytuple.push_back(energy); }
void RunAction::FillProcesstuple    (G4int proc)        {Processtuple.push_back(proc); }
void RunAction::FillHPindextuple    (G4int hpindex)     {HPindextuple.push_back(hpindex); }

//==================================================================================================
// Clear functions

void RunAction::ClearDOIbintuple()      {DOIbintuple.clear(); }
void RunAction::ClearDetIDtuple()       {DetIDtuple.clear(); }
void RunAction::ClearEvtNtuple()        {EvtNtuple.clear(); }
void RunAction::ClearHitNumtuple()      {HitNumtuple.clear(); }
void RunAction::ClearTrackIDtuple()     {TrackIDtuple.clear(); }
void RunAction::ClearEnergytuple()      {Energytuple.clear(); }
void RunAction::ClearProcesstuple()     {Processtuple.clear(); }
void RunAction::ClearHPindextuple()     {HPindextuple.clear(); }

//==================================================================================================
// Get functions

vector<G4int>         RunAction::GetDOIbintuple()      {return DOIbintuple; }
vector<G4int>         RunAction::GetDetIDtuple()       {return DetIDtuple; }
vector<G4int>         RunAction::GetEvtNtuple()        {return EvtNtuple; }
vector<G4int>         RunAction::GetHitNumtuple()      {return HitNumtuple; }
vector<G4int>         RunAction::GetTrackIDtuple()     {return TrackIDtuple; }
vector<G4float>       RunAction::GetEnergytuple()      {return Energytuple; }
vector<G4int>         RunAction::GetProcesstuple()     {return Processtuple; }
vector<G4int>         RunAction::GetHPindextuple()     {return HPindextuple; }

//==================================================================================================

void RunAction::PrintToTextFile(){

    G4String filename = "output/";
    filename.append(GetOutputFilename());
    filename.append(".txt");

    ofstream myfile;
    myfile.open (filename, ios::out | ios::app);    // use append so we can write from mulitple runs

    for (int i = 0; i < int(EvtNtuple.size()); i++){
        myfile
        << EvtNtuple[i] << "\t"
        << HitNumtuple[i] << "\t"
        << TrackIDtuple[i] << "\t"
        << Energytuple[i] << "\t"
        << DetIDtuple[i] << "\t"
        << Processtuple[i] << "\t"
        << DOIbintuple[i] << "\t"
        << HPindextuple[i] << "\t"
        << "\n";
    }
    myfile.close();
}

//==================================================================================================

void RunAction::PrintToBinaryFile(){

    G4String filename = "output/";
    filename.append(GetOutputFilename());
    filename.append(".bin");

    ofstream myfile;
    myfile.open (filename, ios::out | ios::app | ios::binary);   // use append so we can write from mulitple runs

    for (int i = 0; i < int(EvtNtuple.size()); i++){

        myfile.write(reinterpret_cast<char*>(&EvtNtuple[i]),    4);  // 4 Bytes (1 to 1e9)
        myfile.write(reinterpret_cast<char*>(&HitNumtuple[i]),  1);  // 1 Byte  (1 to 5)
        myfile.write(reinterpret_cast<char*>(&TrackIDtuple[i]), 1);  // 1 Byte  (1 to 10)
        myfile.write(reinterpret_cast<char*>(&Energytuple[i]),  4);  // 4 Bytes (float)
        myfile.write(reinterpret_cast<char*>(&DetIDtuple[i]),   1);  // 1 Byte  (1 to 192)
        myfile.write(reinterpret_cast<char*>(&Processtuple[i]), 1);  // 1 Byte  (1 to 4)
        myfile.write(reinterpret_cast<char*>(&DOIbintuple[i]),  1);  // 1 Byte  (1 to 20)
        myfile.write(reinterpret_cast<char*>(&HPindextuple[i]), 2);  // 2 Bytes (1 to 768 or 3072)

    }
    myfile.close();
}

//==================================================================================================

void RunAction::ClearTuples(){

    ClearDOIbintuple();
    ClearDetIDtuple();
    ClearEvtNtuple();
    ClearHitNumtuple();
    ClearTrackIDtuple();
    ClearEnergytuple();
    ClearProcesstuple();
    ClearHPindextuple();

}

//==================================================================================================
