#ifndef runaction_hh
#define runaction_hh

#include "globals.hh"
#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"
#include <vector>

using namespace std;

class G4Run;

class RunAction : public G4UserRunAction
{
public:
	RunAction();
	~RunAction();

	static RunAction* Instance();

public:
	void BeginOfRunAction(const G4Run*);
	void EndOfRunAction(const G4Run*);


    // -----------------------------------------------
    // Analysis tuples

    // Depth of interaction
    vector<G4double> DOItuple;
    void FillDOItuple(G4double);
    void ClearDOItuple();
    vector<G4double> GetDOItuple();

    // Depth of interaction (binned)
    vector<G4int> DOIbintuple;
    void FillDOIbintuple(G4int);
    void ClearDOIbintuple();
    vector<G4int> GetDOIbintuple();

    // Detector ID
    vector<G4int> DetIDtuple;
    void FillDetIDtuple(G4int);
    void ClearDetIDtuple();
    vector<G4int> GetDetIDtuple();

    // Event Number
    vector<G4int> EvtNtuple;
    void FillEvtNtuple(G4int);
    void ClearEvtNtuple();
    vector<G4int> GetEvtNtuple();

    // Hit Number
    vector<G4int> HitNumtuple;
    void FillHitNumtuple(G4int);
    void ClearHitNumtuple();
    vector<G4int> GetHitNumtuple();

    // Track ID
    vector<G4int> TrackIDtuple;
    void FillTrackIDtuple(G4int);
    void ClearTrackIDtuple();
    vector<G4int> GetTrackIDtuple();

    // Energy
    vector<G4float> Energytuple;
    void FillEnergytuple(G4float);
    void ClearEnergytuple();
    vector<G4float> GetEnergytuple();

    // Process
    vector<G4int> Processtuple;
    void FillProcesstuple(G4int);
    void ClearProcesstuple();
    vector<G4int> GetProcesstuple();

    // HP Index
    vector<G4int> HPindextuple;
    void FillHPindextuple(G4int);
    void ClearHPindextuple();
    vector<G4int> GetHPindextuple();

	// Time
    vector<G4float> Timetuple;
    void FillTimetuple(G4float);
    void ClearTimetuple();
    vector<G4float> GetTimetuple();

    // Clear all tuples
    void ClearTuples();

    // -----------------------------------------------

    // Print to file
    void PrintToTextFile();
    G4bool printtext;
    inline void SetPrintText(G4bool ww){printtext = ww;}
    inline G4bool GetPrintText(){return printtext;}

    void PrintToBinaryFile();
    G4bool printbin;
    inline void SetPrintBinary(G4bool qq){printbin = qq;}
    inline G4bool GetPrintBinary(){return printbin;}

    G4String outputfilename;
    inline void SetOutputFilename(G4String fn){outputfilename = fn;}
    inline G4String GetOutputFilename(){return outputfilename;}


private:
	static RunAction* fgInstance;

};

#endif
