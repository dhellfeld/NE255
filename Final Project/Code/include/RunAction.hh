#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "globals.hh"

class G4Run;
class PrimaryGeneratorAction;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    
    int response[12288][1920] = {};  // allow up tp nside=32 and 10 doi bins
    inline int PrintResponse(G4int a, G4int b){return response[a][b];}
    inline void FillResponse(G4int a, G4int b){response[a][b] += 1;}
    
    G4int evt;
    inline void SetEventID(G4int x){evt = x;}
    inline G4int GetEventID(){return evt;}
    
    static RunAction* Instance();
    
private:
    static RunAction* fgInstance;
};


#endif

