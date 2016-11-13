#ifndef stackingaction_hh
#define stackingaction_hh

#include "G4UserStackingAction.hh"
#include "StackingActionMessenger.hh"
#include "globals.hh"

class StackingAction : public G4UserStackingAction {

  public:
    StackingAction();
    virtual ~StackingAction();

    static StackingAction* Instance();

    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);

    // variable to turn on/off the tracking of secondary electrons
    G4bool electrontracking;
    inline void SetElectronTracking(G4bool gt){electrontracking = gt;}
    inline G4bool GetElectronTracking(){return electrontracking;}

    // intialize the messenger pointer
    StackingActionMessenger * stackingactionmessenger;

private:
    static StackingAction* fgInstance;


};

#endif
