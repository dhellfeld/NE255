#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"

class StackingAction : public G4UserStackingAction {
    
  public:
    StackingAction();
    virtual ~StackingAction();
    
    static StackingAction* Instance();
     
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
    
private:
    static StackingAction* fgInstance;

    
};

#endif

