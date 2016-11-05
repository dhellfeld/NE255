#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"


class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step* step);
    
private:
    
public:
    int det;
    int hpi;
    
    int doi;
    int CalculateDOI(const G4Step*);

    
};


#endif
