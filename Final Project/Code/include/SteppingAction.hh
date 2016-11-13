#ifndef steppingaction_hh
#define steppingaction_hh

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class G4Step;

class SteppingAction : public G4UserSteppingAction
{
public:
	SteppingAction();
	~SteppingAction();

	static SteppingAction* Instance();

public:
	void UserSteppingAction(const G4Step*);

private:
	static SteppingAction* fgInstance;


};

#endif
