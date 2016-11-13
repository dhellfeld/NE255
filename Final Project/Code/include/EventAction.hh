#ifndef eventaction_hh
#define eventaction_hh

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "Hit.hh"


class EventAction : public G4UserEventAction
{
public:
	EventAction();
	~EventAction();

	static EventAction* Instance();


public:
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

    void FillTuples(const G4Event*);


private:
	static EventAction* fgInstance;
    G4int HCollID;

};

#endif
