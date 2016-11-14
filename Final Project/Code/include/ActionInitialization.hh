#ifndef actioninitialization_hh
#define actioninitialization_hh

#include "globals.hh"
#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization
{
public:
	ActionInitialization();
	virtual ~ActionInitialization();

public:
    virtual void BuildForMaster() const;
    virtual void Build() const;

private:


};

#endif