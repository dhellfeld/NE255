#ifndef physicslist_hh
#define physicslist_hh 1

#include "globals.hh" 
#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"


class PhysicsList : public G4VUserPhysicsList
{
public:
	PhysicsList();
	~PhysicsList();

	static PhysicsList* Instance();

public:
	void ConstructParticle();
	void ConstructProcess();
    void ConstructEM();
	void SetCuts();


private:
	static PhysicsList* fgInstance;




};

#endif
