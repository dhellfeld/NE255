#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "PrimaryGeneratorActionMessenger.hh"

using namespace std;

class G4ParticleGun;
class G4Event;

struct Angles{
    G4double theta;
    G4double phi;
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();    
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);

  void SetRandomFlag(G4bool value);

private:
  G4ParticleGun*  fParticleGun;

public:
    PrimaryGeneratorActionMessenger * primarygeneratoractionmessenger;
    
    // Healpix
    G4int HP_index;
    inline void SetHP_index(G4int hpi){HP_index = hpi;}
    inline G4int GetHP_index(){return HP_index;}
    
    vector<struct Angles> HPangles;
    void ReadInHEALPixAngles(G4String);
    inline vector<struct Angles> GetHPangles(){return HPangles;}
    
    G4String HPindexing;
    inline void SetHPindexing(G4String _hpindexing){HPindexing = _hpindexing;}
    inline G4String GetHPindexing(){return HPindexing;}
    
    G4int HP_Nside;
    inline void SetHPNside(G4int _hpNside){HP_Nside = _hpNside;}
    inline G4int GetHPNside(){return HP_Nside;}
    
    void SetTheta(G4double);
    void SetPhi(G4double);
    
    G4double theta;
    G4double phi;
    
    //
    
    G4int numpart;
    inline void SetNumPart(G4int x){numpart = x;}
    inline G4int GetNumPart(){return numpart;}

    G4int evt;
    G4int hp;
    
    G4int hpstart, hpfinish;
    inline void SetHPStart(G4int x){hpstart = x;}
    inline void SetHPFinish(G4int x){hpfinish = x;}
    inline G4int GetHPStart(){return hpstart;}
    inline G4int GetHPFinish(){return hpfinish;}
    
    
    //
    
    G4bool doisetting;
    inline G4bool GetDOIsetting(){return doisetting;}
    inline void SetDOIsetting(G4bool a){doisetting = a;}
    
    //
    
    
    static PrimaryGeneratorAction* Instance();

private:
    static PrimaryGeneratorAction* fgInstance;

};

#endif
