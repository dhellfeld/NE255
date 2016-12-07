#ifndef primarygeneratoraction_hh
#define primarygeneratoraction_hh

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorActionMessenger.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

using namespace std;
using namespace CLHEP;

class G4Event;
class G4ParticleGun;
class DetectorConstruction;

struct Angles{
    G4double theta;
    G4double phi;
};

class SourceInfo {
    G4ThreeVector dir;
	G4ThreeVector pos;
  public:
	void SetDir(G4ThreeVector dir_){ dir = dir_;}
	void SetPos(G4ThreeVector pos_){ pos = pos_;}
    G4ThreeVector GetDir() {return dir;}
	G4ThreeVector GetPos() {return pos;}
};


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	PrimaryGeneratorAction();
	~PrimaryGeneratorAction();

	static PrimaryGeneratorAction* Instance();

    PrimaryGeneratorActionMessenger * primarygeneratoractionmessenger;

public:
	void GeneratePrimaries(G4Event*);

    void SetTheta(G4double);
    void SetPhi(G4double);

    SourceInfo FarFieldSource(G4double, G4double);
    SourceInfo NearFieldSource(G4double, G4double, G4double);
    SourceInfo FarFieldRingSource();
    G4bool farfieldsource;
    G4bool nearfieldsource;
    G4double nearfieldsourcedist;
    G4bool farfieldringsource;
    inline void SetFarFieldSource(G4bool a){farfieldsource = a;}
    inline G4bool GetFarFieldSource(){return farfieldsource;}
    inline void SetNearFieldSource(G4bool a ){nearfieldsource = a;}
    inline G4bool GetNearFieldSource(){return nearfieldsource;}
    inline void SetNearFieldSourceDist(G4double a ){nearfieldsourcedist = a;}
    inline G4double GetNearFieldSourceDist(){return nearfieldsourcedist;}
    inline void SetFarFieldRingSource(G4bool a){farfieldringsource = a;}
    inline G4bool GetFarFieldRingSource(){return farfieldringsource;}

    G4int acqtime;
    G4int GetAcquisitionTime(){return acqtime;}
    void SetAcquisitionTime(G4int a){acqtime = a;}

    G4double srcstrength;
    G4double GetSourceStrength(){return srcstrength;}
    void SetSourceStrength(G4double a){srcstrength = a;}


    inline G4double GetTheta(){return theta;}
    inline G4double GetPhi(){return phi;}


    // Healpix
    G4int HP_index;
    void SetHP_index(G4int);
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

    //

protected:
    inline G4ThreeVector GetIsotropicMomentumDirection() const;
    inline G4ThreeVector GetHalfIsotropicMomentumDirection() const;
    inline G4ThreeVector GetConeMomentumDirection(G4double) const;

    inline G4double GetRandTheta();
    inline G4double GetRandPhi();

private:
	static PrimaryGeneratorAction* fgInstance;

    G4double theta;
    G4double phi;

    G4ParticleGun* gun;
    DetectorConstruction* detector;

};

//----------------------------------------------------------------------------//
inline G4ThreeVector PrimaryGeneratorAction::GetIsotropicMomentumDirection() const {

    G4double u = 2*G4UniformRand()-1;
    G4double v = 2.*pi*G4UniformRand();
    return G4ThreeVector(sqrt(1-u*u)*cos(v),
                         sqrt(1-u*u)*sin(v),
                         u);
}
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
inline G4ThreeVector PrimaryGeneratorAction::GetHalfIsotropicMomentumDirection() const {

    G4double u = G4UniformRand()-1;
    G4double v = 2.*pi*G4UniformRand();
    return G4ThreeVector(sqrt(1-u*u)*cos(v),
                         sqrt(1-u*u)*sin(v),
                         u);
}
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
inline G4ThreeVector PrimaryGeneratorAction::GetConeMomentumDirection(G4double ang) const {

    // open a cone of degree theta
    G4double _theta = ang;

    G4double u = (-cos(_theta*pi/180)+1)*G4UniformRand()-1;
    G4double v = 2.*pi*G4UniformRand();
    return G4ThreeVector(sqrt(1-u*u)*cos(v),
                         sqrt(1-u*u)*sin(v),
                         u);
}
//----------------------------------------------------------------------------//

inline G4double PrimaryGeneratorAction::GetRandTheta(){

    return G4UniformRand()*360.;

}

//----------------------------------------------------------------------------//

inline G4double PrimaryGeneratorAction::GetRandPhi(){

    return acos(2.*G4UniformRand()-1.);
}

//----------------------------------------------------------------------------//





#endif
