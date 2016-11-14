#include "Hit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

using namespace std;

G4ThreadLocal G4Allocator<Hit>* HitAllocator = 0;

//==================================================================================================

Hit::Hit()
: G4VHit(),
  fTrackID(-1),
  fHP(0),
  fProc(G4int()),
  fVol(G4int()),
  fEnergy(0.),
  fDOI(0.)
{}

//==================================================================================================

Hit::~Hit() {}

//==================================================================================================

Hit::Hit(const Hit& right): G4VHit(){

    fTrackID   = right.fTrackID;
    fEnergy    = right.fEnergy;
    fVol       = right.fVol;
    fProc      = right.fProc;
    fDOI       = right.fDOI;
    fHP        = right.fHP;
}

//==================================================================================================

const Hit& Hit::operator=(const Hit& right){

    fTrackID   = right.fTrackID;
    fEnergy    = right.fEnergy;
    fVol       = right.fVol;
    fProc      = right.fProc;
    fDOI       = right.fDOI;
    fHP        = right.fHP;

    return *this;
}

//==================================================================================================

G4int Hit::operator==(const Hit& right) const{

    return ( this == &right ) ? 1 : 0;
}

//==================================================================================================
