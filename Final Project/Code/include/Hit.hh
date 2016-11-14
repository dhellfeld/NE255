#ifndef hit_hh
#define hit_hh

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

using namespace std;

class Hit : public G4VHit
{
public:
    Hit();
    Hit(const Hit&);
    virtual ~Hit();

    // operators
    const Hit& operator=(const Hit&);
    G4int operator==(const Hit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);


    // Set methods
    void SetTrackID     (G4int track)          {fTrackID = track;}
    void SetEnergy      (G4float de)           {fEnergy = de;}
    void SetVol         (G4int volname)        {fVol = volname;}
    void SetProcess     (G4int procname)       {fProc = procname;}
    void SetDOI         (G4double doi)         {fDOI = doi;}
    void SetHPindex     (G4int hpindex)        {fHP = hpindex;}

    // Get methods
    G4int GetTrackID() const               {return fTrackID;}
    G4float GetEnergy() const              {return fEnergy;}
    G4int GetVol() const                   {return fVol;}
    G4int GetProcess() const               {return fProc;}
    G4double GetDOI() const                {return fDOI;}
    G4int GetHPindex() const               {return fHP;}

private:

    G4int         fTrackID, fHP, fProc, fVol;
    G4float       fEnergy;
    G4double      fDOI;

};


typedef G4THitsCollection<Hit> HitsCollection;

extern G4ThreadLocal G4Allocator<Hit>* HitAllocator;

inline void* Hit::operator new(size_t)
{
    if(!HitAllocator)
        HitAllocator = new G4Allocator<Hit>;
    return (void *) HitAllocator->MallocSingle();
}


inline void Hit::operator delete(void *hit)
{
    HitAllocator->FreeSingle((Hit*) hit);
}


#endif
