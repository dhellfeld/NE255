#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "DetectorConstructionMessenger.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

using namespace std;

class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    
    static DetectorConstruction* Instance();

  public:
    virtual G4VPhysicalVolume* Construct();
  
private:
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    
    static DetectorConstruction* fgInstance;

    
  public:
    vector<G4ThreeVector> centers;
    inline vector<G4ThreeVector> GetDetCenters(){return centers;}
    
    vector<G4RotationMatrix> rotationmat;
    inline vector<G4RotationMatrix> GetRotationMat(){return rotationmat;}
    
    vector<G4int> GetRandomMask();
    G4String BinToHex(vector<G4int>);
    vector<G4int> HexToBin(G4String);
    void SetMask(vector<G4int>);
    vector<G4int> _mask;
    
    DetectorConstructionMessenger * detectorconstructionmessenger;

    void UpdateGeometry();

};

#endif

