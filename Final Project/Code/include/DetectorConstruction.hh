#ifndef detectorconstruction_hh
#define detectorconstruction_hh

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "SensitiveDetector.hh"
#include "DetectorConstructionMessenger.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4TwoVector.hh"

using namespace std;

class G4Box;
class G4Tubs;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
	DetectorConstruction();
	~DetectorConstruction();

	static DetectorConstruction* Instance();

    DetectorConstructionMessenger * detectorconstructionmessenger;


public:
	G4VPhysicalVolume* Construct();

    inline const G4double& GetWorldDimensions() const { return world_dim; }
    inline const G4ThreeVector& GetTargetPosition() const { return detector_pos; }
    inline const G4ThreeVector& GetTargetDimension() const { return detector_dim; }

    virtual void ConstructSDandField();

    vector<G4int> GetRandomMask();
	vector<G4int> GetFullMask();
    void SetMask(vector<G4int>);
    G4String BinToHex(vector<G4int>);
    vector<G4int> HexToBin(G4String);

    void UpdateGeometry();
    void SetDetDim(G4ThreeVector);

    void CheckOverlapsOn();

    vector<G4ThreeVector> centers;
    inline vector<G4ThreeVector> GetDetCenters(){return centers;}

    vector<G4RotationMatrix> rotationmat;
    inline vector<G4RotationMatrix> GetRotationMat(){return rotationmat;}


    G4String detindexing;
    inline void SetDetIndexing(G4String di){detindexing = di;}
    inline G4String GetDetIndexing(){return detindexing;}


protected:
    virtual G4VPhysicalVolume* ConstructWorld();
    virtual void ConstructMaterials();

protected:
    G4VPhysicalVolume* worldPhys;
    G4Material* mworld;
    G4Material* mdetector;

    G4double world_dim;          // World is a sphere so we need radius
    G4ThreeVector detector_dim;  // detector are cubes, need xyz dimensions
    G4ThreeVector detector_pos;
    G4RotationMatrix detector_rot;

    // bool to check overlapping geometry (can be time consuming... default is false)
    //   can be turned on with CheckOverlapsOn()
    bool _checkoverlaps;


private:
	static DetectorConstruction* fgInstance;

    vector<G4int> _mask;


};

#endif
