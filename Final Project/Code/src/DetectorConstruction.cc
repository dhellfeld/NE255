#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"

#include "G4RunManager.hh"



using namespace std;
using namespace CLHEP;


//==================================================================================================

DetectorConstruction* DetectorConstruction::fgInstance = 0;

//==================================================================================================

DetectorConstruction* DetectorConstruction::Instance() {
    
    return fgInstance;
}

// ====================================================================================================

DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction(){

    //SetMask(GetRandomMask());
    SetMask(HexToBin("21524FA478BD521AB44791322B545C979943A029753854BB"));
    //SetMask(HexToBin("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    
    // Create a new messenger class
    detectorconstructionmessenger = new DetectorConstructionMessenger(this);
    
    fgInstance = this;

}

// ====================================================================================================

DetectorConstruction::~DetectorConstruction(){
    
    // Delete messenger class
    delete detectorconstructionmessenger;
    
    fgInstance = 0;

}

// ====================================================================================================

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

// ====================================================================================================

void DetectorConstruction::DefineMaterials()
{ 
    // Define CZT material
    G4Element* elCd = new G4Element("Cadmium"  ,"Cd", 48., 112.411*g/mole);
    G4Element* elZn = new G4Element("Zinc"     ,"Zn", 30., 65.39  *g/mole);
    G4Element* elTe = new G4Element("Tellurium","Te", 52., 127.6  *g/mole);
    
    G4Material* CZT = new G4Material("CZT", 5.78*g/cm3, 3);
    CZT->AddElement(elCd, 9);
    CZT->AddElement(elZn, 1);
    CZT->AddElement(elTe, 10);
    
    // the NIST manager has simple materials
    G4NistManager* nist = G4NistManager::Instance();
    nist->FindOrBuildMaterial("G4_Galactic");
    
}

// ====================================================================================================

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    
    // Geometry parameters
    G4double world_dim = 30*cm;
    G4ThreeVector detector_dim = G4ThreeVector(0.5*cm, 0.5*cm, 0.5*cm);
    G4ThreeVector detector_pos;
    
    
    // Get Materials
    G4Material* mdetector = G4Material::GetMaterial("CZT");
    G4Material* mworld    = G4Material::GetMaterial("G4_Galactic");
    
    
    // ---------------------------------------
    // Create world volume
    // ---------------------------------------
    
    G4VSolid* worldSolid = new G4Orb("World",  // name
                                     world_dim
                                     );
    
    G4LogicalVolume* worldLog = new G4LogicalVolume(worldSolid, // corresponding solid volume
                                                    mworld,     // material
                                                    "World"     // name
                                                    );
    
    G4VPhysicalVolume* worldPhys = new G4PVPlacement(0,       // no rotations so give it a null pointer
                                       G4ThreeVector(0.),     // the placement of the volume at (0, 0, 0)
                                       "World",               // name
                                       worldLog,              // the corresponding logical volume -- gives
                                                              //   volume the material (and the dimensions
                                                              //   via the solid assigned to the logical volume)
                                       0,                     // the mother physical volume (the world doesn't
                                                              //   have a mother so pass it a null pointer
                                       false,                 // many = false (no copies)
                                       0                      // the replica id number (only > 0 if copies exist)
                                       );
    
    
    
    // ---------------------------------------
    // Create detectors
    // ---------------------------------------
    
    G4VSolid* detectorSolid = new G4Box("DetectorSolid",      // name
                                        detector_dim.x(),     // x half-length
                                        detector_dim.y(),     // y half-length
                                        detector_dim.z()      // z half-length
                                        );
    
     G4LogicalVolume* detectorLog = new G4LogicalVolume(detectorSolid,    // target solid
                                        mdetector,                        // target material
                                        "DetectorLog"                     // name
                                        );

    
    // Pull in detector center verticies from file
    G4String centfilename = "geo/ring/centervertices_Ring.txt";
    ifstream centerfile(centfilename);
    G4double x_, y_, z_;
    string line;
    centers.clear();
    if (centerfile.is_open()){
        while (getline(centerfile,line)){
            istringstream(line) >> x_ >> y_ >> z_;
            centers.push_back(G4ThreeVector(x_*mm, y_*mm ,z_*mm));
        }
        centerfile.close();
    }
    
    // Pull in rotation matrices from file
    G4String rotfilename = "geo/ring/rotationmatrices_Ring.txt";
    ifstream myfile3(rotfilename);
    G4double x1_,x2_,x3_,y1_,y2_,y3_,z1_,z2_,z3_;
    rotationmat.clear();
    if (myfile3.is_open()){
        while (getline(myfile3,line)){
            istringstream(line) >> x1_ >> x2_ >> x3_ >> y1_ >> y2_ >> y3_ >> z1_ >> z2_ >> z3_;
            rotationmat.push_back(G4RotationMatrix(G4ThreeVector(x1_,x2_,x3_), G4ThreeVector(y1_,y2_,y3_), G4ThreeVector(z1_,z2_,z3_)));
        }
        myfile3.close();
    }
    
    
    // ---------------------------------------
    // Place masked detectors (192 slots)
    // ---------------------------------------
    
    // Loop through mask vector
    for (int i = 0; i < 192; i++){
        
        // Detctor ID are just numbers
        G4String detID = to_string(i+1);
        
        // if mask value is 1 then fill, if not leave empty
        if (_mask[i]){
            
            detector_pos = centers[i];
            G4Transform3D transform = G4Transform3D(rotationmat[i],detector_pos);
            
            new G4PVPlacement(transform,                // the placement of the volume in center
                              detID,                    // name
                              detectorLog,              // the corresponding logical volume -- gives
                                                        //   volume the material (and the dimensions
                                                        //   via the solid assigned to the logical volume)
                              worldPhys,                // inside the world so the world is the mother
                                                        //   physical volume,
                              false,                    // many = false (no copies)
                              0                         // the replica id number (only > 0 if copies exist)
                              );

        }
    }
    
    
    // ---------------------------------------
    // Set up visualization of geometry
    // ---------------------------------------
    
    G4VisAttributes* detector_vis_att = new G4VisAttributes(G4Color(1.0,1.0,1.0));
    detector_vis_att->SetForceSolid(true);
    detector_vis_att->SetVisibility(true);
    detectorLog->SetVisAttributes(detector_vis_att);

    
    G4VisAttributes* world_vis_att = new G4VisAttributes(G4Color(1.,0.,0.));
    world_vis_att->SetForceWireframe(true);
    world_vis_att->SetVisibility(true);
    worldLog->SetVisAttributes(world_vis_att);
    
    
   
    // Return the physical world
    return worldPhys;

}

//==================================================================================================


vector<G4int> DetectorConstruction::GetRandomMask(){
    
    vector<G4int> mask;
    G4double rand;
    
    // Get random number, convert to 1 or 0, fill mask array
    for (int q = 0; q < 192; q++){
        rand = G4UniformRand();
        if (rand > 0.5){mask.push_back(1);}
        else {mask.push_back(0);}
    }
    
    return mask;
}

//==================================================================================================

void DetectorConstruction::SetMask(vector<G4int> inputmask){
    
    _mask = inputmask;
}

//==================================================================================================

G4String DetectorConstruction::BinToHex(vector<G4int> mask){
    
    // Digit to hex table
    const char *hex_dig[] = {"0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f"};
    
    G4String hex_;
    
    // Loop through mask arry, convert chucks of 4 into respective hex form, append to string
    for (int q = 0; q < int(mask.size()/4); q++){
        int idx = 8*mask[0+(q*4)] + 4*mask[1+(q*4)] + 2*mask[2+(q*4)] + mask[3+(q*4)];
        hex_.append(hex_dig[idx]);
    }
    
    return hex_;
}

//==================================================================================================

vector<G4int> DetectorConstruction::HexToBin(G4String hex_){
    
    vector<G4int> mask;
    
    // Loop through hex string and create a binary string with hex to binary table
    G4String mask_string = "";
    for (int i = 0; i < int(hex_.length()); i++)
    {
        switch (hex_[i])
        {
            case '0': mask_string.append ("0000"); break;
            case '1': mask_string.append ("0001"); break;
            case '2': mask_string.append ("0010"); break;
            case '3': mask_string.append ("0011"); break;
            case '4': mask_string.append ("0100"); break;
            case '5': mask_string.append ("0101"); break;
            case '6': mask_string.append ("0110"); break;
            case '7': mask_string.append ("0111"); break;
            case '8': mask_string.append ("1000"); break;
            case '9': mask_string.append ("1001"); break;
            case 'a': mask_string.append ("1010"); break;
            case 'b': mask_string.append ("1011"); break;
            case 'c': mask_string.append ("1100"); break;
            case 'd': mask_string.append ("1101"); break;
            case 'e': mask_string.append ("1110"); break;
            case 'f': mask_string.append ("1111"); break;
            case 'A': mask_string.append ("1010"); break;
            case 'B': mask_string.append ("1011"); break;
            case 'C': mask_string.append ("1100"); break;
            case 'D': mask_string.append ("1101"); break;
            case 'E': mask_string.append ("1110"); break;
            case 'F': mask_string.append ("1111"); break;
        }
    }
    
    // Convert string into an int array
    for (int i = 0; i < 192; i++){mask.push_back((G4int)mask_string[i]-48);}
    
    return mask;
}

//==================================================================================================

void DetectorConstruction::UpdateGeometry(){
    
    //G4RunManager::GetRunManager()->DefineWorldVolume(ConstructWorld());
    G4RunManager::GetRunManager()->ReinitializeGeometry();   // pass it "true" so it clears the geometry first...
    //G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//==================================================================================================
