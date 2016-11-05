#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

using namespace CLHEP;
using namespace std;

//==================================================================================================

PrimaryGeneratorAction* PrimaryGeneratorAction::fgInstance = 0;

//==================================================================================================

PrimaryGeneratorAction* PrimaryGeneratorAction::Instance() {
    
    return fgInstance;
    
}

// ====================================================================================================

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   HP_index(1),
   HPindexing("Ring"),
   HP_Nside(16),
   numpart(1),
   hpstart(1),
   hpfinish(3072),
   doisetting(true)
{
    
    fgInstance = this;
    
    G4int nparticles = 1;
    fParticleGun = new G4ParticleGun(nparticles);
    
    G4ParticleDefinition* particleDef= G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    
    fParticleGun->SetParticleDefinition(particleDef);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm, 0.*cm, 0.5*m));
    fParticleGun->SetParticleEnergy(60.*keV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    
    theta = 0.;
    phi = 0.;
    
    // Create a new messenger class
    primarygeneratoractionmessenger = new PrimaryGeneratorActionMessenger(this);
    
    // Pull in HEALPix angles
    ReadInHEALPixAngles("HEALpix/n=16/HEALPix_PhiTheta_Ring.txt");
        
}

// ====================================================================================================

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    
    fgInstance = 0;
    
    delete primarygeneratoractionmessenger;
}

// ====================================================================================================

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    //G4int hp = int(3072.*G4UniformRand());
    
    evt = RunAction::Instance()->GetEventID();
    
    hp = floor(evt/numpart)+hpstart;
    //G4cout << hp << "\n";
    
    SetHP_index(hp);
    
    SetTheta(HPangles[hp-1].theta*(180./pi));
    SetPhi(HPangles[hp-1].phi*(180./pi));
    
    
    // - -
    // Far field source at angle to origin
    // - -
    
    // Set the direction of the rays, given theta and phi (calculate position, take all negative values)
    G4double r = 15.*cm;
    G4ThreeVector dir;
    dir.setX(-r*cos(theta*(pi/180))*sin(phi*(pi/180)));
    dir.setY(-r*sin(theta*(pi/180))*sin(phi*(pi/180)));
    dir.setZ(-r*cos(phi  *(pi/180)));
    
    
    // Uniformly sample a disk that just covers the entire detector (start at z-axis (0,0,1) and then rotate)
    G4double rand_r = G4UniformRand();
    G4double rand_theta = 2. * pi * G4UniformRand();
    G4double sphererad = 8.*cm;
    G4ThreeVector disk_pos(0);
    disk_pos.setX(sphererad * sqrt(rand_r) * cos(rand_theta));
    disk_pos.setY(sphererad * sqrt(rand_r) * sin(rand_theta));
    disk_pos.setZ(r);
    
    // Get position by rotating the z-oriented disk
    G4ThreeVector pos = disk_pos.rotateY(phi*deg).rotateZ(theta*deg);
    
    // Set results to the gun
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(dir);
    
    // Use the particle gun to define a primary particle for the event
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
}

// ====================================================================================================


void PrimaryGeneratorAction::ReadInHEALPixAngles(G4String filename){
    
    HPangles.clear();
    
    ifstream myfile(filename);
    string line;
    if (myfile.is_open()){
        Angles ang;
        while (getline(myfile,line)){
            istringstream(line) >> ang.phi >> ang.theta;
            HPangles.push_back(ang);
        }
        myfile.close();
    }

}

//==================================================================================================

void PrimaryGeneratorAction::SetTheta(G4double theta_){
    
    theta = theta_;
    
}

//==================================================================================================

void PrimaryGeneratorAction::SetPhi(G4double phi_){
    
    phi = phi_;
    
}

//==================================================================================================


