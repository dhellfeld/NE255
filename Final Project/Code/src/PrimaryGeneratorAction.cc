#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "DetectorConstruction.hh"

using namespace CLHEP;
using namespace std;

//==================================================================================================

PrimaryGeneratorAction* PrimaryGeneratorAction::fgInstance = 0;

//==================================================================================================

PrimaryGeneratorAction* PrimaryGeneratorAction::Instance() {

	return fgInstance;

}

//==================================================================================================

PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), HPindexing("Ring"), HP_Nside(8), gun(0), detector(0) {

	fgInstance = this;
    G4int nparticles = 1;
    gun = new G4ParticleGun(nparticles);

    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = ptable->FindParticle("gamma");
    if(!particle) G4Exception("PrimaryGeneratorAction", "null pointer to particle",
                              FatalException, "Find particle did not return a particle");

    gun->SetParticleDefinition(particle);
    gun->SetParticlePosition(G4ThreeVector(0.*cm, 0.*cm, 0.5*m));
    gun->SetParticleEnergy(60.*keV);
    gun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

    theta = 0.;
    phi = 0.;

	farfieldsource = true;
	nearfieldsource = false;
	nearfieldsourcedist = 50*cm;

    // Create a new messenger class
    primarygeneratoractionmessenger = new PrimaryGeneratorActionMessenger(this);

    // Pull in HEALPix angles
    ReadInHEALPixAngles("HEALpix/n=8/HEALPix_PhiTheta_Ring.txt");

}

//==================================================================================================

PrimaryGeneratorAction::~PrimaryGeneratorAction() {

    // delete the particle gun
    delete gun;

    detector = 0;

	fgInstance = 0;

    // Delete messenger class
    delete primarygeneratoractionmessenger;
}

//==================================================================================================

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    if(!detector) detector = DetectorConstruction::Instance();

	SourceInfo sourceinfo;
	if (farfieldsource){
		sourceinfo = FarFieldSource(theta, phi);
	}
	else if (nearfieldsource){
		sourceinfo = NearFieldSource(theta, phi, nearfieldsourcedist);
	}

    // Set results to the gun
    gun->SetParticlePosition(sourceinfo.GetPos());
    gun->SetParticleMomentumDirection(sourceinfo.GetDir());

    // Use the particle gun to define a primary particle for the event
    gun->GeneratePrimaryVertex(anEvent);
}

//==================================================================================================

SourceInfo PrimaryGeneratorAction::FarFieldSource(G4double theta_, G4double phi_){

	// Set the direction of the rays, given theta and phi (calculate position, take all negative values)
	G4double r = 15.*cm;
	G4ThreeVector dir;
	dir.setX(-r*cos(theta_*(pi/180))*sin(phi_*(pi/180)));
	dir.setY(-r*sin(theta_*(pi/180))*sin(phi_*(pi/180)));
	dir.setZ(-r*cos(phi_  *(pi/180)));

	// Uniformly sample a disk that just covers the entire detector (start at z-axis (0,0,1) and then rotate)
	G4double rand_r = G4UniformRand();
	G4double rand_theta = 2. * pi * G4UniformRand();
	G4double sphererad = 8.*cm;
	G4ThreeVector disk_pos(0);
	disk_pos.setX(sphererad * sqrt(rand_r) * cos(rand_theta));
	disk_pos.setY(sphererad * sqrt(rand_r) * sin(rand_theta));
	disk_pos.setZ(r);

	// Get position by rotating the z-oriented disk
	G4ThreeVector pos = disk_pos.rotateY(phi_*deg).rotateZ(theta_*deg);

	SourceInfo sourceinfo_;
	sourceinfo_.SetDir(dir);
	sourceinfo_.SetPos(pos);

	return sourceinfo_;

}

//==================================================================================================

SourceInfo PrimaryGeneratorAction::NearFieldSource(G4double theta_, G4double phi_, G4double dist_){

	// Isotropic rays
	//G4ThreeVector dir = GetIsotropicMomentumDirection();
	G4double ang = atan(9.*cm / dist_) * 180. / CLHEP::pi;
	G4ThreeVector dir = (GetConeMomentumDirection(ang)).rotateY(phi_*deg).rotateZ(theta_*deg);

	// Get position using theta, phi, and distance
	G4ThreeVector pos(dist_*cos(theta_*(pi/180))*sin(phi_*(pi/180)),
					  dist_*sin(theta_*(pi/180))*sin(phi_*(pi/180)),
					  dist_*cos(phi_  *(pi/180)));

	SourceInfo sourceinfo_;
	sourceinfo_.SetDir(dir);
	sourceinfo_.SetPos(pos);

	return sourceinfo_;

}

//==================================================================================================

void PrimaryGeneratorAction::SetTheta(G4double theta_){

    theta = theta_;

}
//==================================================================================================


//==================================================================================================

void PrimaryGeneratorAction::SetPhi(G4double phi_){

    phi = phi_;

}

//==================================================================================================

void PrimaryGeneratorAction::SetHP_index(G4int hpindex){

    HP_index = hpindex;
}

//==================================================================================================

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
