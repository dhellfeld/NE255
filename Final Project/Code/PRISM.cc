#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Analysis.hh"
#include "RunAction.hh"

#include <fstream>

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " examplea [-m macro ] [-u UIsession] [-t nThreads] [-f outputfilename]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}


#include "PhysicsList.hh"
typedef PhysicsList ThePhysicsList;



int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;
  G4String output;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-f" ) output = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }
  
  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
if ( nThreads > 0 ) {
    runManager->SetNumberOfThreads(nThreads);
}
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  DetectorConstruction* detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);


    
  runManager->SetUserInitialization(new ThePhysicsList);
  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  //runManager->SetUserInitialization(physicsList);
    
  ActionInitialization* actionInitialization
     = new ActionInitialization(detConstruction);
  runManager->SetUserInitialization(actionInitialization);
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute mac/init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute mac/gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

    
    
        
    
    /*
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->CloseFile();
     */
    
    //clock_t start;
    //start = clock();
    //G4cout << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";

    G4String filename = output;
    
    ofstream myfile;
    myfile.open (filename, ios::out | ios::binary);
    
    int dets;
    if (PrimaryGeneratorAction::Instance()->GetDOIsetting()){dets = 1920;}
    else {dets = 192;}
    
    int nside = PrimaryGeneratorAction::Instance()->GetHPNside();
    int hpi = 12*nside*nside;
    
    for (int i = 0; i < dets; i++){
        for (int j = 0; j < hpi; j++){
            
            G4int val = RunAction::Instance()->PrintResponse(j,i);
            myfile.write(reinterpret_cast<char*>(&val), 2);   // 2 Bytes (values up to 2^16)

        }
    }
    myfile.close();
    

    /*
    vector<struct Angles> angles = PrimaryGeneratorAction::Instance()->GetHPangles();
    
    for (int j = 0; j < hpi; j++){
    
        G4double theta_ = angles[j].theta;
        G4double phi_ = angles[j].phi;
     
        
        G4double r = 15.*cm;
        //G4cout << j+1 << "\t" << theta_ << ", " << phi_ << "\t\t" << (-r*cos(theta_)*sin(phi_)) << ", " << (-r*sin(theta_)*sin(phi_)) << ", " << (-r*cos(phi_)) << "\n";
        G4cout << (-r*cos(theta_)*sin(phi_)) << ", " << (-r*sin(theta_)*sin(phi_)) << ", " << (-r*cos(phi_)) << "\n";
        
    }
     */

        

    
    //RunAction::Instance()->PrintResponseToBinary();
    
    
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
