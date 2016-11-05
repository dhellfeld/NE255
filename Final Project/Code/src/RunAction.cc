#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>

using namespace std;


//==================================================================================================

RunAction* RunAction::fgInstance = 0;

//==================================================================================================

RunAction* RunAction::Instance() {
    
    return fgInstance;
    
}


// ====================================================================================================

RunAction::RunAction(): G4UserRunAction(){
    
    fgInstance = this;
    
  // set printing event number per each event
  //G4RunManager::GetRunManager()->SetPrintProgress(30720000);

/*
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //
  
  // Creating histograms
  analysisManager->CreateH2("1","Response", 3072, 0.5, 3072.5,192, 0.5, 192.5);
*/

}

// ====================================================================================================

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();
    
    fgInstance = 0;
}

// ====================================================================================================

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  

/*
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // Open an output file
    //
    G4String fileName = "output";
    analysisManager->OpenFile(fileName);
*/


}

// ====================================================================================================

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{


/*
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // save histograms
    //
    analysisManager->Write();
    analysisManager->CloseFile();
*/
    
    
    
    
/*
    for(int i = 0; i < 3072; i++){
       for(int j = 0; j < 192; j++){
        
           G4cout << blah[i][j] << "  ";
       }
        G4cout << "\t";
    }
 */
    
    
    
}

// ====================================================================================================
