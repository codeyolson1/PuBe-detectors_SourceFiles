
#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UserRunAction.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4ParticleHPManager.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "QGSP_BIC_AllHP.hh"
#include "FTFP_BERT_HP.hh"
#include "globals.hh"
#include "G4ThermalNeutrons.hh"

// Execute with:
//
// ./pu-be macro_file.mac det_option
//
// det_option: optional - use to specify either 'he3' or 'bf3' for he3/bf3 detector options, respectively.

int main(int argc, char** argv)
{
  // Check to see if in Batch or interactive mode based on command line args:
  G4UIExecutive* ui = 0;
  if ( argc == 1 ){
    ui = new G4UIExecutive(argc, argv);
  }
  // Set the random number engine and seed. Report to std::output
  G4Random::setTheEngine(new CLHEP::MixMaxRng);
  G4Random::setTheSeed(time(NULL));
  G4cout << "Seed " << G4Random::getTheSeed() << G4endl;
  G4MTRunManager* runManager = new G4MTRunManager;
  // Generate and connect with physics list. Add Thermal neutrons for low energy scattering, default cut value, and verbosity.
  G4VModularPhysicsList* physicsList = new QGSP_BIC_AllHP();
  physicsList->RegisterPhysics( new G4ThermalNeutrons());
  physicsList->SetDefaultCutValue(700*CLHEP::um);
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
  // Get Detector info from argument to determine analysis info, geometry, etc. and initialize DetectorConstruction.
  G4bool isHe3;
  G4String detName = argv[2];
  if (detName == "he3") {
    G4cout << "He3 argument passed." << G4endl;
    isHe3 = true;
  } else if (detName == "bf3") {
    G4cout <<"BF3 argument passed." << G4endl;
    isHe3 = false;
  } else {
    G4cout << "Invalid Detector Argument Passed. Defaulting to He3." << G4endl;
    isHe3 = true;
  }
  runManager->SetUserInitialization(new DetectorConstruction(isHe3));
  runManager->SetUserInitialization(new ActionInitialization(isHe3));
  runManager->SetVerboseLevel(0);
  // Additional options for HP physics:
  G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes( false );
  G4ParticleHPManager::GetInstance()->SetDoNotAdjustFinalState( false );
  G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation( false );
  G4ParticleHPManager::GetInstance()->SetNeglectDoppler( false );
  G4ParticleHPManager::GetInstance()->SetProduceFissionFragments( false );
  //G4ParticleHPManager::GetInstance()->SetUseWendtFissionModel( false );
  G4ParticleHPManager::GetInstance()->SetUseNRESP71Model( false );
  
  // Visulaization manager for interactive mode or static displays.
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  // Interactive session generation:
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( !ui ) {
    // batch mode - Apply macros directly 
    UImanager->ApplyCommand("/control/macroPath ../macros/");
    UImanager->ApplyCommand("/control/execute init.mac");
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  } else {
    // Interactive Mode:
    //UImanager->ApplyCommand("/control/macroPath /reactorBay_sourceFiles");
    UImanager->ApplyCommand("/control/macroPath ../macros/");
    UImanager->ApplyCommand("/control/execute init.mac");
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;

}