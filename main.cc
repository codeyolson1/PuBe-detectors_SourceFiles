
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
#include "PhysicsList.hh"
#include "G4ThermalNeutrons.hh"

int main(int argc, char** argv)
{
  G4UIExecutive* ui = 0;
  if ( argc == 1 ){
    ui = new G4UIExecutive(argc, argv);
  }

  G4Random::setTheEngine(new CLHEP::MixMaxRng);
  G4Random::setTheSeed(time(NULL));
  G4cout << "Seed " << G4Random::getTheSeed() << G4endl;
  G4MTRunManager* runManager = new G4MTRunManager;

 

  G4VModularPhysicsList* physicsList = new QGSP_BIC_AllHP();
  physicsList->RegisterPhysics( new G4ThermalNeutrons());

  runManager->SetUserInitialization(new DetectorConstruction());
  //G4VModularPhysicsList* physicsList = new QGSP_BIC_HP();
  physicsList->SetDefaultCutValue(700*CLHEP::um);
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
  runManager->SetVerboseLevel(0);
  G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes( false );
  G4ParticleHPManager::GetInstance()->SetDoNotAdjustFinalState( false );
  G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation( false );
  G4ParticleHPManager::GetInstance()->SetNeglectDoppler( false );
  G4ParticleHPManager::GetInstance()->SetProduceFissionFragments( false );
  //G4ParticleHPManager::GetInstance()->SetUseWendtFissionModel( false );
  G4ParticleHPManager::GetInstance()->SetUseNRESP71Model( false );
  runManager->SetUserInitialization(new ActionInitialization());

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

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