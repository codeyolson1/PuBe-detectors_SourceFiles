// Source code for RunActionMessenger().
// Created by Codey Olson on August 17, 2021.

/// \file RunActionMessenger.cc
/// \file Source code for RunActionMessenger class.

#include "RunAction.hh"
#include "RunActionMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

RunActionMessenger::RunActionMessenger(RunAction* myRunAction)
: G4UImessenger(), fRunAction(myRunAction)
{
  fRADir = new G4UIdirectory("/RunAction/");
  fRADir->SetGuidance("Control parameters set within RunAction.");

  fFileName = new G4UIcmdWithAString("/RunAction/FileName", this);
  fFileName->SetGuidance("Set the file name for output.");
  fFileName->SetParameterName("choice", false);
  fFileName->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//
//

RunActionMessenger::~RunActionMessenger()
{
  delete fRADir;
  delete fFileName;
}

//
//

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newVal)
{
  if (command == fFileName) {
    fRunAction->SetFileName(newVal);
  }
}