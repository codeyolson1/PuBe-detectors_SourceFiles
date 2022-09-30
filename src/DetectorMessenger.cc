// Source code for DetectorMessenger().
// Created by Codey Olson on August 17, 2021.

/// \file DetectorMessenger.cc
/// \file Source code for DetectorMessenger class.

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* myDet)
: G4UImessenger(), fDetConst(myDet)
{
  fDetDir = new G4UIdirectory("/Detector/");
  fDetDir->SetGuidance("Control of detector characteristics");

  fDetOffset = new G4UIcmdWithADoubleAndUnit("/Detector/SetOffset", this);
  fDetOffset->SetUnitCandidates(G4UIcommand::UnitsList(G4UIcommand::CategoryOf("cm")));
  fDetOffset->SetRange("Offset > 0.");
  fDetOffset->SetGuidance("Set the distance for the detector from edge of table (in cm). Default is 2 ft (60.96 cm).");
  fDetOffset->SetParameterName("Offset", false);
  fDetOffset->AvailableForStates(G4State_PreInit, G4State_Idle);


/*
  G4UIparameter* SizeZParam = new G4UIparameter("SizeZ", 'd', false);
  SizeZParam->SetGuidance("SizeZ");
  SizeZParam->SetParameterRange("SizeZ > 0.");
  fScintThick->SetParameter(SizeZParam);

  G4UIparameter* unitParameter = new G4UIparameter("ThickUnit", 's', false);
  unitParameter->SetGuidance("unit of thickness");
  unitParameter->SetParameterCandidates(G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm")));
  fScintThick->SetParameter(unitParameter);
*/
/*
  fScintThick = new G4UIcmdWithAString("/Detector/SetThickness", this);
  fScintThick->SetGuidance("Set the thickness for the scintillator.");
  fScintThick->SetParameterName("choice", false);
  fScintThick->AvailableForStates(G4State_PreInit, G4State_Idle);
  */
}

//
//

DetectorMessenger::~DetectorMessenger()
{
  delete fDetDir;
  delete fDetOffset;
}

//
//

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newVal)
{

  if (command == fDetOffset) {
    fDetConst->SetDetOffset(fDetOffset->GetNewDoubleValue(newVal));
  }
}
