// Source file for RunAction() class.
// Created by Codey Olson on May 28, 2021.
/// \file RunAction.cc
/// \brief Source code for class RunAction.
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "Analysis.hh"
#include <G4WorkerThread.hh>
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSSphereSurfaceFlux.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4StatAnalysis.hh"
#include "G4RunManager.hh"
#include <iostream>
#include <fstream>
RunAction::RunAction() : G4UserRunAction()
{

  fMessenger = new RunActionMessenger(this);
  outFileName = "He3Full.root";
}
//
//
RunAction::~RunAction()
{}
//
//
G4Run* RunAction::GenerateRun()
{
  return new Run();
}
//
//
void RunAction::BeginOfRunAction(const G4Run*)
{
  Analysis* myAnalysis = Analysis::GetAnalysis();
  myAnalysis->Book(outFileName);
  myAnalysis->OpenFile(outFileName);
} 
//
//
void RunAction::EndOfRunAction(const G4Run* aRun)
{
  Analysis* myAnalysis = Analysis::GetAnalysis();
  if (IsMaster()) {
    G4cout << "End of Global Run" << G4endl;
    myAnalysis->Save();
    myAnalysis->Close();
    myAnalysis->CheckConvergence();
  } else {
    myAnalysis->Save();
  }
}

//
//

void RunAction::SetFileName(G4String name)
{
  outFileName = name;
}
