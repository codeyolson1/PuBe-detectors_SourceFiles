// Source code for class Run().
// Created by Codey Olson on June 2, 2021.

/// \file Run.cc
/// \brief Source code for Run class.

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Threading.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4THitsMap.hh"
#include "G4HCofThisEvent.hh"
#include "G4PrimaryVertex.hh"
#include "G4ThreeVector.hh"
#include "G4WorkerThread.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4StatAnalysis.hh"

#include <atomic>

Run::Run()
{
  EDepPerEvent = {};
}

//
//

Run::~Run()
{}

//
//

void Run::Merge(const G4Run* aRun)
{
  G4Run::Merge(aRun);

  const Run* localRun = static_cast<const Run*>(aRun);
  for (size_t i = 0; i != localRun->EDepPerEvent.size(); i++) {
    EDepPerEvent.push_back(localRun->EDepPerEvent[i]);
  }
}

//
//

void Run::RecordEvent(const G4Event* anEvent)
{
  G4int eventNum = anEvent->GetEventID();
  if (eventNum % 100000 == 0) {
    std::cout << "Event " << eventNum << " started." << std::endl;
  }

  Analysis* myAnalysis = Analysis::GetAnalysis();
  G4SDManager* sdMan = G4SDManager::GetSDMpointer();
  // Get Primary Energy information:
  G4PrimaryVertex* pVertex = anEvent->GetPrimaryVertex();
  G4ThreeVector primPos = pVertex->GetPosition();
  G4PrimaryParticle* primary = pVertex->GetPrimary();
  if (primary->GetG4code()->GetParticleName() == "neutron") {
    G4double primEnergy = primary->GetKineticEnergy();
    myAnalysis->FillPrimaryEne(primEnergy/MeV);
    myAnalysis->FillPrimaryPos(primPos.getX()/cm, primPos.getY()/cm);
  }
  //G4cout << "Primary Energy is: " << energy/MeV << G4endl;
  G4HCofThisEvent* hce = anEvent->GetHCofThisEvent();
  G4int collID = sdMan->GetCollectionID("Helium-3/EnergyDep");
  if (!hce) return;
  G4THitsMap<G4double>* eventMap = 0;
  eventMap = static_cast<G4THitsMap<G4double>*>(hce->GetHC(collID));
  if (eventMap && eventMap->entries() >= 1) {
    G4double val = 0.;
    for (auto itr = eventMap->begin(); itr != eventMap->end(); itr++) {
      val += *itr->second;
    }
    if (val > 0./MeV) {
      myAnalysis->FillEDep(val/MeV);
    }
  }

  G4Run::RecordEvent(anEvent);
}