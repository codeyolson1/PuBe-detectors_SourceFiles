// Definition of world geometry and detectors.
// Created by Codey Olson on May 10, 2021.

/// \file DetectorConstruction.cc
/// \brief Definition of world geometry and detectors.

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSCellFlux.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#define _USE_MATH_DEFINES 
#include <math.h>
#include <iomanip>
#include <iostream>
#include <string>


DetectorConstruction::DetectorConstruction(G4bool he3)
: G4VUserDetectorConstruction()
{
  isHe3 = he3;
  fmats = {};
  ConstructMaterials();
}

//
//

DetectorConstruction::~DetectorConstruction()
{}

//
//

void DetectorConstruction::ConstructMaterials()
{
  // Get instance of nist material manager:
  G4NistManager* nist = G4NistManager::Instance();

  // Create materials and input into dictionary.
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  fmats["air"] = air;
  G4Material* poly = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  fmats["poly"] = poly;

  // BF3 Gas
// Material info from :
// https://gitlab.cern.ch/clemenci/Geant4-srcs/-/blob/92686251452762ac5947193b5f02ba43b77f546b/examples/extended/hadronic/FissionFragment/src/FFDetectorConstruction.cc
    G4double const B10Enrichment = 0.96;
    G4double const B11Enrichment = 0.04;
    G4Isotope* const iB10
        = new G4Isotope("iB10",                         // name
                        5,                              // ZZZ
                        10,                             // AAA
                        10.0129370 * (g / mole));       // molecular weight
    G4Isotope* const iB11
        = new G4Isotope("iB11",                         // name
                        5,                              // ZZZ
                        11,                             // AAA
                        11.0093054 * (g / mole));       // molecular weight
    // Now create the elements and add the isotopes
    G4Element* const B10
        = new G4Element("B10",                          // name
                        "B10",                          // symbol
                        1);                             // number of isotopes
    B10->AddIsotope(iB10,                               // isotope
                     1.0);                              // abundance
    G4Element* const B11
        = new G4Element("B11",                          // name
                        "B11",                          // symbol
                        1);                             // number of isotopes
    B11->AddIsotope(iB11,                               // isotope
                     1.0);                              // abundance
    G4Element* const flouride = nist->FindOrBuildElement("F");
    // Calculate the mass fractions
    const G4double BF3MolecularWeight = B10->GetA() * B10Enrichment
                                        + B11->GetA() * B11Enrichment
                                        + flouride->GetA() * 3;
    const G4double B10MassFraction = (B10->GetA() * B10Enrichment)
                                     / BF3MolecularWeight;
    const G4double B11MassFraction = (B11->GetA() * B11Enrichment)
                                     / BF3MolecularWeight;
    const G4double flourideMassFraction = (flouride->GetA() * 3)
                                          / BF3MolecularWeight;
    // create the material and add the elements
    fmats["enrBF3"] = new G4Material("BF3_96E",                // name
                              2.73E-3 * (g/cm3),          // density
                              3, kStateGas, 293.*kelvin, 1.0*atmosphere);                       // number of components
    fmats["enrBF3"]->AddElement(B10,                           // element
                         B10MassFraction);              // mass fraction
    fmats["enrBF3"]->AddElement(B11,                           // element
                         B11MassFraction);              // mass fraction
    fmats["enrBF3"]->AddElement(flouride,                      // element
                         flourideMassFraction);         // mass fraction

  // He3 Gas       
  G4double atomicMass = 3.01602932197*g/mole;
  G4Isotope* he3 = new G4Isotope("He3", 2, 3, atomicMass);
  G4Element* He3 = new G4Element("Helium3", "He3", 1);
  He3->AddIsotope(he3, 100*perCent);
  //G4double pressure = 4.053*bar;
  G4double pressure = 4.0*atmosphere;
  G4double temperature = 293*kelvin;
  //G4double molar_constant = CLHEP::Avogadro*CLHEP::k_Boltzmann;  //from clhep
  G4double density = 5.39E-4*(g/cm3);
  G4cout << "He3 density: " << density/(g/cm3) << G4endl;
  G4Material* Helium3 = new G4Material("Helium3", density, 1, kStateGas, temperature, pressure);
  Helium3->AddElement(He3, 100*perCent);
  fmats["he3"] = Helium3;

  G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  fmats["steel"] = steel;

  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  fmats["aluminum"] = aluminum;

  G4Material* galactic = nist->FindOrBuildMaterial("G4_Galactic");
  fmats["galactic"] = galactic;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;

  //
  // World:
  // Params:
  G4double worldX, worldY, worldZ;
  worldX = 11.54*cm; 
  worldY = 9.54*cm; 
  worldZ = 11.*cm;
  // Construction:
  G4Box* solidWorld = new G4Box("World", 0.5*worldX, 0.5*worldY,0.5*worldZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, fmats["air"], "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);




  // Source and Shielding Bucket
  G4double shieldTopD = 38.1*cm;
  G4double shieldBottomD = 33.02*cm;
  G4double shieldH = 45.72*cm;
  //G4ThreeVector shieldCenter = G4ThreeVector(-48.26*cm, -22.86*cm, -148.49*cm); // Storage
  G4ThreeVector shieldCenter = G4ThreeVector(-127.635*cm, -22.86*cm, -176.53*cm); // Operation
  // Construction:
  // BeamPort`
  G4Tubs* beamDummy = new G4Tubs("BeamDummy", 0., 2.407*cm, 19.05*0.5*cm, 0., 360.*deg);
  G4Cons* shieldDummy = new G4Cons("ShieldDummy", 0, shieldBottomD*0.5, 0, shieldTopD*0.5, shieldH*0.5, 0, 360.*deg);
  G4Tubs* sourceDummy = new G4Tubs("SourceDummy", 0, 0.5*2.53492*cm, 0.5*5.461*cm, 0, 360.*deg);
  G4VSolid* shieldDummySource = new G4SubtractionSolid("ShieldDummySource", shieldDummy, sourceDummy, 0, G4ThreeVector(0, 0, -6.7675*cm));
  G4RotationMatrix* rotateX = new G4RotationMatrix();
  rotateX->rotateX(90.*deg);  
  G4VSolid* beamIntersection = new G4IntersectionSolid("BeamIntersection", shieldDummySource, beamDummy, rotateX, G4ThreeVector(0, -9.525*cm, -6.7675*cm));
  G4VSolid* shieldSolid = new G4SubtractionSolid("ShieldSolid", shieldDummySource, beamDummy, 0, G4ThreeVector());
  G4Tubs* beamDummyAir = new G4Tubs("BeamDummyAir", 0., 1.907*cm, 19.05*cm, 0., 360.*deg);
  G4VSolid* beamSolid = new G4SubtractionSolid("BeamSolid", beamIntersection, beamDummyAir, rotateX, G4ThreeVector(0, 0, -6.7675*cm));
  G4LogicalVolume* beamLogic = new G4LogicalVolume(beamSolid, fmats["lead"], "BeamSolid");
  new G4PVPlacement(0, G4ThreeVector(-127.635*cm, -22.86*cm, -176.53*cm), beamLogic, "BeamSolid", logicWorld, false, 0, checkOverlaps);
  G4LogicalVolume* shieldLogic = new G4LogicalVolume(shieldDummySource, fmats["BPoly5"], "Shield");
  new G4PVPlacement(0, shieldCenter, shieldLogic, "Shield", logicWorld, false, 0, checkOverlaps);
  // Source:
  // Params:
  G4double stainlessD = 2.53492*cm;
  G4double stainlessH = 5.461*cm;
  G4double tantalumD = 2.37236*cm;
  G4double tantalumH = 4.572*cm;
  G4double PuBeD = 2.06756*cm;
  G4double PuBeH = 3.814*cm;
  G4ThreeVector PuBeCenter = G4ThreeVector(0, 0, 0); // Local coords for shield:
  // Construction:
  G4Tubs* stainlessSolid = new G4Tubs("StainlessShell", 0, 0.5*stainlessD, 0.5*stainlessH, 0, 360.*deg);
  G4LogicalVolume* stainlessLogic = new G4LogicalVolume(stainlessSolid, fmats["steel"], "StainlessShell");
  G4VisAttributes* sourceAttr =  new G4VisAttributes(G4Colour(1.,0.,0.));
  sourceAttr->SetForceSolid(true);
  stainlessLogic->SetVisAttributes(sourceAttr);
  new G4PVPlacement(0, PuBeCenter, stainlessLogic, "StainlessShell", shieldLogic, false, 0, checkOverlaps);
  G4Tubs* tantalumSolid = new G4Tubs("TantalumShell", 0, 0.5*tantalumD, 0.5*tantalumH, 0, 360.*deg);
  G4LogicalVolume* tantalumLogic = new G4LogicalVolume(tantalumSolid, fmats["Tantalum"], "TantalumShell");
  new G4PVPlacement(0, G4ThreeVector(), tantalumLogic, "TantalumShell", stainlessLogic, false, 0, checkOverlaps);
  G4Tubs* PuBeSolid = new G4Tubs("PuBeSource", 0, 0.5*PuBeD, 0.5*PuBeH, 0, 360.*deg);
  G4LogicalVolume* PuBeLogic = new G4LogicalVolume(PuBeSolid, fmats["PuBe"], "PuBeSource");
  new G4PVPlacement(0, G4ThreeVector(), PuBeLogic, "PuBe Source", tantalumLogic, false, 0, checkOverlaps);
  // Detector and Placement:
  if (isHe3) {
    G4double tubeDiam;
    G4double tubeHeight;
    G4double modx, mody, modz;
    // Tube and moderator dimensions:
    tubeDiam = 2.54*cm; // With shell
    tubeHeight = 10.*cm;
    modx = tubeDiam + 4.*cm; mody = tubeDiam + 2.*cm; modz = tubeHeight;

    // Tube Construction
    G4Tubs* ssShellSolid = new G4Tubs("SS Shell", 0, 0.5*(tubeDiam + 0.2*cm), 0.5*(tubeHeight + 0.2*cm), 0, 360.*deg);
    G4LogicalVolume* ssShellLogic = new G4LogicalVolume(ssShellSolid, fmats["steel"], "SS Shell");
    //new G4PVPlacement(0, G4ThreeVector(), ssShellLogic, "SS Shell", logicWorld, false, 0, checkOverlaps); 
    // Visual Stuff for shells
    G4VisAttributes* shellAttr = new G4VisAttributes(G4Colour(192., 192., 192.)); // silver
    shellAttr->SetForceWireframe(true);
    ssShellLogic->SetVisAttributes(shellAttr);
    // helium3 fill gas:
    G4Tubs* he3GasSolid = new G4Tubs("He3 Gas", 0, 0.5*(tubeDiam), 0.5*(tubeHeight), 0, 360.*deg);
    G4LogicalVolume* he3GasLogic = new G4LogicalVolume(he3GasSolid, fmats["he3"], "He3 Gas");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), he3GasLogic, "He3 Gas", logicWorld, false, 0, checkOverlaps); 
    G4cout << "Gas volume: " << he3GasSolid->GetCubicVolume()/cm3 << G4endl;
    // Visual Stuff for gas
    G4VisAttributes* gasAttr = new G4VisAttributes(G4Colour(255., 0., 0.)); // red
    gasAttr->SetForceSolid(true);
    he3GasLogic->SetVisAttributes(gasAttr);
    //Moderator:
    // Dummies for subtraction solid:
    G4Box* moderatorDummy = new G4Box("He3 Moderator Dummy", 0.5*modx, 0.5*mody, 0.5*modz);
    G4Tubs* moderatorVoidDummy = new G4Tubs("He3 Void Dummy", 0, 0.5*(tubeDiam), 0.5*(tubeHeight + 1*cm), 0, 360.*deg);
    // Final solid:
    G4VSolid* he3ModeratorSolid = new G4SubtractionSolid("He3 Moderator", moderatorDummy, moderatorVoidDummy, 0, G4ThreeVector());
    G4LogicalVolume* he3ModeratorLogic = new G4LogicalVolume(he3ModeratorSolid, fmats["poly"], "He3 Moderator");
    new G4PVPlacement(0, G4ThreeVector(), he3ModeratorLogic, "He3 Moderator", logicWorld, false, 0, checkOverlaps);
    G4cout << "Moderator volume: " << he3ModeratorSolid->GetCubicVolume()/cm3 << G4endl;
    // Visual Stuff for moderator
    G4VisAttributes* moderatorAttr = new G4VisAttributes(G4Colour()); // white
    moderatorAttr->SetForceSolid(true);
    he3ModeratorLogic->SetVisAttributes(moderatorAttr);
  } else {
    G4double tubeDiam;
    G4double tubeHeight;
    G4double modx, mody, modz;

    // Tube and moderator dimensions:
    tubeDiam = 4.4*cm;
    tubeHeight = 10.0*cm;
    modx = tubeDiam*2. + 4.5*cm; mody = tubeDiam + 2.*cm; modz = tubeHeight;

    // Construct BF3 Detectors:
    // SS Shells
    G4Tubs* bf3ShellSolid1 = new G4Tubs("BF3 Shell1", 0, 0.5*(tubeDiam + 0.2*cm), 0.5*(tubeHeight + 0.2*cm), 0., 360.*deg);
    G4LogicalVolume* bf3ShellLogic1 = new G4LogicalVolume(bf3ShellSolid1, fmats["steel"], "BF3 Shell1");
    //new G4PVPlacement(0, G4ThreeVector(tubeDiam*0.5 + 0.5*cm, 0, 0), bf3ShellLogic1, "BF3 Shell1", logicWorld, false, 0, checkOverlaps);
    G4Tubs* bf3ShellSolid2 = new G4Tubs("BF3 Shell2", 0, 0.5*(tubeDiam + 0.2*cm), 0.5*(tubeHeight + 0.2*cm), 0., 360.*deg);
    G4LogicalVolume* bf3ShellLogic2 = new G4LogicalVolume(bf3ShellSolid2, fmats["steel"], "BF3 Shell2");
    //new G4PVPlacement(0, G4ThreeVector(-tubeDiam*0.5 - 0.5*cm, 0, 0), bf3ShellLogic2, "BF3 Shell2", logicWorld, false, 0, checkOverlaps);
    // Visual Stuff for shells
    G4VisAttributes* shellAttr = new G4VisAttributes(G4Colour(192., 192., 192.)); // silver
    shellAttr->SetForceWireframe(true);
    bf3ShellLogic1->SetVisAttributes(shellAttr);
    bf3ShellLogic2->SetVisAttributes(shellAttr);
    // BF3 fill gas:
    G4Tubs* bf3GasSolid1 = new G4Tubs("BF3 Gas1", 0, 0.5*(tubeDiam), 0.5*(tubeHeight), 0, 360.*deg);
    G4LogicalVolume* bf3GasLogic1 = new G4LogicalVolume(bf3GasSolid1, fmats["enrBF3"], "BF3 Gas1");
    new G4PVPlacement(0, G4ThreeVector(tubeDiam*0.5 + 0.25*cm, 0, 0), bf3GasLogic1, "BF3 Gas1", logicWorld, false, 0, checkOverlaps);
    G4Tubs* bf3GasSolid2 = new G4Tubs("BF3 Gas2", 0, 0.5*(tubeDiam), 0.5*(tubeHeight), 0, 360.*deg);
    G4LogicalVolume* bf3GasLogic2 = new G4LogicalVolume(bf3GasSolid2, fmats["enrBF3"], "BF3 Gas2");
    new G4PVPlacement(0, G4ThreeVector(-(tubeDiam)*0.5 - 0.25*cm, 0, 0), bf3GasLogic2, "BF3 Gas2", logicWorld, false, 0, checkOverlaps);
    G4cout << "BF3 gas volume: " << bf3GasSolid1->GetCubicVolume()/cm3 + bf3GasSolid2->GetCubicVolume()/cm3 << G4endl;
    // Visual Stuff for gas
    G4VisAttributes* gasAttr = new G4VisAttributes(G4Colour(255., 0., 0.)); // red
    gasAttr->SetForceSolid(true);
    bf3GasLogic1->SetVisAttributes(gasAttr);
    bf3GasLogic2->SetVisAttributes(gasAttr);

    // Moderator:
    G4Box* moderatorDummy1 = new G4Box("BF3 Moderator Dummy", 0.5*modx, 0.5*mody, 0.5*modz);
    G4Tubs* moderatorVoidDummy1 = new G4Tubs("BF3 Moderator Void Dummy", 0, 0.5*(tubeDiam), 0.5*(tubeHeight + 1.*cm), 0, 360.*deg);
    G4VSolid* bf3ModeratorTemp = new G4SubtractionSolid("Mod Temp", moderatorDummy1, moderatorVoidDummy1, 0, G4ThreeVector((tubeDiam)*0.5 + 0.25*cm, 0, 0));
    G4VSolid* bf3ModeratorSolid = new G4SubtractionSolid("BF3 Moderator", bf3ModeratorTemp, moderatorVoidDummy1, 0, G4ThreeVector(-(tubeDiam)*0.5 - 0.25*cm, 0, 0));
    G4LogicalVolume* moderatorBF3Logic = new G4LogicalVolume(bf3ModeratorSolid, fmats["poly"], "ModeratorBF3");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), moderatorBF3Logic, "ModeratorBF3", logicWorld, false, 0, checkOverlaps);
    G4cout << "Moderator volume: " << bf3ModeratorSolid->GetCubicVolume()/cm3 << G4endl;
    // Visual Stuff for moderator
    G4VisAttributes* moderatorAttr = new G4VisAttributes(G4Colour()); // white
    moderatorAttr->SetForceSolid(true);
    moderatorBF3Logic->SetVisAttributes(moderatorAttr);
  }
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  G4SDParticleFilter* nFilter = new G4SDParticleFilter("NeutronFilter");
  nFilter->add("proton");
  nFilter->add("triton");
  nFilter->addIon(2, 3);
  nFilter->add("deuteron");
  nFilter->add("alpha");
  nFilter->add("neutron");

  G4MultiFunctionalDetector* he3Detector = new G4MultiFunctionalDetector("Helium-3");
  G4SDManager::GetSDMpointer()->AddNewDetector(he3Detector);
  G4VPrimitiveScorer* energyDep = new G4PSEnergyDeposit("EnergyDep");
  he3Detector->RegisterPrimitive(energyDep);
  energyDep->SetFilter(nFilter);
  SetSensitiveDetector("He3 Gas", he3Detector);

}