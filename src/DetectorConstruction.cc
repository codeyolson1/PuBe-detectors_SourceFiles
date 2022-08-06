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
#include "G4UnitsTable.hh"
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
  G4Material* concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
  fmats["concrete"] = concrete;
  G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  fmats["steel"] = steel;
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  fmats["aluminum"] = aluminum;
  G4Material* galactic = nist->FindOrBuildMaterial("G4_Galactic");
  fmats["galactic"] = galactic;
  G4Material* Tantalum = nist->FindOrBuildMaterial("G4_Ta");
  fmats["Tantalum"] = Tantalum;
  G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
  fmats["lead"] = lead;

  // Material characteristics from ShieldWerx.
  G4Element* boron = new G4Element("Boron", "B", 2);
  G4Isotope* b10 = new G4Isotope("Boron10", 5, 10, 10.012936862*g/mole);
  G4Isotope* b11 = new G4Isotope("Boron11", 5, 11, 11.0093051662*g/mole);
  boron->AddIsotope(b10, 19.6*perCent);
  boron->AddIsotope(b11, 80.4*perCent);
  G4Material* BPoly5 = new G4Material("5% Borated Polyethylene", 0.96*g/cm3, 2);
  BPoly5->AddElement(boron, 5.0*perCent);
  BPoly5->AddMaterial(poly, 95.0*perCent);
  fmats["BPoly5"] = BPoly5;
  // BF3 Gas// Material info from :
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

  G4Material* Berillyum = nist->FindOrBuildMaterial("G4_Be");
  G4Element* Plutonium = new G4Element("Plutonium", "Pu", 4);
  G4Isotope* Pu239 = new G4Isotope("Plutonium239", 94, 239, 239.0521617*g/mole);
  G4Isotope* Pu240 = new G4Isotope("Plutonium240", 94, 240, 240.0538118*g/mole);
  G4Isotope* Pu241 = new G4Isotope("Plutonium241", 94, 241, 241.0568497*g/mole);
  G4Isotope* Pu242 = new G4Isotope("Plutonium242", 94, 242, 242.0587410*g/mole);
  // Decayed percentages:
  Plutonium->AddIsotope(Pu239, 91.66472859*perCent);
  Plutonium->AddIsotope(Pu240, 8.25386571*perCent);
  Plutonium->AddIsotope(Pu241, 0.0498354615*perCent);
  Plutonium->AddIsotope(Pu242, 0.0315702471*perCent);
  G4Material* PuBe = new G4Material("Pu-Be Mixture", 3.764135004*g/cm3, 2); 
  PuBe->AddElement(Plutonium, 66.39667705*perCent);
  PuBe->AddMaterial(Berillyum, 33.60332295*perCent);
  fmats["PuBe"] = PuBe;
  
  // From MCNP Material Compendium: Wood (Southern Pine)
  /*
  G4Material* wood = new G4Material("Wood", 0.64*g/cm3, 8);
  G4Element* Ca = nist->FindElement(20);
  G4Element* K = nist->FindElement(19);
  G4Element* S = nist->FindElement(16);
  G4Element* Mg = nist->FindElement(12);
  G4Element* O = nist->FindElement(8);
  G4Element* N = nist->FindElement(7);
  G4Element* C = nist->FindElement(6);
  G4Element* H = nist->FindElement(1);
  wood->AddElementByMassFraction(Ca, 0.001988);
  wood->AddElementByMassFraction(K, 0.001988);
  wood->AddElementByMassFraction(S, 0.004970);
  wood->AddElementByMassFraction(Mg, 0.001988);
  wood->AddElementByMassFraction(O, 0.427435);
  wood->AddElementByMassFraction(N, 0.004970);
  wood->AddElementByMassFraction(C, 0.497018);
  wood->AddElementByMassFraction(H, 0.059642);
  //fmats["wood"] = wood;
  */

}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;

  //
  // World:
  // Params:
  G4double worldX, worldY, worldZ;
  worldX = 381.*cm;
  worldY = 482.6*cm; 
  worldZ = 317.5*cm;
  // Construction:
  G4Box* solidWorld = new G4Box("World", 0.5*worldX, 0.5*worldY,0.5*worldZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, fmats["air"], "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

  //
  // Floor
  // Params:
  G4double floorX = 378.46*cm;
  G4double floorY = 477.52*cm;
  G4double floorZ = 30.48*cm;
  G4ThreeVector floorCenter = G4ThreeVector(0, 0, -worldZ*0.5 + floorZ*0.5 + 5.08*cm);
  // Construction
  G4Box* floorSolid = new G4Box("Floor", 0.5*floorX, 0.5*floorY, 0.5*floorZ);
  G4LogicalVolume* floorLogic = new G4LogicalVolume(floorSolid, fmats["concrete"], "Floor");
  new G4PVPlacement(0, floorCenter, floorLogic, "Floor", logicWorld, false, 0, checkOverlaps);
  G4VisAttributes* concAttr =  new G4VisAttributes(G4Colour::Grey());
  concAttr->SetForceSolid(true);
  floorLogic->SetVisAttributes(concAttr);

  //
  // South Wall
  // Params:
  G4double sWallX = 30.48*cm;
  G4double sWallY = 477.52*cm;
  G4double sWallZ = 276.86*cm;
  G4ThreeVector sWallCenter = G4ThreeVector(worldX*0.5 - sWallX*0.5 - 1.27*cm, worldY*0.5 - sWallY*0.5 - 2.54*cm, worldZ*0.5 - sWallZ*0.5 - 5.08*cm);
  // Construction:
  G4Box* sWallSolid = new G4Box("SouthWall", 0.5*sWallX, 0.5*sWallY, 0.5*sWallZ);
  G4LogicalVolume* sWallLogic = new G4LogicalVolume(sWallSolid, fmats["concrete"], "SouthWall");
  new G4PVPlacement(0, sWallCenter, sWallLogic, "SouthWall", logicWorld, false, 0, checkOverlaps);
  sWallLogic->SetVisAttributes(concAttr);

  //
  // Pillar
  // Params:
  G4double pillarX = 134.62*cm;
  G4double pillarY = 106.68*cm;
  G4double pillarZ = 246.38*cm;
  G4ThreeVector pillarCenter = G4ThreeVector(worldX*0.5 - pillarX*0.5 - sWallX - 1.27*cm, worldY*0.5 - pillarY*0.5 - 2.54*cm, -worldZ*0.5 + pillarZ*0.5 + floorZ + 5.08*cm);
  // Construction:
  G4Box* pillarSolid = new G4Box("Pillar", 0.5*pillarX, 0.5*pillarY, 0.5*pillarZ);
  G4LogicalVolume* pillarLogic = new G4LogicalVolume(pillarSolid, fmats["concrete"], "Pillar");
  new G4PVPlacement(0, pillarCenter, pillarLogic, "Pillar", logicWorld, false, 0, checkOverlaps);
  pillarLogic->SetVisAttributes(concAttr);

  //
  // Lower Pillar
  // Params
  G4double lowPillarX, lowPillarY, lowPillarZ;
  lowPillarX = pillarX;
  lowPillarY = 96.52*cm;
  lowPillarZ = 132.08*cm;
  G4ThreeVector lowPillarCenter = G4ThreeVector(worldX*0.5 - lowPillarX*0.5 - sWallX - 1.27*cm, worldY*0.5 - lowPillarY*0.5 - pillarY - 2.54*cm, -worldZ*0.5 + lowPillarZ*0.5 + floorZ + 5.08*cm);
  // Construction
  G4Box* lowPillarSolid = new G4Box("LowPillar", 0.5*lowPillarX, 0.5*lowPillarY, 0.5*lowPillarZ);
  G4LogicalVolume* lowPillarLogic = new G4LogicalVolume(lowPillarSolid, fmats["concrete"], "LowPillar");
  new G4PVPlacement(0, lowPillarCenter, lowPillarLogic, "LowPillar", logicWorld, false, 0, checkOverlaps);
  lowPillarLogic->SetVisAttributes(concAttr);
  
  //
  // Mesanine
  // Params:
  G4double mesX = 213.36*cm;
  G4double mesY = 162.56*cm;
  G4double mesZ = 198.12*cm;
  G4ThreeVector mesCenter = G4ThreeVector(-worldX*0.5 + mesX*0.5 + 1.27*cm, worldY*0.5 - mesY*0.5 - 2.54*cm, -worldZ*0.5 + mesZ*0.5 + floorZ + 5.08*cm);
  // Construction:
  G4Box* mesanineSolid = new G4Box("Mesanine", 0.5*mesX, 0.5*mesY, 0.5*mesZ);
  G4LogicalVolume* mesanineLogic = new G4LogicalVolume(mesanineSolid, fmats["concrete"], "Mesanine");
  new G4PVPlacement(0, mesCenter, mesanineLogic, "Mesanine", logicWorld, false, 0, checkOverlaps);
  mesanineLogic->SetVisAttributes(concAttr);

  //
  // Table
  // Params
  G4double tableX, tableY, tableZ;
  tableX = 121.92*cm;
  tableY = 182.88*cm;
  tableZ = 3.175*cm;
  G4ThreeVector tableCenter = G4ThreeVector(worldX*0.5 - tableX*0.5 - sWallX - 71.12*cm - 1.27*cm, -worldY*0.5 + tableY*0.5 + 30.48*cm + 2.54*cm, -worldZ*0.5 - tableZ*0.5 + floorZ + 76.2*cm + 5.08*cm);
  // Construction
  G4Box* tableTopSolid = new G4Box("TableTop", 0.5*tableX, 0.5*tableY, 0.5*tableZ);
  G4LogicalVolume* tableTopLogic = new G4LogicalVolume(tableTopSolid, fmats["steel"], "TableTop");
  new G4PVPlacement(0, tableCenter, tableTopLogic, "TableTop", logicWorld, false, 0, checkOverlaps);

  
  // Source and Shielding Bucket
  G4double shieldTopD = 38.1*cm;
  G4double shieldBottomD = 33.02*cm;
  G4double shieldH = 45.72*cm;
  G4ThreeVector shieldCenter = G4ThreeVector(tableCenter.x(), tableCenter.y() - tableY*0.5 + shieldBottomD*0.5, tableCenter.z() + tableZ*0.5 + shieldH*0.5);
  // Construction:
  // BeamPort
  G4Tubs* beamDummy = new G4Tubs("BeamDummy", 0., 2.407*cm, 19.05*0.5*cm, 0., 360.*deg);
  G4Cons* shieldDummy = new G4Cons("ShieldDummy", 0, shieldBottomD*0.5, 0, shieldTopD*0.5, shieldH*0.5, 0, 360.*deg);
  G4Tubs* sourceDummy = new G4Tubs("SourceDummy", 0, 0.5*2.53492*cm, 0.5*5.461*cm, 0, 360.*deg);
  G4VSolid* shieldDummySource = new G4SubtractionSolid("ShieldDummySource", shieldDummy, sourceDummy, 0, G4ThreeVector(0, 0, -6.7675*cm));
  G4RotationMatrix* rotateX = new G4RotationMatrix();
  rotateX->rotateX(90.*deg);  
  G4VSolid* beamIntersection = new G4IntersectionSolid("BeamIntersection", shieldDummySource, beamDummy, rotateX, G4ThreeVector(0, 9.525*cm, -6.7675*cm));
  G4VSolid* shieldSolid = new G4SubtractionSolid("ShieldSolid", shieldDummySource, beamDummy, 0, G4ThreeVector(0, 9.525*cm, -6.7675*cm));
  G4LogicalVolume* shieldLogic = new G4LogicalVolume(shieldSolid, fmats["BPoly5"], "Shield");
  new G4PVPlacement(0, shieldCenter, shieldLogic, "Shield", logicWorld, false, 0, checkOverlaps);
  G4Tubs* beamDummyAir = new G4Tubs("BeamDummyAir", 0., 1.907*cm, 19.05*cm, 0., 360.*deg);
  G4VSolid* beamSolid = new G4SubtractionSolid("BeamSolid", beamIntersection, beamDummyAir, rotateX, G4ThreeVector(0, 0, -6.7675*cm));
  G4LogicalVolume* beamLogic = new G4LogicalVolume(beamSolid, fmats["lead"], "BeamSolid");  
  G4RotationMatrix* rotatebeam = new G4RotationMatrix();
  rotatebeam->rotateX(-180.*deg);  
  rotatebeam->rotateY(180.*deg);
  //new G4PVPlacement(rotatebeam, G4ThreeVector(shieldCenter.x(), shieldCenter.y(), shieldCenter.z()), beamLogic, "BeamSolid", logicWorld, false, 0, checkOverlaps);
  

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
    //new G4PVPlacement(0, G4ThreeVector(0, 0, 0), he3GasLogic, "He3 Gas", logicWorld, false, 0, checkOverlaps); 
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
    //new G4PVPlacement(0, G4ThreeVector(), he3ModeratorLogic, "He3 Moderator", logicWorld, false, 0, checkOverlaps);
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