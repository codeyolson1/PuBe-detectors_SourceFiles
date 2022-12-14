// Class Definition of DetectorConstruction().
// Created by Codey Olson on May 8, 2021.

/// \file DetectorConstruction.hh
/// \file Definition of DetectorConstruction class.

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "globals.hh"
#include "G4Region.hh"
#include "DetectorMessenger.hh"
#include "G4PVPlacement.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

// Define detector/world geomteries and materials.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(G4bool);
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
    DetectorMessenger* fMessenger;

  private:
    std::map<std::string, G4Material*> fmats;
    G4bool isHe3;
    G4Region* scatteringRegion;
    G4PVPlacement* fFluxScorerPlacement;
    G4Box* fFluxScorerSolid;
    G4ThreeVector tableCenter;
    G4double fdetOffset;

  public:
    void ConstructMaterials();
    inline G4Region* GetScatteringRegion() const { return scatteringRegion; }
    void SetDetOffset(G4double);
    

};

#endif