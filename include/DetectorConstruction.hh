// Class Definition of DetectorConstruction().
// Created by Codey Olson on May 8, 2021.

/// \file DetectorConstruction.hh
/// \file Definition of DetectorConstruction class.

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "globals.hh"
#include "G4Region.hh"
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

  private:
    std::map<std::string, G4Material*> fmats;
    G4bool isHe3;
    G4Region* scatteringRegion;

  public:
    void ConstructMaterials();
    inline G4Region* GetScatteringRegion() const { return scatteringRegion; }

    

};

#endif