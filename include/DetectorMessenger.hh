// Header file for DetectorMessenger().
// Created by Codey Olson on August 17, 2021.

/// \file DetectorMessenger.hh
/// \file HEader file for DetectorMessenger class.

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
    virtual ~DetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    DetectorConstruction* fDetConst;
    G4UIdirectory* fDetDir;
    G4UIcmdWithADoubleAndUnit* fDetOffset;
};


#endif