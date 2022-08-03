// Header file for RunActionMessenger().
// Created by Codey Olson on August 17, 2021.

/// \file RunActionMessenger.hh
/// \file HEader file for RunActionMessenger class.

#ifndef RunActionMessenger_h
#define RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;

class RunActionMessenger: public G4UImessenger
{
  public:
    RunActionMessenger(RunAction*);
    virtual ~RunActionMessenger();
    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    RunAction* fRunAction;
    G4UIdirectory* fRADir;
    G4UIcmdWithAString* fFileName;
};
#endif