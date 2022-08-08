// Source file for ActionInitialization().
// Created by Codey Olson on May 28, 2021.

/// \file ActionInitialization.cc
/// \brief Source code for ActionInitialization class.

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "G4Threading.hh"

ActionInitialization::ActionInitialization(G4bool He3) : G4VUserActionInitialization()
{
  isHe3 = He3;
}

//
//

ActionInitialization::~ActionInitialization()
{}

//
//

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction(isHe3));
}

//
//

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new RunAction(isHe3));
}