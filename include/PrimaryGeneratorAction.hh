// Class definition of PrimaryGeneratorAction().
// Created by Codey Olson on May 28, 2021;

/// \file PrimaryGeneratorAction.hh
/// \brief Class definition of PrimaryGeneratorAction.

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* );

    const G4GeneralParticleSource* GetParticleGun() const { return fParticleGun; }
  
  private:
    G4GeneralParticleSource* fParticleGun;
};

#endif