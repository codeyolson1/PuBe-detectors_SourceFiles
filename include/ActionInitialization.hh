// Class definition for ActionInitialization()
// Created by Codey Olson on May 28, 2021

/// \file ActionInitialization.hh
/// \brief Definition of ActionInitialization Class


#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

// Initialize user actions (Event, Run, Step, etc.)

class ActionInitialization : public G4VUserActionInitialization
{
    public:
        ActionInitialization();
        virtual ~ActionInitialization();

        virtual void BuildForMaster() const;
        virtual void Build() const;

    private:
};

#endif