// Class definition of Run().
// Created by Codey Olson on May 28, 2021;

/// \file Run.hh
/// \brief Class definition of Run.

#ifndef Run_h
#define Run_h 1

#include <G4Run.hh>

//Dummy class to show how to use MPI merging of
//user defined G4Run. The dummycounter it's a unique number among threads.
//Merge actually sums the dummyCounter so for the global run that is  Sum(i, 0<i<Nt)
class Run : public G4Run {

    public:
      Run();
      virtual ~Run();
      void Merge(const G4Run*);
      void RecordEvent(const G4Event* anEvent);

      G4int GetNumEntries() const { return EDepPerEvent.size(); }
      G4double GetEDepAtEvent(G4int i) const { return EDepPerEvent[i]; }
    private:
      std::vector<G4double> EDepPerEvent;
};

#endif 
