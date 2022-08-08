// Header file for analysis class.
// Created by Codey Olson on May 28, 2021.

/// \file Analysis.hh
/// \brief Selection of analysis classes

#ifndef Analysis_h
#define Analysis_h 1

#include <tools/histo/h1d>
#include <tools/histo/h2d>


#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

class Analysis {
  public:
    Analysis(G4bool);
    ~Analysis();

    static Analysis* GetAnalysis(G4bool);

    void Book(G4String);
    void EndOfRun();

    void OpenFile(const G4String& fname);
    void Save();
    void Close(G4bool reset = true);

    void FillEDep(G4double eDep, G4int);
    void FillPrimaryEne(G4double);
    void FillPrimaryPos(G4double, G4double);
    void CheckConvergence();

  private:
    Analysis();
    DISALLOW_COPY_AND_ASSIGN(Analysis);

    G4int eDepHist;
    G4int eDepHist1;
    G4int eDepHist2;
    G4int eDepHistTot;
    G4int primEneHist;
    G4int primPosHist;
    G4String convergenceName;
    G4bool isHe3;
};

#endif