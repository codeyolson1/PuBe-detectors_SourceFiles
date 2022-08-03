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
    ~Analysis();

    static Analysis* GetAnalysis();

    void Book(G4String);
    void EndOfRun();

    void OpenFile(const G4String& fname);
    void Save();
    void Close(G4bool reset = true);

    void FillEDep(G4double eDep);
    void FillPrimaryEne(G4double);
    void FillPrimaryPos(G4double, G4double);
    void CheckConvergence();

  private:
    Analysis();
    DISALLOW_COPY_AND_ASSIGN(Analysis);

    G4int eDepHist;
    G4int primEneHist;
    G4int primPosHist;
    G4String convergenceName;
};

#endif