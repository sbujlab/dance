#ifndef __TaEigenRegression_hh__
#define __TaEigenRegression_hh__

#include "TTree.h"
#include "TMatrixD.h"
#include "TCut.h"

#include "TaAccumulator.hh"
#include "TaConfig.hh"
#include "TaInput.hh"
#include "VAnalysisModule.hh"
#include <vector>

using namespace std;

class TaEigenRegression: public VAnalysisModule{
public:
  TaEigenRegression(){};
  TaEigenRegression(Int_t ana_index, TaConfig *aConfig);
  virtual ~TaEigenRegression(){};
  void Process(TaOutput* fOutput);
  vector<vector<Double_t> > Solve(TMatrixD, TMatrixD);
  void CorrectTree();
  void WriteSummary();
  virtual  TMatrixD GetDetMonCovMatrix(Int_t imini);
  TMatrixD GetMonMonCovMatrix(Int_t imini);
  Double_t GetCovariance(TaChannel*, TaChannel*,Int_t);
  
  vector<Double_t> GetColumnVector(TMatrixD, Int_t icol);
  vector<Double_t> GetRowVector(TMatrixD, Int_t irow);
  
private:
  vector<TaChannel*> fEigenVar; // uncorrelated BPMs
  ClassDef(TaEigenRegression,0);
};

#endif
