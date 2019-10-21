#ifndef __TaRegression_hh__
#define __TaRegression_hh__

#include "TTree.h"
#include "TMatrixD.h"
#include "TCut.h"

#include "TaAccumulator.hh"
#include "TaConfig.hh"
#include "TaInput.hh"
#include "VAnalysisModule.hh"
#include <vector>

using namespace std;

class TaRegression: public VAnalysisModule{
public:
  TaRegression(TString type, TString name, TaConfig *aConfig);
  virtual ~TaRegression(){};
  void Process(TaOutput* fOutput);
  TMatrixD Solve(TMatrixD, TMatrixD);
  void CorrectTree();
  void WriteSummary();
  TMatrixD GetDetMonCovMatrix(Int_t imini);
  TMatrixD GetMonMonCovMatrix(Int_t imini);
  Double_t GetCovariance(TaChannel*, TaChannel*,Int_t);

  vector<Double_t> GetColumnVector(TMatrixD, Int_t icol);
  vector<Double_t> GetRowVector(TMatrixD, Int_t irow);
  
private:
  ClassDef(TaRegression,0);
};

#endif
