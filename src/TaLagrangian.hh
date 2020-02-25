#ifndef __TaLagrangian_hh__
#define __TaLagrangian_hh__

#include <vector>

#include "TTree.h"
#include "TMatrixD.h"
#include "TCut.h"

#include "TaAccumulator.hh"
#include "TaConfig.hh"
#include "TaInput.hh"
#include "TaRegression.hh"

using namespace std;

class TaLagrangian: public TaRegression{
public:
  TaLagrangian(){};
  TaLagrangian(Int_t ana_index, TaConfig *aConfig);
  virtual ~TaLagrangian(){};

  void LoadConstraint(Int_t, TaConfig*);
  vector< vector<Double_t> > Solve(TMatrixD CovDM, TMatrixD CovMM);
  TMatrixD GetDetMonCovMatrix(Int_t mini);
  TString GetBaseName(TString);
  void LoadRawDVList(TString);
  void LoadRawDVList(vector<TString>);
private:
  Int_t nCoil;
  TMatrixD detConstraints;
  TMatrixD monConstraints;
  vector<TString> fRawDVlist;
  Int_t FindRawDVIndexFromList(TString raw_name);

  ClassDef(TaLagrangian,0);
};
#endif
