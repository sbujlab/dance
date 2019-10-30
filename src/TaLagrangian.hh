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
  TMatrixD Solve(TMatrixD CovDM, TMatrixD CovMM);
  TString GetBaseName(TString);
private:
  Int_t nCoil;
  TMatrixD detConstraints;
  TMatrixD monConstraints;

  ClassDef(TaLagrangian,0);
};
#endif
