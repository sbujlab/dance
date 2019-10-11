#ifndef __TaAnalysis_hh__
#define __TaAnalysis_hh__
#include "Rtypes.h"
#include "TString.h"
#include <vector>
#include "TaConfig.hh"
#include "TaInput.hh"
#include "TaDithering.hh"
#include "TaLagrangian.hh"

class TaAnalysis {
public:
  TaAnalysis(TaInput *aInput);
  virtual ~TaAnalysis();
  Bool_t LoadConfig(TaConfig *aConfig);
  Bool_t isArrayMatched(vector<TString>,vector<TString>);
  Bool_t Process();
  void End();
private:
  TaDithering *fDithering;
  TaLagrangian *fLagrangian;
  TaConfig *fConfig;
  TaInput *fInput;
  TFile *output_rootfile;
  TString analysisType;
  ClassDef(TaAnalysis,0)
};
#endif
