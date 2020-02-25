#ifndef __TaCorrection_hh__
#define __TaCorrection_hh__

#include "TaConfig.hh"
#include "TaInput.hh"
#include "VAnalysisModule.hh"
#include "TMatrixD.h"

class TaCorrection: public VAnalysisModule{
public:
  TaCorrection(){};
  TaCorrection(Int_t ana_index, TaConfig *aConfig);
  virtual ~TaCorrection(){};

  void Process(TaOutput* fOutput);
  Bool_t LoadSlopeMatrix(Int_t ana_index,TaConfig *aConfig);
  TString GetBaseName(TString);
  Int_t FindRawDVIndexFromList(TString);
  void LoadRawDVList(vector<TString>);
  void LoadRawDVList(TString );
  vector<Double_t> GetSlopeVector(TString);
private:
  vector<TString> fRawDVList;
  TMatrixD slope_matrix;
  ClassDef(TaCorrection,0);
};
#endif
