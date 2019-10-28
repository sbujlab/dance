#ifndef __TaSuperCycle_hh__
#define __TaSuperCycle_hh__

#include "Rtypes.h"
#include "TObject.h"
#include "TString.h"

#include <vector>

#include "TaConfig.hh"
#include "TaAccumulator.hh"
#include "TaDataElement.hh"
#include "TaPrinter.hh"

using namespace std;

class TaSuperCycle: public TObject{
public:
  TaSuperCycle(){};
  virtual ~TaSuperCycle(){};

  void RegisterDependentVarArray(vector<TaDataElement*>);
  void LoadDetectorList( vector<TString> );
  void RegisterCoilArray(vector<TaDataElement*>);

  void ConstructTreeBranches(TaOutput *aOutput);

  void InitAccumulators();
  void UpdateSamples(Int_t );
  void CalcSensitivities();
  void PrintSensitivities();
  void WriteToPrinter(TaPrinter *fPrinter);
  inline void SetCycleID(Int_t id){ cycID = id;};
  inline Int_t GetCycleID(){ return cycID;};

  Int_t GetNumberOfValues(){ return nDependentVar*nCoil;};
  inline Int_t GetIndex(pair<TString,TString> in_pr){ return fSensitivityMap[in_pr];};
  // void WriteSensTree(TOutput *aOutput);
  inline Double_t GetSensitivity(Int_t index){ return fSensitivity[index];};
  inline Double_t GetErrorBar(Int_t index){ return fSensitivity_err[index];};
  inline Double_t GetNSamples(Int_t index){ return fSamples[index];};
  inline Int_t GetSize(){return fSamples.size();};
private:

  Int_t nDependentVar;
  Int_t nCoil;
  vector<TString> fDetectorList;
  
  vector<TaDataElement*> fDependentVarArray;
  vector<TaDataElement*> fCoilArray;
  vector<Bool_t> isDetectorFlag;

  vector< vector<TaAccumulator> > fCovarianceArray; 
  vector< vector<TaAccumulator> > fDepVarianceArray;
  vector<TaAccumulator> fCoilVarianceArray;

  vector<Double_t> fSamples;
  vector<Double_t> fSensitivity;
  vector<Double_t> fSensitivity_err;

  map< pair<TString, TString>, Int_t > fSensitivityMap;

  Int_t cycID;
  ClassDef(TaSuperCycle,0);
};

#endif
