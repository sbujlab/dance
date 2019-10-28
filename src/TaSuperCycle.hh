#ifndef __TaSuperCycle_hh__
#define __TaSuperCycle_hh__

#include "Rtypes.h"
#include "TObject.h"
#include "TString.h"

#include "TaConfig.hh"
#include "TaAccumulator.hh"
#include "TaDataElement.hh"
#include <vector>

using namespace std;

class TaSuperCycle: public TObject{
public:
  TaSuperCycle(){};
  virtual ~TaSuperCycle(){};

  void RegisterDependentVarArray(vector<TaDataElement*>);
  void LoadDetectorList( vector<TString> );
  void RegisterCoilArray(vector<TaDataElement*>);

  void InitAccumulators();
  void UpdateSamples(Int_t );
  void CalcSensitivities();

  inline void SetCycleID(Int_t id){ cycID = id;};
  inline Int_t GetCycleID(){ return cycID;};
  inline Int_t GetNSamples(Int_t i){ return fSamples[i];};

  // void WriteSensTree(TOutput *aOutput);
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


  vector<Int_t> fSamples;
  vector<Double_t> fSensitivity;
  vector<Double_t> fSensitivity_err;
  map< pair<TString, TString>, Int_t > fSensitivityMap;

  Int_t cycID;
  ClassDef(TaSuperCycle,0);
};

#endif
