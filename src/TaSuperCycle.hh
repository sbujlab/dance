#ifndef __TaSuperCycle_hh__
#define __TaSuperCycle_hh__

#include "Rtypes.h"
#include "TObject.h"
#include "TString.h"
#include "TMatrixD.h"

#include <vector>

#include "TaOutput.hh"
#include "TaConfig.hh"
#include "TaAccumulator.hh"
#include "TaDataElement.hh"
#include "TaPrinter.hh"

using namespace std;

class TaSuperCycle: public TObject{
public:
  TaSuperCycle();
  virtual ~TaSuperCycle(){};

  void EnableDeviceErrorCut(){ kDeviceErrorCut=kTRUE;};
  void RegisterDependentVarArray(vector<TaDataElement*>);
  void LoadDetectorList( vector<TString> );
  void RegisterCoilArray(vector<TaDataElement*>);
  void ConfigSlopesCalculation(TaConfig*);
  void ConstructTreeBranches(TaOutput *aOutput);

  void InitAccumulators();
  void UpdateSamples(Int_t );
  void CalcSensitivities();

  void PrintSensitivities();
  void WriteToPrinter(TaPrinter *fPrinter);

  void CalcSlopes(); 
  Bool_t MakeMatrixFromList(vector<TString>, vector<TString>, TMatrixD&);
  Bool_t MakeMatrixByName(TString, vector<TString>, TMatrixD&);
  Bool_t GetMatrixSolution(TMatrixD lhs, TMatrixD rhs, TMatrixD &sol);
  
  inline void SetCycleID(Int_t id){ cycID = id;};
  inline Int_t GetCycleID(){ return cycID;};

  inline Int_t GetNumberOfValues(){ return nDependentVar*nCoil;};
  inline Int_t GetIndex(pair<TString,TString> in_pr){ return fSensitivityMap[in_pr];};
  // void WriteSensTree(TOutput *aOutput);
  inline Double_t GetSensitivity(Int_t index){ return fSensitivity[index];};
  inline Double_t GetErrorBar(Int_t index){ return fSensitivity_err[index];};
  inline Double_t GetNSamples(Int_t index){ return fSamples[index];};
  inline Int_t GetSize(){return fSamples.size();};

private:
  Int_t cycID;
  Bool_t kDeviceErrorCut;
  Int_t nDependentVar;
  Int_t nCoil;
  vector<TString> fDetectorList;
  
  vector<TaDataElement*> fDependentVarArray;
  vector<TaDataElement*> fCoilArray;
  vector<Bool_t> isDetectorFlag;

  vector< vector<TaAccumulator> > fCovarianceArray; 
  vector< vector<TaAccumulator> > fDepVarianceArray;
  vector< vector<TaAccumulator> >  fCoilVarianceArray;

  vector<Double_t> fSamples;
  vector<Double_t> fSensitivity;
  vector<Double_t> fSensitivity_err;

  map< pair<TString, TString>, Int_t > fSensitivityMap;
  map< TString, Int_t> fDVIndexMapByName;

  // for slopes calculation
  vector<vector<TString> >  fcoil_list;
  vector<vector<TString> > fmonitor_list;

  vector<TString> slope_tree_name;

  vector< vector<Double_t> >  fSlopeContainer;
  vector< vector<Double_t> > fFlagContainer;
  vector< vector<pair<TString,TString> > > fKeyContainer;
  
public:
  inline TString GetTreeName(Int_t index){return slope_tree_name[index];};
  Int_t GetNumberOfSlopeMode(){return slope_tree_name.size();};
  void ConstructSlopeTreeBranch(TaOutput *aOuput, 
				Int_t ana_index,
				vector<Double_t> &fBranchValues,
				vector<Double_t> &fFlagValues);
  void FillSlopeTree(TaOutput *aOuput,
		     Int_t ana_index, 
		     vector<Double_t> &fBranchValues,
		     vector<Double_t> &fFlagValues);

  ClassDef(TaSuperCycle,0);
};

#endif
