#ifndef __TaDitAna_hh__
#define __TaDitAna_hh__
#include "TTree.h"
#include "TMatrixD.h"
#include "TCut.h"

#include <vector>

#include "TaSuperCycle.hh"
#include "TaConfig.hh"
#include "TaInput.hh"
#include "TaDataElement.hh"

using namespace std;

class TaDitAna: public TObject{
public:
  TaDitAna(TaConfig *aConfig);
  ~TaDitAna();

  Bool_t LoadModulationData(TaInput *aInput);
  // Bool_t CalcSensitivities();
  // void WriteOutput(TaOutput* aOutput);
  // void PrintSummary();

private:

  void RegisterRawDataElements(vector<TString> device_array);
  void ProcessDefinitions(vector<pair<TString,TString> > fDefinitions);
  vector<TaDataElement*> BuildDataElementArray( vector<TString> device_array);
  void RegisterBranchAddress(TTree* );
  TCut bmod_cut;
  TString tree_name;
  
  TaSuperCycle protoCycle;
  vector<TaSuperCycle> fSuperCycleArray;

  vector<TaDataElement*> fCoilArray;
  vector<TaDataElement*> fDependentVarArray;
  vector<TaDataElement*> fRawDataElementArray;
  
  map<TString, TaDataElement*> fDataElementMap;
  // map<TString, TaDataElement* > fRawDataElementMap;
  // map< TString, TaDataElement* > fDefinedElementMap;

  ClassDef(TaDitAna,0);
};

#endif
