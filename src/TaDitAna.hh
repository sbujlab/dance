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
class TaDefinition;
class TaDitAna: public TObject{
public:
  TaDitAna(TaConfig *aConfig);
  ~TaDitAna(){};

  Bool_t LoadModulationData(TaInput *aInput);
  void Process();
  void WriteToTree(TaOutput* aOutput);
  void PrintSummary(TaOutput* aOutput);

private:
  void RegisterDataElements(vector<TaDefinition*> device_array);
  void RegisterDataElements(vector<TString> device_array);
  void ConnectDataElements();
  vector<TaDataElement*> BuildDataElementArray( vector<TaDefinition*> device_array);
  vector<TaDataElement*> BuildDataElementArray( vector<TString> device_array);
  void RegisterBranchAddress(TTree* );
  TCut bmod_cut;
  TString tree_name;
  
  TaSuperCycle templateCycle;
  vector<TaSuperCycle> fSuperCycleArray;

  vector<TaDataElement*> fCoilArray;
  vector<TaDataElement*> fDependentVarArray;
  
  map<TString, TaDataElement*> fDataElementMap;
  
  ClassDef(TaDitAna,0);
};

#endif
