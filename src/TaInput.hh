#ifndef __TaInput_hh__
#define __TaInput_hh__

#include "TaConfig.hh"

#include <vector>
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;
class TaConfig;
class TaInput{
public:
  TaInput(TaConfig *aconfig);
  virtual ~TaInput();

  Bool_t LoadROOTFile();
  inline TTree* GetEvtTree() const { return evt_tree;};
  inline TTree* GetMulTree() const { return mul_tree;};
  inline TTree* GetSlowTree() const { return slow_tree;};
  inline Int_t GetRunNumber() const { return run_number;};
  inline Int_t GetSegNumber() const { return seg_number;};
  inline void SetRunNumber(Int_t i){run_number=i;};
  inline void SetExtFileName(TString str){
    ext_filename=str;
    isExternalConstraint=kTRUE;
  };
  inline TString GetExtFileName(){return ext_filename;};
  inline Bool_t UseExternalConstraint() const {return isExternalConstraint;};
private:
  TFile* input_file;
  TaConfig *fConfig;
  TTree* evt_tree;
  TTree* mul_tree;
  TTree* slow_tree;
  Int_t run_number;
  Int_t seg_number;
  TString ext_filename;
  Bool_t isExternalConstraint;
  ClassDef(TaInput,0);
};

#endif
