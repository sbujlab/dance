#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TaConfig.hh"
#include "TaAccumulator.hh"
#include "TaPrinter.hh"

#ifndef __TaOutput_hh__
#define __TaOutput_hh__

using namespace std;
typedef struct{Double_t ppm, ppb, um,nm;} UNIT;
typedef struct{Double_t hw_sum,block0,block1,block2,block3;} ROOTDATA;

class TaConfig;
class TaOutput{
public:
  TaOutput();
  TaOutput(TaConfig* );
  virtual ~TaOutput();
  
  void ConstructTreeBranch(TString treeName, TString branchName,Double_t &value);
  void ConstructTreeBranch(TString treeName, TString branchName,TString desc, void* value);
  void ConstructTreeBranch(TString treeName, TString branchName,TString desc, ROOTDATA &value);
  void ConstructStatTreeBranch(TString treeName, TString branchName,STAT &value);
  void FillTree(TString treeName);

  void Write();
  void Close();
  inline Int_t GetRunNumber(){return run_number;};
  inline TaPrinter* GetPrinter(){return fPrinter;};
private:
  Int_t nBranches;
  Int_t run_number;
  vector<TBranch*> fBranchArray;
  map< pair<TString,TString>, Int_t> fBranchIndex;
  map<TString, TTree*> fTreeArrayByName;
  vector<TTree*> fTreeArray;
  TFile *outputFile;
  TaPrinter* fPrinter;
  UNIT parity_scale;
  ClassDef(TaOutput,0);
};

#endif
