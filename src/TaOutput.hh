#ifndef __TaOutput_hh__
#define __TaOutput_hh__

#include "TTree.h"

#include "TaConfig.hh"
#include "TaAccumulator.hh"
#include "TaPrinter.hh"

using namespace std;
typedef struct{Double_t ppm, ppb, um,nm;} UNIT;
class TaConfig;
class TaOutput{
public:
  TaOutput();
  TaOutput(TaConfig* );
  virtual ~TaOutput();
  
  void ConstructTreeBranch(TString treeName, TString branchName,Double_t &value);
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
