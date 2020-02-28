#ifndef __TaChannel_hh__
#define __TaChannel_hh__

#include "TaOutput.hh"
#include "TaAccumulator.hh"
#include "TaConfig.hh"
#include "TaDefinition.hh"
using namespace std;
class TaOutput;
class TaAccumulator;

class TaChannel {
public:
  TaChannel();
  TaChannel(TString tree, TaDefinition* adef);
  TaChannel(TString tree, TString name);
  virtual ~TaChannel();
  void ConnectChannels(vector<TaChannel*>, vector<Double_t> );
  void DefineSubtraction(TaChannel*,TaChannel*);
  void CalcCombination();
  void FillDataArray();
  Double_t FillOutputValue();

  void AccumulateRunSum();
  void AccumulateMiniSum();

  void UpdateRunStat(){fRunAccumulator.UpdateStat(run_stat);};
  void UpdateMiniStat(){fMiniAccumulator.UpdateStat(mini_stat);};
  void ResetMiniAccumulator(){fMiniAccumulator.Zero();};
  void ResetRunAccumulator(){fRunAccumulator.Zero();};

  void ConstructTreeBranch(TaOutput*);
  void ConstructMiniTreeBranch(TaOutput*,TString);
  void ConstructSlopeBranch(TaOutput*,TString);
  void ConstructSumTreeBranch(TaOutput*,TString);

  TString GetChannelBaseName();
  Double_t GetEntry(Int_t ie);
  Double_t fOutputValue;  
  Double_t fBranchValue; // input from raw JAPAN output
  
  inline void SetDefUsage(Bool_t flag){kUseDefinition = flag;};
  inline Bool_t IsUsingDefinition(){return kUseDefinition;};
  vector<TString> GetRawChannelList(){
    return myDefinition->GetRawChannelList();
  };
  vector<Double_t> GetFactorArray(){
    return myDefinition->GetFactorArray();
  };
  Bool_t HasUserDefinition(){
    return myDefinition->HasUserDefinition();
  };
  
private:
  TString fTreeName;
  TaDefinition* myDefinition;
  Bool_t kUseDefinition;
  TString fChannelName;
  vector<TaChannel*> fChannels;
  vector<Double_t> fPrefactors;
  vector<Double_t> fDataArray;

  TaAccumulator fRunAccumulator;
  TaAccumulator fMiniAccumulator;
  STAT mini_stat;
  STAT run_stat;
  ClassDef(TaChannel,0);
};

#endif
