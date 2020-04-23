#include "TaOutput.hh"
#include "TaAccumulator.hh"
#include "TaConfig.hh"
#include "TaDefinition.hh"
#include "TTree.h"
#include "TLeaf.h"

#ifndef __TaChannel_hh__
#define __TaChannel_hh__

using namespace std;
class TaOutput;
class TaAccumulator;
class TaChannel {
public:
  TaChannel();
  TaChannel(TString tree, TaDefinition* adef);
  TaChannel(TString tree, TString name);
  virtual ~TaChannel();
  void RegisterBranchAddress(TBranch *aBranch);
  void ConnectChannels(vector<TaChannel*>, vector<Double_t> );
  void DefineSubtraction(TaChannel*,TaChannel*);
  void CalcCombination();
  void FillDataArray();
  Double_t FillOutputValue();

  void AccumulateRunSum();
  void AccumulateMiniSum();

  void UpdateRunStat(){
    fRunAccumulator.UpdateStat(run_stat);
    for(int i=0;i<4;i++)
      fRunAccumulator_block[i].UpdateStat(run_stat_block[i]);
  };
  void UpdateMiniStat(){
    fMiniAccumulator.UpdateStat(mini_stat);
    for(int i=0;i<4;i++)
      fMiniAccumulator_block[i].UpdateStat(mini_stat_block[i]);
  };
  
  void ResetMiniAccumulator(){
    fMiniAccumulator.Zero();
    for(int i=0;i<4;i++)
      fMiniAccumulator_block[i].Zero();
  };
  void ResetRunAccumulator(){
    fRunAccumulator.Zero();
    for(int i=0;i<4;i++)
      fRunAccumulator_block[i].Zero();
  };

  void ConstructTreeBranch(TaOutput*);
  void ConstructTreeBranch(TaOutput*,TString leaflist);
  void ConstructMiniTreeBranch(TaOutput*,TString);
  void ConstructSlopeBranch(TaOutput*,TString);
  void ConstructSumTreeBranch(TaOutput*,TString);

  TString GetChannelBaseName();
  Double_t GetEntry(Int_t ie);
  ROOTDATA fOutputValue;  
  
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
  vector<ROOTDATA> fDataArray;

  TaAccumulator fRunAccumulator;
  TaAccumulator fMiniAccumulator;
  STAT mini_stat;
  STAT run_stat;
  
  TaAccumulator fRunAccumulator_block[4];
  TaAccumulator fMiniAccumulator_block[4];
  STAT mini_stat_block[4];
  STAT run_stat_block[4];

  ClassDef(TaChannel,0);
};

#endif
