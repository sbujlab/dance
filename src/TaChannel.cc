#include "TaChannel.hh"
ClassImp(TaChannel);

TaChannel::TaChannel():fOutputValue(0.0){
  kUseDefinition = kFALSE;
}

TaChannel::TaChannel(TString tree, TaDefinition* adef)
  :fOutputValue(0.0){
  fTreeName = tree;
  fChannelName = adef->GetName();
  myDefinition = adef;
  kUseDefinition = kFALSE;
}

TaChannel::TaChannel(TString tree, TString channel)
  :fOutputValue(0.0){
  fTreeName = tree;
  fChannelName = channel;
  kUseDefinition = kFALSE;
}

TaChannel::~TaChannel(){}

void TaChannel::ConnectChannels(vector<TaChannel*> in_channels, 
				vector<Double_t> in_prefactors){
  SetDefUsage(kTRUE);
  fChannels = in_channels;
  fPrefactors = in_prefactors;
}

void TaChannel::DefineSubtraction(TaChannel* raw, TaChannel* correction){
  fChannels.clear();
  fPrefactors.clear();
  vector<TaChannel*> ch_pair= {raw,correction};
  vector<Double_t> factor_pair= {1,-1};
  ConnectChannels(ch_pair,factor_pair);
}

void TaChannel::CalcCombination(){
  Int_t nch = fChannels.size();
  fOutputValue=0.0;
  for(int i=0;i<nch;i++){
    fOutputValue+=fPrefactors[i]*(fChannels[i]->fOutputValue);
  }
}

void TaChannel::FillDataArray(){
  FillOutputValue();
  fDataArray.push_back(fOutputValue);
}

Double_t TaChannel::GetEntry(Int_t ie){
  fOutputValue = fDataArray[ie];
  return fDataArray[ie];
}

void TaChannel::ConstructTreeBranch(TaOutput* fOutput){
  fOutput->ConstructTreeBranch(fTreeName,fChannelName,fOutputValue);
}

void TaChannel::ConstructMiniTreeBranch(TaOutput* fOutput,TString treename){
  fOutput->ConstructStatTreeBranch(treename,fChannelName, mini_stat);
}

void TaChannel::ConstructSlopeBranch(TaOutput* fOutput,TString treename){
  Int_t nIV = fPrefactors.size();
  TString dv_base = GetChannelBaseName();
  for(int iiv=0;iiv<nIV;iiv++){
    TString iv_base = fChannels[iiv]->GetChannelBaseName();
    TString branch_name = dv_base+"_"+iv_base;
    fOutput->ConstructTreeBranch(treename,branch_name,fPrefactors[iiv]);
  }
}

void TaChannel::ConstructSumTreeBranch(TaOutput* fOutput,TString treename){
  fOutput->ConstructStatTreeBranch(treename,fChannelName, run_stat);
}

Double_t TaChannel::FillOutputValue(){
  if(kUseDefinition){
    CalcCombination();
  }
  return fOutputValue;
}
void TaChannel::AccumulateRunSum(){
  fRunAccumulator.Update(fOutputValue);
}
void TaChannel::AccumulateMiniSum(){
  fMiniAccumulator.Update(fOutputValue);
}

TString TaChannel::GetChannelBaseName(){
  TString base_name = fChannelName;
  base_name.ReplaceAll("cor_","");
  base_name.ReplaceAll("asym_","");
  base_name.ReplaceAll("diff_","");
  return base_name;
}
