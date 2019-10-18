#include "TaChannel.hh"
ClassImp(TaChannel);

TaChannel::TaChannel():fOutputValue(0.0){
}

TaChannel::TaChannel(TString tree, TString channel, Int_t index)
  :fOutputValue(0.0){
  fTreeName = tree;
  fChannelName = channel;
  fChannelIndex = index;
}

TaChannel::~TaChannel(){}

void TaChannel::ConnectChannels(vector<TaChannel*> in_channels, 
				vector<Double_t> in_prefactors){
  fChannels = in_channels;
  fPrefactors = in_prefactors;
}

void TaChannel::ProcessCombination(){
  Int_t nch = fChannels.size();
  for(int i=0;i<nch;i++){
    fOutputValue+=fPrefactors[i]*(fChannels[i]->fOutputValue);
  }
}

void TaChannel::FillDataArray(){
  fDataArray.push_back(fOutputValue);
}

void TaChannel::ApplyCorrection(){
  ProcessCombination();
}

Double_t TaChannel::GetEntry(Int_t ie){
  return fDataArray[ie];
}

void TaChannel::ConstructTreeBranch(TaOutput* fOutput){
  fOutput->ConstructTreeBranch(fTreeName,fChannelName,fOutputValue);
}

void TaChannel::ConstructMiniTreeBranch(TaOutput* fOutput,TString treename){
  fOutput->ConstructStatTreeBranch(treename,fChannelName, mini_stat);
}
void TaChannel::ConstructSumTreeBranch(TaOutput* fOutput,TString treename){
  fOutput->ConstructStatTreeBranch(treename,fChannelName, run_stat);
}

void TaChannel::FillOutputValue(){
  fOutputValue=fBranchValue;
}

void TaChannel::AccumulateRunSum(){
  fRunAccumulator.Update(fOutputValue);
}
void TaChannel::AccumulateMiniSum(){
  fMiniAccumulator.Update(fOutputValue);
}
