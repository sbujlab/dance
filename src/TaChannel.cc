#include "TaChannel.hh"
ClassImp(TaChannel);

TaChannel::TaChannel(){
  kUseDefinition = kFALSE;
  fOutputValue.hw_sum=0.0;
}

TaChannel::TaChannel(TString tree, TaDefinition* adef) {
  fTreeName = tree;
  fChannelName = adef->GetName();
  myDefinition = adef;
  kUseDefinition = kFALSE;
  fOutputValue.hw_sum=0.0;
}

TaChannel::TaChannel(TString tree, TString channel){
  fTreeName = tree;
  fChannelName = channel;
  kUseDefinition = kFALSE;
  fOutputValue.hw_sum=0.0;
}

TaChannel::~TaChannel(){}

void TaChannel::RegisterBranchAddress(TBranch* aBranch){
  TLeaf* hw_sum = aBranch->GetLeaf("hw_sum");
  if(hw_sum!=NULL){
    hw_sum->SetAddress(&(fOutputValue.hw_sum));
    aBranch->GetLeaf("block0")->SetAddress(&(fOutputValue.block0));
    aBranch->GetLeaf("block1")->SetAddress(&(fOutputValue.block1));
    aBranch->GetLeaf("block2")->SetAddress(&(fOutputValue.block2));
    aBranch->GetLeaf("block3")->SetAddress(&(fOutputValue.block3));
  }else
    aBranch->SetAddress(&(fOutputValue.hw_sum));
  
}

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
  fOutputValue.hw_sum=0.0;
  fOutputValue.block0=0.0;
  fOutputValue.block1=0.0;
  fOutputValue.block2=0.0;
  fOutputValue.block3=0.0;
  for(int i=0;i<nch;i++){
    // fChannels[i]->FillOutputValue(); // FIXME: think about a better recursive approach
    fOutputValue.hw_sum +=fPrefactors[i]*(fChannels[i]->fOutputValue.hw_sum);
    fOutputValue.block0 +=fPrefactors[i]*(fChannels[i]->fOutputValue.block0);
    fOutputValue.block1 +=fPrefactors[i]*(fChannels[i]->fOutputValue.block1);
    fOutputValue.block2 +=fPrefactors[i]*(fChannels[i]->fOutputValue.block2);
    fOutputValue.block3 +=fPrefactors[i]*(fChannels[i]->fOutputValue.block3);
  }
}

void TaChannel::FillDataArray(){
  FillOutputValue(); // What's the purpose ... Load values from combination or input 
  fDataArray.push_back(fOutputValue);
}

Double_t TaChannel::GetEntry(Int_t ie){
  fOutputValue = fDataArray[ie];
  return fOutputValue.hw_sum;
}

void TaChannel::ConstructTreeBranch(TaOutput* fOutput,TString leaflist){
  fOutput->ConstructTreeBranch(fTreeName,fChannelName,leaflist,fOutputValue);
}

void TaChannel::ConstructTreeBranch(TaOutput* fOutput){
  fOutput->ConstructTreeBranch(fTreeName,fChannelName,fOutputValue.hw_sum);
}

void TaChannel::ConstructMiniTreeBranch(TaOutput* fOutput,TString treename){
  fOutput->ConstructStatTreeBranch(treename,fChannelName, mini_stat);
  for(int iblk=0;iblk<4;iblk++){
    TString suffix = Form("_block%d",iblk);
    fOutput->ConstructStatTreeBranch(treename,fChannelName+suffix, mini_stat_block[iblk]);
  }
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
  for(int iblk=0;iblk<4;iblk++){
    TString suffix = Form("_block%d",iblk);
    fOutput->ConstructStatTreeBranch(treename,fChannelName+suffix, run_stat_block[iblk]);
  }
}

Double_t TaChannel::FillOutputValue(){
  if(kUseDefinition)
    CalcCombination();
  return fOutputValue.hw_sum;
}
void TaChannel::AccumulateRunSum(){
  fRunAccumulator.Update(fOutputValue.hw_sum);
  fRunAccumulator_block[0].Update(fOutputValue.block0);
  fRunAccumulator_block[1].Update(fOutputValue.block1);
  fRunAccumulator_block[2].Update(fOutputValue.block2);
  fRunAccumulator_block[3].Update(fOutputValue.block3);
}
void TaChannel::AccumulateMiniSum(){
  fMiniAccumulator.Update(fOutputValue.hw_sum);
  fMiniAccumulator_block[0].Update(fOutputValue.block0);
  fMiniAccumulator_block[1].Update(fOutputValue.block1);
  fMiniAccumulator_block[2].Update(fOutputValue.block2);
  fMiniAccumulator_block[3].Update(fOutputValue.block3);
}

TString TaChannel::GetChannelBaseName(){
  TString base_name = fChannelName;
  base_name.ReplaceAll("cor_","");
  base_name.ReplaceAll("asym_","");
  base_name.ReplaceAll("diff_","");
  return base_name;
}
