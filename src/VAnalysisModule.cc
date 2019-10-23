#include "VAnalysisModule.hh"

ClassImp(VAnalysisModule);
VAnalysisModule::VAnalysisModule(){
  
}

void VAnalysisModule::Init(TString type, TString name,TaConfig *aConfig){
  Int_t ana_index = aConfig->GetAnalysisIndex(type,name);
  tree_name = aConfig->GetAnalysisParameter(ana_index,"tree_name");
  branch_prefix = aConfig->GetAnalysisParameter(ana_index,"branch_prefix");
  sDVlist = aConfig->GetDVlist(type,name);
  sIVlist = aConfig->GetIVlist(type,name);
  Int_t nDV = sDVlist.size();
  vector<Double_t> subtraction{1,-1};
  for(int ich=0;ich<nDV;ich++){
    TaChannel *aOutputChannel = new TaChannel(tree_name,branch_prefix+sDVlist[ich]);
    TaChannel *aCorrection = new TaChannel(tree_name,"cor_"+sDVlist[ich]);
    fOutputChannels.push_back(aOutputChannel);
    fCorrections.push_back(aCorrection);
  }
}

void VAnalysisModule::LoadInput(TaInput *aInput){

  fChannelCutFlag = aInput->GetChannelCutFlag();
  minirun_range = aInput->GetMiniRange();
  Int_t nDV = sDVlist.size();
  for(int idv=0;idv<nDV;idv++){
    TaChannel *aChannel = aInput->GetChannel(sDVlist[idv]);
    fDependentVar.push_back(aChannel);
  }

  for(int ich=0;ich<nDV;ich++)
    fOutputChannels[ich]->DefineSubtraction(fDependentVar[ich],fCorrections[ich]);

  Int_t nIV = sIVlist.size();
  for(int iiv=0;iiv<nIV;iiv++){
    TaChannel *aChannel = aInput->GetChannel(sIVlist[iiv]);
    fIndependentVar.push_back(aChannel);
  }

}

void VAnalysisModule::ConstructOutputs(TaOutput* fOutput){
  Int_t nDV = sDVlist.size();
  for(int ich=0;ich<nDV;ich++){
    fOutputChannels[ich]->ConstructTreeBranch(fOutput);
    fOutputChannels[ich]->ConstructMiniTreeBranch(fOutput,"mini_"+tree_name);
    fOutputChannels[ich]->ConstructSumTreeBranch(fOutput,"sum_"+tree_name);

    fCorrections[ich]->ConstructTreeBranch(fOutput);
    fCorrections[ich]->ConstructMiniTreeBranch(fOutput,"mini_"+tree_name);
    fCorrections[ich]->ConstructSumTreeBranch(fOutput,"sum_"+tree_name);
  }
  
}

void VAnalysisModule::GetEntry(Int_t ievt){
  fChannelCutFlag->GetEntry(ievt);
  Int_t nDV = sDVlist.size();
  for(int idv=0;idv<nDV;idv++)
    fDependentVar[idv]->GetEntry(ievt);
  Int_t nIV = sIVlist.size();
  for(int iiv=0;iiv<nIV;iiv++)
    fIndependentVar[iiv]->GetEntry(ievt);
}

void VAnalysisModule::CalcCombination(){
  Int_t nDV = sDVlist.size();
  for(int ich=0;ich<nDV;ich++){
    fCorrections[ich]->CalcCombination();
    fOutputChannels[ich]->CalcCombination();
  }
}

void VAnalysisModule::AccumulateMiniSum(){
  Int_t nDV = sDVlist.size();
  if(fChannelCutFlag->fOutputValue){
    for(int ich=0;ich<nDV;ich++){
      fCorrections[ich]->AccumulateMiniSum();
      fOutputChannels[ich]->AccumulateMiniSum();
    }
  }
}

void VAnalysisModule::AccumulateRunSum(){
  Int_t nDV = sDVlist.size();
  if(fChannelCutFlag->fOutputValue){
    for(int ich=0;ich<nDV;ich++){
      fCorrections[ich]->AccumulateRunSum();
      fOutputChannels[ich]->AccumulateRunSum();
    }
  }
}

void VAnalysisModule::ResetMiniAccumulator(){
  Int_t nDV = sDVlist.size();
  for(int ich=0;ich<nDV;ich++){
    fCorrections[ich]->ResetMiniAccumulator();
    fOutputChannels[ich]->ResetMiniAccumulator();
  }
}

void VAnalysisModule::ResetRunAccumulator(){
  Int_t nDV = sDVlist.size();
  for(int ich=0;ich<nDV;ich++){
    fCorrections[ich]->ResetRunAccumulator();
    fOutputChannels[ich]->ResetRunAccumulator();
  }
}

void VAnalysisModule::UpdateMiniStat(){
  Int_t nDV = sDVlist.size();
  for(int ich=0;ich<nDV;ich++){
    fCorrections[ich]->UpdateMiniStat();
    fOutputChannels[ich]->UpdateMiniStat();
  }
}

void VAnalysisModule::UpdateRunStat(){
  Int_t nDV = sDVlist.size();
  for(int ich=0;ich<nDV;ich++){
    fCorrections[ich]->UpdateRunStat();
    fOutputChannels[ich]->UpdateRunStat();
  }
}
