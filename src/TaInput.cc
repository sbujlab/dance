#include "TaInput.hh"
#include "TaOutput.hh"

#include <iostream>
#include "TLeaf.h"
#include "TBranch.h"
#include "TString.h"
#include "TEventList.h"
#include "TCut.h"

ClassImp(TaInput);
using namespace std;

TaInput::TaInput(TaConfig *aConfig){
  run_number= aConfig->GetRunNumber();
  seg_number=0;
  isExternalConstraint=kFALSE;
  InitChannels(aConfig);
  input_name = aConfig->GetInputName();
  if(run_number==0){
    Ssiz_t pos_head = input_name.Last('_');
    Ssiz_t pos_end = input_name.Last('.');
    Ssiz_t length = pos_end-pos_head-1;
    TString run_dot_seg = TString(input_name(pos_head+1,length));
    Ssiz_t pos_dot = run_dot_seg.Last('.');
    run_number = TString(run_dot_seg(0,pos_dot)).Atoi();
    run_dot_seg.Remove(0,pos_dot+1);
    seg_number = run_dot_seg.Atoi();
    aConfig->SetRunNumber(run_number);
  }
  input_path=aConfig->GetConfigParameter("input_path");
  input_prefix=aConfig->GetConfigParameter("input_prefix");
  minirun_size = (aConfig->GetConfigParameter("minirun_size")).Atof();
}
TaInput::~TaInput(){}

Bool_t TaInput::LoadROOTFile(){
#ifdef NOISY
  cout <<  __FUNCTION__ << endl;
#endif  
  if(run_number!=0){
    if(input_path=="") 
      input_path = TString(getenv("QW_ROOTFILES")); 
    input_name = input_path+"/"+input_prefix
      +Form("%d.%03d.root",(int)run_number,(int)seg_number);
    cout << " -- Opening "
	 << input_name << endl;
    input_file = TFile::Open(input_name);
  }
  else{
    cout << " -- Opening "
	 << input_name << endl;
    input_file = TFile::Open(input_name);
  }

  if(input_file==NULL){
    cerr << __PRETTY_FUNCTION__
	 << " Error: Input file is not found !! " << endl;
    return kFALSE;
  }

  evt_tree = (TTree*)input_file->Get("evt");
  mul_tree = (TTree*)input_file->Get("mul");
  mulc_tree = (TTree*)input_file->Get("mulc");
  if(mulc_tree!=NULL)
    mul_tree->AddFriend(mulc_tree);
  
  return kTRUE;
}

void TaInput::InitChannels(TaConfig *aConfig){
#ifdef NOISY
  cout <<  __FUNCTION__ << endl;
#endif  
  vector<TString> device_list = aConfig->GetDeviceList();
  int nDevice = device_list.size();
  for(int i=0;i<nDevice;i++){
    TaChannel* aChannel = new TaChannel("mul",device_list[i]);
    fChannelArray.push_back(aChannel);
    fChannelNames.push_back(device_list[i]);
    fChannelMap[device_list[i]]=i;
  }

  fChannelErrorFlag = new TaChannel("mul","ErrorFlag");
  fChannelArray.push_back(fChannelErrorFlag);
  fChannelNames.push_back("ErrorFlag");
  fChannelMap["ErrorFlag"]=nDevice;

  fChannelCutFlag = new TaChannel("mul","ok_cut");
}

void TaInput::WriteRawChannels(TaOutput *aOutput){
#ifdef NOISY
  cout <<  __PRETTY_FUNCTION__ << endl;
#endif  

  TCut default_cut("ErrorFlag==0");
  TEventList *elist_mul = new TEventList("elist_mul");  
  Int_t nGoodPatterns = mul_tree->Draw(">>+elist_mul", default_cut, "goff");
  cout << " -- nGoodPatterns in Mul Tree : " << nGoodPatterns << endl;
  
  mul_tree->SetBranchStatus("*",0);
  int nch = fChannelNames.size();
  for(int i=0;i<nch;i++){
    mul_tree->SetBranchStatus(fChannelNames[i],1);
    Double_t* fValue_ptr = &(fChannelArray[i]->fBranchValue);
    TBranch *aBranch = mul_tree->GetBranch(fChannelNames[i]);
    if(aBranch!=NULL){
      TLeaf* aLeaf = aBranch->GetLeaf("hw_sum");
      if(aLeaf!=NULL)
	aLeaf->SetAddress(fValue_ptr);
      else
	aBranch->SetAddress(fValue_ptr);
    }
    else
      cout << "TBranch " <<fChannelNames[i] << " not found " << endl;
  }
  
  for(int ich=0;ich<nch;ich++){
    fChannelArray[ich]->ConstructTreeBranch(aOutput);
    fChannelArray[ich]->ConstructMiniTreeBranch(aOutput,"mini");
    fChannelArray[ich]->ConstructSumTreeBranch(aOutput,"sum");
  }
  fChannelCutFlag->ConstructTreeBranch(aOutput);
  
  Double_t mini_id=0;
  aOutput->ConstructTreeBranch("sum","run",run_number);
  aOutput->ConstructTreeBranch("mini","mini",mini_id);
  aOutput->ConstructTreeBranch("mul","mini",mini_id);
  Int_t ievt =0;
  Double_t goodCounts=0;
  Int_t mini_start =0;
  Int_t mini_end =0;
  Bool_t isGoodPattern=kFALSE;
  while((mul_tree->GetEntry(ievt))>0){
    if(elist_mul->GetIndex(ievt)!=-1){
      isGoodPattern=kTRUE;
      goodCounts++;
      fChannelCutFlag->fOutputValue = 1;
      fChannelCutFlag->FillDataArray();
    }
    else{
      isGoodPattern=kFALSE;
      fChannelCutFlag->fOutputValue = 0;
      fChannelCutFlag->FillDataArray();
    }
    for(int ich=0;ich<nch;ich++){
      fChannelArray[ich]->FillOutputValue();
      fChannelArray[ich]->FillDataArray();
      if(isGoodPattern){
	fChannelArray[ich]->AccumulateRunSum();
	fChannelArray[ich]->AccumulateMiniSum();
      }
    }

    aOutput->FillTree("mul");
    ievt++;

    if(goodCounts==minirun_size){
      mini_end = ievt-1;
      minirun_range.push_back(make_pair(mini_start,mini_end));
      cout << " -- Mini-run ends at event: " << mini_end << endl;
      mini_start = mini_end+1;
      cout << " -- Next Mini-run starts at event: " << mini_start << endl;
      goodCounts = 0;
      int nMinirun = minirun_range.size();
      Bool_t is_last_minrun = kFALSE;
      if(nGoodPatterns-minirun_size*nMinirun<minirun_size){
	is_last_minrun = kTRUE;
	mini_start = minirun_range[nMinirun-1].first;
	minirun_range.pop_back();
	
	cout << " -- Meeting last mini-run, " << endl;
	cout << " -- the rest will be merged into this mini-run  "  << endl;
      }
      if(!is_last_minrun){
	for(int ich=0;ich<nch;ich++){
	  fChannelArray[ich]->UpdateMiniStat();
	  fChannelArray[ich]->ResetMiniAccumulator();
	}
	aOutput->FillTree("mini");
	mini_id++;
      }
    }
  }
  cout << " -- last mini-run ends at event: " << ievt-1 << endl;
  minirun_range.push_back(make_pair(mini_start,ievt-1));

  for(int ich=0;ich<nch;ich++){
    fChannelArray[ich]->UpdateMiniStat();
    fChannelArray[ich]->UpdateRunStat();
    fChannelArray[ich]->ResetMiniAccumulator();
  }
  aOutput->FillTree("mini");
  aOutput->FillTree("sum");
}


void TaInput::Close(){
  input_file->Close();
}


TaChannel* TaInput::GetChannel(TString name){
  Int_t index = fChannelMap[name];
  return fChannelArray[index];
}
