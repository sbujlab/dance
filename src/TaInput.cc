#include "TaInput.hh"
#include <iostream>

ClassImp(TaInput);
using namespace std;

TaInput::TaInput(TaConfig *aconfig){
  fConfig=aconfig;
  run_number=0;
  seg_number=0;
  isExternalConstraint=kFALSE;
}
TaInput::~TaInput(){}
Bool_t TaInput::LoadROOTFile(){
#ifdef NOISY
  cout <<  __FUNCTION__ << endl;
#endif  
  TString input_name;
  if(run_number!=0){
    TString input_path=fConfig->GetInputPath();
    TString input_prefix=fConfig->GetInputPrefix();
    if(input_path=="") 
      input_path = TString(getenv("QW_ROOTFILES")); 
    input_name = input_path+"/"+input_prefix
      +Form("%d.%03d.root",run_number,seg_number);
    cout << " -- Opening "
	 << input_name << endl;
    input_file = TFile::Open(input_name);
  }
  else{
    cout << " -- Opening "
	 << fConfig->GetInputName() << endl;
    input_file = TFile::Open(fConfig->GetInputName());
    input_name = fConfig->GetInputName();
    Ssiz_t pos_head = input_name.Last('_');
    Ssiz_t pos_end = input_name.Last('.');
    Ssiz_t length = pos_end-pos_head-1;
    TString run_dot_seg = TString(input_name(pos_head+1,length));
    Ssiz_t pos_dot = run_dot_seg.Last('.');
    run_number = TString(run_dot_seg(0,pos_dot)).Atoi();
    run_dot_seg.Remove(0,pos_dot+1);
    seg_number = run_dot_seg.Atoi();
  }

  if(input_file==NULL){
    cerr << __PRETTY_FUNCTION__
	 << " Error: Input file is not found !! " << endl;
    return kFALSE;
  }

  if(fConfig->GetAnalysisType()=="regression"){
    mul_tree = (TTree*)input_file->Get("mul");
  }
  else if(fConfig->GetAnalysisType()=="sens"){
    evt_tree = (TTree*)input_file->Get("evt");
  }
  else if(fConfig->GetAnalysisType()=="lagrangian"){
    evt_tree = (TTree*)input_file->Get("evt");
    mul_tree = (TTree*)input_file->Get("mul");
  }
  else{
    cout << " Error: Analysis type is not defined" << endl;
    cout << " Error: Analysis Aborted" << endl;
    return kFALSE;
  }
  return kTRUE;
}
