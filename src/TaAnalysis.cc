#include "TaAnalysis.hh"

ClassImp(TaAnalysis);

TaAnalysis::TaAnalysis(TaInput *aInput){
  fInput=aInput;
}
TaAnalysis::~TaAnalysis(){}

Bool_t TaAnalysis::LoadConfig(TaConfig* aConfig){
#ifdef NOISY
  cout <<  __FUNCTION__ << endl;
#endif  
  fConfig = aConfig;
  analysisType = fConfig->GetAnalysisType();
  TString output_path = fConfig->GetOutputPath();
  TString output_prefix = fConfig->GetOutputPrefix();
  if(output_prefix=="")
    output_prefix = "prexPrompt_"+analysisType+"_";
  Int_t run_number = fInput->GetRunNumber();
  Int_t seg_number = fInput->GetSegNumber();
  TString run_seg = Form("%d_%03d",run_number,seg_number);
  TString output_name =output_path+output_prefix+run_seg+".root";
  cout << " -- Writing Output to " 
       <<output_name << endl;
  output_rootfile = TFile::Open(output_name,"RECREATE");
  if(output_rootfile==NULL){
    cerr << " ** Error : Can't open " 
	 << output_name << endl
      	 << " Analysis aborted " 
	 << endl;
    return kFALSE;
  }
  if(analysisType=="sens"){
    fDithering = new TaDithering(fInput);
    return fDithering->LoadConfig(fConfig);
  }
  else if(analysisType=="lagrangian"){
    if(!fInput->UseExternalConstraint()){
      fDithering = new TaDithering(fInput);
      if(!fDithering->LoadConfig(fConfig)){
	return kFALSE;
      }
    }
    fLagrangian = new TaLagrangian(fInput);
    return fLagrangian->LoadConfig(fConfig);
  }
  else if (analysisType=="regression"){
    fLagrangian = new TaLagrangian(fInput);
    return fLagrangian->LoadConfig(fConfig);
  }
  else
    return kFALSE;
}

Bool_t TaAnalysis::Process(){
#ifdef NOISY
  cout <<  __FUNCTION__ << endl;
#endif  
  if(analysisType=="sens"){
    if(fDithering->LoadModulationData()){
      fDithering->ComputeSensitivities();
      fDithering->PrintSummary();
      output_rootfile->cd();
      (fDithering->GetSensTree())->Write();
      TMatrixD detMatrix = fDithering->GetDetSensMatrix();
      TMatrixD monMatrix = fDithering->GetMonSensMatrix();
      output_rootfile->WriteObject(&detMatrix,"DetSens");
      output_rootfile->WriteObject(&monMatrix,"MonSens");
      return kTRUE;
    }
    else
      return kFALSE;
  }
  else if(analysisType=="regression"){
    fLagrangian->ComputeCorrelation();
    fLagrangian->ComputeSlopes();
    fLagrangian->CorrectTree();
    fLagrangian->WriteSummary();
    output_rootfile->cd();
    fLagrangian->GetMiniTree()->Write();
    fLagrangian->GetSlopeTree()->Write();
    fLagrangian->GetCorrectionTree()->Write();
    return kTRUE;
  }
  else if(analysisType=="lagrangian"){
    TMatrixD detMatrix,monMatrix;
    if(fInput->UseExternalConstraint()){
      TString ext_filename= fInput->GetExtFileName();
      TFile* ext_rootfile = TFile::Open(ext_filename);
      if (ext_rootfile == NULL){
	cerr << " ** Error: "
	     << ext_filename 
	     << " doesn't exist !! " << endl;
	cerr << " Analysis aborted " << endl;
	return kFALSE;
      }
      vector<TString>* src_det_array = (vector<TString>*)ext_rootfile->Get("det_array");
      vector<TString>* src_mon_array = (vector<TString>*)ext_rootfile->Get("mon_array");
      vector<TString>* src_coil_array = (vector<TString>*)ext_rootfile->Get("coil_array");
      if( src_coil_array==NULL ||
	  src_mon_array==NULL ||
	  src_det_array==NULL ) {
	cerr << " ** Error: detector/monitor/coil array not found" << endl;
	cerr << "Analysis Aborted " << endl;
	return kFALSE;
      }
      vector<TString> det_array=fConfig->GetDetArray();
      vector<TString> mon_array=fConfig->GetMonArray();
      vector<TString> coil_array=fConfig->GetCoilArray();
      if(!isArrayMatched(det_array,*src_det_array))
	 return kFALSE;
      if(!isArrayMatched(mon_array,*src_mon_array))
	 return kFALSE;
      if(!isArrayMatched(coil_array,*src_coil_array))
	 return kFALSE;

      TMatrixD* detMatrix_ptr = (TMatrixD*)ext_rootfile->Get("DetSens");
      TMatrixD* monMatrix_ptr = (TMatrixD*)ext_rootfile->Get("MonSens");
      if(detMatrix_ptr!=NULL && monMatrix_ptr!=NULL){
	cout << " -- Loading sensitivity constraint from "
	     << ext_filename << endl;
	int nrow_det = detMatrix_ptr->GetNrows();
	int ncol_det = detMatrix_ptr->GetNcols();
	for(int irow=0;irow<nrow_det;irow++){
	  for(int icol=0;icol<ncol_det;icol++){
	    if((*detMatrix_ptr)[irow][icol]==0){
	      cerr << " ** Error: detector constraint matrix has at least on missing coil "<< endl;
	      cerr << " ** Analysis aborted " << endl;
	      return kFALSE;
	    }
	  }
	}
	
	int nrow_mon = monMatrix_ptr->GetNrows();
	int ncol_mon = monMatrix_ptr->GetNcols();
	for(int irow=0;irow<nrow_mon;irow++){
	  for(int icol=0;icol<ncol_mon;icol++){
	    if((*monMatrix_ptr)[irow][icol]==0){
	      cerr << " ** Error: monitor constraint matrix has at least on missing coil "<< endl;
	      cerr << " ** Analysis aborted " << endl;
	      return kFALSE;
	    }
	  }
	}

	detMatrix.ResizeTo(*detMatrix_ptr);
	detMatrix=*detMatrix_ptr;
	monMatrix.ResizeTo(*monMatrix_ptr);
	monMatrix=*monMatrix_ptr;
      }
      else{
	cerr << " ** Error: constraint matrices not found in "
	     << ext_filename << endl;
	return kFALSE;
      }
      output_rootfile->WriteObject(&detMatrix,"DetSens");
      output_rootfile->WriteObject(&monMatrix,"MonSens");
      ext_rootfile->Close();
    } // end of if UseExternalConstraint
    else{
      if(fDithering->LoadModulationData()){
	fDithering->ComputeSensitivities();
	fDithering->PrintSummary();
	output_rootfile->cd();
	(fDithering->GetSensTree())->Write();
	detMatrix.ResizeTo(fDithering->GetDetSensMatrix());
	monMatrix.ResizeTo(fDithering->GetMonSensMatrix());
	detMatrix = fDithering->GetDetSensMatrix();
	monMatrix = fDithering->GetMonSensMatrix();
	output_rootfile->WriteObject(&detMatrix,"DetSens");
	output_rootfile->WriteObject(&monMatrix,"MonSens");
      }
      else
	return kFALSE;
    }
    fLagrangian->LoadConstraint(detMatrix,monMatrix);
    fLagrangian->ComputeCorrelation();
    fLagrangian->ComputeSlopes();
    fLagrangian->CorrectTree();
    fLagrangian->WriteSummary();
    output_rootfile->cd();
    fLagrangian->GetMiniTree()->Write();
    fLagrangian->GetSlopeTree()->Write();
    fLagrangian->GetCorrectionTree()->Write();
    return kTRUE;
  }
  else
    return kFALSE;
}
void TaAnalysis::End(){
  vector<TString> det_array = fConfig->GetDetArray();
  vector<TString> mon_array = fConfig->GetMonArray();
  vector<TString> coil_array = fConfig->GetCoilArray();
  output_rootfile->WriteObject(&det_array,"det_array");
  output_rootfile->WriteObject(&mon_array,"mon_array");
  output_rootfile->WriteObject(&coil_array,"coil_array");
  output_rootfile->Close();
}

Bool_t TaAnalysis::isArrayMatched(vector<TString> a, vector<TString> b){
  if(a.size()!=b.size())
    return kFALSE;
  else{
    int nsize=a.size();
    for(int i=0;i<nsize;i++){
      if(a[i]!=b[i]){
	cerr << " ** Error : Array element mismatched where "
	     << a[i] 
	     << " != "
	     << b[i] << endl;
	cerr << " ** Analysis aborted ! "<< endl;
	return kFALSE;
      }
    }
  }
  return kTRUE;
}
