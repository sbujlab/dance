#include "MyFunctions.C"
#include "../src/TaAccumulator.cc"
void MakeSensMatrix(TString slug_name);
void MakeSensMatrix(Int_t slug_id){
  TString slug_name = Form("./crex-runlist/slug%d.list",slug_id);
  MakeSensMatrix(slug_name);
}
void MakeSensMatrix(TString slug_name){
  FILE *runlist = fopen(slug_name,"r");
  int run_idx;
  vector<int> fRunNumberList;
  while(!feof(runlist)){
    fscanf(runlist,"%d\n",&run_idx);
    fRunNumberList.push_back(run_idx);
  }
  fclose(runlist);

  Int_t coilID[]={1,3,4,6,7};
  vector<TString> dv_array = {"usl","usr","dsl","dsr",
			      "sam1","sam2","sam3","sam4",
			      "sam5","sam6","sam7","sam8",
			      "bpm11X","bpm11Y","bpm12X","bpm12Y","bpm16X","bpm16Y",
			      "bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm1X","bpm1Y"};
  vector<TString> coil_array;
  Int_t nCoil = sizeof(coilID)/sizeof(*coilID);
  for(int icoil=0;icoil<nCoil;icoil++)
    coil_array.push_back(Form("bmod_trim%d",coilID[icoil]));

  Int_t nDV = dv_array.size();
  TMatrixD sens_matrix(nDV,nCoil);
  TChain *sens_tree = new TChain("sens");
  Int_t nrun = fRunNumberList.size();
  for(int i=0;i<nrun;i++){
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.root",
			  fRunNumberList[i]);
    sens_tree->Add(filename);
  }

  Int_t nEntries = sens_tree ->GetEntries();
  Double_t cycNumber;
  
  Double_t *fSens = new Double_t[nCoil*nDV];
  Double_t *fError = new Double_t[nCoil*nDV];
  for(int idv=0;idv<nDV;idv++)
    for(int icoil=0;icoil<nCoil;icoil++){
      sens_tree->SetBranchAddress(Form("%s_coil%d",dv_array[idv].Data(),coilID[icoil]), 
				  &(fSens[idv*nCoil+icoil]));
      sens_tree->SetBranchAddress(Form("%s_coil%d_err",dv_array[idv].Data(),coilID[icoil]), 
				  &(fError[idv*nCoil+icoil]));
    }

  sens_tree->SetBranchAddress("cycID",&cycNumber);
  vector<TaAccumulator> fAccumulator(nCoil*nDV);
  for(int ievt=0;ievt<nEntries;ievt++){
    sens_tree->GetEntry(ievt);

    for(int idv=0;idv<nDV;idv++)
      for(int icoil=0;icoil<nCoil;icoil++){
	if(fError[idv*nCoil+icoil]>0){
	  fAccumulator[idv*nCoil+icoil].Update(fSens[idv*nCoil+icoil]);
	}
      }

  }
  for(int idv=0;idv<nDV;idv++)
    for(int icoil=0;icoil<nCoil;icoil++)
      sens_matrix[idv][icoil] = fAccumulator[idv*nCoil+icoil].GetMean1();

  TString out_filename = Form("./dit-coeffs/prex_sens_matrix.%d-%d.root",
			      fRunNumberList[0],fRunNumberList[nrun-1]);
  TFile *matrix_output = TFile::Open(out_filename,"RECREATE");
  matrix_output->WriteObject(&sens_matrix,"sens_matrix");
  matrix_output->WriteObject(&dv_array,"dv_array");
  matrix_output->WriteObject(&coil_array,"coil_array");
  matrix_output->Close();

}

