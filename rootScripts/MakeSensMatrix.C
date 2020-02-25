#include "MyFunctions.C"

void MakeSensMatrix(){
  Int_t coilID[]={1,3,4,6,7};
  vector<TString> dv_array{"usl","usr","dsl","dsr",
			   "bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm11X","bpm12X"};
  vector<TString> coil_array;
  Int_t nCoil = sizeof(coilID)/sizeof(*coilID);
  for(int icoil=0;icoil<nCoil;icoil++){
    coil_array.push_back(Form("bmod_trim%d",coilID[icoil]));
  }
  Int_t nDV = dv_array.size();
  TMatrixD sens_matrix(nDV,nCoil);

  TString filename="./dit-coeffs/prexPrompt_ditcoeffs_4978.root";
  TFile* input = TFile::Open(filename);
  TTree* sens_tree = (TTree*)input->Get("sens");
  Int_t nEntries = sens_tree ->GetEntries();
  Double_t cycNumber;
  Double_t *fSens = new Double_t[nCoil*nDV];

  for(int idv=0;idv<nDV;idv++)
    for(int icoil=0;icoil<nCoil;icoil++)
      sens_tree->SetBranchAddress(Form("%s_coil%d",dv_array[idv].Data(),coilID[icoil]), &(fSens[idv*nCoil+icoil]));

  sens_tree->SetBranchAddress("cycID",&cycNumber);
  sens_tree->GetEntry(0);
  
  for(int idv=0;idv<nDV;idv++)
    for(int icoil=0;icoil<nCoil;icoil++)
      sens_matrix[idv][icoil] = fSens[idv*nCoil+icoil];

  
  TFile *matrix_output = TFile::Open("./dit-coeffs/prex_sens_matrix.root","RECREATE");
  matrix_output->WriteObject(&sens_matrix,"sens_matrix");
  matrix_output->WriteObject(&dv_array,"dv_array");
  matrix_output->WriteObject(&coil_array,"coil_array");
  matrix_output->Close();
  input->Close();    
}

