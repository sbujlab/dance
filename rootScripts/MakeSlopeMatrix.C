#include "MyFunctions.C"

void MakeSlopeMatrix(){
  vector<TString> mon_array{"bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm11X12X"};
  vector<TString> det_array{"usl","usr","dsl","dsr"};
  Int_t nMon = mon_array.size();
  Int_t nDet = det_array.size();
  TString filename="./dit-coeffs/prexPrompt_ditcoeffs_4978.root";
  TFile* input = TFile::Open(filename);
  TTree* slope_tree = (TTree*)input->Get("dit_slope1");
  Int_t nEntries = slope_tree ->GetEntries();
  Double_t cycNumber;
  slope_tree->SetBranchAddress("cycID",&cycNumber);
  Double_t *fSlope = new Double_t[nMon*nDet];
  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){
      slope_tree->SetBranchAddress(Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data()),
				   &(fSlope[idet*nMon+imon]));
    }
  }
  slope_tree->GetEntry(0);
  TMatrixD slope_matrix(nDet,nMon);
  for(int idet=0;idet<nDet;idet++)
    for(int imon=0;imon<nMon;imon++)
      slope_matrix[idet][imon] = fSlope[idet*nMon+imon];
  
  TFile *matrix_output = TFile::Open("./dit-coeffs/prex_slope_matrix.root","RECREATE");
  matrix_output->WriteObject(&slope_matrix,"slope_matrix");
  matrix_output->WriteObject(&det_array,"dv_array");
  matrix_output->WriteObject(&mon_array,"iv_array");
  matrix_output->Close();
  input->Close();    
}

