#include "../src/TaAccumulator.cc"
#include "MyFunctions.C"
void MakeSlopeMatrix(Int_t slug_id, TString tree_name);
void MakeSlopeMatrix(TString slug_name, TString tree_name);
void MakeSlopeMatrix(Int_t slug_id,TString tree_name){
  TString slug_name = Form("./crex-runlist/slug%d.list",slug_id);
  MakeSlopeMatrix(slug_name,tree_name);
}
void MakeSlopeMatrix(TString slug_name,TString tree_name){
  FILE *runlist = fopen(slug_name,"r");
  int run_idx;
  vector<int> fRunNumberList;
  while(!feof(runlist)){
    fscanf(runlist,"%d\n",&run_idx);
    fRunNumberList.push_back(run_idx);
  }
  fclose(runlist);

  vector<TString> det_array = {"usl","usr","dsl","dsr",
			       "sam1","sam2","sam3","sam4",
			       "sam5","sam6","sam7","sam8"};

  vector<TString> mon_set1={"bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm12X"};
  vector<TString> mon_set2={"bpm1X","bpm4aY","bpm4eX","bpm4eY","bpm12X"};
  vector<TString> mon_array;
  if(tree_name=="dit_slope1" || tree_name=="dit_slope2"|| tree_name=="dit_slope3")
    mon_array = mon_set1;
  else
    mon_array = mon_set2;
  Int_t nMon = mon_array.size();
  Int_t nDet = det_array.size();

  TChain *slope_tree = new TChain(tree_name);
  Int_t nrun = fRunNumberList.size();
  for(int i=0;i<nrun;i++){
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.root",
			  fRunNumberList[i]);
    slope_tree->Add(filename);
  }

  Int_t nEntries = slope_tree ->GetEntries();
  Double_t cycNumber;
  slope_tree->SetBranchAddress("cycID",&cycNumber);
  Double_t *fSlope = new Double_t[nMon*nDet];
  Double_t *kFlag = new Double_t[nMon*nDet];
  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){
      slope_tree->SetBranchAddress(Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data()),
				   &(fSlope[idet*nMon+imon]));
      slope_tree->SetBranchAddress(Form("%s_%s_flag",det_array[idet].Data(),mon_array[imon].Data()),
				   &(kFlag[idet*nMon+imon]));

    }
  }
  vector<TaAccumulator> fAccumulator(nMon*nDet);
  TMatrixD slope_matrix(nDet,nMon);

  for(int ievt=0;ievt<nEntries;ievt++){
    slope_tree->GetEntry(ievt);
    for(int idet=0;idet<nDet;idet++)
      for(int imon=0;imon<nMon;imon++){
	if(kFlag[idet*nMon+imon])
	  fAccumulator[idet*nMon+imon].Update(fSlope[idet*nMon+imon]);
      }
  }
  for(int idet=0;idet<nDet;idet++)
    for(int imon=0;imon<nMon;imon++)
      slope_matrix[idet][imon] = fAccumulator[idet*nMon+imon].GetMean1();

  TString out_filename = Form("./dit-coeffs/prex_%s_matrix.%d-%d.root",
			      tree_name.Data(),
			      fRunNumberList[0],fRunNumberList[nrun-1]);

  TFile *matrix_output = TFile::Open(out_filename,"RECREATE");
  matrix_output->WriteObject(&slope_matrix,"slope_matrix");
  matrix_output->WriteObject(&det_array,"dv_array");
  matrix_output->WriteObject(&mon_array,"iv_array");
  matrix_output->Close();

}

