#include "../src/TaAccumulator.cc"
#include "MyFunctions.C"
void AverageSlope(Int_t slug_id,TString tree_name){
  // FIXME : Get Run list and split list
  // Load Runlist
  // Load RunInfo
  // Load badcycle cuts
  // Load Split Definition
  
  vector<TString> det_array = {"usl","usr","dsl","dsr"};
			       // "sam1","sam2","sam3","sam4",
			       // "sam5","sam6","sam7","sam8"};
  vector<TString> at_array={"atl1","atl2","atr1","atr2"};

  vector<TString> mon_array={"bpm4aX","bpm4aY","bpm4eX","bpm4eY"};

  if(slug_id<=3)
    mon_array.push_back("bpm12X");
  else if(slug_id<=94)
    mon_array.push_back("bpm11X12X");
  else
    mon_array.push_back("bpm12X");

  if(slug_id>=26)
    det_array.insert(det_array.end(),at_array.begin(),at_array.end());

  Int_t nMon = mon_array.size();
  Int_t nDet = det_array.size();

  TChain *slope_tree = new TChain(tree_name);
  Int_t nrun = fRunList.size();
  for(int i=0;i<nrun;i++){
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.root",
			  fRunList[i]);
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
  
  // FIXME: Cut on Arm flag , bad cycle, bad coil 

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


  // Output :
  // TMatrixD to rootfiles: dv_array definition
  // JAPAN mapfiles 
  // Plots cycle wise slope vs averaged slopes
  
  TString out_filename = Form("./dit-coeffs/prex_%s_matrix.%d-%d.root",
			      tree_name.Data(),
			      fRunList[0],fRunList[nrun-1]);

  TFile *matrix_output = TFile::Open(out_filename,"RECREATE");
  matrix_output->WriteObject(&slope_matrix,"slope_matrix");
  matrix_output->WriteObject(&det_array,"dv_array");
  matrix_output->WriteObject(&mon_array,"iv_array");
  matrix_output->Close();

}

