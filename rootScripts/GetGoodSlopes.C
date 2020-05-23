#include "utilities.cc"
#include "plot_util.cc"
void GetGoodSlopes(Int_t slug_id);
void GetGoodSlopes(){
  for(int i=1;i<=94;i++)
    GetGoodSlopes(i);
}
void GetGoodSlopes(Int_t slug_id){
  TString tree_name;
  if(slug_id<=76 || slug_id==501)
    tree_name ="dit_slope1";
  else if(slug_id<=94 )
    tree_name ="dit_slope3";
  else if(slug_id>=143)
    tree_name ="dit1x_slope1";
  else 
    tree_name ="dit_slope1";

  Bool_t kCrex = kFALSE;
  if(slug_id>=100 && slug_id<500)
    kCrex = kTRUE;
  map<Int_t, vector<Int_t> > fBadCycleMap = LoadBadCycleList(kCrex);

  vector<Int_t> fRunList = LoadRunListBySlug(slug_id);
  map<Int_t, Int_t> fArmMap = LoadArmMapBySlug(slug_id);
  vector<TString> det_array = {"usl","usr","dsl","dsr"};
  vector<TString> at_array={"atl1","atl2","atr1","atr2"};
  vector<TString> sam_array={"sam1","sam2","sam3","sam4",
			     "sam5","sam6","sam7","sam8"};
  vector<TString> mon_array={"bpm4aX","bpm4eX","bpm4aY","bpm4eY"};
  vector<TString> mon1x_array={"bpm1X","bpm4eX","bpm4aY","bpm4eY","bpm12X"};
  if(slug_id<=2)
    mon_array.push_back("bpm12X");
  else if(slug_id<=94 || slug_id==501)
    mon_array.push_back("bpm11X12X");
  else
    mon_array.push_back("bpm12X");

  if(slug_id>=143 && kCrex && slug_id!=501)
    mon_array=mon1x_array;

  if(slug_id>=26)
    det_array.insert(det_array.end(),at_array.begin(),at_array.end());
  if(slug_id>=100 && slug_id!=501)
    det_array.insert(det_array.end(),sam_array.begin(),sam_array.end());

  Int_t nMon = mon_array.size();
  Int_t nDet = det_array.size();

  TChain *slope_tree = new TChain(tree_name);
  Int_t nrun = fRunList.size();
  for(int i=0;i<nrun;i++){
    Int_t seg_number =0;
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.%03d.root",
			  fRunList[i],seg_number);
    while(gSystem->AccessPathName(filename)==0){
      slope_tree->Add(filename);
      seg_number++;
      filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.%03d.root",
		    fRunList[i],seg_number);
    }

  }

  Int_t nCycles = slope_tree ->GetEntries();
  if(nCycles==0)
    return;

  Double_t cycNumber;
  slope_tree->SetBranchAddress("cycID",&cycNumber);
  Double_t runNumber;
  slope_tree->SetBranchAddress("run",&runNumber);

  vector<Double_t> fSlope(nMon*nDet);
  vector<Double_t> kFlag(nMon*nDet);
  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){
      slope_tree->SetBranchAddress(Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data()),
				   &(fSlope[idet*nMon+imon]));
      slope_tree->SetBranchAddress(Form("%s_%s_flag",det_array[idet].Data(),mon_array[imon].Data()),
				   &(kFlag[idet*nMon+imon]));

    }
  }

  // ++++++++++
  TString rootfile_name = Form("./slopes/slug%d_dit_slope_by_cycle.root",slug_id);
  TFile* slope_output = TFile::Open(rootfile_name,"RECREATE");
  slope_output->cd();
  TTree* dit_tree = new TTree("dit","dit");
  vector<Double_t> fSlope_val(nDet*nMon);
  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){  
      TString branch_name =Form("%s_%s",
				det_array[idet].Data(),
				mon_array[imon].Data());
      dit_tree->Branch(branch_name,&fSlope_val[idet*nMon+imon]);
      
      // if(slug_id<=2){
      // 	if(mon_array[imon]=="bpm12X")
      // 	  dit_tree->Branch(Form("%s_bpmE",det_array[idet].Data()),
      // 			   &fSlope_val[idet*nMon+imon]);
      // }else{
      // 	if(mon_array[imon]=="bpm11X12X")
      // 	  dit_tree->Branch(Form("%s_bpmE",det_array[idet].Data()),
      // 			   &fSlope_val[idet*nMon+imon]);
      // }
    }
  }
  Int_t fRun,fCycID,fArmFlag;
  dit_tree->Branch("run",&fRun);
  dit_tree->Branch("cycID",&fCycID);
  dit_tree->Branch("arm_flag",&fArmFlag);

  for(int ievt=0;ievt<nCycles;ievt++){
    slope_tree->GetEntry(ievt);
    fRun  = runNumber;
    fCycID = cycNumber;
    Int_t arm_flag = fArmMap[runNumber];
    fArmFlag = arm_flag;

    if(IsBadCycle(fBadCycleMap,cycNumber))
      continue;

    Bool_t kGoodCycle = kTRUE;
    for(int idet=0;idet<nDet;idet++){
      // Not checking slope flag for bad arm data
      if(arm_flag==1 && det_array[idet].Contains("l"))
	continue;
      if(arm_flag==2 && det_array[idet].Contains("r"))
	continue;
      
      for(int imon=0;imon<nMon;imon++){
	if(!kFlag[idet*nMon+imon])
	  kGoodCycle = kFALSE;
      }// end of mon loop
    } // end of det loop

    if(kGoodCycle){
      for(int idet=0;idet<nDet;idet++)
	for(int imon=0;imon<nMon;imon++)
	  if(kFlag[idet*nMon+imon])
	    fSlope_val[idet*nMon+imon] = fSlope[idet*nMon+imon];

      dit_tree->Fill();
    }
  } // end of ievt loop
  
  dit_tree->Write();
  slope_output->Close();
  cout << " Writing " << rootfile_name << endl;

}

