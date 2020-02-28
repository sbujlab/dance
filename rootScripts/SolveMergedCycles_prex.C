#include "TSolver.cc"
#include "utilities.cc"
void SolveMergedCycles_prex(Int_t);
void SolveMergedCycles_prex(TString,Int_t);
void SolveMergedCycles_prex(Int_t slug_number=6){
  TString runlist = Form("./prex-runlist/simple_list/slug%d.list",slug_number);
  SolveMergedCycles_prex(runlist,slug_number);
}
void SolveMergedCycles_prex(TString runlist, Int_t slug_number=9999){
  vector<Int_t> fRunList = ParseRunList(runlist);
  map< Int_t, vector<Int_t> > fblmap = GetCycleBlackList();
  gStyle->SetOptStat(0);
  TFile *merged_rf = TFile::Open(Form("./dit-coeffs/slug%d.root",slug_number));
  TTree *sens = (TTree*)merged_rf->Get("sens");
  Int_t coil_index[] = {1,3,4,6,7,2,5};
  vector<Int_t> required_index;
  vector<Int_t> candidate_set1= {1,3,4,6,7}; // run 3130-4734;
  vector<Int_t> candidate_set2= {3,5,4,6,7}; // run 4735-4980;
  if(slug_number<=76){
    required_index= candidate_set1;
    sens->AddFriend("dit_slope1");
  }else if(slug_number<=94){
    required_index=candidate_set2;
    sens->AddFriend("dit_slope3");
  }else{
    required_index=candidate_set1;
    sens->AddFriend("dit_slope1");
  }
  vector<TString> mon_set1={"bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm11X12X"};  // run 3404 - 4980
  vector<TString> mon_set2={"bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm12X"}; // run 3130-3403
  vector<TString> mon_array;
  if(slug_number<=3)
    mon_array = mon_set2;
  else if(slug_number<=94)
    mon_array = mon_set1;
  else
    mon_array = mon_set2;
  TString det_array[]={"usl","usr"};
  
  const Int_t ncoil = sizeof(coil_index)/sizeof(*coil_index);
  const Int_t nmon = mon_array.size();
  const Int_t ndet = sizeof(det_array)/sizeof(*det_array);

  vector<Double_t> fdummy_det(ndet,0.0);
  vector<Double_t> fdummy_mon(nmon,0.0);

  vector<vector<Double_t> > det_val(ncoil,fdummy_det);
  vector<vector<Double_t> > det_err(ncoil,fdummy_det);
  vector<vector<Double_t> > mon_val(ncoil,fdummy_mon);
  vector<vector<Double_t> > mon_err(ncoil,fdummy_mon);

  Double_t cycID;
  Double_t run;
  sens->SetBranchAddress("cycID",&cycID);
  sens->SetBranchAddress("run",&run);
  TString branch_name;
  for(int icoil=0;icoil<ncoil;icoil++){
    for(int imon=0;imon<nmon;imon++){
      branch_name = Form("%s_coil%d_err",
			 mon_array[imon].Data(),coil_index[icoil]);
				 
      sens->SetBranchAddress(branch_name, &mon_err[icoil][imon]);

      branch_name = Form("%s_coil%d",
			 mon_array[imon].Data(),coil_index[icoil]);
      
      sens->SetBranchAddress(branch_name, &mon_val[icoil][imon]);
    }
    for(int idet=0;idet<ndet;idet++){
      branch_name = Form("%s_coil%d_err",
			 det_array[idet].Data(),coil_index[icoil]);
				 
      sens->SetBranchAddress(branch_name, &det_err[icoil][idet]);

      branch_name = Form("%s_coil%d",
			 det_array[idet].Data(),coil_index[icoil]);

      sens->SetBranchAddress(branch_name, &det_val[icoil][idet]);
    }
  }

  Double_t slope1_val[ndet][nmon];
  Double_t slope1_flag[ndet][nmon];
  for(int idet=0;idet<ndet;idet++){
    for(int imon=0;imon<nmon;imon++){
      TString ch_name = Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data());
      TString flag_name = Form("%s_%s_flag",det_array[idet].Data(),mon_array[imon].Data());
      
      sens->SetBranchAddress(ch_name,&slope1_val[idet][imon]);
      sens->SetBranchAddress(flag_name,&slope1_flag[idet][imon]);
    }
  }
  vector< Double_t> fproto_vec;
  vector< vector<Double_t > > fSlope1(ndet*nmon,fproto_vec);
  vector< vector<Double_t > > fSlope1Xcord(ndet*nmon,fproto_vec);
  
  Int_t nevt = sens->GetEntries();
  Bool_t kMatch;
  Int_t ncycles=0;
  Int_t last_run = 0;
  vector<TSolver> fSolver;
  TSolver* aSolver;
  vector<Double_t> fCycleIDBuff;
  vector<Int_t> fCycleArray;
  TFile* output = TFile::Open(Form("./rootfiles/slug%d_dit_slope_run.root",slug_number),"RECREATE");
  TTree *slope_tree  = new TTree("slope","slope from merged cycles");
  Double_t fCycID, fRun;
  vector< vector<Double_t> > fSlopes(ndet,fdummy_mon);
  vector< vector<Double_t> > fSlopesXcord; // [isol][icyc]

  Int_t fErrorFlag;
  // slope_tree->Branch("cycID",&fCycID);
  slope_tree->Branch("run",&fRun);
  slope_tree->Branch("ErrorFlag",&fErrorFlag,"ErrorFlag/I");
  for(int idet=0;idet<ndet;idet++){
    for(int imon=0;imon<nmon;imon++){
      TString ch_name=Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data());
      slope_tree->Branch(ch_name,&fSlopes[idet][imon]);
    }
  }
  TMultiGraph* fMgArray[ndet][nmon];
  for(int idet=0;idet<ndet;idet++)
    for(int imon=0;imon<nmon;imon++)
      fMgArray[idet][imon]=new TMultiGraph();
  
  for(int ievt=0;ievt<nevt;ievt++){
    sens->GetEntry(ievt);
    if(find(fRunList.begin(),fRunList.end(),(Int_t)run)==fRunList.end())
      kMatch = kFALSE;
    else
      kMatch = kTRUE;
    if(kMatch){
      fCycleArray.push_back(cycID);
      ncycles++;
      if(run!=last_run ){
	cout << "loading solver : " << run << endl;
	if(last_run!=0){
	  fSolver.push_back(*aSolver);
	  fSlopesXcord.push_back(fCycleIDBuff);
	}
	aSolver = new TSolver(run);
	aSolver->SetCondition(required_index);
	last_run = run;
	fCycleIDBuff.clear();
      }
      // ==== Loading Sensitivity data ====
      fCycleIDBuff.push_back(ncycles-1);
      for(int icoil=0;icoil<ncoil;icoil++){
	if(IsGoodCoil(det_err[icoil])
	   && !IsBadCycle(fblmap,cycID,coil_index[icoil] ))
	  aSolver->LoadSensivity(coil_index[icoil],det_val[icoil],mon_val[icoil]);
      }
      // ==== END of Loading  ===
      
      // ==== Loading Slope data ====
      for(int idet=0;idet<ndet;idet++)
	for(int imon=0;imon<nmon;imon++)
	  if(slope1_flag[idet][imon] && !IsBadCycle(fblmap,cycID) ){
	    fSlope1[idet*nmon+imon].push_back(slope1_val[idet][imon]*1e3);
	    fSlope1Xcord[idet*nmon+imon].push_back(ncycles-1);
	  }
      // ==== END of Loading  ===
      
    }// end of if kMatch
  } // end of event loop
  fSolver.push_back(*aSolver);
  fSlopesXcord.push_back(fCycleIDBuff);

  vector<Int_t> fRunLabel;
  vector<Double_t> fRunLabelXcord;

  auto itsol = fSolver.begin();
  while(itsol!=fSolver.end()){
    Int_t count = itsol-fSolver.begin();
    fRun = (*itsol).GetID();
    fRunLabel.push_back(fRun);
    fRunLabelXcord.push_back(fSlopesXcord[count][0]);

    cout << " -- Solving run " << fRun << endl;
    if((*itsol).SolveMatrix()){
      (*itsol).WriteSolution(fSlopes);
      // auto icyc = fCycleIDBuff.begin();
      // while(icyc!=fCycleIDBuff.end()){
      //   fCycID = (*icyc);
      slope_tree->Fill();
      //   icyc++;
      // }
      for(int idet=0;idet<ndet;idet++){
	for(int imon=0;imon<nmon;imon++){
	  TGraph* g1= GraphAverageSlope(fSlopes[idet][imon],fSlopesXcord[count]);
	  fMgArray[idet][imon]->Add(g1,"l");
	}
      }
    }
    itsol++;
  }
  
  map<Int_t, vector<vector<Double_t> > > fSlugAvgSlopeMap;
  fSlugAvgSlopeMap = GetSlugAveragedSlopesMap(slug_number,det_array,ndet,
					     mon_array,nmon);
  cout << " -- loading averaged slopes " << endl;
  auto it_run = fRunLabel.begin();
  Int_t counter =0;
  while(it_run!=fRunLabel.end()){
    if(fSlugAvgSlopeMap.find(*it_run)==fSlugAvgSlopeMap.end()){
      counter++;
      it_run++;
      continue;
    }
    vector<vector<Double_t > > fSlugAvgSlope = fSlugAvgSlopeMap[*it_run];
    for(int idet=0;idet<ndet;idet++){
      for(int imon=0;imon<nmon;imon++){
  	TGraph* g1= GraphAverageSlope(fSlugAvgSlope[idet][imon],
				      fSlopesXcord[counter],0.5);
  	g1->SetLineColor(kBlue);
  	fMgArray[idet][imon]->Add(g1,"l");
      }
    }
    counter++;
    it_run++;
  }
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->cd();
  // c2->SetGridx();
  c2->SetRightMargin(0.015);
  for(int idet=0;idet<ndet;idet++){
    for(int imon=0;imon<nmon;imon++){
      TString title = Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data());
      TGraph *gtest = GraphVector(fSlope1[idet*nmon+imon],fSlope1Xcord[idet*nmon+imon]);
      fMgArray[idet][imon]->Add(gtest,"p");
      fMgArray[idet][imon]->Draw("A");
      fMgArray[idet][imon]->SetTitle(title);
      fMgArray[idet][imon]->GetYaxis()->SetTitle("ppm/um");
      TH1F *hfmg = fMgArray[idet][imon]->GetHistogram();
      hfmg->GetXaxis()->Set(ncycles,-0.5,ncycles-0.5);
      for(int i=0;i<ncycles;i++)
	hfmg->GetXaxis()->SetBinLabel(i+1,Form("%d",fCycleArray[i]));
      hfmg->SetMarkerColor(kWhite);
      hfmg->Draw("same p");

      Double_t ymin = fMgArray[idet][imon]->GetYaxis()->GetXmin();
      Double_t ymax = fMgArray[idet][imon]->GetYaxis()->GetXmax();
      Int_t ntext = fRunLabelXcord.size();
      for(int i=0;i<ntext;i++){
	TText *t1 = new TText(fRunLabelXcord[i]-0.5,ymax,Form("%d",fRunLabel[i]));
	TLine *line1 = new TLine(fRunLabelXcord[i]-0.5,ymin,
				 fRunLabelXcord[i]-0.5,ymax);
	line1->SetLineColor(kMagenta);
	line1->SetLineWidth(1);
	line1->SetLineStyle(7);
	t1->SetTextAngle(27);
	t1->SetTextSize(0.03);
	line1->Draw("same");
	t1->Draw("same");
      }
      c2->SaveAs(Form("./pdf/slug%d_dit_slope_%s.pdf",slug_number,title.Data()));
    }
  }
  gSystem->Exec(Form("pdfunite $(ls -rt ./pdf/slug%d_dit_slope_*.pdf) ./pdf/slug%d_dit_slope.pdf",
		     slug_number,slug_number));
  gSystem->Exec(Form("rm -f ./pdf/slug%d_dit_slope_*.pdf ",slug_number));
  output->cd();
  slope_tree->Write();
  output->Close();
}


		  
