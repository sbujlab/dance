#include "../src/TaAccumulator.cc"
#include "TSolver.cc"
#include "utilities.cc"
#include "plot_util.cc"
void SolveMergedCycles(Int_t slug_id, Bool_t kMatrixOutput=kFALSE);
void SolveMergedCycles(){
  for(int i=1;i<=94;i++)
    SolveMergedCycles(i,kTRUE);
}
void SolveMergedCycles(Int_t slug_id,Bool_t kMatrixOutput){

  Bool_t kCrex = kFALSE;
  if(slug_id>=100 && slug_id<500)
    kCrex = kTRUE;
  map<Int_t, vector<Int_t> > fBadCycleMap = LoadBadCycleList(kCrex);
  vector<Int_t> fRunList = LoadRunListBySlug(slug_id);
  map<Int_t, Int_t> fArmMap = LoadArmMapBySlug(slug_id);
  vector< vector<Int_t> > fRunListArray = LoadSplitListBySlug(slug_id);
  Int_t nSplits = fRunListArray.size();
  map<Int_t, Int_t > fSplitMap;
  vector< Int_t> range_low;
  vector< Int_t> range_up;
  for(int isplit=0;isplit<nSplits;isplit++){
    vector<Int_t> this_list = fRunListArray[isplit];
    Int_t nrun = this_list.size();
    for(int i=0;i<nrun;i++)
      fSplitMap[ this_list[i] ] = isplit;
    range_low.push_back(this_list[0]);
    range_up.push_back(this_list[nrun-1]);
  }
  vector<TString> det_array = {"usl","usr","dsl","dsr"};
  vector<TString> sam_array={"sam1","sam2","sam3","sam4",
			     "sam5","sam6","sam7","sam8"};
  vector<TString> at_array={"atl1","atl2","atr1","atr2"};
  vector<TString> mon_array={"bpm4aX","bpm4eX","bpm4aY","bpm4eY"};
  vector<TString> mon1x_array={"bpm1X","bpm4eX","bpm4aY","bpm4eY","bpm12X"};
  if(slug_id<=2)
    mon_array.push_back("bpm12X");
  else if(slug_id<=94 || slug_id==501)
    mon_array.push_back("bpm11X12X");
  else
    mon_array.push_back("bpm12X");

  if(slug_id>=143 && slug_id!=501)
    mon_array=mon1x_array;

  if(slug_id>=26)
    det_array.insert(det_array.end(),at_array.begin(),at_array.end());
  if(slug_id>=100 && slug_id!=501)
    det_array.insert(det_array.end(),sam_array.begin(),sam_array.end());

  Int_t nMon = mon_array.size();
  Int_t nDet = det_array.size();
  
  Int_t coil_index[] = {1,3,4,6,7,2,5};
  vector<Int_t> required_index;
  vector<Int_t> candidate_set1= {1,3,4,6,7}; // run 3130-4734;
  vector<Int_t> candidate_set2= {3,5,4,6,7}; // run 4735-4980;
  
  if(slug_id<=76)
    required_index= candidate_set1;
  else if(slug_id<=94)
    required_index=candidate_set2;
  else
    required_index=candidate_set1;
  Int_t nCoil = sizeof(coil_index)/sizeof(*coil_index);
  TChain *sens_tree = new TChain("sens");
  Int_t nrun = fRunList.size();
  for(int i=0;i<nrun;i++){
    Int_t seg_number =0;    
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.%03d.root",
			  fRunList[i],seg_number);
    while(gSystem->AccessPathName(filename)==0){
      sens_tree->Add(filename);
      seg_number++;
      filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.%03d.root",
		    fRunList[i],seg_number);
    }
  }

  Int_t nCycles = sens_tree ->GetEntries();
  if(nCycles==0)
    return;
  Double_t cycNumber;
  sens_tree->SetBranchAddress("cycID",&cycNumber);
  Double_t runNumber;
  sens_tree->SetBranchAddress("run",&runNumber);
  vector<Double_t> fdummy_det(nDet);
  vector<Double_t> fdummy_mon(nMon);
  vector<vector<Double_t> > det_val(nCoil,fdummy_det);
  vector<vector<Double_t> > det_err(nCoil,fdummy_det);
  vector<vector<Double_t> > mon_val(nCoil,fdummy_mon);
  vector<vector<Double_t> > mon_err(nCoil,fdummy_mon);

  TString branch_name;
  for(int icoil=0;icoil<nCoil;icoil++){
    for(int imon=0;imon<nMon;imon++){
      branch_name = Form("%s_coil%d_err",
			 mon_array[imon].Data(),coil_index[icoil]);
				 
      sens_tree->SetBranchAddress(branch_name, &mon_err[icoil][imon]);

      branch_name = Form("%s_coil%d",
			 mon_array[imon].Data(),coil_index[icoil]);
      
      sens_tree->SetBranchAddress(branch_name, &mon_val[icoil][imon]);
    }
    for(int idet=0;idet<nDet;idet++){
      branch_name = Form("%s_coil%d_err",
			 det_array[idet].Data(),coil_index[icoil]);
				 
      sens_tree->SetBranchAddress(branch_name, &det_err[icoil][idet]);

      branch_name = Form("%s_coil%d",
			 det_array[idet].Data(),coil_index[icoil]);

      sens_tree->SetBranchAddress(branch_name, &det_val[icoil][idet]);
    }
  }

  TSolver fProtoSolver;
  fProtoSolver.SetCondition(required_index);
  vector<TSolver> fSolver(nDet,fProtoSolver);
  vector< vector<TSolver> > fSolverArrayByRange(nSplits,fSolver);
  vector< vector<TSolver> > fSolverArrayByRun;//(fRunList.size(),fSolver);
  vector<Int_t> fRunLabel;
  vector<Int_t> fCycleNumber;
  vector<Double_t> fRunLabelXcord;
  vector< vector<Double_t> > fSplitXcord(nSplits);
  vector< vector<vector<Double_t> > > fRunXcord(nDet);
  vector< vector<Double_t> > fRunCycleList(fRunList.size());
  vector< vector<Double_t> > fAveragedSlopeByRange(nDet*nMon);//[idet*nMon+imon][split]
  vector< vector<Double_t> > fAveragedSlopeByRun(nDet*nMon);//[idet*nMon+imon][run]

  Int_t last_runnumber=0;
  Int_t run_count = -1;
  for(int ievt=0;ievt<nCycles;ievt++){
    sens_tree->GetEntry(ievt);
    fCycleNumber.push_back(cycNumber);
    if(last_runnumber!=runNumber){
      cout << "reading run " <<runNumber <<endl;
      last_runnumber = runNumber;
      fRunLabel.push_back(runNumber);
      fRunLabelXcord.push_back(ievt);
      run_count++;
      fSolverArrayByRun.push_back(fSolver);
    }
    Int_t split_id = fSplitMap[runNumber];
    fSplitXcord[split_id].push_back(ievt);
    fRunCycleList[run_count].push_back(ievt);
    if(IsBadCycle(fBadCycleMap,cycNumber))
      continue;
    Int_t arm_flag = fArmMap[runNumber];
    for(int icoil=0;icoil<nCoil;icoil++){
      if(!IsGoodCoil(mon_err[icoil]))
	 continue;
      for(int idet=0;idet<nDet;idet++){
	if(arm_flag==1 && det_array[idet].Contains("l"))
	  continue;
	if(arm_flag==2 && det_array[idet].Contains("r"))
	  continue;

	fSolverArrayByRange[split_id][idet].LoadSensitivity(coil_index[icoil],
							    det_val[icoil][idet],
							    mon_val[icoil]);
	fSolverArrayByRun[run_count][idet].LoadSensitivity(coil_index[icoil],
							   det_val[icoil][idet],
							   mon_val[icoil]);

      } // end of det loop
    } // end of coil loop
  } // end of ievt loop
  
  // ++++++++++ Solving
  for(int irun=0;irun<=run_count;irun++){
    cout << "run " <<  fRunLabel[irun]<< endl;
    for(int idet=0;idet<nDet;idet++){
      vector<Double_t> fSlope_buff(nMon,0.0);
      if(fSolverArrayByRun[irun][idet].SolveMatrix())
	fSolverArrayByRun[irun][idet].WriteSolution(fSlope_buff);
      for(int imon=0;imon<nMon;imon++)
	fAveragedSlopeByRun[idet*nMon+imon].push_back( fSlope_buff[imon] );
      fRunXcord[idet].push_back(fRunCycleList[irun]);
    }
  }

  for(int isplit=0;isplit<nSplits;isplit++){
    ofstream mapfile;
    TString fullpath="./mapfiles_ovcn/";
    TString range_tag;
    int low = range_low[isplit];
    int up = range_up[isplit];
    if(low==up)
      range_tag = Form("%d",low);
    else
      range_tag = Form("%d-%d",low,up);

    TString mapname = "prex_combiner_dit."+range_tag+".map";
    fullpath+=mapname;
    mapfile.open(fullpath.Data());
    cout << "Writing map " << fullpath << endl;
    for(int idet=0;idet<nDet;idet++){
      cout << det_array[idet] << endl;
      vector<Double_t> fSlope_buff(nMon,0.0);
      if(fSolverArrayByRange[isplit][idet].SolveMatrix())
	fSolverArrayByRange[isplit][idet].WriteSolution(fSlope_buff);
      mapfile << "[asym:@dit_asym_"<<det_array[idet]<<"]" << endl;
      mapfile << "asym:"<<det_array[idet] << ",1.0"<< endl;
      for(int imon=0;imon<nMon;imon++){
	TString iv_name = mon_array[imon];
	if(iv_name=="bpm11X12X")
	  iv_name = "diff_bpmE";
	mapfile << "diff:"<<iv_name << ","
		<< -fSlope_buff[imon]<< endl;
	// Yeah, has a minus sign.
	fAveragedSlopeByRange[idet*nMon+imon].push_back(fSlope_buff[imon]);
      }
      mapfile << endl;
    }
    mapfile.close();
  }

  // ++++++++++ Writing TTree
  TString rootfile_name = Form("./slopes/slug%d_dit_slope_merged_cycle_ovcn.root",slug_id);
  cout << " Creating " << rootfile_name << endl;
  TFile* avg_output = TFile::Open(rootfile_name,"RECREATE");
  avg_output->cd();
  TTree* dit_range = new TTree("dit1","dit slopes by ranges");
  TTree* dit_run = new TTree("dit2","dit slopes by run");
  vector<Double_t> fSlope_val(nDet*nMon);
  Int_t fRun,fCounter;  

  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){  
      TString branch_name =Form("%s_%s",
				det_array[idet].Data(),
				mon_array[imon].Data());
      dit_range->Branch(branch_name,&fSlope_val[idet*nMon+imon]);
    }
  }

  dit_range->Branch("run",&fRun);
  dit_range->Branch("range",&fCounter);
  cout << " Filling run-range averaged slopes" << endl;
  for(int isplit=0;isplit<nSplits;isplit++){
    vector<Int_t> this_list = fRunListArray[isplit];
    Int_t nrun = this_list.size();
    fCounter = isplit;
    for(int idet=0;idet<nDet;idet++)
      for(int imon=0;imon<nMon;imon++)
	fSlope_val[idet*nMon+imon] = fAveragedSlopeByRange[idet*nMon+imon][isplit];

    for(int i=0;i<nrun;i++){
      fRun = this_list[i];
      dit_range->Fill();
    }
  }
  dit_range->Write();
  /// --- slope by run

  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){  
      TString branch_name =Form("%s_%s",
				det_array[idet].Data(),
				mon_array[imon].Data());
      dit_run->Branch(branch_name,&fSlope_val[idet*nMon+imon]);
    }
  }

  dit_run->Branch("run",&fRun);

  cout << " Filling  run-by-run slopes " << endl;
  for(int i=0;i<=run_count;i++){
    fRun = fRunLabel[i];
    for(int idet=0;idet<nDet;idet++)
      for(int imon=0;imon<nMon;imon++)
	fSlope_val[idet*nMon+imon] = fAveragedSlopeByRun[idet*nMon+imon][i];
    dit_run->Fill();
  }
  dit_run->Write();
  
  TDirectory *graph_dir = avg_output->GetDirectory("/")->mkdir("graph");
  TDirectory *canvas_dir = avg_output->GetDirectory("/")->mkdir("canvas");
  
  TString input_name = Form("./slopes/slug%d_dit_slope_cyclewise_average.root",slug_id);
  TFile *cycle_input = TFile::Open(input_name);
  TDirectory *input_dir = cycle_input->GetDirectory("graph");

  // TString input_name_merged = Form("./slopes/slug%d_dit_slope_merged_cycle_5coils.root",slug_id);
  // TFile *merged_input = TFile::Open(input_name_merged);
  // TDirectory *input_dir_merged = merged_input->GetDirectory("graph");

  // ++++++++++ Plots
  cout << " Making Plots " << endl;
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->cd();
  c2->SetRightMargin(0.015);
  Int_t plot_counter=0;
  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){
      c2->Clear("D");
      TLegend leg(0.985,0.90,0.75,0.65);
      Bool_t kSlugDone=kFALSE;
      Bool_t kRunDone=kFALSE;
      Bool_t k5x5SlugDone=kFALSE;
      TMultiGraph *fmg = new TMultiGraph();
      TString title = Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data());
      fmg->SetName(title);
      for(int isplit=0;isplit<nSplits;isplit++){
	TGraph *g_avg = GraphAverageSlope(fAveragedSlopeByRange[idet*nMon+imon][isplit],
					  fSplitXcord[isplit],0.5);
	if(g_avg==NULL) continue;
	g_avg->SetLineColor(kRed);
	g_avg->SetLineStyle(7);
	fmg->Add(g_avg,"l");
	if(!kSlugDone){
	  leg.AddEntry(g_avg,"over constraint slug avg","l");
	  kSlugDone=kTRUE;
	}
      }
      
      for(int irun=0;irun<fRunXcord[idet].size();irun++){
	if(fRunXcord[idet][irun].size()==0)
	  continue;
      	TGraph *g_avg = GraphAverageSlope(fAveragedSlopeByRun[idet*nMon+imon][irun],
      					  fRunXcord[idet][irun],0.5);
	if(g_avg==NULL) continue;
      	g_avg->SetLineColor(kRed);
      	fmg->Add(g_avg,"l");
	if(!kRunDone){
	  leg.AddEntry(g_avg,"over constraint run avg","l");
	  kRunDone=kTRUE;
	}
      }
      TMultiGraph* cycle_mg = (TMultiGraph*)input_dir->Get(title);
      fmg->Add(cycle_mg);
      TIter next(cycle_mg->GetListOfGraphs());
      TGraph *gint;
      
      while ( (gint=(TGraph*)next()) ){
      	if((gint->GetMarkerStyle())==47 && gint->GetMarkerColor()==kBlue)
	  leg.AddEntry(gint,"5x5 cyclewise","p");
	if((gint->GetLineStyle())==1 && gint->GetLineColor()==kBlue){
	  if(!k5x5SlugDone){
	    leg.AddEntry(gint,"5x5 slug avg.","l");
	    k5x5SlugDone=kTRUE;
	  }	  
	}
      }
      
      // TMultiGraph* merged_mg = (TMultiGraph*)input_dir_merged->Get(title);
      // TIter next_merged(merged_mg->GetListOfGraphs());
      // Bool_t k5coilSlugMergedDone=kFALSE;
      // while ( (gint=(TGraph*)next_merged()) ){
      // 	if((gint->GetLineStyle())==7 && gint->GetLineColor()==kRed){
      // 	  gint->SetLineStyle(1);
      // 	  gint->SetLineColor(kBlue);
      // 	  fmg->Add(gint);
      // 	  if(!k5coilSlugMergedDone){
      // 	    leg.AddEntry(gint,"5x5 slug merged","l");
      // 	    k5coilSlugMergedDone=kTRUE;
      // 	  }
      // 	}
      // }
      
      fmg->Draw("A");
      fmg->SetTitle(title);
      fmg->GetYaxis()->SetTitle("ppm/um");
      TH1F *hfmg = fmg->GetHistogram();
      hfmg->GetXaxis()->Set(nCycles,-0.5,nCycles-0.5);
      for(int i=0;i<nCycles;i++)
	hfmg->GetXaxis()->SetBinLabel(i+1,Form("%d",fCycleNumber[i]));
      hfmg->SetMarkerColor(kWhite);
      hfmg->Draw("same p");
 
      Double_t ymin = fmg->GetYaxis()->GetXmin();
      Double_t ymax = fmg->GetYaxis()->GetXmax();
      fmg->GetYaxis()->SetRangeUser(ymin, ymax+0.5*(ymax-ymin));
      ymax = ymax+0.5*(ymax-ymin);
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
      leg.Draw("same");

      graph_dir->cd();
      fmg->Write();
      canvas_dir->cd();
      c2->SetName(title);
      c2->Write();
      c2->SaveAs(Form("./plots/slug%d_dit_slope_buff_%003d.pdf",slug_id,plot_counter++));
    }
  }
  gSystem->Exec(Form("pdfunite `ls ./plots/slug%d_dit_slope_buff_*.pdf` ./plots/slug%d_dit_slope.pdf",
		     slug_id,slug_id));
  gSystem->Exec(Form("rm -f ./plots/slug%d_dit_slope_buff_*.pdf ",slug_id));

  canvas_dir->Write();
  graph_dir->Write();
  avg_output->Close();
  cout << " Writing " << rootfile_name << endl;

  //++++++++++++++
  if(kMatrixOutput){
    for(int isplit=0;isplit<nSplits;isplit++){
      TString range_tag;
      int low = range_low[isplit];
      int up = range_up[isplit];
      if(low==up)
	range_tag = Form("%d",low);
      else
	range_tag = Form("%d-%d",low,up);

      TMatrixD slope_matrix(nDet,nMon);  
      for(int idet=0;idet<nDet;idet++)
	for(int imon=0;imon<nMon;imon++)
	  slope_matrix[idet][imon] = fAveragedSlopeByRange[idet*nMon+imon][isplit];
  
      TString out_filename = Form("./matrices/prex_ovcn_slope_matrix.%s.root",
				  range_tag.Data());

      TFile *matrix_output = TFile::Open(out_filename,"RECREATE");
      matrix_output->WriteObject(&slope_matrix,"slope_matrix");
      matrix_output->WriteObject(&det_array,"dv_array");
      matrix_output->WriteObject(&mon_array,"iv_array");
      matrix_output->Close();
    }
  }
}

