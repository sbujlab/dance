#include "../src/TaAccumulator.cc"
#include "utilities.cc"
#include "plot_util.cc"
void AverageSlope(Int_t slug_id,TString tree_name);
void AverageSlope(){
  TString tree_name;
  for(int i=1;i<=94;i++){
    if(i<=76)
      tree_name ="dit_slope1";
    else if(i<=94)
      tree_name ="dit_slope3";
    else
      tree_name ="dit_slope1";

    AverageSlope(i,tree_name);
  }
}
void AverageSlope(Int_t slug_id,TString tree_name){
  map<Int_t, vector<Int_t> > fBadCycleMap = LoadBadCycleList();
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
			       // "sam1","sam2","sam3","sam4",
			       // "sam5","sam6","sam7","sam8"};
  vector<TString> at_array={"atl1","atl2","atr1","atr2"};
  vector<TString> mon_array={"bpm4aX","bpm4eX","bpm4aY","bpm4eY"};
  vector<TString> mon1x_array={"bpm1X","bpm4eX","bpm4aY","bpm4eY","bpm12X"};
  if(slug_id<=3)
    mon_array.push_back("bpm12X");
  else if(slug_id<=94)
    mon_array.push_back("bpm11X12X");
  else
    mon_array.push_back("bpm12X");

  if(slug_id>=143)
    mon_array=mon1x_array;

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

  vector<TaAccumulator> fAccumulator(nMon*nDet);
  vector< vector<TaAccumulator> > fAccumulatorArray(nSplits,fAccumulator);
  
  vector<Int_t> fRunLabel;
  vector<Int_t> fCycleNumber;
  vector<Double_t> fRunLabelXcord;
  vector< vector<Double_t> > fSplitXcord(nSplits);
  vector< vector<Double_t> > fAveragedSlope(nDet*nMon);//[idet*nMon+imon][split]
  vector< vector<Double_t> > fSlopeXcord(nDet*nMon); //[idet*nMon+imon][cycle]
  vector< vector<Double_t> > fSlopeValues(nDet*nMon); //[idet*nMon+imon][cycle]
  Int_t last_runnumber=0;
  for(int ievt=0;ievt<nCycles;ievt++){
    slope_tree->GetEntry(ievt);
    fCycleNumber.push_back(cycNumber);
    if(last_runnumber!=runNumber){
      last_runnumber = runNumber;
      fRunLabel.push_back(runNumber);
      fRunLabelXcord.push_back(ievt); 
    }
    Int_t split_id = fSplitMap[runNumber];
    fSplitXcord[split_id].push_back(ievt);

    if(IsBadCycle(fBadCycleMap,cycNumber))
      continue;
    Int_t arm_flag = fArmMap[runNumber];
    for(int idet=0;idet<nDet;idet++){
      if(arm_flag==1 && det_array[idet].Contains("l"))
	continue;
      if(arm_flag==2 && det_array[idet].Contains("r"))
	continue;
      for(int imon=0;imon<nMon;imon++){
	if(kFlag[idet*nMon+imon]){
	  fAccumulatorArray[split_id][idet*nMon+imon].Update(fSlope[idet*nMon+imon]);
	  fSlopeValues[idet*nMon+imon].push_back(fSlope[idet*nMon+imon]);
	  fSlopeXcord[idet*nMon+imon].push_back(ievt);
	}
      }// end of mon loop
    } // end of det loop
  } // end of ievt loop
  
  // ++++++++++

  for(int isplit=0;isplit<nSplits;isplit++){
    ofstream mapfile;
    TString fullpath="./mapfiles/";
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
      mapfile << "[asym:@dit_asym_"<<det_array[idet]<<"]" << endl;
      mapfile << "asym:"<<det_array[idet] << ",1.0"<< endl;
      for(int imon=0;imon<nMon;imon++){
	TString iv_name = mon_array[imon];
	if(iv_name=="bpm11X12X")
	  iv_name = "diff_bpmE";
	mapfile << "diff:"<<iv_name << ","
		<< -fAccumulatorArray[isplit][idet*nMon+imon].GetMean1()<< endl;
	// Yeah, has a minus sign.
	fAveragedSlope[idet*nMon+imon].push_back(fAccumulatorArray[isplit][idet*nMon+imon].GetMean1());

      }
      mapfile << endl;
    }
    mapfile.close();
  }
  

  // ++++++++++
  TString rootfile_name = Form("./rootfiles/slug%d_dit_slope_cyclewise_average.root",slug_id);
  TFile* avg_output = TFile::Open(rootfile_name,"RECREATE");
  avg_output->cd();
  TTree* dit_tree = new TTree("dit","dit");
  vector<Double_t> fSlope_val(nDet*nMon);
  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){  
      TString branch_name =Form("%s_%s",
				det_array[idet].Data(),
				mon_array[imon].Data());
      dit_tree->Branch(branch_name,&fSlope_val[idet*nMon+imon]);
    }
  }
  Int_t fRun,fCounter;
  dit_tree->Branch("run",&fRun);
  dit_tree->Branch("range",&fCounter);
  
  for(int isplit=0;isplit<nSplits;isplit++){
    vector<Int_t> this_list = fRunListArray[isplit];
    Int_t nrun = this_list.size();
    fCounter = isplit;
    for(int i=0;i<nrun;i++){
      fRun = this_list[i];
      for(int idet=0;idet<nDet;idet++)
	for(int imon=0;imon<nMon;imon++)
	  fSlope_val[idet*nMon+imon] = fAveragedSlope[idet*nMon+imon][isplit];
      dit_tree->Fill();
    }
  }
  dit_tree->Write();
  
  TDirectory *graph_dir = avg_output->GetDirectory("/")->mkdir("graph");
  TDirectory *canvas_dir = avg_output->GetDirectory("/")->mkdir("canvas");
  // ++++++++++
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->cd();
  c2->SetRightMargin(0.015);
  Int_t plot_counter=0;
  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){
      TMultiGraph *fmg = new TMultiGraph();
      TString title = Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data());
      fmg->SetName(title);
      TGraph *g_cyclewise = GraphVector(fSlopeValues[idet*nMon+imon],fSlopeXcord[idet*nMon+imon]);

      for(int isplit=0;isplit<nSplits;isplit++){
	TGraph *g_avg = GraphAverageSlope(fAveragedSlope[idet*nMon+imon][isplit],
					  fSplitXcord[isplit],0.5);
	if(g_avg==NULL) continue;
	g_avg->SetLineColor(kBlue);
	fmg->Add(g_avg,"l");
      }
      fmg->Add(g_cyclewise,"p");
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
      graph_dir->cd();
      fmg->Write();
      canvas_dir->cd();
      c2->SetName(title);
      c2->Write();
      c2->SaveAs(Form("./plots/slug%d_dit_slope_buff_%003d.pdf",slug_id,plot_counter++));
    }
  }
  gSystem->Exec(Form("pdfunite `ls ./plots/slug%d_dit_slope_buff_*.pdf` ./plots/slug%d_dit_slope_cyclewise_average.pdf",
		     slug_id,slug_id));
  gSystem->Exec(Form("rm -f ./plots/slug%d_dit_slope_buff_*.pdf ",slug_id));

  canvas_dir->Write();
  graph_dir->Write();
  avg_output->Close();
  cout << " Writing " << rootfile_name << endl;
  // TMatrixD slope_matrix(nDet,nMon);  
  // for(int idet=0;idet<nDet;idet++)
  //   for(int imon=0;imon<nMon;imon++)
  //     slope_matrix[idet][imon] = fAccumulator[idet*nMon+imon].GetMean1();


  
  // TString out_filename = Form("./dit-coeffs/prex_%s_matrix.%d-%d.root",
  // 			      tree_name.Data(),
  // 			      fRunList[0],fRunList[nrun-1]);

  // TFile *matrix_output = TFile::Open(out_filename,"RECREATE");
  // matrix_output->WriteObject(&slope_matrix,"slope_matrix");
  // matrix_output->WriteObject(&det_array,"dv_array");
  // matrix_output->WriteObject(&mon_array,"iv_array");
  // matrix_output->Close();

}

