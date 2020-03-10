#include "utilities.cc"
#include "plot_util.cc"
void ResidualSensByCycle(Int_t slug_number,Bool_t kOverConstraint);

void ResidualSensByCycle(){
  for(int i=1;i<=94;i++)
    ResidualSensByCycle(i,kFALSE);
}

void ResidualSensByCycle(Int_t slug_number ,Bool_t kOverConstraint){

  map<Int_t,Int_t> fArmMap = LoadArmMapBySlug(slug_number);
  vector<Int_t> fRunList = LoadRunListBySlug(slug_number);
  map< Int_t, vector<Int_t> > fblmap = LoadBadCycleList();

  TString tree_name;
  if(slug_number<=76)
    tree_name ="dit_slope1";
  else if(slug_number<=94)
    tree_name ="dit_slope3";
  else
    tree_name ="dit_slope1";
  if(kOverConstraint)
    tree_name="dit_slope_lsq";
  
  TChain *slope_tree = new TChain(tree_name);
  Int_t nrun = fRunList.size();
  for(int i=0;i<nrun;i++){
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.root",
			  fRunList[i]);
    slope_tree->Add(filename);
  }

  TChain *sens = new TChain("sens");
  for(int i=0;i<nrun;i++){
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.root",
			  fRunList[i]);
    sens->Add(filename);
  }
  sens->AddFriend(slope_tree);

  vector<Int_t> coil_index;
  vector<Int_t> coil_set1={5,2,1,3,4,6,7};
  vector<Int_t> coil_set2={1,2,3,5,4,6,7};

  if(slug_number<=76)
    coil_index = coil_set1;
  else if(slug_number<=94)
    coil_index = coil_set2;
  else
    coil_index = coil_set1;

  vector<TString> mon_set1={"bpm4aX","bpm4eX","bpm4aY","bpm4eY","bpm11X12X"};  // run 3404 - 4980
  vector<TString> mon_set2={"bpm4aX","bpm4eX","bpm4aY","bpm4eY","bpm12X"}; // run 3130-3403
  vector<TString> mon_array;

  if(slug_number<=3)
    mon_array = mon_set2;
  else if(slug_number<=94)
    mon_array = mon_set1;
  else
    mon_array = mon_set2;

  vector<TString> det_array = {"usl","usr","dsl","dsr"};
  // "sam1","sam2","sam3","sam4",
  // "sam5","sam6","sam7","sam8"};
  vector<TString> at_array={"atl1","atl2","atr1","atr2"};
  if(slug_number>=26)
    det_array.insert(det_array.end(),at_array.begin(),at_array.end());

  const Int_t ncoil = coil_index.size();
  const Int_t nmon = mon_array.size();
  const Int_t ndet = det_array.size();
  vector<Double_t> fdummy_det(ndet,0.0);
  vector<Double_t> fdummy_mon(nmon,0.0);

  vector<vector<Double_t> > det_val(ncoil,fdummy_det);
  vector<vector<Double_t> > det_err(ncoil,fdummy_det);
  vector<vector<Double_t> > mon_val(ncoil,fdummy_mon);
  vector<vector<Double_t> > mon_err(ncoil,fdummy_mon);
  vector<vector<Double_t> > slope_val(ndet,fdummy_mon);
  vector<vector<Double_t> > slope_flag(ndet,fdummy_mon);
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

  for(int idet=0;idet<ndet;idet++){
    for(int imon=0;imon<nmon;imon++){
      TString ch_name = Form("%s_%s",det_array[idet].Data(),mon_array[imon].Data());
      TString flag_name = Form("%s_%s_flag",det_array[idet].Data(),mon_array[imon].Data());
      sens->SetBranchAddress(ch_name,&slope_val[idet][imon]);
      sens->SetBranchAddress(flag_name,&slope_flag[idet][imon]);
    }
  }

  TString pdf_label;
  if(kOverConstraint)
    pdf_label="ovcn_by_cycle";
  else
    pdf_label="5x5_by_cycle";

  TFile *output = TFile::Open(Form("./residuals/slug%d_%s.root",slug_number,pdf_label.Data()),"RECREATE");
  TTree *residual_tree = new TTree("res","Residual Sensitivity Tree");
  vector<Double_t> fRes_ptr(ndet);
  vector<Double_t> fResSq_ptr(ndet);
  for(int idet=0;idet<ndet;idet++){
    for(int icoil=0;icoil<ncoil;icoil++){
      residual_tree->Branch(det_array[idet],&fRes_ptr[idet]);
      residual_tree->Branch(det_array[idet]+"_sq",&fResSq_ptr[idet]);
    }
  }
  residual_tree->Branch("run",&run);
  residual_tree->Branch("cycID",&cycID);
  Int_t coil_id;
  residual_tree->Branch("coil",&coil_id);
  Int_t ncycles=0;
  vector< Double_t > fdmy_vec;
  vector< vector<Double_t > > fRes_cyc(ndet*ncoil,fdmy_vec);
  vector< vector<Double_t > > fResXcord_cyc(ndet*ncoil,fdmy_vec);
  vector< vector<Double_t > > fRes(ndet*ncoil,fdmy_vec);
  vector< vector<Double_t > > fResRaw(ndet*ncoil,fdmy_vec);
  vector< vector<Double_t > > fResXcord(ndet*ncoil,fdmy_vec);
  vector< Int_t > fCycleArray;
  TH1D hist1("hist1","",30,-1.5,1.5);
  vector<TH1D> hPull(ndet*ncoil,hist1);
  Int_t nevt = sens->GetEntries();
  Bool_t kMatch;
  vector<Int_t> fRunLabel;
  vector<Double_t> fRunLabelXcord;
  Int_t last_run = 0.0;
  for(int ievt=0;ievt<nevt;ievt++){
    sens->GetEntry(ievt);
    if(find(fRunList.begin(),fRunList.end(),(Int_t)run)==fRunList.end())
      kMatch = kFALSE;
    else
      kMatch = kTRUE;
    if(kMatch){
      ncycles++;
      fCycleArray.push_back(cycID);
      if(run!=last_run){
	fRunLabel.push_back(run);
	fRunLabelXcord.push_back(ncycles-1);
	last_run=run;
      }
      
      Int_t arm_flag = fArmMap[(Int_t)run];
      for(int icoil=0;icoil<ncoil;icoil++){
	coil_id = coil_index[icoil];
	if(IsGoodCoil(mon_err[icoil])
	   && !IsBadCycle(fblmap,cycID,coil_index[icoil] )){
	  for(int idet=0;idet<ndet;idet++){
	    if( (arm_flag==1 && det_array[idet].Contains("l")) ||
		(arm_flag==2 && det_array[idet].Contains("r")) ){
	      fRes_ptr[idet] = 0.0;
	      fResSq_ptr[idet] = -1;
	      continue;
	    }
	    fResRaw[idet*ncoil+icoil].push_back(det_val[icoil][idet]);
	    fResXcord[idet*ncoil+icoil].push_back(ncycles-1);
	    if(slope_flag[idet][0]>0){
	      double residual = compute_residual(det_val[icoil][idet],mon_val[icoil],slope_val[idet]);
	      fRes_ptr[idet] = residual;
	      fResSq_ptr[idet] = residual*residual;
	      fRes_cyc[idet*ncoil+icoil].push_back(residual);
 	      fResXcord_cyc[idet*ncoil+icoil].push_back(ncycles-1);
	      hPull[idet*ncoil+icoil].Fill(residual*1e6);
	    }
	  } // end of det loop
	} // end of if its good coil
	residual_tree->Fill();
      } // end of coil loop
    }// end of if kMatch
  } // end of event loop
  output->cd();
  residual_tree->Write();
  output->Close();
  
  gStyle->SetTitleSize(0.3);
  gStyle->SetLabelSize(0.07);
  gStyle->SetOptStat("emrou");
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.35);
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->cd();
  TPad *p1 = new TPad("p1","",0,0.5,0.7,1);
  p1->Draw();
  TPad *p2 = new TPad("p2","",0,0,0.7,0.5);
  p2->Draw();
  TPad *p3 = new TPad("p3","",0.7,0,1,1);
  p3->Draw();
  
  p1->SetRightMargin(0.015);
  p1->SetTopMargin(0.15);
  p1->SetBottomMargin(0);
  p1->SetGridx();
  p2->SetTopMargin(0);
  p2->SetRightMargin(0.015);
  p2->SetGridx();

  Int_t ntext = fRunLabelXcord.size();
  Int_t page_counter=0;
  for(int idet=0;idet<ndet;idet++){
    for(int icoil=0;icoil<ncoil;icoil++){
      TString title = Form("%s_coil%d",det_array[idet].Data(),coil_index[icoil]);
      TGraph *g_res_cyc = GraphVector(fRes_cyc[idet*ncoil+icoil],fResXcord_cyc[idet*ncoil+icoil],1e6);
      TGraph *g_raw = GraphVector(fResRaw[idet*ncoil+icoil],fResXcord[idet*ncoil+icoil],1e6);
      p1->cd();
      g_res_cyc->SetMarkerColor(kBlue);
      g_res_cyc->SetMarkerStyle(20);
      TMultiGraph *mg_res = new TMultiGraph();
      mg_res->Add(g_res_cyc,"lp");
      mg_res->Draw("A");
      mg_res->SetTitle("Residual Sensitivity in " + title);
      mg_res->GetYaxis()->SetTitle("ppm/count");
      mg_res->GetYaxis()->SetTitleSize(0.06);
      mg_res->GetYaxis()->SetTitleOffset(0.6);
      mg_res->GetYaxis()->SetLabelSize(0.08);
      TH1F *hres = mg_res->GetHistogram();
      hres->GetXaxis()->Set(ncycles,-0.5,ncycles-0.5);
      mg_res->GetXaxis()->SetNdivisions(ncycles);

      Double_t ymin1 = mg_res->GetYaxis()->GetXmin();
      Double_t ymax1 = mg_res->GetYaxis()->GetXmax();
      for(int i=0;i<ntext;i++){
	TText *t1 = new TText(fRunLabelXcord[i],ymax1,Form("%d",fRunLabel[i]));
	TLine *line1 = new TLine(fRunLabelXcord[i],ymin1,
				 fRunLabelXcord[i],ymax1);
	line1->SetLineColor(kBlack);
	line1->SetLineWidth(1);
	line1->SetLineStyle(7);
	t1->SetTextAngle(27);
	t1->SetTextSize(0.05);
	line1->Draw("same");
	t1->Draw("same");
      }

      p2->cd();
      g_raw->SetMarkerStyle(20);
      g_raw->SetMarkerColor(kBlack);
      g_raw->Draw("ALP");
      g_raw->SetTitle(" ");
      g_raw->GetYaxis()->SetTitle("Raw Sensitivity (ppm/count)");
      g_raw->GetYaxis()->SetTitleSize(0.06);
      g_raw->GetYaxis()->SetTitleOffset(0.8);
      g_raw->GetYaxis()->SetLabelSize(0.08);
      TH1F *hraw = g_raw->GetHistogram();
      hraw->GetXaxis()->Set(ncycles,-0.5,ncycles-0.5);
      for(int i=0;i<ncycles;i++)
	hraw->GetXaxis()->SetBinLabel(i+1,Form("%d",fCycleArray[i]));
      hraw->SetMarkerColor(kWhite);
      hraw->Draw("same p");

      Double_t ymin = g_raw->GetYaxis()->GetXmin();
      Double_t ymax = g_raw->GetYaxis()->GetXmax();
      for(int i=0;i<ntext;i++){
	TText *t1 = new TText(fRunLabelXcord[i]-0.5,ymax,Form("%d",fRunLabel[i]));
	TLine *line1 = new TLine(fRunLabelXcord[i],ymin,
				 fRunLabelXcord[i],ymax);
	line1->SetLineColor(kBlack);
	line1->SetLineWidth(1);
	line1->SetLineStyle(7);
	t1->SetTextAngle(27);
	t1->SetTextSize(0.03);
	line1->Draw("same");
	t1->Draw("same");
      }

      p3->cd();
      hPull[idet*ncoil+icoil].Draw();
      hPull[idet*ncoil+icoil].GetXaxis()->SetTitle("residual (ppm/count)");

      c2->SaveAs(Form("./plots/slug%d_dit_res_buff_%003d.pdf",slug_number,page_counter++));
    } // end of coil loop
  } // end of det loop
  
  
  gSystem->Exec(Form("pdfunite $(ls ./plots/slug%d_dit_res_buff_*.pdf) ./plots/slug%d_dit_res_%s.pdf",slug_number,slug_number,pdf_label.Data()));

  gSystem->Exec(Form("rm -f ./plots/slug%d_dit_res_buff_*.pdf ",slug_number));

}
