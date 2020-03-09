#include "utilities.cc"
void ProbeCycle_prex(Int_t slug_number);
void ProbeCycle_prex(){
  for(int i=1;i<=94;i++)
    ProbeCycle_prex(i);
}
void ProbeCycle_prex(Int_t slug_number){
  vector<Int_t> fRunList = LoadRunListBySlug(slug_number);

  Bool_t kCrex = kFALSE;
  if(slug_number>=100 && slug_number<500)
    kCrex = kTRUE;
  map< Int_t, vector<Int_t> > fblmap = LoadBadCycleList(kCrex);
  gStyle->SetOptStat(0);
  TChain *sens = new TChain("sens");
  Int_t nrun = fRunList.size();
  for(int i=0;i<nrun;i++){
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.root",
			  fRunList[i]);
    sens->Add(filename);
  }

  vector<Int_t> coil_index;
  vector<Int_t> coil_set1={1,3,4,6,7,2,5};
  vector<Int_t> coil_set2={3,5,4,6,7,1,2};
  if(slug_number<=76)
    coil_index= coil_set1;
  else if(slug_number<=94)
    coil_index=coil_set2;
  else
    coil_index= coil_set1;

  const Int_t ncoil = coil_index.size();
  TString dev_array[]={"usl","usr"};
  const Int_t ndev = sizeof(dev_array)/sizeof(*dev_array);

  Double_t val_array[ndev][ncoil];
  Double_t err_array[ndev][ncoil];
  Double_t cycID;
  Double_t run;
  sens->SetBranchAddress("cycID",&cycID);
  sens->SetBranchAddress("run",&run);
  for(int icoil=0;icoil<ncoil;icoil++){
    for(int idev=0;idev<ndev;idev++){
      TString branch_name = Form("%s_coil%d_err",
				 dev_array[idev].Data(),coil_index[icoil]);
				 
      sens->SetBranchAddress(branch_name, &err_array[idev][icoil]);

      branch_name = Form("%s_coil%d",
			 dev_array[idev].Data(),coil_index[icoil]);

      sens->SetBranchAddress(branch_name, &val_array[idev][icoil]);
    }
  }
  
  Int_t nevt = sens->GetEntries();
  Bool_t kMatch;
  Int_t ncycles=0;
  vector<Int_t> fCycleID;
  vector< vector<Int_t> > fCoilFlags;
  vector<TText*> fTextRun;
  vector<TLine*> fLine;
  Int_t last_run = 0;
  for(int ievt=0;ievt<nevt;ievt++){
    sens->GetEntry(ievt);
    if(find(fRunList.begin(),fRunList.end(),(Int_t)run)==fRunList.end())
      kMatch = kFALSE;
    else
      kMatch = kTRUE;
    vector<Int_t> fCoilBuff;
    if(kMatch){
      ncycles++;
      if(run!=last_run){
	last_run = run;
	TText *t1 = new TText(ncycles-1,ncoil,Form("%d",last_run));
	TLine *line1 = new TLine(ncycles-1,0,ncycles-1,ncoil);
	fTextRun.push_back(t1);
	fLine.push_back(line1);
      }
      
      vector<Int_t> fCoilBuff;
      for(int icoil=0;icoil<ncoil;icoil++){
	if(err_array[0][icoil]>0 || err_array[1][icoil]>0)
	  fCoilBuff.push_back(1);
	else
	  fCoilBuff.push_back(0);
      }
      fCycleID.push_back(cycID);
      fCoilFlags.push_back(fCoilBuff);
      
    }// end of if kMatch
  } // end of event loop
  
  TH2D *h2d = new TH2D("h2d","",
		       ncycles,0,ncycles,ncoil,0,ncoil);
  for(int icyc=0;icyc<ncycles;icyc++){
    vector<Int_t> fCoilData = fCoilFlags[icyc];
    for(int icoil=0;icoil<ncoil;icoil++){
      Int_t weight = (icoil>4 ? 10:20);
      if(IsBadCycle(fblmap,fCycleID[icyc],coil_index[icoil]))
	weight=1;
      h2d->Fill(Form("%d",fCycleID[icyc]),
		Form("%d",coil_index[icoil]),
		weight*fCoilData[icoil]);
    }
  }
  h2d->GetYaxis()->SetLabelSize(0.07);
  h2d->GetXaxis()->SetLabelSize(0.04);
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetRightMargin(0.015);
  h2d->Draw("COL");
  
  vector<TText*>::iterator it_text = fTextRun.begin();
  while(it_text!=fTextRun.end()){
    (*it_text)->SetTextAngle(27);
    (*it_text)->SetTextSize(0.03);
    (*it_text)->Draw("same");
    it_text++;
  }
  
  vector<TLine*>::iterator it_line = fLine.begin();
  while(it_line!=fLine.end()){
    (*it_line)->SetLineColor(kRed);
    (*it_line)->SetLineWidth(3);
    (*it_line)->Draw("same");
    it_line++;
  }
  c1->SaveAs(Form("./plots/slug%d_dit_cyc.pdf",slug_number));
}

