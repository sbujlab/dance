#include "../src/TaAccumulator.cc"
#include "utilities.cc"
#include "plot_util.cc"
void AverageSensitivity(Int_t slug_id, Bool_t kMatrixOutput=kFALSE);
void AverageSensitivity(){
  for(int i=1;i<=94;i++)
    AverageSensitivity(i,kFALSE);
}
void AverageSensitivity(Int_t slug_id, Bool_t kMatrixOutput){
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

  vector<TString> device_array;
  vector<TString> prex_set1={"bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm12X","bpm8X"};
  vector<TString> prex_set2={"bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm11X12X","bpm11X","bpm12X","bpm1X","bpm16X"};
  vector<TString> crex_set={"bpm1X","bpm4aY","bpm4eX","bpm4eY","bpm12X","bpm11X","bpm4aX"};
  if(slug_id<=3)
    device_array = prex_set1;
  else if(slug_id<=94 || slug_id==501)
    device_array = prex_set2;
  else
    device_array = crex_set;
  
  vector<TString> det_array={"usl","usr","dsl","dsr"};
  vector<TString> at_array={"atl1","atl2","atr1","atr2"};
  vector<TString> sam_array={"sam1","sam2","sam3","sam4",
			     "sam5","sam6","sam7","sam8"};

  device_array.insert(device_array.end(),det_array.begin(),det_array.end());  
  if(slug_id>=26)
    device_array.insert(device_array.end(),at_array.begin(),at_array.end());
  if(slug_id>=100)
    device_array.insert(device_array.end(),sam_array.begin(),sam_array.end());

  vector<TString> coil_array;

  Int_t coil_index[]={1,3,4,6,7,2,5};
  Int_t nCoil = sizeof(coil_index)/sizeof(*coil_index); 
  for(int icoil=0;icoil<nCoil;icoil++)
    coil_array.push_back(Form("bmod_trim%d",coil_index[icoil]));

  Int_t nDev = device_array.size();

  TChain *sens_tree = new TChain("sens");
  Int_t nrun = fRunList.size();
  for(int i=0;i<nrun;i++){
    TString filename=Form("./dit-coeffs/prexPrompt_ditcoeffs_%d.root",
			  fRunList[i]);
    sens_tree->Add(filename);
  }
  Int_t nCycles = sens_tree ->GetEntries();
  if(nCycles==0)
    return;
  Double_t cycNumber;
  sens_tree->SetBranchAddress("cycID",&cycNumber);
  Double_t runNumber;
  sens_tree->SetBranchAddress("run",&runNumber);

  vector<Double_t> sens_val(nDev*nCoil);
  vector<Double_t> sens_err(nDev*nCoil);

  TString branch_name;
  for(int icoil=0;icoil<nCoil;icoil++){
    for(int idev=0;idev<nDev;idev++){
      branch_name = Form("%s_coil%d_err",
			 device_array[idev].Data(),coil_index[icoil]);
      sens_tree->SetBranchAddress(branch_name, &sens_err[idev*nCoil+icoil]);
      branch_name = Form("%s_coil%d",
			 device_array[idev].Data(),coil_index[icoil]);
      sens_tree->SetBranchAddress(branch_name, &sens_val[idev*nCoil+icoil]);
    }
  }
  vector<Int_t> fCycleNumber;
  vector<TaAccumulator> fAccumulator(nDev*nCoil);
  vector< vector<TaAccumulator> > fAccumulatorArray(nSplits,fAccumulator);
  vector< vector<Double_t> > fSplitXcord(nSplits); //[isplit][cycle]
  vector< vector<Double_t> > fSensValue(nDev*nCoil); //[idev*nMon+imon][cycle]
  vector< vector<Double_t> > fSensError(nDev*nCoil); //[idev*nMon+imon][cycle] 
  vector< vector<Double_t> > fSensXcord(nDev*nCoil); //[idev*nMon+imon][cycle]
  for(int ievt=0;ievt<nCycles;ievt++){
    sens_tree->GetEntry(ievt);
    fCycleNumber.push_back(cycNumber);
    Int_t split_id = fSplitMap[runNumber];
    fSplitXcord[split_id].push_back(cycNumber);
    if(IsBadCycle(fBadCycleMap,cycNumber))
      continue;
    Int_t arm_flag = fArmMap[runNumber];
      for(int idev=0;idev<nDev;idev++){

	if(arm_flag==1 && device_array[idev].Contains("l"))
	  continue;
	if(arm_flag==2 && device_array[idev].Contains("r"))
	  continue;
	
	for(int icoil=0;icoil<nCoil;icoil++){
	if(sens_err[idev*nCoil+icoil]>0){
	  fSensValue[idev*nCoil+icoil].push_back(sens_val[idev*nCoil+icoil]);
	  fSensError[idev*nCoil+icoil].push_back(sens_err[idev*nCoil+icoil]);
	  fSensXcord[idev*nCoil+icoil].push_back(cycNumber);
	  fAccumulatorArray[split_id][idev*nCoil+icoil].Update(sens_val[idev*nCoil+icoil]);
	}
      }// end of coil loop
    } // end of device loop
  } // end of ievt loop

  // ++++++++++
  TCanvas *c1  = new TCanvas("c1","c1",1400,700);
  c1->SetGridx();
  c1->SetGridy();
  c1->Divide(4,2);
  c1->Print(Form("./plots/slug%d_dit_sens.pdf[",slug_id),"pdf");
  for(int idev=0;idev<nDev;idev++){
    for(int icoil=0;icoil<nCoil;icoil++){
      c1->cd(icoil+1);
      double_t scale =1e6;
      if(device_array[idev].Contains("bpm"))
	scale =1e3;
      TMultiGraph *mg = new TMultiGraph();
      TGraphErrors *g_cyclewise = GraphErrorsVector(fSensValue[idev*nCoil+icoil],fSensError[idev*nCoil+icoil],fSensXcord[idev*nCoil+icoil],scale);
      mg->Add(g_cyclewise,"P");
      for(int isplit=0;isplit<nSplits;isplit++){
	TGraph *g_avg = GraphAverageSlope(fAccumulatorArray[isplit][idev*nCoil+icoil].GetMean1(),fSplitXcord[isplit],0,scale);
	if(g_avg==NULL) continue;
	g_avg->SetLineColor(kRed);
	mg->Add(g_avg,"l");
      }
      mg->Draw("A");
      TString title = device_array[idev]+Form("_coil%d",coil_index[icoil]);
      if(scale==1e6)
	title+="(ppm/count)";
      else if(scale==1e3)
	title+="(um/count)";
      mg->SetTitle(title+";cycle ID;sensitivity");

    }
    c1->Print(Form("./plots/slug%d_dit_sens.pdf",slug_id));
  }
  c1->Print(Form("./plots/slug%d_dit_sens.pdf]",slug_id),"pdf");

  // ++++++++++
  if(kMatrixOutput){
    for(int isplit=0;isplit<nSplits;isplit++){
      TMatrixD sensMatrix(nDev,nCoil);

      for(int idev=0;idev<nDev;idev++)
	for(int icoil=0;icoil<nCoil;icoil++)
	  sensMatrix[idev][icoil]=fAccumulatorArray[isplit][idev*nCoil+icoil].GetMean1();

      TString range_tag;
      int low = range_low[isplit];
      int up = range_up[isplit];
      if(low==up)
	range_tag = Form("%d",low);
      else
	range_tag = Form("%d-%d",low,up);

      TString output_filename = Form("./matrices/prex_sens_matrix.%s.root",range_tag.Data());
      cout << "Writing output: " << output_filename << endl;
      TFile *matrix_output = TFile::Open(output_filename,"RECREATE");
      matrix_output->WriteObject(&sensMatrix,"sens_matrix");
      matrix_output->WriteObject(&device_array,"dv_array");
      matrix_output->WriteObject(&coil_array,"coil_array");
      matrix_output->Close();
    }
  }
}

