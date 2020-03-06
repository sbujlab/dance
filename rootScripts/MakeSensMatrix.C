#include "../src/TaAccumulator.cc"
#include "utilities.cc"
#include "plot_util.cc"

void MakeSensMatrix(Int_t slug_id){

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

  vector<TString> device_array{"bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm1X","bpm12X","bpm11X"};

  vector<TString> det_array{"usl","usr","dsl","dsr"};
  device_array.insert(device_array.end(),det_array.begin(),det_array.end());  
  vector<TString> at_array={"atl1","atl2","atr1","atr2"};
  if(slug_id>=26)
    device_array.insert(device_array.end(),at_array.begin(),at_array.end());

  vector<TString> coil_array;
  Int_t nCoil = 7; 
  for(int icoil=1;icoil<=nCoil;icoil++)
    coil_array.push_back(Form("bmod_trim%d",icoil));

  Int_t nDev = device_array.size();
  TMatrixD sensMatrix(nDev,nCoil);

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

  Int_t cycNumber;
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


  
  TCanvas *c1  = new TCanvas("c1","c1",1400,700);
  c1->SetGridx();
  c1->SetGridy();
  c1->Divide(3,2);
  gStyle->SetOptFit(1);
  TF1 *fpol0 = new TF1("fpol0","[0]",0,1e6);
  for(int imon=0;imon<nMon;imon++){
    for(int icoil=0;icoil<nCoil;icoil++){
      c1->cd(icoil+1);
      TMultiGraph *mge = new TMultiGraph();
      TString draw_cmd = mon_array[imon]+"_coil"+coilID[icoil];
      TString good_cut =draw_cmd+"_err>0";
      if(fCycleCutMap.find(coilID[icoil])!=fCycleCutMap.end())
	good_cut+= Form(" && (%s)",fCycleCutMap[coilID[icoil]].Data());
      Int_t npt =sens_tree->Draw(draw_cmd+":cycID:"+draw_cmd+"_err",
				 good_cut,"goff");
      
      Double_t *sens = sens_tree->GetV1();
      Double_t *cycID = sens_tree->GetV2();
      Double_t *sens_err = sens_tree->GetV3();
      TGraphErrors *ge = new TGraphErrors(npt,cycID,sens,0,sens_err);
      ge->SetMarkerStyle(20);
      mge->Add(ge);
      ge->Fit("fpol0","QR");
      monCMatrix(imon,icoil)=fpol0->GetParameter(0);

      Int_t npt_bad =sens_tree->Draw(draw_cmd+":cycID:"+draw_cmd+"_err",
				     "!("+good_cut+")","goff");
      Double_t *sens_bad = sens_tree->GetV1();
      Double_t *cycID_bad = sens_tree->GetV2();
      Double_t *sens_err_bad = sens_tree->GetV3();
      TGraphErrors *ge_bad = new TGraphErrors(npt_bad,cycID_bad,sens_bad,0,0);
      ge_bad->SetMarkerColor(kRed);
      ge_bad->SetMarkerStyle(47);
      ge_bad->SetMarkerSize(1.5);
      mge->Add(ge_bad);
      mge->Draw("AP");
      mge->SetTitle(Form("%s;cycleID",draw_cmd.Data()));
    }
    if(imon==0)
      c1->Print(Form("./pdf/run%s_sens.pdf[",key.Data()),"pdf");
    c1->Print(Form("./pdf/run%s_sens.pdf",key.Data()));
  }
  for(int idet=0;idet<nDet;idet++){
    for(int icoil=0;icoil<nCoil;icoil++){
      c1->cd(icoil+1);
      TMultiGraph *mge = new TMultiGraph();
      TString draw_cmd = det_array[idet]+"_coil"+coilID[icoil];
      TString good_cut = draw_cmd+"_err>0";
      if(fCycleCutMap.find(coilID[icoil])!=fCycleCutMap.end())
	good_cut+= Form(" && (%s)",fCycleCutMap[coilID[icoil]].Data());
      Int_t npt =sens_tree->Draw(draw_cmd+":cycID:"+draw_cmd+"_err",
				 good_cut,"goff");
      Double_t *sens = sens_tree->GetV1();
      Double_t *cycID = sens_tree->GetV2();
      Double_t *sens_err = sens_tree->GetV3();
      TGraphErrors *ge = new TGraphErrors(npt,cycID,sens,0,sens_err);
      ge->SetMarkerStyle(20);
      mge->Add(ge);
      ge->Fit("fpol0","QR");
      detCMatrix(idet,icoil)=fpol0->GetParameter(0);
      
      Int_t npt_bad =sens_tree->Draw(draw_cmd+":cycID:"+draw_cmd+"_err",
				     "!("+good_cut+")","goff");
      Double_t *sens_bad = sens_tree->GetV1();
      Double_t *cycID_bad = sens_tree->GetV2();
      Double_t *sens_err_bad = sens_tree->GetV3();
      TGraphErrors *ge_bad = new TGraphErrors(npt_bad,cycID_bad,sens_bad,0,0);
      ge_bad->SetMarkerColor(kRed);
      ge_bad->SetMarkerStyle(47);
      ge_bad->SetMarkerSize(1.5);
      mge->Add(ge_bad);
      mge->Draw("AP");
      mge->SetTitle(Form("%s;cycleID",draw_cmd.Data()));

    }
    c1->Print(Form("./pdf/run%s_sens.pdf",key.Data()));
  }
  c1->Print(Form("./pdf/run%s_sens.pdf]",key.Data()),"pdf");
  
  TFile *matrix_output = TFile::Open(Form("./dit-coeffs/prex_sens_matrix.%s.root",key.Data())
				     ,"RECREATE");

  matrix_output->WriteObject(&detCMatrix,"DetSens");
  matrix_output->WriteObject(&monCMatrix,"MonSens");
  matrix_output->WriteObject(&det_array,"det_array");
  matrix_output->WriteObject(&mon_array,"mon_array");
  matrix_output->WriteObject(&coil_array,"coil_array");
  matrix_output->Close();
  input->Close();    
}

