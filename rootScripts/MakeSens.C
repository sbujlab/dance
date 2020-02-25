#include "MyFunctions.C"

void MakeSensMatrix(){
  Int_t coilID[]={1,3,4,6,7};
  vector<TString> mon_array{"bpm4aX","bpm4aY","bpm4eX","bpm4eY","bpm11X","bpm12X"};
  vector<TString> det_array{"usl","usr","dsl","dsr"};
  vector<TString> coil_array;
  Int_t nCoil = sizeof(coilID)/sizeof(*coilID);
  for(int icoil=0;icoil<nCoil;icoil++){
    coil_array.push_back(Form("bmod_trim%d",coilID[icoil]));
  }
  Int_t nMon = mon_array.size();
  Int_t nDet = det_array.size();
  TMatrixD detCMatrix(nDet,nCoil);
  TMatrixD monCMatrix(nMon,nCoil);
  TString filename="./dit-coeffs/prex_sens_matrix.root"
  TFile* input = TFile::Open(filename);
  TTree* sens_tree = (TTree*)input->Get("sens");
  Int_t nEntries = sens_tree ->GetEntries();
  Int_t cycNumber;
  sens_tree->SetBranchAddress("cycID",&cycNumber);
  sens_tree->GetEntry(0);
  Int_t cycNumber_min = cycNumber;
  Int_t cycNumber_max = cycNumber;
  for(int ievt=1;ievt<nEntries;ievt++){
    sens_tree->GetEntry(ievt);
    if(cycNumber>cycNumber_max)
      cycNumber_max=cycNumber;
    if(cycNumber<cycNumber_min)
      cycNumber_min=cycNumber;
  }
  std::map<Int_t, TCut> fCycleCutMap = GetBadCycleCut(cycNumber_min,
						      cycNumber_max);
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
      TCut good_cut = TCut(draw_cmd+"_err>0"+"&&"+user_cut);
      if(fCycleCutMap.find(coilID[icoil])!=fCycleCutMap.end())
	good_cut+=fCycleCutMap[coilID[icoil]];
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
				     !good_cut,"goff");
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
      c1->Print(Form("%s_sens.pdf[",key.Data()),"pdf");
    c1->Print(Form("%s_sens.pdf",key.Data()));
  }
  for(int idet=0;idet<nDet;idet++){
    for(int icoil=0;icoil<nCoil;icoil++){
      c1->cd(icoil+1);
      TMultiGraph *mge = new TMultiGraph();
      TString draw_cmd = det_array[idet]+"_coil"+coilID[icoil];
      TCut good_cut = TCut(draw_cmd+"_err>0"+"&&"+user_cut);
      if(fCycleCutMap.find(coilID[icoil])!=fCycleCutMap.end())
	good_cut+=fCycleCutMap[coilID[icoil]];
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
				     !good_cut,"goff");
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
    c1->Print(Form("%s_sens.pdf",key.Data()));
  }
  c1->Print(Form("%s_sens.pdf]",key.Data()),"pdf");

  TFile *matrix_output = TFile::Open(Form("./lagrange/sens/%s_sens_matrix.root",key.Data()),"RECREATE");
  matrix_output->WriteObject(&detCMatrix,"DetSens");
  matrix_output->WriteObject(&monCMatrix,"MonSens");
  matrix_output->WriteObject(&det_array,"det_array");
  matrix_output->WriteObject(&mon_array,"mon_array");
  matrix_output->WriteObject(&coil_array,"coil_array");
  matrix_output->Close();
  input->Close();    
}

