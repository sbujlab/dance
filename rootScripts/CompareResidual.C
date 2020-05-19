#include "plot_util.cc"
void CompareResidual(){
  
  // vector<TString> graph_legend={"13467 cyclewise","35467 cyclewise","overconstraint cyclewise"};
  // vector<TString> file_tag={"dit_slope1","dit_slope3","dit_slope_lsq"};
  vector<TString> graph_legend={"5x5 slug avg","overconstraint slug avg","overconstraint run avg", "overconstraint cyclewise"};
  vector<TString> file_tag={"slug_5x5_merged","slug_ovcn","run_ovcn","dit_slope_lsq"};
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Print("compare_residual.pdf[");
  for(int icoil=1;icoil<=7;icoil++){
    TMultiGraph *mg_usl = new TMultiGraph();
    TMultiGraph *mg_usr = new TMultiGraph();
  
    TString coil_cut=Form("&& coil==%d",icoil);
    TLegend *leg = new TLegend(0.95,0.95,0.65,0.65);
    TLegend *leg2 = new TLegend(0.95,0.95,0.65,0.65);
    Int_t ntype = file_tag.size();
    for(int itype=0;itype<ntype;itype++){
      cout << file_tag[itype]  << endl;
      vector<Double_t> fres_usl;
      vector<Double_t> fx_usl;
      vector<Double_t> fres_usr;
      vector<Double_t> fx_usr;
      for(int i=1;i<=94;i++){
	TFile *input = TFile::Open(Form("./residuals/slug%d_%s.root",i,file_tag[itype].Data()));
	if(input==NULL)
	  continue;
	TTree *res = (TTree*)input->Get("res");
	if(res==NULL)
	  continue;
	int npt;
	double res_sq_sum;
	npt=res->Draw("usl_sq","usl_sq>0"+coil_cut,"goff");
	double *y_ptr = res->GetV1();
	res_sq_sum=0;
	for(int ipt=0;ipt<npt;ipt++){
	  res_sq_sum+=y_ptr[ipt];
	}
	if(npt!=0){
	  fres_usl.push_back( sqrt(res_sq_sum/npt) );
	  fx_usl.push_back(i);
	  cout << "slug " << i
	       << " res_sq_sum:" << sqrt(res_sq_sum)*1e6 << "\t"
	       << " npt:"<<npt << "\t"
	       << " res_norm:" << sqrt(res_sq_sum/npt)*1e6 << endl;
	}else
	  cout << "slug " << i << " has zero pt" << endl;

	npt=res->Draw("usr_sq","usr_sq>0"+coil_cut,"goff");
	y_ptr = res->GetV1();
	res_sq_sum=0;
	for(int ipt=0;ipt<npt;ipt++){
	  res_sq_sum+=y_ptr[ipt];
	}
	if(npt!=0){
	  fres_usr.push_back( sqrt(res_sq_sum/npt) );
	  fx_usr.push_back(i);
	}else
	  cout << "slug " << i << " has zero pt" << endl;

	input->Close();
      }// end of slug loop
      TGraph* gusl = GraphVector(fres_usl,fx_usl,1e6);
      gusl->SetMarkerColor(itype+1);
      gusl->SetLineColor(itype+1);
      gusl->SetLineWidth(2);
      gusl->SetMarkerStyle(20);
      gusl->SetMarkerSize(0.7);
      mg_usl->Add(gusl,"lp");
      TGraph* gusr = GraphVector(fres_usr,fx_usr,1e6);
      gusr->SetMarkerColor(itype+1);
      gusr->SetLineColor(itype+1);
      gusr->SetLineWidth(2);
      gusr->SetMarkerStyle(20);
      gusr->SetMarkerSize(0.7);
      mg_usr->Add(gusr,"lp");

      leg->AddEntry(gusl,graph_legend[itype],"lp");
      leg2->AddEntry(gusr,graph_legend[itype],"lp");
    } // end of types loop
    c1->cd();
    mg_usl->Draw("A");
    mg_usl->SetTitle(Form("usl residual Coil%d sensitivity;Slug; ppm/count",icoil));
    double ymax= mg_usl->GetYaxis()->GetXmax();
    double ymin=mg_usl->GetYaxis()->GetXmin();
    mg_usl->GetYaxis()->SetRangeUser(ymin,ymax+0.4*(ymax-ymin));
    leg->Draw("same");
    c1->Print("compare_residual.pdf");
    c1->Clear("D");
    mg_usr->Draw("A");
    mg_usr->SetTitle(Form("usr residual Coil%d sensitivity;Slug; ppm/count",icoil));
    ymax= mg_usr->GetYaxis()->GetXmax();
    ymin=mg_usr->GetYaxis()->GetXmin();
    mg_usr->GetYaxis()->SetRangeUser(ymin,ymax+0.4*(ymax-ymin));
    leg2->Draw("same");
    c1->Print("compare_residual.pdf");
  } // end of coil loop
  c1->Print("compare_residual.pdf]");
}
