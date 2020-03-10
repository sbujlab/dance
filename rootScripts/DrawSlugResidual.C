#include "plot_util.cc"
void DrawSlugResidual(Bool_t kCyclewise){

  
  TMultiGraph *mg_usl = new TMultiGraph();
  TMultiGraph *mg_usr = new TMultiGraph();
  
  vector<TString> graph_legend;
  vector<TString> file_tag;

  vector<TString> graph_legend_cyclewise={"5x5 cyclewise slope","overconstraint cyclewise slope"};
  vector<TString> file_tag_cyclewise={"5x5_by_cycle","ovcn_by_cycle"};

  vector<TString> graph_legend_avg={"5x5 slug averaged slope",
					  "overconstraint run avg. slope",
					  "overconstraint slug avg. slope"};
  vector<TString> file_tag_avg={"slug_5x5_avg",
				"run_ovcn","slug_ovcn"};

  if(kCyclewise){
    graph_legend=graph_legend_cyclewise;
    file_tag=file_tag_cyclewise;
  }else{
    graph_legend=graph_legend_avg;
    file_tag=file_tag_avg;
  }

  TLegend *leg = new TLegend(0.95,0.95,0.65,0.65);
  Int_t ntype = file_tag.size();
  for(int itype=0;itype<ntype;itype++){
    cout << file_tag[itype]  << endl;
    vector<Double_t> fres_usl;
    vector<Double_t> fx_usl;
    vector<Double_t> fres_usr;
    vector<Double_t> fx_usr;
    for(int i=1;i<=94;i++){
      TFile *input = TFile::Open(Form("./residuals/slug%d_%s.root",i,file_tag[itype].Data()));
      TTree *res = (TTree*)input->Get("res");
      if(res==NULL)
	continue;
      int npt;
      double res_sq_sum;
      TString redundant_cut="";
      // if(itype==0){
      // 	if(i<=76)
      // 	  redundant_cut = "&&(coil==2 || coil==5)";
      // 	else if(i<=94)
      // 	  redundant_cut = "&&(coil==1 || coil==2)";
      // }
      npt=res->Draw("usl_sq","usl_sq>0"+redundant_cut,"goff");
      double *y_ptr = res->GetV1();
      res_sq_sum=0;
      for(int ipt=0;ipt<npt;ipt++){
	res_sq_sum+=y_ptr[ipt];
      }
      if(npt!=0){
	fres_usl.push_back( sqrt(res_sq_sum/npt) );
	fx_usl.push_back(i);
      }else
	cout << "slug " << i << "has zero pt" << endl;

      npt=res->Draw("usr_sq","usr_sq>0"+redundant_cut,"goff");
      y_ptr = res->GetV1();
      res_sq_sum=0;
      for(int ipt=0;ipt<npt;ipt++){
	res_sq_sum+=y_ptr[ipt];
      }
      if(npt!=0){
	fres_usr.push_back( sqrt(res_sq_sum/npt) );
	fx_usr.push_back(i);
      }
      input->Close();
    }
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
  }
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->cd();
  mg_usl->Draw("A");
  mg_usl->SetTitle("usl residual sensitivity;Slug; ppm/count");
  double ymax= mg_usl->GetYaxis()->GetXmax();
  double ymin=mg_usl->GetYaxis()->GetXmin();
  mg_usl->GetYaxis()->SetRangeUser(ymin,ymax+0.2*(ymax-ymin));
  leg->Draw("same");
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->cd();
  mg_usr->Draw("A");
  mg_usr->SetTitle("usr residual sensitivity;Slug; ppm/count");
  ymax= mg_usr->GetYaxis()->GetXmax();
  ymin=mg_usr->GetYaxis()->GetXmin();
  mg_usr->GetYaxis()->SetRangeUser(ymin,ymax+0.2*(ymax-ymin));


  leg->Draw("same");
}
