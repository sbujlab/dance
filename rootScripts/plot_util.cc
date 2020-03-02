TGraph* GraphAverageSlope(Double_t  fyval,
			  vector<Double_t> fxval, Double_t ext=0.0){
  Double_t scale = 1e3; // ppm/um
  const Int_t npt = fxval.size()+2;
  
  Double_t xarr[npt];
  Double_t yarr[npt];
  for(int i=0;i<npt;i++){
    if(i==0)
      xarr[i] = fxval[0]-ext;
    else if(i==npt-1)
      xarr[i] = fxval[npt-3]+ext;
    else
      xarr[i] = fxval[i-1];
    
    yarr[i] = fyval*scale;
  }

  TGraph *g1 = new TGraph(npt,xarr,yarr);
  g1->SetMarkerStyle(47);
  g1->SetMarkerSize(2);
  g1->SetMarkerColor(kRed);
  g1->SetLineColor(kRed);
  g1->SetLineWidth(2);
  return g1;
}

TGraph* GraphVector(vector<Double_t> fyval,
		    vector<Double_t> fxval){
  Double_t scale = 1e3; // ppm/um
  const Int_t npt = fxval.size();
  Double_t xarr[npt];
  Double_t yarr[npt];
  for(int i=0;i<npt;i++){
    xarr[i] = fxval[i];
    yarr[i] = fyval[i]*scale;
  }

  TGraph *g1 = new TGraph(npt,xarr,yarr);
  g1->SetMarkerStyle(47);
  g1->SetMarkerSize(1.5);
  g1->SetMarkerColor(kBlue);
  return g1;
}
