#include "TaLagrangian.hh"
#include "TBranch.h"
#include "TLeaf.h"
#include "TEventList.h"
#include "TStopwatch.h"
ClassImp(TaLagrangian);

TaLagrangian::TaLagrangian(TaInput *aInput){
  fInput = aInput;
  mul_tree=fInput->GetMulTree();
  run_number = fInput->GetRunNumber();
  seg_number = fInput->GetSegNumber();
  correct_tree = new TTree("redi","Correction Tree");
  mini_tree = new TTree("mini","Run Statistics Tree");
  mini_tree->SetMarkerStyle(20);
  slope_tree = new TTree("slopes","Slopes Tree");
  slope_tree->SetMarkerStyle(20);
}
TaLagrangian::~TaLagrangian(){}

Bool_t TaLagrangian::LoadConfig(TaConfig *fConfig){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  det_array = fConfig->GetDetArray();
  mon_array = fConfig->GetMonArray();
  nDet = det_array.size();
  nMon = mon_array.size(); 
  nCoil = (fConfig->GetCoilArray()).size();
  if(det_array.size()==0 ||
     mon_array.size()==0){
    cerr << " ** Error: Empty Detector or Monitor Array " << endl;
    cerr << " ** Error: Nothing I can do... Bye-bye! " << endl;
    return kFALSE;
  }
  TObjArray* branch_array = mul_tree->GetListOfBranches();
  for(int idet=0;idet<nDet;idet++){
    TObject *obj_ptr = branch_array->FindObject("asym_"+det_array[idet]);
    if(obj_ptr==NULL){
      cerr<< " ** Error : Branch " 
	  << "asym_"+det_array[idet] 
	  << " is not found " << endl;
      return kFALSE;
    }
  }

  for(int imon=0;imon<nMon;imon++){
    TObject *obj_ptr = branch_array->FindObject("diff_"+mon_array[imon]);
    if(obj_ptr==NULL){
      cerr<< " ** Error : Branch " 
	  << "diff_"+mon_array[imon] 
	  << " is not found " << endl;
      return kFALSE;
    }
  }
  cout << " -- Configuring Alias  " << endl;
  vector<pair<TString,ElementsVector> > aliasConfig = fConfig->GetAliasArray();
  vector<pair<TString,ElementsVector> >::iterator iter = aliasConfig.begin();
  while(iter!=aliasConfig.end()){
    TString alias_name = (*iter).first;
    Bool_t kMatched= kTRUE;
    ElementsVector::iterator itelem=((*iter).second).begin();
    AliasMap this_map;
    while(itelem!=((*iter).second).end()){
      TString elem_name = (*itelem).second;
      for(int idet=0;idet<nDet;idet++){
  	if(elem_name==det_array[idet]){
  	  pair<Double_t,Int_t> this_pair = make_pair((*itelem).first,idet);
  	  this_map.push_back(this_pair);
  	  kMatched =kTRUE;
  	  break;
  	}
  	kMatched=kFALSE;
      }
      if(!kMatched){
  	cout << " ** Warning: alias element "
  	     << elem_name 
  	     << " is not found in detector list " << endl;
  	break;
      }
      itelem++;
    }
    if(kMatched){
      cout << " -- Loading alias "
  	   << alias_name  << endl;
#ifdef DEBUG
      AliasMap::iterator it_map = this_map.begin();
      while(it_map!=this_map.end()){
	cout << " factor: " << (*it_map).first  << " "
	     << " det index: " << (*it_map).second << " ";
	it_map++;
      }
      cout << endl;
#endif
      alias_array.push_back(alias_name);
      alias_map_array.push_back(this_map);
    }
    iter++;
  }
  nCombo = alias_array.size();
  if(fConfig->GetAnalysisType()=="regression")
    nCoil=0;
  mini_size = fConfig->GetMiniSize();
  cout << " -- Mini-Run Size " << mini_size << endl;
  custom_cut = fConfig->GetCustomCut();
  output_path = fConfig->GetOutputPath();
  if(custom_cut!="")
    cout << " -- Using Mul Tree custom cut "
	 << custom_cut.GetTitle()
	 << endl;
  return kTRUE;
}

Bool_t TaLagrangian::ComputeCorrelation(){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  if(mul_tree==NULL){
    cerr << " ** Error : Mul Tree doesn't exisit ! "
	 << " Analysis aborted " << endl;
    return kFALSE;
  }
  TCut default_cut("ErrorFlag==0");
  TEventList *elist_mul = new TEventList("elist_mul");  
  nGoodPatterns = mul_tree->Draw(">>+elist_mul",
				 default_cut&&custom_cut,
				 "goff");
  Int_t ntotal = mul_tree->GetEntries();
  cout << " -- nGoodPatterns in Mul Tree : " << nGoodPatterns << endl;
  
  vector<AccVector> DetMonAccumulator(nDet); //[idet][imon]
  vector<AccVector> MonAccumulator(nMon); // [imon][jmon]
  AccVector DetAccumulator(nDet); // [idet]
  vector<Double_t> asym_det(nDet);
  vector<Double_t> diff_mon(nMon);
  for(int idet=0;idet<nDet;idet++){
    TLeaf* this_leaf=mul_tree->GetBranch("asym_"+det_array[idet])->GetLeaf("hw_sum");
    this_leaf->SetAddress(&asym_det[idet]);
  }
  for(int imon=0;imon<nMon;imon++){
    TLeaf* this_leaf=mul_tree->GetBranch("diff_"+mon_array[imon])->GetLeaf("hw_sum");
    this_leaf->SetAddress(&diff_mon[imon]);
  }
  Double_t kErrorFlag;
  mul_tree->SetBranchAddress("ErrorFlag",&kErrorFlag);

  Vec2D CovDetMon(nDet); //[idet][imon]
  Vec2D CovMonMon(nMon); //[imon][jmon] symmetric

  for(int idet=0;idet<nDet;idet++){
    DetMonAccumulator[idet].resize(nMon);
    CovDetMon[idet].resize(nMon);
  }
  for(int imon=0;imon<nMon;imon++){
    MonAccumulator[imon].resize(nMon);
    CovMonMon[imon].resize(nMon);
  }

  Int_t goodCounts=0;
  Int_t start = 0;
  Int_t end = 0;
  double ppm = 1e-6;
  double um = 1e-3;
  TStopwatch tsw;
  dataRawDet.resize(nDet);
  dataMon.resize(nMon);
  for(int ievt=0;ievt<ntotal;ievt++){
    mul_tree->GetEntry(ievt);
    for(int idet=0;idet<nDet;idet++)
      dataRawDet[idet].push_back(asym_det[idet]);
    for(int imon=0;imon<nMon;imon++)
      dataMon[imon].push_back(diff_mon[imon]);
    dataErrorFlag.push_back(kErrorFlag);
    if(!elist_mul->Contains(ievt)) {
      cutFlag.push_back(kFALSE);
      continue;
    }
    cutFlag.push_back(kTRUE);
    goodCounts++;
    for(int idet=0;idet<nDet;idet++){
      DetAccumulator[idet].Update(asym_det[idet]);
      for(int imon=0;imon<nMon;imon++){
	DetMonAccumulator[idet][imon].Update(asym_det[idet],diff_mon[imon]);
      }
    }
    for(int imon=0;imon<nMon;imon++){
      for(int jmon=imon;jmon<nMon;jmon++){ // just upper right 
	MonAccumulator[imon][jmon].Update(diff_mon[imon],diff_mon[jmon]);
      }
    }
    if(goodCounts==mini_size){
      end = ievt;
      mini_range.push_back(make_pair(start,end));
      cout << " -- Mini-run ends at event: " << end << endl;
      start = end+1;
      cout << " -- Next Mini-run starts at event: " << start << endl;
      goodCounts = 0;

      int nMini = mini_range.size();
      Bool_t is_last_minrun = kFALSE;
      if(nGoodPatterns-mini_size*nMini<mini_size){
	is_last_minrun = kTRUE;
	Int_t last_start = (mini_range.back()).first;
	mini_range.pop_back();
	mini_range.push_back(make_pair(last_start,ntotal-1));
	cout << " -- Meeting last mini-run, " << endl;
	cout << " -- the rest will be merged into this mini-run  "  << endl;
      }
      for(int idet=0;idet<nDet;idet++){
	for(int imon=0;imon<nMon;imon++){
	  CovDetMon[idet][imon] = DetMonAccumulator[idet][imon].GetM2()/ppm/um;
	  CovDetMon[idet][imon] = CovDetMon[idet][imon]/DetMonAccumulator[idet][imon].GetN();
	  if(!is_last_minrun)
	    DetMonAccumulator[idet][imon].Zero();
	}
      }
      for(int imon=0;imon<nMon;imon++){
	for(int jmon=imon;jmon<nMon;jmon++){
	  CovMonMon[imon][jmon] = MonAccumulator[imon][jmon].GetM2()/um/um;
	  CovMonMon[imon][jmon] = CovMonMon[imon][jmon]/MonAccumulator[imon][jmon].GetN();
	  if(!is_last_minrun)
	    MonAccumulator[imon][jmon].Zero();
	}
      }
      if(!is_last_minrun){
	miniCovDetMon.push_back(CovDetMon);
	miniCovMonMon.push_back(CovMonMon);
      }
    }
  }// End of Mul Tree Event Loop
  cout << " -- last mini-run ends at event: " << ntotal-1 << endl;
  if(nGoodPatterns<mini_size)
    mini_range.push_back(make_pair(start,ntotal-1));

  for(int idet=0;idet<nDet;idet++){
    for(int imon=0;imon<nMon;imon++){
      CovDetMon[idet][imon] = DetMonAccumulator[idet][imon].GetM2()/ppm/um;
      CovDetMon[idet][imon] = CovDetMon[idet][imon]/DetMonAccumulator[idet][imon].GetN();
    }
  }
  for(int imon=0;imon<nMon;imon++){
    for(int jmon=imon;jmon<nMon;jmon++){
      CovMonMon[imon][jmon] = MonAccumulator[imon][jmon].GetM2()/um/um;
      CovMonMon[imon][jmon] = CovMonMon[imon][jmon]/MonAccumulator[imon][jmon].GetN();
    }
  }
  miniCovDetMon.push_back(CovDetMon);
  miniCovMonMon.push_back(CovMonMon);

  cout << " -- Mul Event loop Done " << endl;
  tsw.Print();
  return kTRUE;
}

void TaLagrangian::LoadConstraint(TMatrixD detMatrix,TMatrixD monMatrix){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  kUseConstraint = kTRUE;
  if(detMatrix.GetNrows()!=nDet){
    cerr << " ** Warning :"
	 << " Dimension of Detector constraint matrix mismatched "  << endl
	 << "   Source matrix Row Upper bound:" <<detMatrix.GetNrows()	 
	 << endl
	 << " ** Warning : no constraints in use " << endl;
    detCmatrix.ResizeTo(0,0);
    nCoil=0;
    kUseConstraint=kFALSE;
  }
  else{
    detCmatrix.ResizeTo(nDet,nCoil);
    detCmatrix=detMatrix;
  }
  if(monMatrix.GetNrows()!=nMon){
    cerr << " ** Warning :"
	 << " Dimension of Monitor constraint matrix mismatched " << endl
	 << "   Source matrix Row Upper bound:" <<monMatrix.GetNrows()	 
	 << endl
	 << " ** Warning : no constraints in use " << endl;
    monCmatrix.ResizeTo(0,0);
    nCoil=0;
    kUseConstraint=kFALSE;
  }
  else{
    monCmatrix.ResizeTo(nMon,nCoil);
    monCmatrix=monMatrix;
  }
}

void TaLagrangian::ComputeSlopes(){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  if(detCmatrix.GetNrows()==0 || monCmatrix.GetNrows()==0)
    nCoil=0;
  TMatrixD lhsM_template(nMon+nCoil,nMon+nCoil);    
  TMatrixD rhsM_template(nMon+nCoil,nDet);    
  TMatrixD monCtrans(nCoil,nMon);
  TMatrixD detCtrans(nCoil,nDet);
  if(kUseConstraint){
    monCtrans.Transpose(monCmatrix);
    detCtrans.Transpose(detCmatrix);
    lhsM_template.SetSub(0,nMon,monCmatrix);
    lhsM_template.SetSub(nMon,0,monCtrans);
    rhsM_template.SetSub(nMon,0,detCtrans);
  }
  Int_t nMini = mini_range.size();
  for(int imini=0;imini<nMini;imini++){
    Vec2D CovMonMon = miniCovMonMon[imini];
    Vec2D CovDetMon = miniCovDetMon[imini];
    TMatrixD lhsM(lhsM_template);
    TMatrixD rhsM(rhsM_template);
    for(int irow=0;irow<nMon;irow++)
      for(int icol=irow;icol<nMon;icol++){
	lhsM[irow][icol] = CovMonMon[irow][icol];
	if(icol!=irow){
	  lhsM[icol][irow] = CovMonMon[irow][icol];
	}
      }
    for(int irow=0;irow<nMon;irow++)
      for(int icol=0;icol<nDet;icol++)
	rhsM[irow][icol] = CovDetMon[icol][irow];
    
    TMatrixD invlhsM = lhsM.Invert();
    TMatrixD solutionM = invlhsM*rhsM;
#ifdef NOISY	
    cout << " -- Slopes Matrix (nMon x nDet) " << endl;
    TMatrixD slopeM(nMon,nDet);
    slopeM=solutionM.GetSub(0,nMon-1,0,nDet-1);
    slopeM.Print();
#endif
    Vec2D DetMonSlopes(nDet);
    for(int idet=0;idet<nDet;idet++){
      DetMonSlopes[idet].resize(nMon);
      for(int imon=0; imon<nMon; imon++)
	DetMonSlopes[idet][imon] = solutionM[imon][idet];
    }
    miniDetMonSlopes.push_back(DetMonSlopes);
#ifdef NOISY
    if(kUseConstraint){
      cout << " -- Check Slope-Sensitivity Consistency " << endl;
      TMatrixD ZeroM = monCtrans*slopeM- detCtrans;
      ZeroM.Print();
    }
#endif 
  } // end of mini loop
}

void TaLagrangian::CorrectTree(){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  vector<Double_t> raw_det(nDet);
  vector<Double_t> cor_det(nDet);
  vector<Double_t> raw_combo(nCombo);
  vector<Double_t> cor_combo(nCombo);
  vector<Double_t> diff_mon(nMon);

  statMon.resize(nMon);
  statRawDet.resize(nDet);
  statCorDet.resize(nDet);
  statRawCombo.resize(nCombo);
  statCorCombo.resize(nCombo);

  AccVector accMon(nMon);
  AccVector accRawDet(nDet);
  AccVector accCorDet(nDet);
  AccVector accRawCombo(nCombo);
  AccVector accCorCombo(nCombo);

  accMon_avg.resize(nMon); 
  accRawDet_avg.resize(nDet);
  accCorDet_avg.resize(nDet);
  accRawCombo_avg.resize(nCombo);
  accCorCombo_avg.resize(nCombo);

  typedef struct {Double_t ppm,ppb,um,nm;} UNIT;
  UNIT parity_scale;
  parity_scale.ppm = 1e-6;
  parity_scale.ppb = 1e-9;
  parity_scale.um = 1e-3;
  parity_scale.nm = 1e-6;

  for(int imon=0;imon<nMon;imon++)
    mini_tree->Branch("diff_"+mon_array[imon],&statMon[imon],
		     "mean/D:error:rms:nsamples");

  for(int idet=0;idet<nDet;idet++){
    mini_tree->Branch("asym_"+det_array[idet],&statRawDet[idet],
		     "mean/D:error:rms:nsamples");
    mini_tree->Branch("redi_asym_"+det_array[idet],&statCorDet[idet],
		     "mean/D:error:rms:nsamples");
  }
  for(int icom=0;icom<nCombo;icom++){
    mini_tree->Branch("asym_"+alias_array[icom],&statRawCombo[icom],
		     "mean/D:error:rms:nsamples");
    mini_tree->Branch("redi_asym_"+alias_array[icom],&statCorCombo[icom],
		     "mean/D:error:rms:nsamples");
  }
  mini_tree->Branch("unit",&parity_scale,"ppm/D:ppb:um:nm");
  mini_tree->Branch("run",&run_number);
  mini_tree->Branch("seg",&seg_number);

  Vec2D ComboSlopes(nCombo);
  Vec2D DetMonSlopes(nDet);;

  for(int icom=0;icom<nCombo;icom++)
    ComboSlopes[icom].resize(nMon);

  for(int idet=0;idet<nDet;idet++)
    DetMonSlopes[idet].resize(nMon);

  for(int idet=0;idet<nDet;idet++)
    for(int imon=0;imon<nMon;imon++)
      slope_tree->Branch(det_array[idet]+"_"+mon_array[imon],&DetMonSlopes[idet][imon]);
  for(int icom=0;icom<nCombo;icom++)
    for(int imon=0;imon<nMon;imon++)
      slope_tree->Branch(alias_array[icom]+"_"+mon_array[imon],&ComboSlopes[icom][imon]);

  for(int imon=0;imon<nMon;imon++)
    correct_tree->Branch("diff_"+mon_array[imon],&diff_mon[imon]);
  for(int idet=0;idet<nDet;idet++){
    correct_tree->Branch("asym_"+det_array[idet],&raw_det[idet]);
    correct_tree->Branch("redi_asym_"+det_array[idet],&cor_det[idet]);
  }
  for(int icom=0;icom<nCombo;icom++){
    correct_tree->Branch("asym_"+alias_array[icom],&raw_combo[icom]);
    correct_tree->Branch("redi_asym_"+alias_array[icom],&cor_combo[icom]);
  }
  correct_tree->Branch("unit",&parity_scale,"ppm/D:ppb:um:nm");
  Bool_t kCutFlag;
  correct_tree->Branch("ok_cut",&kCutFlag);
  Double_t kErrorFlag;
  correct_tree->Branch("ErrorFlag",&kErrorFlag);
  Int_t mini_id;
  correct_tree->Branch("mini",&mini_id);

  mini_tree->Branch("mini",&mini_id);
  slope_tree->Branch("mini",&mini_id);

  Int_t nMini = mini_range.size();
  for(int imini=0;imini<nMini;imini++){
    mini_id = imini;
    int start_evt = mini_range[imini].first;
    int end_evt = mini_range[imini].second;
    cout << " -- Correction Mini-run: " << imini << endl;
    cout << " -- Event Number : " 
	 << start_evt
	 << " - "
	 << end_evt << endl;
    DetMonSlopes = miniDetMonSlopes[imini];
    for(int icom=0;icom<nCombo;icom++){
      AliasMap myMap = alias_map_array[icom];
      int nelem = myMap.size();
      for(int imon=0;imon<nMon;imon++){
	double this_slope =0.0;
	for(int ielem =0;ielem<nelem;ielem++){
	  int det_index = myMap[ielem].second;
	  this_slope += myMap[ielem].first * DetMonSlopes[det_index][imon];
	}
	ComboSlopes[icom][imon] = this_slope;
      }
    }
    miniComboSlopes.push_back(ComboSlopes);
    for(int ievt =start_evt ;ievt<=end_evt;ievt++){
      for(int idet=0;idet<nDet;idet++){
	raw_det[idet]=dataRawDet[idet][ievt];
	cor_det[idet]=dataRawDet[idet][ievt];
      }
      for(int icom=0;icom<nCombo;icom++){
	AliasMap myMap = alias_map_array[icom];
	int nelem = myMap.size();
	raw_combo[icom] = 0.0;
	for(int ielem =0;ielem<nelem;ielem++){
	  int det_index = myMap[ielem].second;
	  raw_combo[icom]+= myMap[ielem].first * raw_det[det_index];
	}
	cor_combo[icom] = raw_combo[icom];
      }
      for(int imon=0;imon<nMon;imon++)
	diff_mon[imon]=dataMon[imon][ievt];
      kErrorFlag=dataErrorFlag[ievt];
      kCutFlag=cutFlag[ievt];
      if(kCutFlag==kTRUE){
	for(int idet=0;idet<nDet;idet++)
	  for(int imon=0;imon<nMon;imon++)
	    cor_det[idet]-=DetMonSlopes[idet][imon]/1e3*diff_mon[imon]; 

	for(int icom=0;icom<nCombo;icom++)
	  for(int imon=0;imon<nMon;imon++)
	    cor_combo[icom]-=ComboSlopes[icom][imon]*diff_mon[imon]/1e3;

	for(int imon=0;imon<nMon;imon++)
	  accMon[imon].Update(diff_mon[imon]);
	for(int idet=0;idet<nDet;idet++){
	  accCorDet[idet].Update(cor_det[idet]);
	  accRawDet[idet].Update(raw_det[idet]);
	}
	for(int icom=0;icom<nCombo;icom++){
	  accCorCombo[icom].Update(cor_combo[icom]);
	  accRawCombo[icom].Update(raw_combo[icom]);
	}
      } // End of if kCutFlag==true
      correct_tree->Fill();
    } // end of event loop
    // Fill Slope Tree for each mini run;
    slope_tree->Fill();
    // Fill Run Stat Tree for each mini run;
    for(int imon=0;imon<nMon;imon++){
      UpdateRunStat(statMon[imon],accMon[imon]);
      accMon_avg[imon].Merge(accMon[imon]);
      accMon[imon].Zero();
    }
    for(int idet=0;idet<nDet;idet++){
      UpdateRunStat(statRawDet[idet],accRawDet[idet]);
      UpdateRunStat(statCorDet[idet],accCorDet[idet]);
      accRawDet_avg[idet].Merge(accRawDet[idet]);
      accCorDet_avg[idet].Merge(accCorDet[idet]);
      accRawDet[idet].Zero();
      accCorDet[idet].Zero();
    }
    for(int icom=0;icom<nCombo;icom++){
      UpdateRunStat(statRawCombo[icom],accRawCombo[icom]);
      UpdateRunStat(statCorCombo[icom],accCorCombo[icom]);
      accRawCombo_avg[icom].Merge(accRawCombo[icom]);
      accCorCombo_avg[icom].Merge(accCorCombo[icom]);
      accRawCombo[icom].Zero();
      accCorCombo[icom].Zero();
    }
    mini_tree->Fill(); 
  } // end of mini loop
}

void TaLagrangian::WriteSummary(){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  TString txt_filename = output_path+"lagrange_";
  txt_filename +=Form("%d_%03d.txt",run_number,seg_number);
  FILE *summary_txt = fopen(txt_filename.Data(),"w");
  if(summary_txt==NULL){
    cerr << " Error : Can not Open" 
	 << txt_filename
	 << endl
	 << " -- No Text Summary "
	 << endl;
    return;
  }else{
    cout << " -- Writing text summary to "
	 << txt_filename
	 << endl;
    double ppm = 1e-6;
    double um = 1e-3;

    vector<RUN_STATS> statMon_avg(nMon);
    vector<RUN_STATS> statRawDet_avg(nDet);
    vector<RUN_STATS> statCorDet_avg(nDet);
    vector<RUN_STATS> statRawCombo_avg(nCombo);
    vector<RUN_STATS> statCorCombo_avg(nCombo);

    for(int imon=0;imon<nMon;imon++){
      UpdateRunStat(statMon_avg[imon],accMon_avg[imon]);
    }
    for(int idet=0;idet<nDet;idet++){
      UpdateRunStat(statRawDet_avg[idet],accRawDet_avg[idet]);
      UpdateRunStat(statCorDet_avg[idet],accCorDet_avg[idet]);
    }
    for(int icom=0;icom<nCombo;icom++){
      UpdateRunStat(statRawCombo_avg[icom],accRawCombo_avg[icom]);
      UpdateRunStat(statCorCombo_avg[icom],accCorCombo_avg[icom]);
    }

    fprintf(summary_txt, " == Run Average Statistics == \n");
    fprintf(summary_txt, 
	    "%-12s\t%-12s\t%-12s\t%-20s\t%-15s\t%-15s\t%-12s\n",
	    "",
	    "Raw Mean+/-Err(ppm)","RMS(ppm)",
	    "Regressed Mean+/-Err(ppm)","RMS (ppm)", 
	    "Correction on Mean(ppm)","RMS of BPM Coherent Sum(ppm)");
    for(int idet=0;idet<nDet;idet++){
      double CoherentNoise = sqrt( pow(statRawDet_avg[idet].rms,2)-
				   pow(statCorDet_avg[idet].rms,2));
      fprintf(summary_txt, 
	      "%-12s\t%-3.2f+/-%-9.2f\t%-12.2f\t%-3.2f+/-%-17.2f\t%-17.2f\t%-15.2f\t%-12.2f\n",
	      det_array[idet].Data(), 
	      statRawDet_avg[idet].mean/ppm,statRawDet_avg[idet].error/ppm,statRawDet_avg[idet].rms/ppm, 
	      statCorDet_avg[idet].mean/ppm,statCorDet_avg[idet].error/ppm,statCorDet_avg[idet].rms/ppm, 
	      (statRawDet_avg[idet].mean-statCorDet_avg[idet].mean)/ppm,CoherentNoise/ppm);
    }

    for(int icom=0;icom<nCombo;icom++){
      double CoherentNoise = sqrt( pow(statRawCombo_avg[icom].rms,2)-
				   pow(statCorCombo_avg[icom].rms,2));
      fprintf(summary_txt, 
	      "%-12s\t%-3.2f+/-%-9.2f\t%-12.2f\t%-3.2f+/-%-17.2f\t%-17.2f\t%-15.2f\t%-12.2f\n",
	      alias_array[icom].Data(), 
	      statRawCombo_avg[icom].mean/ppm,statRawCombo_avg[icom].error/ppm,statRawCombo_avg[icom].rms/ppm, 
	      statCorCombo_avg[icom].mean/ppm,statCorCombo_avg[icom].error/ppm,statCorCombo_avg[icom].rms/ppm, 
	      (statRawCombo_avg[icom].mean-statCorCombo_avg[icom].mean)/ppm,CoherentNoise/ppm);
    }

    fprintf(summary_txt, " \n --Monitors: \n" );
    fprintf(summary_txt,"%-13s\t", "");
    for(int imon=0; imon<nMon; imon++)
      fprintf(summary_txt,"%-13s \t", mon_array[imon].Data());
    fprintf(summary_txt," \n");

    fprintf(summary_txt,"%-13s \t", "Mean+/-Err(um)");
    for(int imon=0; imon<nMon; imon++)
      fprintf(summary_txt,"%-3.2f+/-%-4.2f \t", 
	      statMon_avg[imon].mean/um,statMon_avg[imon].error/um);
    fprintf(summary_txt," \n");

    fprintf(summary_txt,"%-13s \t", "RMS(um)");
    for(int imon=0; imon<nMon; imon++)
      fprintf(summary_txt,"%-13.2f \t", statMon_avg[imon].rms/um);
    fprintf(summary_txt," \n");

    for(int imon=0;imon<nMon;imon++)
      mini_tree->SetBranchAddress("diff_"+mon_array[imon],&statMon[imon]);
    for(int idet=0;idet<nDet;idet++){
      mini_tree->SetBranchAddress("asym_"+det_array[idet],&statRawDet[idet]);
      mini_tree->SetBranchAddress("redi_asym_"+det_array[idet],&statCorDet[idet]);
    }
    for(int icom=0;icom<nCombo;icom++){
      mini_tree->SetBranchAddress("asym_"+alias_array[icom],&statRawCombo[icom]);
      mini_tree->SetBranchAddress("redi_asym_"+alias_array[icom],&statCorCombo[icom]);
    }
    Int_t nMinirun = mini_tree->GetEntries();
    for(int kMinirun=0; kMinirun<nMinirun;kMinirun++){
      mini_tree->GetEntry(kMinirun);

      fprintf(summary_txt,"\n\n\n\n"); 
      fprintf(summary_txt, " == Mini-run: %d == \n", kMinirun);
      fprintf(summary_txt, 
	      "%-12s\t%-12s\t%-12s\t%-20s\t%-15s\t%-15s\t%-12s\n",
	      "",
	      "Raw Mean+/-Err(ppm)","RMS(ppm)",
	      "Regressed Mean+/-Err(ppm)","RMS (ppm)", 
	      "Correction on Mean(ppm)","RMS of BPM Coherent Sum(ppm)");

      for(int idet=0;idet<nDet;idet++){
	double CoherentNoise = sqrt( pow(statRawDet[idet].rms,2)-
				     pow(statCorDet[idet].rms,2));
	fprintf(summary_txt, 
		"%-12s\t%-3.2f+/-%-9.2f\t%-12.2f\t%-3.2f+/-%-17.2f\t%-17.2f\t%-15.2f\t%-12.2f\n",
		det_array[idet].Data(), 
		statRawDet[idet].mean/ppm,statRawDet[idet].error/ppm,statRawDet[idet].rms/ppm, 
		statCorDet[idet].mean/ppm,statCorDet[idet].error/ppm,statCorDet[idet].rms/ppm, 
		(statRawDet[idet].mean-statCorDet[idet].mean)/ppm,CoherentNoise/ppm);
      }

      for(int icom=0;icom<nCombo;icom++){
	double CoherentNoise = sqrt( pow(statRawCombo[icom].rms,2)-
				     pow(statCorCombo[icom].rms,2));
	fprintf(summary_txt, 
		"%-12s\t%-3.2f+/-%-9.2f\t%-12.2f\t%-3.2f+/-%-17.2f\t%-17.2f\t%-15.2f\t%-12.2f\n",
		alias_array[icom].Data(), 
		statRawCombo[icom].mean/ppm,statRawCombo[icom].error/ppm,statRawCombo[icom].rms/ppm, 
		statCorCombo[icom].mean/ppm,statCorCombo[icom].error/ppm,statCorCombo[icom].rms/ppm, 
		(statRawCombo[icom].mean-statCorCombo[icom].mean)/ppm,CoherentNoise/ppm);
      }

      fprintf(summary_txt, " \n --Monitors: \n" );
      fprintf(summary_txt,"%-13s\t", "");
      for(int imon=0; imon<nMon; imon++)
	fprintf(summary_txt,"%-13s \t", mon_array[imon].Data());
      fprintf(summary_txt," \n");

      fprintf(summary_txt,"%-13s \t", "Mean+/-Err(um)");
      for(int imon=0; imon<nMon; imon++)
	fprintf(summary_txt,"%-3.2f+/-%-4.2f \t", 
		statMon[imon].mean/um,statMon[imon].error/um);
      fprintf(summary_txt," \n");

      fprintf(summary_txt,"%-13s \t", "RMS(um)");
      for(int imon=0; imon<nMon; imon++)
	fprintf(summary_txt,"%-13.2f \t", statMon[imon].rms/um);
      fprintf(summary_txt," \n");

      fprintf(summary_txt, " \n --Slope (ppm/um): \n");
      fprintf(summary_txt, " %-13s\t ","");
      for(int imon=0; imon<nMon; imon++)
	fprintf(summary_txt,"%-13s \t", mon_array[imon].Data());
      fprintf(summary_txt," \n");

      for(int idet=0;idet<nDet;idet++){
	fprintf(summary_txt,"%-13s \t", det_array[idet].Data());
	for(int imon=0; imon<nMon;imon++)
	  fprintf(summary_txt,"%-12.1f\t", 
		  miniDetMonSlopes[kMinirun][idet][imon]);
	fprintf(summary_txt," \n");
      }

      for(int icom=0;icom<nCombo;icom++){
      	fprintf(summary_txt,"%-13s \t", alias_array[icom].Data());
      	for(int imon=0; imon<nMon;imon++)
      	  fprintf(summary_txt,"%-12.1f\t",
      		  miniComboSlopes[kMinirun][icom][imon]);
      	fprintf(summary_txt," \n");
      }
    }  // End of Mini Run loop
    fclose(summary_txt);
    cout << "Done summary  text" << endl;
  }
}

void TaLagrangian::UpdateRunStat(RUN_STATS& runstat , 
				 TaAccumulator accumulator){
  runstat.mean = accumulator.GetMean1();
  runstat.nsamples = accumulator.GetN();
  runstat.rms = sqrt(accumulator.GetM2()/accumulator.GetN());
  runstat.error = runstat.rms/sqrt(runstat.nsamples);
}
