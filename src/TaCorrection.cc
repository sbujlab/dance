#include "TaCorrection.hh"
ClassImp(TaCorrection);

TaCorrection::TaCorrection(Int_t ana_index,TaConfig *fConfig){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  Init(ana_index,fConfig);
  vector<TaDefinition*> fDVDefList = fConfig->GetDefList("dv");
  auto iter_def = fDVDefList.begin();
  while(iter_def!=fDVDefList.end()){
    TString myName = (*iter_def)->GetName();
    if( find(sDVlist.begin(),sDVlist.end(),myName)!=sDVlist.end()){
      if( (*iter_def)->HasUserDefinition())
	LoadRawDVList( (*iter_def)->GetRawChannelList());
      else
	LoadRawDVList(myName);
    }
    iter_def++;
  }
  LoadSlopeMatrix(ana_index,fConfig);
}

void TaCorrection::Process(TaOutput *fOutput){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  ConstructOutputs(fOutput);
  Int_t nMini =  minirun_range.size();
  Int_t nDV = sDVlist.size();
  for(int ich=0;ich<nDV;ich++){
    vector<Double_t> fprefactors = GetSlopeVector( sDVlist[ich]);
    fCorrections[ich]->ConnectChannels(fIndependentVar,fprefactors);
    fCorrections[ich]->ConstructSlopeBranch(fOutput,"mini_"+tree_name);
  }
  
  for(int imini=0; imini<nMini;imini++){
    int istart = minirun_range[imini].first;
    int iend = minirun_range[imini].second;
    for(int ievt=istart;ievt<=iend;ievt++){
      GetEntry(ievt);
      CalcCombination();
      AccumulateMiniSum();
      AccumulateRunSum();
      if(!kOutputMiniOnly)
	fOutput->FillTree(tree_name);
    } // end of  event loop
    UpdateMiniStat();
    fOutput->FillTree("mini_"+tree_name);
    ResetMiniAccumulator();
  } // end of mini run loop
  UpdateRunStat();
  fOutput->FillTree("sum_"+tree_name);
}

vector<Double_t> TaCorrection::GetSlopeVector(TString dv_name){
  TaChannel* ch_ptr = fDVMaps[dv_name];
  Int_t nIV  = sIVlist.size();
  vector<Double_t> fret(nIV,0);
  if( !(ch_ptr->HasUserDefinition())){
    Int_t idx = FindRawDVIndexFromList(dv_name);
    for(int icol=0;icol<nIV;icol++)
      fret[icol]= slope_matrix[idx][icol];
  }else{
    vector<TString> fRawElementList = ch_ptr->GetRawChannelList();
    vector<Double_t> fPrefactors = ch_ptr->GetFactorArray();
    Int_t nelem = fPrefactors.size();
    for(int i=0;i<nelem;i++){
      Int_t idx = FindRawDVIndexFromList(fRawElementList[i]);
      for(int imon=0;imon<nIV;imon++)
	fret[imon] += fPrefactors[i]*slope_matrix[idx][imon];
    }
  }
  return fret;
}


Bool_t TaCorrection::LoadSlopeMatrix(Int_t ana_index, TaConfig* aConfig){
  
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  TString ext_format = aConfig->GetAnalysisParameter(ana_index,"slope_input");
  TString ext_filename = aConfig->FindExtRootfile(ext_format);
  TFile ext_file(ext_filename);
  if(!ext_file.IsOpen()){
    cout << ext_filename <<" doesn't exist!" << endl;
    return kFALSE;
  }

  vector<TString> *ext_dv_array =
    (vector<TString>*)ext_file.Get("dv_array");
  vector<TString> *ext_iv_array =
    (vector<TString>*)ext_file.Get("iv_array");
  TMatrixD* input_matrix = (TMatrixD*) ext_file.Get("slope_matrix");
  Int_t nDV  = fRawDVList.size();
  Int_t nIV  = sIVlist.size();
  slope_matrix.ResizeTo(nDV,nIV);
  for(int i=0;i<nDV;i++){
    auto iter_dv = find((*ext_dv_array).begin(),(*ext_dv_array).end(),
			  GetBaseName(fRawDVList[i]));
    Int_t index_dv = iter_dv-(*ext_dv_array).begin();
    for(int j=0;j<nIV;j++){
      auto iter_iv = find((*ext_iv_array).begin(),(*ext_iv_array).end(),
			  GetBaseName(sIVlist[j]));
      Int_t index_iv = iter_iv-(*ext_iv_array).begin();
      slope_matrix[i][j]= (*input_matrix)[index_dv][index_iv];
    }
  }
  return kTRUE;
}

TString TaCorrection::GetBaseName(TString in){
  in.ReplaceAll("diff_","");
  in.ReplaceAll("asym_","");
  return in;
}

Int_t TaCorrection::FindRawDVIndexFromList(TString raw_name){
  Int_t index = -1;
  auto iter_find = find(fRawDVList.begin(),fRawDVList.end(),
			raw_name);
  if(iter_find!=fRawDVList.end())
    index = iter_find - fRawDVList.begin();

  return index;
}

void TaCorrection::LoadRawDVList(vector<TString> fRawElementList){
  auto iter_ele = fRawElementList.begin();
  while(iter_ele!=fRawElementList.end()){
    if(FindRawDVIndexFromList(*iter_ele)==-1)
      fRawDVList.push_back(*iter_ele);
    iter_ele++;
  }
}

void TaCorrection::LoadRawDVList(TString raw_name){
  if(FindRawDVIndexFromList(raw_name)==-1)
    fRawDVList.push_back(raw_name);
}

