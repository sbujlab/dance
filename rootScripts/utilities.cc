vector<Int_t> LoadRunListBySlug(Int_t slug_id){
  vector<Int_t> fRet;
  TString fullpath =Form("prex-runlist/simple_list/slug%d.list",slug_id);
  FILE *alist = fopen(fullpath.Data(),"r");
  while(!feof(alist)){
    Int_t run_number=0;
    fscanf(alist,"%d\n",&run_number);
    if(run_number!=0)
      fRet.push_back(run_number);
  }

  fclose(alist);
  return fRet;
}

map<Int_t,Int_t> LoadArmMapBySlug(Int_t slug_id){
  map<Int_t,Int_t> fret;
  TString info_filename =Form("prex-runlist/slug%d_info.list",slug_id);
  ifstream slug_info;
  slug_info.open(info_filename.Data());
  TString sline;
  while(sline.ReadLine(slug_info)){
    TObjArray *token = sline.Tokenize(',');
    Int_t run_number = (((TObjString*)(token->At(0)))->GetString()).Atoi();
    Int_t arm_flag = (((TObjString*)(token->At(6)))->GetString()).Atoi();
    fret[run_number]=arm_flag;
  }
  slug_info.close();
  return fret;
}

vector< vector<Int_t> > LoadSplitListBySlug(Int_t slug_id){
  vector< vector< Int_t> > fRet;
  TString filename = "rootScripts/splits.map";
  FILE *input = fopen(filename,"r");
  char list_char[256];
  Int_t myslug;
  while(!feof(input)){
    fscanf(input,"%d:%s\n",&myslug,list_char);
    if(slug_id!=myslug)
      continue;
    TString list_string = TString(list_char);
    TObjArray *token = list_string.Tokenize(',');
    Int_t nEntries = token->GetEntries();
    vector<Int_t> fRunList;
    for(int i=0;i<nEntries;i++){
      Int_t run_id = (((TObjString*)(token->At(i)))->GetString()).Atoi();
      fRunList.push_back(run_id);
    }
    fRet.push_back(fRunList);
  }

  if(fRet.size()==0)
    fRet.push_back(LoadRunListBySlug(slug_id));
  return fRet;
}

map<Int_t,vector<Int_t> > LoadBadCycleList(){
  map<Int_t,vector<Int_t > > fRet;
  cout << " -- Loading bad supercycle list " << endl;
  FILE *black_list = fopen("rootScripts/badcycle.list","r");
  if(black_list==NULL){
    cerr << " -- Error: bad cycle list is not found " << endl;
    return fRet;
  }
  
  while(!feof(black_list)){
    Int_t cycle=0;
    Int_t coil=0;
    fscanf(black_list,"%d,%d\n",&cycle,&coil);
    if(cycle!=0){
      if(coil==-1)
	for(int i=1;i<=7;i++)
	  fRet[cycle].push_back(i);
      else
	fRet[cycle].push_back(coil);
    }
  }
  fclose(black_list);
  return fRet;
}

Bool_t IsBadCycle(map<Int_t,vector<Int_t> > fMap, Int_t cycID, Int_t coil_idx=-1){

  if(fMap.find(cycID)==fMap.end())
    return kFALSE;
  
  if(coil_idx==-1) // only checking cycle number
    return kTRUE;
  
  vector<Int_t> coil_list = fMap[cycID];
  if(find(coil_list.begin(),coil_list.end(),coil_idx) == coil_list.end())
    return kFALSE;
  else
    return kTRUE;
}
