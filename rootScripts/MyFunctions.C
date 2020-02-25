std::map<Int_t, TString> GetBadCycleCut(Int_t,Int_t);
std::vector<Int_t> ParseLine(char*);

std::map<Int_t, TString>  GetBadCycleCut(Int_t low, Int_t up){
  FILE *blacklist = fopen("./scripts/black.list","r");
  char cyc_c[256];
  char coilString[256];

  std::map< Int_t, vector<pair<Int_t,Int_t> > > fBadCycleMap;
  while(!feof(blacklist)){

    fscanf(blacklist,"%s\t%s\n",cyc_c,coilString);
    TString cycString = TString((const char*)cyc_c);
    Int_t cyc_low, cyc_up;
    if(cycString.Contains('-')){
      Ssiz_t pos_t = cycString.First('-');
      Ssiz_t length_t = cycString.Length();
      TString start_str = cycString(0,pos_t);
      TString end_str = cycString(pos_t+1,length_t-pos_t-1);
      cyc_low =start_str.Atoi();
      cyc_up =end_str.Atoi();
    }else{
      cyc_low =cycString.Atoi();
      cyc_up =cycString.Atoi();
    }

    if( cyc_low<=up && cyc_up>=low ){
      std::vector<Int_t> coilArray  = ParseLine(coilString);
      auto coil_iter = coilArray.begin();
      while(coil_iter!=coilArray.end()){
	fBadCycleMap[*coil_iter].push_back(make_pair(cyc_low,cyc_up));
	coil_iter++;
      }
    }

  }
  fclose(blacklist);

  std::map< Int_t, TString>  fBadCycleCut;
  auto iter_map = fBadCycleMap.begin();
  while(iter_map!=fBadCycleMap.end()){
    Int_t coilID = (*iter_map).first;
    auto cycpair_iter = ((*iter_map).second).begin();
    TString myCut;
    while(cycpair_iter!=((*iter_map).second).end()){

      if(cycpair_iter==((*iter_map).second).begin())
	myCut = Form("(!(cycID>=%d && cycID<=%d))",
		     (*cycpair_iter).first,
		     (*cycpair_iter).second);
      else
	myCut += Form("&& (!(cycID>=%d && cycID<=%d))",
		      (*cycpair_iter).first,
		      (*cycpair_iter).second);

      cycpair_iter++;
    }
    fBadCycleCut.insert(std::make_pair(coilID,myCut));
    iter_map++;
  }
  return fBadCycleCut;
}

std::vector<Int_t> ParseLine( char* input_char){
  TString input_string = TString((const char*)input_char);
  std::vector<Int_t> ret_array;
  TString delim=",";
  if(input_string.Contains("-1")){
    for(int i=1;i<=7;i++){
      ret_array.push_back(i);
    }
    return ret_array;
  }
  
  while(input_string.Contains(delim)){
    Ssiz_t end_pt = input_string.First(delim);
    TString new_entry = input_string(0,end_pt);
    ret_array.push_back(new_entry.Atoi());
    input_string.Remove(0,end_pt+1);
  }
  ret_array.push_back(input_string.Atoi());
  return ret_array;
}
