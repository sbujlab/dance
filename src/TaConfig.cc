#include "TaConfig.hh"

ClassImp(TaConfig);
using namespace std;
  
TaConfig::TaConfig(){
}
TaConfig::~TaConfig(){}
Bool_t TaConfig::ParseFile(TString fileName){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  configName=fileName;
  ifstream configFile;
  cout << " -- Opening " << configName << endl;
  configFile.open(configName.Data());
  if(configFile==NULL){
    cerr << __PRETTY_FUNCTION__ 
	 << " Error: failed to open config file "
	 << configName << endl;
    return kFALSE;
  }
  TString comment = "#";
  TString sline;
  vector<TString > vecLine;
  
  while(sline.ReadLine(configFile) ){
    if(sline.Contains(comment))
      continue;
    vecLine.push_back(sline);
  }
  std::vector<TString>::iterator iter = vecLine.begin();
  vector<TString> vecStr;
  while (iter!=vecLine.end()){
    vecStr = ParseLine(*iter, " ");
    iter++;
    if(vecStr[0].Contains("input_path")){
      input_path = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("input_prefix")){
      input_prefix = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("output_path")){
      output_path = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("output_prefix")){
      output_prefix = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("analysisType")){
      analysisType = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("minirun_size")){
      minirun_size = vecStr[1].Atoi();
      continue;
    }
    if(vecStr[0].Contains("event_cut")){
      event_cut = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("custom_cut")){
      custom_cut = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("det")){
      det_array.push_back(vecStr[1]);
      continue;
    }
    if(vecStr[0].Contains("mon")){
      mon_array.push_back(vecStr[1]);
      continue;
    }
    if(vecStr[0].Contains("coil")){
      coil_array.push_back(vecStr[1]);
      continue;
    }
    if(vecStr[0].Contains("alias")){
      alias_array.push_back(ParseAlias(vecStr[1]));
#ifdef DEBUG
      pair<TString,ElementsVector> ret = ParseAlias(vecStr[1]);
      cout << "-- Getting Alias " << ret.first << "= " ;
      ElementsVector::iterator iter = (ret.second).begin();
      while(iter!=(ret.second).end()){
	cout << (*iter).first << "*";
	cout << (*iter).second << "\t";
	iter++;
      }
      cout << endl;
#endif
    }
    else{
      cerr << __FILE__ << ":"
  	   << __FUNCTION__ << ":"
  	   << "unknown config parameters "
  	   << vecStr[0]
  	   << " will be ignored. " << endl;
    }
  }

  configFile.close();
  return kTRUE;
}

vector<TString> TaConfig::ParseLine(TString sline, TString delim){
  vector<TString> ret;
  if(sline.Contains(delim)){
    Int_t index = sline.First(delim);
    TString buff = sline(0,index);
    ret.push_back(buff);
    sline.Remove(0,index+1);
    sline.ReplaceAll(" ","");
    ret.push_back(sline);
  }
  return ret;
}

pair<TString,ElementsVector> TaConfig::ParseAlias(TString sline){
  pair<TString,ElementsVector> ret;
  if(sline.Contains("=")){
    Int_t pos_delim=sline.First("=");
    TString aliasName = sline(0,pos_delim);
    sline.Remove(0,pos_delim+1);
    ElementsVector elemVect;
    while(sline.Contains("*")){
      Int_t pos_delim = sline.First("*");
      Double_t factor= TString(sline(0,pos_delim)).Atof();
      sline.Remove(0,pos_delim+1);
      if(sline.Contains("+")){
	Int_t next_sign = sline.First("+");
	TString element_name = sline(0,next_sign);
	elemVect.push_back(make_pair(factor,element_name));
	sline.Remove(0,next_sign);
      }else if(sline.Contains("-")){
	Int_t next_sign = sline.First("-");
	TString element_name = sline(0,next_sign);
	elemVect.push_back(make_pair(factor,element_name));
	sline.Remove(0,next_sign);
      }else{ // it is all the way to the end
	TString element_name = sline;
	elemVect.push_back(make_pair(factor,element_name));
      }
    }
    ret = make_pair(aliasName,elemVect);
  }
  return ret;
}
