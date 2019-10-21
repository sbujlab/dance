#include "TaConfig.hh"
#include "TaRegression.hh"
#include "VAnalysisModule.hh"

ClassImp(TaConfig);
using namespace std;
TaConfig::TaConfig(){

}

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
  TString a_single_space =" ";  
  Bool_t isInModulde = kFALSE;
  vector<TString> vecStr;
  pair<Int_t,TString> myIDName;
  Int_t myIndex;
  Int_t index_ana = 0;
  Int_t channel_count=0;
  while(sline.ReadLine(configFile) ){
    if(sline.Contains(comment))
      continue;

    if(sline.Contains("[")){
      isInModulde = kTRUE;
      pair<TString,TString> aTypeName = GetAnalysisTypeName(sline);
      fAnalysisTypeNames.push_back(aTypeName);
      myIndex = index_ana++;
      fAnalysisMap[aTypeName] = myIndex;
      continue;
    }

    if(!isInModulde){
      vecStr = ParseLine(sline, a_single_space);
      fConfigParameters[vecStr[0]] = vecStr[1];
    }

    if(isInModulde){
      vecStr = ParseLine(sline,a_single_space);
      if(vecStr[0]=="dv"){
	(fDVMap[myIndex]).push_back(vecStr[1]);
	if(device_map.find(vecStr[1])==device_map.end()){
	  device_map[vecStr[1]]=channel_count++;
	  device_list.push_back(vecStr[1]);
	}
      }else if(vecStr[0]=="iv"){
	(fIVMap[myIndex]).push_back(vecStr[1]);
	if(device_map.find(vecStr[1])==device_map.end()){
	  device_map[vecStr[1]]=channel_count++;
	  device_list.push_back(vecStr[1]);
	}
      }else{
	pair<Int_t,TString> thiskey=make_pair(myIndex,vecStr[0]);
	fAnalysisParameters[thiskey]=vecStr[1];
      }
    } // end of if it inside an module section
  } // end of line loop
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

pair<TString,TString> TaConfig::GetAnalysisTypeName(TString sline){
#ifdef DEBUG
  cout << __PRETTY_FUNCTION__ << endl;
#endif

  Ssiz_t start_pt=sline.First('[')+1;
  Ssiz_t length = sline.First(':');
  length = length-start_pt;
  TString type = sline(start_pt,length);

  start_pt = sline.First(':')+1;
  length = sline.Last(']')-1;
  length = length - start_pt+1;
  TString name = sline(start_pt,length);
  cout << type << "\t " << name << endl;
  return make_pair(type,name);
}


vector<TString> TaConfig::GetDVlist(TString type,TString name){
  pair<TString,TString> aTypeName = make_pair(type,name);
  Int_t myIndex = fAnalysisMap[aTypeName];
  return fDVMap[myIndex];
}

vector<TString> TaConfig::GetIVlist(TString type,TString name){
  pair<TString,TString> aTypeName = make_pair(type,name);
  Int_t myIndex = fAnalysisMap[aTypeName];
  return fIVMap[myIndex];
}


TString TaConfig::GetConfigParameter(TString key){
  return fConfigParameters[key];
}

TString TaConfig::GetAnalysisParameter(Int_t index, TString key){
  return fAnalysisParameters[make_pair(index,key)];
}

Int_t TaConfig::GetAnalysisIndex(TString type, TString name){
  return fAnalysisMap[make_pair(type,name)];
}

vector<VAnalysisModule*> TaConfig::GetAnalysisArray(){
  Int_t nMod = fAnalysisTypeNames.size();
  vector<VAnalysisModule*> fAnalysisArray;
  for(int i=0; i<nMod;i++){
    TString type = fAnalysisTypeNames[i].first;
    TString name = fAnalysisTypeNames[i].second;
    VAnalysisModule* anAnalysis;
    if(type=="regression"){
      cout << "type == regression " << endl;
      anAnalysis = new TaRegression(type,name,this);
      fAnalysisArray.push_back(anAnalysis);
    }
  }
  return fAnalysisArray;
}
