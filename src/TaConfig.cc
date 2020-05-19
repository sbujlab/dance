#include "TaConfig.hh"
#include "TaRegression.hh"
#include "TaLagrangian.hh"
#include "TaCorrection.hh"
#include "VAnalysisModule.hh"
#include "TSystemDirectory.h"
ClassImp(TaConfig);
using namespace std;
TaConfig::TaConfig(){}

Bool_t TaConfig::ParseFile(TString fileName){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  configName=fileName;
  ifstream configFile;
  cout << " -- Opening " << configName << endl;
  configFile.open(configName.Data());
  if(!configFile.is_open()){
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
  Int_t myIndex=-1;
  Int_t index_ana = 0;
  while(sline.ReadLine(configFile) ){
    if(sline.Contains(comment)) // FIXME
      continue;
    
    if(sline.Contains("[")){
      isInModulde = kTRUE;
      TString myType = ParseAnalysisType(sline);
      fAnalysisTypeArray.push_back(myType);
      myIndex = index_ana++;
      map<TString, vector<TaDefinition*> >aNewChMap;
      map<TString, TString> aNewParmMap;
      fLocalDeviceMap.push_back(aNewChMap);
      fAnalysisParameters.push_back(aNewParmMap);
      continue;
    }
    
    vecStr = ParseLine(sline, a_single_space);
    if(isKeyWord(vecStr[0])){ // is Channel List
      TaDefinition* aDef =ParseChannelDefinition(vecStr[1]);
      if(isInModulde)
	UpdateDeviceList(fLocalDeviceMap[myIndex][vecStr[0]],aDef);
      else
	UpdateDeviceList(fGlobalDeviceMap[vecStr[0]],aDef);
      
    }else{  // is analysis parameter
      
      if(isInModulde)
	fAnalysisParameters[myIndex][vecStr[0]]=vecStr[1];
      else
	fConfigParameters[vecStr[0]] = vecStr[1];
      
    } // end of parsing channel list / definition
    
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

TString TaConfig::ParseAnalysisType(TString sline){
#ifdef DEBUG
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  Ssiz_t start_pt=sline.First('[')+1;
  Ssiz_t length = sline.First(']');
  length = length-start_pt;
  TString this_type = sline(start_pt,length);
#ifdef DEBUG
  cout << "Parsing type:" << this_type << endl;
#endif
  return this_type;
}

TString TaConfig::GetConfigParameter(TString key){
  return fConfigParameters[key];
}

TString TaConfig::GetAnalysisParameter(Int_t index, TString key){
  return fAnalysisParameters[index][key];
}

vector<VAnalysisModule*> TaConfig::GetAnalysisArray(){
  Int_t nMod = fAnalysisTypeArray.size();
  vector<VAnalysisModule*> fAnalysisArray;
  for(int i=0; i<nMod;i++){
    TString type = fAnalysisTypeArray[i];
    VAnalysisModule* anAnalysis;
    if(type=="regression"){
      cout << "type == regression " << endl;
      anAnalysis = new TaRegression(i,this);
      fAnalysisArray.push_back(anAnalysis);
    }
    if(type=="lagrangian"){
      cout << "type == lagrangian " << endl;
      anAnalysis = new TaLagrangian(i,this);
      fAnalysisArray.push_back(anAnalysis);
    }
    if(type=="dithering"){
      cout << "type == dithering " << endl;
      anAnalysis = new TaCorrection(i,this);
      fAnalysisArray.push_back(anAnalysis);
    }
  }
  return fAnalysisArray;
}

TaDefinition* TaConfig::ParseChannelDefinition(TString input){
  input.ReplaceAll(" ","");
  TaDefinition* aDef;
  if(!input.Contains("=")){
    aDef = new TaDefinition(input);
    aDef->SetDefFlag(kFALSE);
  } else{
    Ssiz_t equal_pos =  input.First('=');
    TString channel_name = input(0,equal_pos);
    aDef = new  TaDefinition(channel_name);
    aDef->SetDefFlag(kTRUE);
    input.Remove(0,equal_pos+1);
    TString def_formula = input;
    while(def_formula.Length()>0){
      Ssiz_t next_plus = def_formula.Last('+');
      Ssiz_t next_minus = def_formula.Last('-');
      Ssiz_t length = def_formula.Length();
      Ssiz_t head =0;
      if(next_minus>next_plus)
	head = next_minus;
      else if(next_plus>next_minus)
	head = next_plus;

      TString extracted = def_formula(head,length-head);
#ifdef DEBUG
      cout << extracted << ":";
#endif
      if(extracted.Contains("*")){
	Ssiz_t aster_pos = extracted.First('*');
	Ssiz_t form_length = extracted.Length();
	TString elementName = extracted(aster_pos+1,form_length-(aster_pos+1));
	TString coeff = extracted(0,aster_pos);
	aDef->AddElement(coeff.Atof(),elementName);
	if(device_map.find(elementName)==device_map.end()){
	  TaDefinition* aRawElement = new TaDefinition(elementName);
	  device_list.push_back(aRawElement);
	  Int_t index = device_list.size()-1;
	  device_map[elementName]=index;      
	}
#ifdef DEBUG
	cout << coeff<< "  " << elementName << endl;
#endif
      }else {
	Double_t coeff;
	if(extracted.Contains("-"))
	  coeff=-1.0;
	else
	  coeff=1.0;
	extracted.ReplaceAll("+","");
	extracted.ReplaceAll("-","");
	aDef->AddElement(coeff,extracted);
	if(device_map.find(extracted)==device_map.end()){
	  TaDefinition *aRawElement = new TaDefinition(extracted);
	  device_list.push_back(aRawElement);
	  Int_t index = device_list.size()-1;
	  device_map[extracted]=index;      
	}
#ifdef DEBUG
	cout << coeff<< "  " << extracted << endl;
#endif
      }
      def_formula.Remove(head,length-head);
    } // end of detecting formula
  }  // end  of if HasUser Def_Formula
  
  TString ch_name = aDef->GetName();
  if(device_map.find(ch_name)==device_map.end()){
    device_list.push_back(aDef);
    Int_t index = device_list.size()-1;
    device_map[ch_name]=index;      
  }
  return aDef;
}

Bool_t TaConfig::isKeyWord(TString input){
  if(input=="det" ||
     input=="mon" ||
     input=="dv" ||
     input=="iv" ||
     input=="coil" ||
     input=="def")
    return kTRUE;
  else
    return kFALSE;
}

void TaConfig::UpdateDeviceList(vector<TaDefinition*> &alist,
			       TaDefinition* aDef){
  TString myName =aDef->GetName();
  auto iter = alist.begin();
  while(iter!=alist.end()){
    if(myName == (*iter)->GetName())
      return;
    iter++;
  }
  alist.push_back(aDef);
}

Bool_t TaConfig::CheckRunRange(TString input){
  TString base_name = input;
  Ssiz_t last_slash = base_name.Last('/');
  base_name.Remove(0,last_slash+1);

  Ssiz_t first_dot = base_name.First('.');
  Ssiz_t last_dot = base_name.Last('.');

  Bool_t kFound = kFALSE;
  TString range = base_name(first_dot+1,last_dot-first_dot);
  if(range.Contains('-')){
    Ssiz_t dash_pos = range.First('-');
    TString lower = range(0,dash_pos);
    TString upper = range(dash_pos+1,range.Length()-(dash_pos+1));
    if(run_number<=upper.Atoi() && run_number>=lower.Atoi())
      kFound = kTRUE;
  } else{
    if(run_number == range.Atoi())
      kFound = kTRUE;
  }
  return kFound;
}

TString TaConfig::FindExtRootfile(TString ext_format){
  TString ext_filename;
  TString target_dir = ext_format;
  Ssiz_t length = target_dir.Length();
  Ssiz_t last_slash = target_dir.Last('/');
  Ssiz_t last_dot = target_dir.Last('.');
  TString ext_prefix = target_dir(last_slash+1,last_dot-last_slash);
  target_dir.Remove(last_slash+1, length-(last_slash+1) );
  const char* dir_name = "";
  const char* path= target_dir.Data();
  TSystemDirectory *sysDir = new TSystemDirectory(dir_name,path);
  TList* fileList = sysDir->GetListOfFiles();
  if(fileList){
    TIter next(fileList);
    TSystemFile* sysfile;
    while( (sysfile=(TSystemFile*)next()) ){
      if(sysfile->IsDirectory())
	continue;
      TString name_buff = sysfile->GetName();
      if(name_buff.Contains(ext_prefix)){
	if(CheckRunRange(name_buff)){
	  ext_filename = target_dir+name_buff;
	  cout << " -- Found run range specified matrix: \n -- " 
	       << ext_filename << endl;
	  return ext_filename;
	}
      }

    }// end of fileList Loop;
  }
  cout << " -- Using default input  matrix: \n -- " 
       << ext_format << endl;

  return ext_format;
}
