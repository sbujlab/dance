#include "TaConfig.hh"
#include "TaRegression.hh"
#include "VAnalysisModule.hh"

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
      if(vecStr[0]=="dv"){
	fDVlist.push_back(vecStr[1]);
	if(device_map.find(vecStr[1])==device_map.end()){
	  device_map[vecStr[1]]=channel_count++;
	  device_list.push_back(vecStr[1]);
	}
      }
      else if(vecStr[0]=="det" || vecStr[0]=="mon"){
	if(vecStr[1].Contains("=")){
	  vector<TString> elements_array=ParseChannelDefinition(vecStr[1]);
	  vector<TString>::iterator iter_ele = elements_array.begin();
	  while(iter_ele!=elements_array.end()){
	    TString myName = *iter_ele;
	    if(device_map.find(myName)==device_map.end()){
	      device_map[myName]=channel_count++;
	      device_list.push_back(myName);
	      fDependentVarArray.push_back(myName);
#ifdef DEBUG
	      cout << myName << endl;
#endif
	      if(iter_ele!=elements_array.begin())
		fRawElementArray.push_back(myName);

	      if(vecStr[0]=="det")
		fDetectorArray.push_back(myName);
	    }
	    iter_ele++;
	  }

	}else{
	  if(device_map.find(vecStr[1])==device_map.end()){
	    device_map[vecStr[1]]=channel_count++;
	    device_list.push_back(vecStr[1]);
#ifdef DEBUG
	      cout << vecStr[1] << endl;
#endif
	    fDependentVarArray.push_back(vecStr[1]);
	    fRawElementArray.push_back(vecStr[1]);

	    if(vecStr[0]=="det")
	      fDetectorArray.push_back(vecStr[1]);
	  } 
	}// end of if it is a raw channel
      }
      else
	fConfigParameters[vecStr[0]] = vecStr[1];
    }
    else{ // else if is InModule
      vecStr = ParseLine(sline, a_single_space);
      if(vecStr[0]=="iv"){
	(fIVMap[myIndex]).push_back(vecStr[1]);
	if(device_map.find(vecStr[1])==device_map.end()){
	  device_map[vecStr[1]]=channel_count++;
	  device_list.push_back(vecStr[1]);
	}
      }else if(vecStr[0]=="mon"){
	(fMonitorMap[myIndex]).push_back(vecStr[1]);

	if(vecStr[1].Contains("=")){
	  vector<TString> elements_array=ParseChannelDefinition(vecStr[1]);
	  vector<TString>::iterator iter_ele = elements_array.begin();
	  while(iter_ele!=elements_array.end()){
	    TString myName = *iter_ele;
	    if(device_map.find(myName)==device_map.end()){
	      device_map[myName]=channel_count++;
	      device_list.push_back(myName);
	      fDependentVarArray.push_back(myName);
#ifdef DEBUG
	      cout << myName << endl;
#endif
	      if(iter_ele!=elements_array.begin())
		fRawElementArray.push_back(myName);
	    }
	    iter_ele++;
	  }

	}else{
	  if(device_map.find(vecStr[1])==device_map.end()){
	    device_map[vecStr[1]]=channel_count++;
	    device_list.push_back(vecStr[1]);

	    fDependentVarArray.push_back(vecStr[1]);
	    fRawElementArray.push_back(vecStr[1]);
	  } 
	}
	
      }else if(vecStr[0]=="coil"){
	(fCoilMap[myIndex]).push_back(vecStr[1]);
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

  TString name;
  if(sline.Contains(':')){
    start_pt = sline.First(':')+1;
    length = sline.Last(']')-1;
    length = length - start_pt+1;
    name = sline(start_pt,length);
  }
  else
    name="";

#ifdef DEBUG
  cout << type << "\t " << name << endl;
#endif

  return make_pair(type,name);
}


// vector<TString> TaConfig::GetDVlist(TString type,TString name){
//   pair<TString,TString> aTypeName = make_pair(type,name);
//   Int_t myIndex = fAnalysisMap[aTypeName];
//   return fDVMap[myIndex];
// }

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

vector<TString> TaConfig::ParseChannelDefinition(TString input){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  vector<TString> elements_array;
  Ssiz_t equal_pos =  input.First('=');
  TString combo_name = input(0,equal_pos);
  elements_array.push_back(combo_name);
  input.Remove(0,equal_pos+1);
  TString def_formula = input;
  def_formula.ReplaceAll(" ","");
  auto iter_combo = find(fDefinedElementArray.begin(),
			 fDefinedElementArray.end(),
			 combo_name);
  if(iter_combo==fDefinedElementArray.end()){
    fDefinedElementArray.push_back(combo_name);
    fDataElementDefinitions.push_back(make_pair(combo_name,def_formula));
  }
#ifdef DEBUG
  cout << "combo_name:" << combo_name << endl ;
  cout << "def_formula:" << def_formula << endl ;
#endif
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
  cout << extracted << endl;
#endif
    if(extracted.Contains("*")){
      Ssiz_t aster_pos = extracted.First('*');
      Ssiz_t form_length = extracted.Length();
      TString elementName = extracted(aster_pos+1,form_length-(aster_pos+1));
      elements_array.push_back(elementName);
#ifdef DEBUG
  cout << elementName << endl;
#endif
    }else {
      extracted.ReplaceAll("+","");
      extracted.ReplaceAll("-","");
      elements_array.push_back(extracted);
#ifdef DEBUG
      cout << extracted << endl;
#endif
    }
    def_formula.Remove(head,length-head);
  }
  
  return elements_array;
}
