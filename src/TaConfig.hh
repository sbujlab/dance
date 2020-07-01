#include "TString.h"
#include "Rtypes.h"

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

#ifndef __TaConfig_hh__
#define __TaConfig_hh__

#include "TaInput.hh"
#include "TaDefinition.hh"
class VAnalysisModule;
class TaRegression;
class TaDefinition;
class TaConfig{
public:
  TaConfig();
  virtual ~TaConfig(){};

  Bool_t ParseFile(TString fileName);
  TString GetConfigParameter(TString key);
  TString GetAnalysisParameter(Int_t index, TString key);
  inline vector<TString> GetAnalysisTypeArray() const {return fAnalysisTypeArray;};


  inline vector<TaDefinition*> GetDefList(TString key){
    return fGlobalDeviceMap[key];};
  inline vector<TaDefinition*> GetDefListByIndex(Int_t ana_idx,TString key){
    return fLocalDeviceMap[ana_idx][key];};
  inline vector<TaDefinition*> GetDeviceDefList() const {return device_list;};

  vector<TString> GetNameList(TString key){
    vector<TString> fRet;
    auto iter =  fGlobalDeviceMap[key].begin();
    while(iter!=fGlobalDeviceMap[key].end()){
      fRet.push_back( (*iter)->GetName());
      iter++;
    }
    return fRet;
  };

  vector<TString> GetNameListByIndex(Int_t ana_idx,TString key){
    vector<TString> fRet;
    auto iter =  fLocalDeviceMap[ana_idx][key].begin();
    while(iter!=fLocalDeviceMap[ana_idx][key].end()){
      fRet.push_back( (*iter)->GetName());
      iter++;
    }
    return fRet;
  };
  
  vector<TString> GetDeviceNameList() {
    vector<TString> fRet;
    auto iter =  device_list.begin();
    while(iter!=device_list.end()){
      fRet.push_back( (*iter)->GetName());
      iter++;
    }
    return fRet;
  };
  inline Int_t GetSegNumber() const { return seg_number;};
  inline Int_t GetRunNumber() const { return run_number;};
  inline TString GetInputName() const { return input_name;};
  void SetInputName(TString str){
    input_name=str;
    ParseRunNumber();
  };
  void SetRunNumber(Int_t i){
    if(i!=0)
      run_number=i;
  };
  inline void SetSegNumber(Int_t i){seg_number=i;};

  vector<VAnalysisModule*> GetAnalysisArray();
  TaDefinition* ParseChannelDefinition(TString );
  Bool_t isKeyWord(TString input );
  void UpdateDeviceList( vector<TaDefinition*>  &alist,
			TaDefinition* aDef);
  
  Bool_t CheckRunRange(TString input);
  TString FindExtRootfile(TString format);
  TString FindConfigByRange(TString format);
  void ParseRunNumber();
private:

  TString configName;
  Int_t run_number;
  Int_t seg_number;
  TString input_name;
  Int_t kUpperBound;
  Int_t kLowerBound;
  map< TString, TString> fConfigParameters;
  vector< map< TString, TString> > fAnalysisParameters;

  vector< TaDefinition*> device_list;
  map< TString, Int_t> device_map;
  map< TString, vector<TaDefinition*> > fGlobalDeviceMap;
  vector< map< TString, vector< TaDefinition*> > >  fLocalDeviceMap;
  vector< TString > fAnalysisTypeArray;

  vector<TString>  ParseLine(TString, TString);
  TString ParseAnalysisType(TString);

  ClassDef(TaConfig,0);
};

#endif
