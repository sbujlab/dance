#ifndef __TaConfig_hh__
#define __TaConfig_hh__

#include "TaInput.hh"
#include "TaOutput.hh"

#include "TString.h"
#include "Rtypes.h"

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

using namespace std;
class TaOutput;
class VAnalysisModule;
class TaRegression;

class TaDefinition{
public:
  TaDefinition(TString input){ nick_name = input;};
  virtual ~TaDefinition();
  TString GetName(){ return nick_name;};
  void AddElement(Double_t factor,TString chname){
    fRawElementArray.push_back(make_pair(factor,chname));
  };
private:
  TString nick_name;
  vector< pair<Double_t, TString> > fRawElementArray;
};

class TaConfig{
public:
  TaConfig();
  virtual ~TaConfig(){};

  Bool_t ParseFile(TString fileName);
  TString GetConfigParameter(TString key);
  TString GetAnalysisParameter(Int_t index, TString key);
  inline vector<TString> GetAnalysisTypeArray() const {return fAnalysisTypeArray;};

  inline vector<TaDefinition> GetNameList(TString key) const {
    return fGlobalDeviceMap[key];};
  
  inline vector<TaDefinition> GetNameListByIndex(Int_t ana_idx,TString key) const {
    return fLocalDeviceMap[ana_idx][key];};
  
  inline vector<TaDefinition> GetDeviceList() const {return device_list;};
  
  inline Int_t GetRunNumber() const { return run_number;};
  inline TString GetInputName() const { return input_name;};
  inline void SetInputName(TString str){input_name=str;};
  inline void SetRunNumber(Int_t i){run_number=i;};
  
  vector<VAnalysisModule*> GetAnalysisArray();
  TaDefinition ParseChannelDefinition(TString );
  void LoadDVChannel(TString);
  Bool_t isKeyWord(TString input );
private:

  TString configName;
  Int_t run_number;
  TString input_name;

  map< TString, TString> fConfigParameters;
  vector< map< TString, TString> > fAnalysisParameters;

  vector< TaDefinition> device_list;
  map< TString, Int_t> device_map;
  map< TString, vector<TaDefinition> > fGlobalDeviceMap;
  vector< map< TString, vector< TaDefinition> > >  fLocalDeviceMap;
  vector< TString > fAnalysisTypeArray;

  vector<TString>  ParseLine(TString, TString);
  TString ParseAnalysisType(TString);

  ClassDef(TaConfig,0);
};

#endif
