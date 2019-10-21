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
class TaConfig{
public:
  TaConfig();
  virtual ~TaConfig(){};

  Bool_t ParseFile(TString fileName);

  vector<TString> GetDVlist(TString type, TString name);
  vector<TString> GetIVlist(TString type, TString name);
  TString GetConfigParameter(TString key);
  TString GetAnalysisParameter(Int_t index, TString key);
  Int_t GetAnalysisIndex(TString type, TString name);

  inline Int_t GetRunNumber() const { return run_number;};
  inline vector<TString> GetDeviceList() const {return device_list;};
  inline TString GetInputName() const { return input_name;};

  inline void SetInputName(TString str){input_name=str;};
  inline void SetRunNumber(Int_t i){run_number=i;};
  
  vector<VAnalysisModule*> GetAnalysisArray();

private:
  TString configName;
  Int_t run_number;
  TString input_name;

  map<TString, TString> fConfigParameters;
  map< pair<Int_t,TString> , TString> fAnalysisParameters;
  map< Int_t , vector<TString> > fDVMap;
  map< Int_t , vector<TString> > fIVMap;
  map<TString,Int_t> device_map;
  vector<TString> device_list;
  vector< pair<TString, TString> > fAnalysisTypeNames;
  map< pair<TString, TString>, Int_t > fAnalysisMap;

  vector<TString>  ParseLine(TString, TString);
  pair<TString,TString> GetAnalysisTypeName(TString);


  ClassDef(TaConfig,0);
};

#endif
