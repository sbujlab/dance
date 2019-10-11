#ifndef __TaConfig_hh__
#define __TaConfig_hh__

#include "TaInput.hh"
#include "TString.h"
#include "Rtypes.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;
typedef pair<Double_t, TString> AliasElement;
typedef vector<AliasElement> ElementsVector;
class TaConfig{
 private:
  TString input_path;
  TString input_prefix;
  TString output_path;
  TString output_prefix;
  Int_t run_num;
  Int_t seg_num;
  Int_t minirun_size;
  TString output_name;
  TString input_name;
  TString analysisType;
  TString event_cut;
  TString custom_cut;

  vector<TString> det_array;
  vector<TString> mon_array;
  vector<TString> coil_array;
  vector< pair<TString,ElementsVector> > alias_array;

  TString configName;
  
  vector<TString>  ParseLine(TString, TString);
  pair<TString,ElementsVector> ParseAlias(TString);
 public:
  TaConfig();
  virtual ~TaConfig();
  Bool_t ParseFile(TString fileName);

  inline TString GetInputName() const { return input_name;};
  inline TString GetOutputName() const { return output_name;};
  inline Int_t GetRunNumber() const { return run_num;};
  inline Int_t GetMiniSize() const { return minirun_size;};
  inline TString GetAnalysisType() const{ return analysisType;};
  inline TString GetInputPath() const {return input_path;};
  inline TString GetInputPrefix() const {return input_prefix;};
  inline TString GetOutputPath() const {return output_path;};
  inline TString GetOutputPrefix() const {return output_prefix;};
  inline TString GetEventCut() const {return event_cut;};  
  inline TString GetCustomCut() const {return custom_cut;};  
  inline vector<TString> GetDetArray() const {return det_array;};
  inline vector<TString> GetMonArray() const {return mon_array;};
  inline vector<TString> GetCoilArray() const {return coil_array;};
  inline vector<pair<TString,ElementsVector> > GetAliasArray() const {return alias_array;};

  // inline void SetOutputName(TString str){output_name=str;};
  inline void SetInputName(TString str){input_name=str;};
  // inline void SetAnalysisType(Int_t i){analysisType=i;};
  // inline void SetRunNumber(Int_t i){run_num=i;};
  // inline void SetConfigFile(TString str){configName=str;};

  ClassDef(TaConfig,0);
};

#endif
