#ifndef __TaInput_hh__
#define __TaInput_hh__

#include "TaConfig.hh"
#include "TaChannel.hh"
#include "TaOutput.hh"

#include <vector>
#include <map>
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;
class TaConfig;
class TaChannel;
class TaDefinition;
class TaOutput;
class TaInput{
public:
  TaInput(TaConfig *aConfig);
  virtual ~TaInput();

  inline void SetRunNumber(Int_t i){run_number=i;};

  inline TTree* GetEvtTree() const { return evt_tree;};
  inline TTree* GetMulTree() const { return mul_tree;};
  inline Int_t GetRunNumber() const { return run_number;};
  inline Int_t GetSegNumber() const { return seg_number;};

  inline TaChannel* GetChannelCutFlag(){return fChannelCutFlag;};
  inline vector< pair<Int_t,Int_t> > GetMiniRange(){return minirun_range;};

  void InitChannels(TaConfig*);
  void ConnectChannels();
  void WriteRawChannels(TaOutput*);
  Bool_t LoadROOTFile();
  TaChannel* GetChannel(TString name);

  void Close();

private:
  TString input_name;
  TString input_path;
  TString input_prefix;

  TFile* input_file;
  TTree* evt_tree;
  TTree* mul_tree;
  TTree* mulc_tree;
  Double_t run_number;
  Double_t seg_number;
  TString ext_filename;
  Bool_t isExternalConstraint;

  Int_t minirun_size;
  vector<pair<Int_t, Int_t> > minirun_range;
  Int_t nEntries;
  
  map<TString, TaChannel*> fChannelMap;
  vector<TaChannel*> fChannelArray;
  vector<TString> fChannelNames;
  
  TaChannel* fChannelErrorFlag;
  TaChannel* fChannelCutFlag;
  ClassDef(TaInput,0);
};

#endif
