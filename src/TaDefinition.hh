#ifndef __TaDefinition_hh__
#define __TaDefinition_hh__
#include "TString.h"
#include "Rtypes.h"

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

class TaDefinition{
public:
  TaDefinition(){ kUserDefinition=kFALSE;};
  TaDefinition(TString input){ nick_name = input; kUserDefinition=kFALSE;};
  virtual ~TaDefinition();
  TString GetName(){ return nick_name;};
  void AddElement(Double_t factor,TString chname);
  vector<Double_t> GetFactorArray();
  vector<TString> GetRawChannelList();
  Bool_t HasUserDefinition(){return kUserDefinition;};
  void SetDefFlag(Bool_t input){ kUserDefinition=input;};

private:
  TString nick_name;
  vector< pair<Double_t, TString> > fRawElementArray;
  Bool_t kUserDefinition;
  ClassDef(TaDefinition,0);
};

#endif
