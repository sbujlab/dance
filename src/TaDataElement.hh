#ifndef __TADATAELEMENT_HH__
#define __TADATAELEMENT_HH__

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TaConfig.hh"
#include "TaDefinition.hh"
using namespace std;
class TaDataElement{
public:
  TaDataElement(TaDefinition* adef){
    myName = adef->GetName();
    myDef = adef;
    kUserDefined = adef->HasUserDefinition();
    if(kUserDefined)
      cout << myName << " has user definition." << endl;
  };
  TaDataElement(TString name){
    myName = name;
    kUserDefined = kFALSE;
  };
  TaDataElement(TString name, TBranch* fBranch_ptr){
    myName= name;
    RegisterBranchAddress(fBranch_ptr);
    kUserDefined = kFALSE;
  };
  virtual ~TaDataElement(){};

  void RegisterBranchAddress(TBranch*);
  Double_t GetHwSum();
  Double_t TestDeviceErrorCode(Int_t ErrorMask=0xFF);
  Double_t GetDeviceErrorCode();
  void AddElement(Double_t, TaDataElement*);
  TString GetName(){return myName;};
  pair<Double_t ,TaDataElement*> GetElement(Int_t i) const {return fElementArray[i];};
  Int_t GetNumberOfElements(){
    if(fElementArray.size()<=1)
      return 1;
    else
      return fElementArray.size();
  };
  Bool_t HasUserDefinition(){
    return kUserDefined;
  };
  vector<Double_t> GetFactorArray(){
    return myDef->GetFactorArray();
  };
  vector<TString> GetRawNameList(){
    return myDef->GetRawChannelList();
  };
private:
  TString myName;
  TaDefinition* myDef;
  vector< pair<Double_t,TaDataElement* > > fElementArray;
  Double_t fhw_sum;
  Double_t fDevice_Error_Code;
  Bool_t kUserDefined;
  ClassDef(TaDataElement,0);
};

#endif
