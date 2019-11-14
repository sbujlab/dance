#ifndef __TADATAELEMENT_HH__
#define __TADATAELEMENT_HH__

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
using namespace std;
class TaDataElement{
public:
  TaDataElement(TString name){
    myName = name;
  }
  TaDataElement(TString name, TBranch* fBranch_ptr){
    myName= name;
    RegisterBranchAddress(fBranch_ptr);
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
  }
private:
  TString myName;
  vector< pair<Double_t,TaDataElement* > > fElementArray;
  Double_t fhw_sum;
  Double_t fDevice_Error_Code;

  ClassDef(TaDataElement,0);
};

#endif
