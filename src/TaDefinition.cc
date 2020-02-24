#include "TaDefinition.hh"

ClassImp(TaDefinition)

TaDefinition::~TaDefinition(){}

void TaDefinition::AddElement(Double_t factor,TString chname){
  fRawElementArray.push_back(make_pair(factor,chname));
}
vector<Double_t> TaDefinition::GetFactorArray(){
  vector<Double_t> fRet;
  auto iter_ele=fRawElementArray.begin();
  while(iter_ele!=fRawElementArray.end()){
    fRet.push_back( (*iter_ele).first);
    iter_ele++;
  }
  return fRet;
}
vector<TString> TaDefinition::GetRawChannelList(){
  vector<TString> fRet;
  auto iter_ele=fRawElementArray.begin();
  while(iter_ele!=fRawElementArray.end()){
    fRet.push_back( (*iter_ele).second);
    iter_ele++;
  }
  return fRet;
}

