#include "TaDataElement.hh"
#include <iostream>
#include "TMath.h"
ClassImp(TaDataElement)

using namespace std;

void TaDataElement::RegisterBranchAddress(TBranch* fBranch_ptr){

  if(fBranch_ptr->GetLeaf("hw_sum")!=NULL){ // a VQWK channel
#ifdef DEBUG
    cout << "Register branch address(VQWK)  " << fBranch_ptr->GetName() << endl;
#endif
    fBranch_ptr->GetLeaf("hw_sum")->SetAddress(&fhw_sum);
  }
  else if (fBranch_ptr->GetLeaf("value")!=NULL){  // a scaler channel
#ifdef DEBUG
    cout << "Register branch address(scaler)  " << fBranch_ptr->GetName() << endl;
#endif
    fBranch_ptr->GetLeaf("value")->SetAddress(&fhw_sum);
  }
  fBranch_ptr->GetLeaf("Device_Error_Code")->SetAddress(&fDevice_Error_Code);
}

Double_t TaDataElement::GetHwSum(){
  Int_t nele = fElementArray.size();
  if(nele>=1){
    fhw_sum =0.0;
    for(int i=0;i<nele;i++){
      double factor = fElementArray[i].first;
      double value = fElementArray[i].second->GetHwSum();
      fhw_sum+=factor*value;
    }
  }

  if(myName=="bmod_trim7")
    return 232.1*sin((fhw_sum+61.94)/1361*2*TMath::Pi())+1693;
  else
    return fhw_sum;
}

Double_t TaDataElement::TestDeviceErrorCode(Int_t ErrorMask){
  Int_t myDEC = (Int_t)GetDeviceErrorCode();
  Double_t ret = (Double_t) ( myDEC&ErrorMask);
  return ret;
}

Double_t TaDataElement::GetDeviceErrorCode(){
  Int_t nele = fElementArray.size();
  if(nele>=1){
    fDevice_Error_Code=0.0;
    Int_t ret_error_code=0;
    for(int i=0;i<nele;i++){
      ret_error_code|= (Int_t)(fElementArray[i].second->GetDeviceErrorCode());
    }
    fDevice_Error_Code = ret_error_code;
  }
  return fDevice_Error_Code;
}

void TaDataElement::AddElement(Double_t weight, TaDataElement* felement){
  fElementArray.push_back(make_pair(weight,felement));
}



