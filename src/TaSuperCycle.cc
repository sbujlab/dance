#include "TaAccumulator.hh"
#include "TaSuperCycle.hh"
#include "TMath.h"
#include "TaPrinter.hh"

ClassImp(TaSuperCycle);
using namespace std;
void TaSuperCycle::LoadDetectorList(vector<TString> input){
  fDetectorList = input;
}
void TaSuperCycle::RegisterDependentVarArray(vector<TaDataElement*> input_array){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  fDependentVarArray = input_array;
  nDependentVar = fDependentVarArray.size();
  for(int idv=0;idv<nDependentVar;idv++){
    auto iter=find(fDetectorList.begin(),fDetectorList.end(),
		   fDependentVarArray[idv]->GetName());
    if(iter!=fDetectorList.end())
      isDetectorFlag.push_back(kTRUE);
    else
      isDetectorFlag.push_back(kFALSE);
  }
}

void TaSuperCycle::RegisterCoilArray(vector<TaDataElement*> input_array){
  fCoilArray = input_array;
  nCoil = fCoilArray.size();
}
void TaSuperCycle::InitAccumulators(){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  fCoilVarianceArray.resize(nCoil);
  fCovarianceArray.resize(nDependentVar,fCoilVarianceArray);
  fDepVarianceArray.resize(nDependentVar,fCoilVarianceArray);
}

void TaSuperCycle::UpdateSamples(Int_t cur_index){

  fCoilVarianceArray[cur_index].Update(fCoilArray[cur_index]->GetHwSum() );
  for(int idv=0;idv<nDependentVar;idv++)
    fDepVarianceArray[idv][cur_index].Update( fDependentVarArray[idv]->GetHwSum() );

  for(int idv=0;idv<nDependentVar;idv++)
    fCovarianceArray[idv][cur_index].Update( fDependentVarArray[idv]->GetHwSum(),
					     fCoilArray[cur_index]->GetHwSum());

}
void TaSuperCycle::CalcSensitivities(){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  Int_t icount=0;
  for(int idv=0;idv<nDependentVar;idv++){
    TString dv_name = fDependentVarArray[idv]->GetName();
    for(int icoil=0;icoil<nCoil;icoil++){
      TString coil_name = fCoilArray[icoil]->GetName();
      fSensitivityMap[make_pair(dv_name,coil_name)] = icount++;
      fSamples.push_back(fCovarianceArray[idv][icoil].GetN());
      if(fCovarianceArray[idv][icoil].GetN()<=50 ||
	 fCoilVarianceArray[icoil].GetM2()/fCoilVarianceArray[icoil].GetN()<10 ){
	fSensitivity.push_back(0.0);
	fSensitivity_err.push_back(-1.0);
	continue;
      }
      double numerator= fCovarianceArray[idv][icoil].GetM2();
      double denominator = fCoilVarianceArray[icoil].GetM2();
      double mySensitiviy =  numerator/denominator;
      if(isDetectorFlag[idv]){
	double detector_norm = fDepVarianceArray[idv][icoil].GetMean1();
	if(detector_norm<=0){
	  fSensitivity.push_back(0.0);
	  fSensitivity_err.push_back(-1.0);
	}else{
	  mySensitiviy /= detector_norm;
	  Double_t a = fDepVarianceArray[idv][icoil].GetM2()- TMath::Power(fCovarianceArray[idv][icoil].GetM2(),2)/fCoilVarianceArray[icoil].GetM2();
	  Double_t b = fCoilVarianceArray[icoil].GetM2();
	  Double_t myError = TMath::Sqrt((a/b)/(fCoilVarianceArray[icoil].GetN()-2));
	  myError /= detector_norm;

	  fSensitivity.push_back(mySensitiviy);
	  fSensitivity_err.push_back(myError);
	}
      } else{
	Double_t a = fDepVarianceArray[idv][icoil].GetM2()- TMath::Power(fCovarianceArray[idv][icoil].GetM2(),2)/fCoilVarianceArray[icoil].GetM2();
	Double_t b = fCoilVarianceArray[icoil].GetM2();
	Double_t myError = TMath::Sqrt((a/b)/(fCoilVarianceArray[icoil].GetN()-2));
	fSensitivity.push_back(mySensitiviy);
	fSensitivity_err.push_back(myError);
      }
    } // end of coil loop
  } // end of dependent variables loop

      cout << cycID << ":" << fSamples.size() << endl;
}

void TaSuperCycle::PrintSensitivities(){
  cout << "-- cycle ID : " << cycID << endl;
  cout << "-- Coil # : " ;
  for(int ic=1;ic<=7;ic++)
    cout << ic << "\t";
  cout << endl;
  cout << "-- nSamples : " ;
  for(int ic=0;ic<7;ic++)
    cout << fCoilVarianceArray[ic].GetN() << "\t";
  cout << endl;  
  for(int idv=0;idv<nDependentVar;idv++){
    Double_t scale = 1;
    if(isDetectorFlag[idv])
      scale = 1e6;
    else
      scale = 1e3;
    TString dv_name = fDependentVarArray[idv]->GetName();
    cout << dv_name << ":";
    for(int ic=1;ic<=7;ic++){
      Int_t index = fSensitivityMap[make_pair(dv_name,Form("bmod_trim%d",ic))];
      cout <<fSensitivity[index]*scale
	   <<" +/- "
	   << fSensitivity_err[index]*scale << "\t" ;
    }
    cout << endl;
  }
}

void TaSuperCycle::WriteToPrinter(TaPrinter* aPrinter){
  aPrinter->InsertHorizontalLine();
  vector<TString> header;
  header.push_back(" ");
  for(int ic=1;ic<=7;ic++)
    header.push_back(Form("bmod_trim%d",ic));
  aPrinter->AddHeader(header);
  for(int idv=0;idv<nDependentVar;idv++){
    Double_t scale = 1;
    if(isDetectorFlag[idv])
      scale = 1e6;
    else
      scale = 1e3;
    TString dv_name = fDependentVarArray[idv]->GetName();
    aPrinter->NewLine();
    aPrinter->AddStringEntry(dv_name);
    for(int ic=1;ic<=7;ic++){
      Int_t index = fSensitivityMap[make_pair(dv_name,Form("bmod_trim%d",ic))];
      aPrinter->AddFloatEntryWithError(fSensitivity[index]*scale,
				       fSensitivity_err[index]*scale);
    }
  }
  aPrinter->InsertHorizontalLine();
}

