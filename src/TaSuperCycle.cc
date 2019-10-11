#include "TaAccumulator.hh"
#include "TaSuperCycle.hh"
#include "TMath.h"

ClassImp(TaSuperCycle);
using namespace std;

TaSuperCycle::TaSuperCycle(){
}

TaSuperCycle::~TaSuperCycle(){
}

// void TaSuperCycle::Init(TaConfig *aConfig){}
void TaSuperCycle::LoadCoilData(Int_t icoil,Double_t trimcard_value){
  if(icoil<nCoil)
    CoilM2[icoil].Update(trimcard_value);
  else{
    cerr << __FUNCTION__ 
	 << " Error:  coil id over limit " << endl;
  }
}

void TaSuperCycle::LoadDetData(Int_t idet, Int_t icoil,
			       Double_t det_value, Double_t trimcard_value){
  
  if(icoil<nCoil && idet< nDet){
    trimcard_value = trimcard_value;
    CovDetCoil[idet][icoil].Update(det_value,trimcard_value);
    DetM2[idet][icoil].Update(det_value);
  }
  else{
    cerr << __FUNCTION__ 
	 << " Error:  coil or detector id over limit " << endl;
  }
}

void TaSuperCycle::LoadMonData(Int_t imon, Int_t icoil,
			       Double_t mon_value, Double_t trimcard_value){
  if(icoil<nCoil && imon< nMon){
    trimcard_value = trimcard_value;
    CovMonCoil[imon][icoil].Update(mon_value,trimcard_value);
    MonM2[imon][icoil].Update(mon_value);
  }
  else{
    cerr << __FUNCTION__ 
	 << " Error:  coil or monitor id over limit " << endl;
  }

}

void TaSuperCycle::Resize(Int_t ndet, Int_t nmon, Int_t ncoil){
  nDet = ndet;
  nMon = nmon;
  nCoil= ncoil;
  AccVector dummy_vs_coil;
  dummy_vs_coil.resize(ncoil);
  CoilM2.resize(ncoil);
  vector<Double_t> dummy_err;
  dummy_err.resize(ncoil,-1.0);
  vector<Double_t> dummy_zero;
  dummy_zero.resize(ncoil,0.0);

  for(int idet=0;idet<ndet;idet++){
    CovDetCoil.push_back(dummy_vs_coil);
    DetM2.push_back(dummy_vs_coil);
    detSens.push_back(dummy_zero);
    detSens_err.push_back(dummy_err);
  }
  for(int imon=0;imon<nmon;imon++){
    CovMonCoil.push_back(dummy_vs_coil);
    MonM2.push_back(dummy_vs_coil);
    monSens.push_back(dummy_zero);
    monSens_err.push_back(dummy_err);
  }
}

void TaSuperCycle::ComputeSensitivities(){
  
  for(int icoil=0;icoil<nCoil;icoil++){
    // Check if it has enough statistics and
    nSamples.push_back(CoilM2[icoil].GetN());
    if( CoilM2[icoil].GetN()<200)
      continue;
    // Check if coil is activated 
    if(CoilM2[icoil].GetM2()/CoilM2[icoil].GetN()<5)
      continue;
    for(int idet=0;idet<nDet;idet++){
      double numerator= CovDetCoil[idet][icoil].GetM2();
      double denominator =CoilM2[icoil].GetM2();
      detSens[idet][icoil] = 1e6*(numerator/denominator)/TMath::Abs(DetM2[idet][icoil].GetMean1()); 
      // Unit (ppm / trimcard_counts)
      double a = DetM2[idet][icoil].GetM2()-TMath::Power(CovDetCoil[idet][icoil].GetM2(),2)/CoilM2[icoil].GetM2();
      double b = CoilM2[icoil].GetM2();
      detSens_err[idet][icoil] = TMath::Sqrt((a/b)/(CoilM2[icoil].GetN()-2));
      detSens_err[idet][icoil] = 1e6*detSens_err[idet][icoil]/TMath::Abs(DetM2[idet][icoil].GetMean1());
    }
    for(int imon=0;imon<nMon;imon++){
      double numerator= CovMonCoil[imon][icoil].GetM2();
      double denominator =CoilM2[icoil].GetM2();
      monSens[imon][icoil] = 1e3*numerator/denominator;
      // Unit (um / trimcard_counts)
      double a = MonM2[imon][icoil].GetM2()-TMath::Power(CovMonCoil[imon][icoil].GetM2(),2)/CoilM2[icoil].GetM2();
      double b = CoilM2[icoil].GetM2();
      monSens_err[imon][icoil] = 1e3*TMath::Sqrt((a/b)/(CoilM2[icoil].GetN()-2));

    }
  }
}


