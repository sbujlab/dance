#include "TaAccumulator.hh"
#include "TMath.h"
ClassImp(TaAccumulator);

TaAccumulator::TaAccumulator(){
  M2=0;
  mu1=0;
  mu2=0;
  n=0;
}

TaAccumulator::~TaAccumulator(){
}

void TaAccumulator::Update(Double_t s, Double_t t){
  n++;
  M2+= (n-1)/n*(s-mu1)*(t-mu2);
  mu1 += (s-mu1)/n;
  mu2 += (t-mu2)/n;
}

void TaAccumulator::Update(Double_t f1){
  Update(f1,f1);
}

void TaAccumulator::Merge(TaAccumulator other){
  Double_t delta_mu1 = other.GetMean1()-mu1;
  Double_t delta_mu2 = other.GetMean2()-mu2;
  Double_t n1 = n ;
  Double_t n2 = other.GetN() ;
  n = n1+n2;
  M2 += other.GetM2() + n1*n2/n*delta_mu1*delta_mu2;
  mu1 += n2*delta_mu1/n;
  mu2 += n2*delta_mu2/n;
}

void TaAccumulator::Zero(){
  M2=0;
  mu1=0;
  mu2=0;
  n=0;
}

void TaAccumulator::UpdateStat(STAT &input){
  input.mean = mu1;
  input.err = TMath::Sqrt(M2)/n;
  input.rms = TMath::Sqrt(M2/n);
  input.m2 = M2;
  input.num_samples = n;
}

