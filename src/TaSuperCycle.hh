#ifndef __TaSuperCycle_hh__
#define __TaSuperCycle_hh__

#include "Rtypes.h"
#include "TObject.h"
#include "TaConfig.hh"
#include "TaAccumulator.hh"
#include <vector>
using namespace std;

class TaSuperCycle: public TObject{
public:
  TaSuperCycle();
  ~TaSuperCycle();
  // void Init(TaConfig* aConfig);
  void LoadCoilData(Int_t, Double_t);
  void LoadDetData(Int_t,Int_t,Double_t,Double_t);
  void LoadMonData(Int_t,Int_t,Double_t,Double_t);
  void Resize(Int_t ndet, Int_t nmon, Int_t ncoil);
  void ComputeSensitivities();
  inline void SetCycleID(Int_t id){ cycID = id;};
  inline Int_t GetCycleID(){ return cycID;};
  inline Double_t GetMonSens(Int_t imon,Int_t ic) const {
    return monSens[imon][ic];
  };
  inline Double_t GetDetSens(Int_t idet,Int_t ic) const {
    return detSens[idet][ic];
  };
  inline Double_t GetMonSens_err(Int_t imon,Int_t ic) const {
    return monSens_err[imon][ic];
  };
  inline Double_t GetDetSens_err(Int_t idet,Int_t ic) const {
    return detSens_err[idet][ic];
  };
  inline Int_t GetNSamples(Int_t i){ return nSamples[i];};
private:
  Int_t nDet;
  Int_t nMon;
  Int_t nCoil;
  vector<Int_t> nSamples;
  // vector<Int_t> coil_index;
  vector<vector<Double_t> > detSens; //[idet][icoil]
  vector<vector<Double_t> > detSens_err;
  vector<vector<Double_t> > monSens; //[imon][icoil]
  vector<vector<Double_t> > monSens_err;

  vector<AccVector> DetM2; // [idet][icoil]
  vector<AccVector> MonM2; // [imon][icoil]
  vector<AccVector> CovDetCoil;
  vector<AccVector> CovMonCoil;
  AccVector CoilM2; //[icoil]

  Int_t cycID;
  ClassDef(TaSuperCycle,0);
};

#endif
