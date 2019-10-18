#ifndef __TaAccumulator_hh__
#define __TaAccumulator_hh__
#include "Rtypes.h"

typedef struct{Double_t mean,err,rms,m2,num_samples;} STAT;

class TaAccumulator{
public:
  TaAccumulator();
  virtual ~TaAccumulator();
  void Update(Double_t s, Double_t t);
  void Update(Double_t ss);
  void UpdateStat(STAT &input);
  void Merge( TaAccumulator);
  void Zero();
  inline Double_t GetM2() const {return M2;};
  inline Double_t GetMean1() const {return mu1;};
  inline Double_t GetMean2() const {return mu2;};
  inline Double_t GetN() const {return n;};
private:
  Double_t M2;
  Double_t mu1;
  Double_t mu2;
  Double_t n;

  ClassDef(TaAccumulator,0);
};
typedef std::vector<TaAccumulator> AccVector;
#endif
