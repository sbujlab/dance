#ifndef __VANALYSIS_MODULE__
#define __VANALYSIS_MODULE__

#include "TaChannel.hh"
#include "TaConfig.hh"
#include "TaInput.hh"
#include "TaOutput.hh"
#include "TaAccumulator.hh"

using namespace std;
class TaInput;
class VAnalysisModule{
public:
  VAnalysisModule();
  virtual ~VAnalysisModule(){};
  virtual void LoadInput(TaInput *aInput);
  virtual void Process(TaOutput *aOutput)=0;
  virtual void End(){};
  void Init(Int_t index, TaConfig *aConfig);
  void ConstructOutputs(TaOutput *aOutput);
  void GetEntry(Int_t ievt);
  void CalcCombination();

  void AccumulateMiniSum();
  void AccumulateRunSum();
  void ResetMiniAccumulator();
  void ResetRunAccumulator();
  void UpdateMiniStat();
  void UpdateRunStat();

protected:

  TaOutput *fOutput;
  vector<pair<Int_t,Int_t> > minirun_range;
  vector<TaChannel*> fDependentVar;
  vector<TaChannel*> fIndependentVar;
  vector<TaChannel*> fOutputChannels;
  vector<TaChannel*> fCorrections;
  TaChannel* fChannelCutFlag;

  vector<TString> sDVlist;
  vector<TString> sIVlist;
  TString tree_name;
  TString branch_prefix;
  TString myType;

  ClassDef(VAnalysisModule,0);
};

#endif
