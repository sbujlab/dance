#ifndef __TaLagrangian_hh__
#define __TaLagrangian_hh__

#include "TTree.h"
#include "TMatrixD.h"
#include "TCut.h"

#include "TaAccumulator.hh"
#include "TaConfig.hh"
#include "TaInput.hh"

#include <vector>

using namespace std;

class TaLagrangian{
  typedef std::vector<std::vector<Double_t> > Vec2D;
  typedef struct{Double_t mean,error,rms,nsamples;} RUN_STATS;
  typedef std::vector< std::pair<Double_t,Int_t> > AliasMap;
public:
  TaLagrangian(TaInput *aInput);
  virtual ~TaLagrangian();

  void LoadConstraint(TMatrixD,TMatrixD);
  Bool_t ComputeCorrelation();
  void ComputeSlopes();
  void CorrectTree();
  void WriteSummary();
  void UpdateRunStat(RUN_STATS&,TaAccumulator);

  Bool_t LoadConfig(TaConfig *aConfig);  

  inline TTree* GetCorrectionTree() const {return correct_tree;};
  inline TTree* GetSlopeTree() const {return slope_tree;};
  inline TTree* GetMiniTree() const {return mini_tree;};

private:
  TaInput* fInput;
  TTree *mul_tree;
  TCut custom_cut;
  TString output_path;
  Int_t run_number;
  Int_t seg_number;
  Int_t mini_size;
  vector< pair< Int_t ,Int_t> > mini_range;

  vector<TString> det_array;
  vector<TString> mon_array;
  vector<TString> alias_array;
  vector<AliasMap> alias_map_array;
  
  Int_t nDet;
  Int_t nMon;
  Int_t nCoil;
  Int_t nCombo;

  Bool_t kUseConstraint;
  TMatrixD detCmatrix;
  TMatrixD monCmatrix;

  vector<Vec2D> miniCovDetMon; 
  vector<Vec2D> miniCovMonMon; 

  Int_t nGoodPatterns;
  vector<Vec2D> miniDetMonSlopes;
  vector<Vec2D> miniComboSlopes; 

  Vec2D dataRawDet;
  Vec2D dataMon;
  vector<Bool_t> cutFlag;
  vector<Double_t> dataErrorFlag;

  vector<RUN_STATS> statMon;
  vector<RUN_STATS> statRawDet;
  vector<RUN_STATS> statCorDet;
  vector<RUN_STATS> statRawCombo;
  vector<RUN_STATS> statCorCombo;

  AccVector accMon_avg;
  AccVector accRawDet_avg;
  AccVector accCorDet_avg;
  AccVector accRawCombo_avg;
  AccVector accCorCombo_avg;

  TTree *correct_tree;
  TTree *slope_tree;
  TTree *mini_tree;
  ClassDef(TaLagrangian,0);
};
#endif
