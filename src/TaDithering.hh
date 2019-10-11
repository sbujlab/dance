#ifndef __TaDithering_hh__
#define __TaDithering_hh__
#include "TTree.h"
#include "TMatrixD.h"
#include "TCut.h"
#include <vector>
#include "TaSuperCycle.hh"
#include "TaConfig.hh"
#include "TaInput.hh"
#include <vector>
using namespace std;
typedef vector<vector<Double_t> > Vec2D;
class TaDithering: public TObject{
 public:
  TaDithering(TaInput *aInput);
  ~TaDithering();

  Bool_t LoadConfig(TaConfig *fConfig);
  Bool_t LoadModulationData();
  Bool_t ComputeSensitivities();
  void PrintSummary();
  TMatrixD GetDetSensMatrix();
  TMatrixD GetMonSensMatrix();

  inline TTree* GetSensTree() const {return sens_tree;};

 private:
  TTree* sens_tree;
  TaInput *fInput;
  TCut bmod_cut;
  vector< TaSuperCycle> fSuperCycleArray;
  vector<TString> det_array;
  vector<TString> mon_array;
  vector<TString> coil_array;
  vector<TString> alias_array;
  vector<Int_t> coil_index;

  Int_t nDet;
  Int_t nMon;
  Int_t nCoil;
  Int_t nCombo;

  Vec2D DetSens; //[idet][icoil]
  Vec2D MonSens; //[imon][icoil]
  Vec2D DetSens_err; //[idet][icoil]
  Vec2D MonSens_err; //[imon][icoil]
  
  ClassDef(TaDithering,0);
};

#endif
