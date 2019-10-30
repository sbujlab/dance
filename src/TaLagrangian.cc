#include "TaLagrangian.hh"
#include "TBranch.h"
#include "TLeaf.h"
#include "TEventList.h"
#include "TStopwatch.h"
ClassImp(TaLagrangian);

TaLagrangian::TaLagrangian(Int_t ana_index,TaConfig *aConfig){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  Init(ana_index,aConfig);
  LoadConstraint(ana_index,aConfig);
}

void TaLagrangian::LoadConstraint(Int_t ana_index, TaConfig *aConfig){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  TFile ext_file(aConfig->GetAnalysisParameter(ana_index,"sens_input"));
  vector<TString> sCoilList = aConfig->GetCoilList(ana_index);
  
  vector<TString> *ext_dv_array =
    (vector<TString>*)ext_file.Get("dv_array");
  vector<TString> *ext_coil_array =
    (vector<TString>*)ext_file.Get("coil_array");
  
  TMatrixD* sens_matrix = (TMatrixD*) ext_file.Get("sens_matrix"); 
  
  nCoil = sCoilList.size();
  Int_t nDet = sDVlist.size();
  Int_t nMon = sIVlist.size();

  detConstraints.ResizeTo(nCoil,nDet);
  monConstraints.ResizeTo(nCoil,nMon);

  for(int icoil=0;icoil<nCoil;icoil++){
    auto iter_coil = find((*ext_coil_array).begin(),(*ext_coil_array).end(),
			  GetBaseName(sCoilList[icoil]));
    Int_t index_coil = iter_coil-(*ext_coil_array).begin();
    for(int idet=0;idet<nDet;idet++){
      auto iter_det = find((*ext_dv_array).begin(),(*ext_dv_array).end(),
			   GetBaseName(sDVlist[idet]));
      Int_t index_det = iter_det-(*ext_dv_array).begin();
      detConstraints(icoil,idet)=(*sens_matrix)[index_det][index_coil];
    }

    for(int imon=0;imon<nMon;imon++){
      auto iter_mon = find((*ext_dv_array).begin(),(*ext_dv_array).end(),
			   GetBaseName(sIVlist[imon]));
      Int_t index_mon = iter_mon-(*ext_dv_array).begin();
      monConstraints(icoil,imon)=(*sens_matrix)[index_mon][index_coil];
    }
  }
#ifdef DEBUG 
  monConstraints.Print();
  detConstraints.Print();
#endif  
  ext_file.Close();
}

TMatrixD TaLagrangian::Solve(TMatrixD CovDM, TMatrixD CovMM){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  Int_t nMon = CovMM.GetNrows();
  Int_t nDet = CovDM.GetNcols();
  TMatrixD lhsM(nMon+nCoil,nMon+nCoil);    
  TMatrixD rhsM(nMon+nCoil,nDet);    
  TMatrixD monConstraintsTran;
  monConstraintsTran.Transpose(monConstraints);
  lhsM.SetSub(0,nMon,monConstraintsTran);
  lhsM.SetSub(nMon,0,monConstraints);
  rhsM.SetSub(nMon,0,detConstraints);
  lhsM.SetSub(0,0,CovMM);
  rhsM.SetSub(0,0,CovDM);

  TMatrixD invlhsM = lhsM.Invert();
  TMatrixD solutionM = invlhsM*rhsM;
#ifdef NOISY	
  cout << " -- Slopes Matrix (nMon x nDet) " << endl;
  TMatrixD slopeM(nMon,nDet);
  slopeM=solutionM.GetSub(0,nMon-1,0,nDet-1);
  slopeM.Print();
#endif

  // #ifdef NOISY
  //       cout << " -- Check Slope-Sensitivity Consistency " << endl;
  //       TMatrixD ZeroM = monCtrans*slopeM- detCtrans;
  //       ZeroM.Print();
  // #endif 
  return slopeM;
}

TString TaLagrangian::GetBaseName(TString in){
  in.ReplaceAll("diff_","");
  in.ReplaceAll("asym_","");
}
