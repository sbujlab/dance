#include "TaAccumulator.hh"
#include "TaSuperCycle.hh"
#include "TaPrinter.hh"

#include "TMath.h"
#include "TDecompLU.h"

ClassImp(TaSuperCycle);
using namespace std;
TaSuperCycle::TaSuperCycle(){
  kDeviceErrorCut = kFALSE;
}
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
    fDVIndexMapByName[ (fDependentVarArray[idv]->GetName()) ]=idv;
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

void TaSuperCycle::ConfigSlopesCalculation(TaConfig *fConfig){
  vector<TString> fAnalysisTypeArray = fConfig->GetAnalysisTypeArray();
  auto iter_anatype = fAnalysisTypeArray.begin();
  while(iter_anatype!=fAnalysisTypeArray.end()){
    if((*iter_anatype)=="slope"){
      Int_t ana_index = iter_anatype-fAnalysisTypeArray.begin();
      fcoil_list.push_back(fConfig->GetCoilList(ana_index));
      fmonitor_list.push_back(fConfig->GetMonitorList(ana_index));
      slope_tree_name.push_back(fConfig->GetAnalysisParameter(ana_index,"tree_name"));
    }
    iter_anatype++;
  }
}

void TaSuperCycle::CalcSlopes(){
  Int_t nMod = slope_tree_name.size();
  for(int imod=0;imod<nMod;imod++){
    TMatrixD mrhs;
    Bool_t kComplete = MakeMatrixFromList( fmonitor_list[imod],fcoil_list[imod], mrhs );
    auto iter_det = fDetectorList.begin();
    vector<Double_t> fSlopeVector;
    vector<Double_t> fFlagVector;
    vector<pair<TString,TString> > fKeyVector;
    while(iter_det!=fDetectorList.end()){
      TMatrixD mlhs;
      Bool_t kComplete_det=MakeMatrixByName( *iter_det, fcoil_list[imod], mlhs );
      Int_t nMon = fmonitor_list[imod].size();
      TMatrixD sol(1, nMon);
      Bool_t isGood = kFALSE;
      if(kComplete && kComplete_det)
	isGood = GetMatrixSolution( mlhs , mrhs, sol);
      
      for(int imon=0;imon<nMon;imon++){
	fSlopeVector.push_back(sol[0][imon]);
	fFlagVector.push_back((Double_t)isGood);
	fKeyVector.push_back(make_pair(*iter_det,fmonitor_list[imod][imon]));
      }
      iter_det++;
    } // end of detector loop
    fSlopeContainer.push_back(fSlopeVector);
    fFlagContainer.push_back(fFlagVector);
    fKeyContainer.push_back(fKeyVector);
  } // Slope Modes loops
}

Bool_t TaSuperCycle::MakeMatrixByName(TString channel_in_row,
				      vector<TString> col_list,
				      TMatrixD &input){
  Bool_t isComplete = kTRUE;
  Int_t nrow = 1;
  Int_t ncol = col_list.size();
  input.ResizeTo(nrow,ncol);

  for(int icol=0;icol<ncol;icol++){
    Int_t sens_index = fSensitivityMap[make_pair(channel_in_row,col_list[icol])];
    if(fSensitivity_err[sens_index]>0)
      input(0,icol) = fSensitivity[sens_index];
    else{
      input(0,icol) = 0.0;
      isComplete = kFALSE;
    }
  }
  return isComplete;
}

Bool_t TaSuperCycle::MakeMatrixFromList(vector<TString> row_list,
					vector<TString> col_list,
					TMatrixD &input){
  Bool_t isComplete = kTRUE;
  Int_t nrow = row_list.size();
  Int_t ncol = col_list.size();
  input.ResizeTo(nrow,ncol);
  for(int irow=0;irow<nrow;irow++){
    for(int icol=0;icol<ncol;icol++){
      Int_t sens_index = fSensitivityMap[make_pair(row_list[irow],col_list[icol])];
      if(fSensitivity_err[sens_index]>0)
	input(irow,icol) = fSensitivity[sens_index];
      else{
	input(irow,icol) = 0.0;
	isComplete = kFALSE;
      }
    }
  }
  return isComplete;
}

Bool_t TaSuperCycle::GetMatrixSolution(TMatrixD lhs, TMatrixD rhs,
				       TMatrixD &sol){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__<< endl;
#endif
  TString isGood = kTRUE;   // FIXME : Test singular

  Int_t nrow = rhs.GetNrows();
  Int_t ncol = rhs.GetNcols();
  if(nrow==ncol){
    TMatrixD inv_rhs = rhs.Invert();
    sol = lhs*inv_rhs;
  }else{
    TMatrixD mX(rhs);
    TMatrixD mXT = rhs.T();
    TMatrixD mXsq = mX*mXT;
    TMatrixD inv = (mXsq).Invert();
    sol = lhs*mXT*inv;
  }

#ifdef DEBUG  
  sol.Print();
#endif
  return isGood;
}

void TaSuperCycle::InitAccumulators(){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  vector<TaAccumulator> fDummy(nCoil);
  fCoilVarianceArray.resize(nDependentVar,fDummy);
  fCovarianceArray.resize(nDependentVar,fDummy);
  fDepVarianceArray.resize(nDependentVar,fDummy);
}

void TaSuperCycle::UpdateSamples(Int_t cur_index){

  for(int idv=0;idv<nDependentVar;idv++){

    if(kDeviceErrorCut &&
       ( fDependentVarArray[idv]->TestDeviceErrorCode()
	 || fCoilArray[cur_index]->TestDeviceErrorCode() ) )
      continue;
    
    fDepVarianceArray[idv][cur_index].Update(fDependentVarArray[idv]->GetHwSum() );

    fCoilVarianceArray[idv][cur_index].Update(fCoilArray[cur_index]->GetHwSum() );

    fCovarianceArray[idv][cur_index].Update(fDependentVarArray[idv]->GetHwSum(),
					    fCoilArray[cur_index]->GetHwSum());
  }

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
	 fCoilVarianceArray[idv][icoil].GetM2()/fCoilVarianceArray[idv][icoil].GetN()<10 ){
	fSensitivity.push_back(0.0);
	fSensitivity_err.push_back(-1.0);
	continue;
      }
      double numerator= fCovarianceArray[idv][icoil].GetM2();
      double denominator = fCoilVarianceArray[idv][icoil].GetM2();
      double mySensitiviy =  numerator/denominator;
      if(isDetectorFlag[idv]){
	double detector_norm = fDepVarianceArray[idv][icoil].GetMean1();
	if(detector_norm<=0){
	  fSensitivity.push_back(0.0);
	  fSensitivity_err.push_back(-1.0);
	}else{
	  mySensitiviy /= detector_norm;
	  Double_t a = fDepVarianceArray[idv][icoil].GetM2()- TMath::Power(fCovarianceArray[idv][icoil].GetM2(),2)/fCoilVarianceArray[idv][icoil].GetM2();
	  Double_t b = fCoilVarianceArray[idv][icoil].GetM2();
	  Double_t myError = TMath::Sqrt((a/b)/(fCoilVarianceArray[idv][icoil].GetN()-2));
	  myError /= detector_norm;

	  fSensitivity.push_back(mySensitiviy);
	  fSensitivity_err.push_back(myError);
	}
      } else{
	Double_t a = fDepVarianceArray[idv][icoil].GetM2()- TMath::Power(fCovarianceArray[idv][icoil].GetM2(),2)/fCoilVarianceArray[idv][icoil].GetM2();
	Double_t b = fCoilVarianceArray[idv][icoil].GetM2();
	Double_t myError = TMath::Sqrt((a/b)/(fCoilVarianceArray[idv][icoil].GetN()-2));
	fSensitivity.push_back(mySensitiviy);
	fSensitivity_err.push_back(myError);
      }
    } // end of coil loop
  } // end of dependent variables loop
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

void TaSuperCycle::ConstructSlopeTreeBranch(TaOutput *aOutput,Int_t index,
					    vector<Double_t> &fBranchValues,
					    vector<Double_t> &fFlagValues){
  Int_t ndim = fSlopeContainer[index].size();
  fBranchValues.resize(ndim,0.0);
  fFlagValues.resize(ndim,0.0);
  
  for(int i=0;i<ndim;i++){
    TString branch_name = Form("%s_%s",
			       (fKeyContainer[index][i].first).Data(),
			       (fKeyContainer[index][i].second).Data());
    aOutput->ConstructTreeBranch(slope_tree_name[index],branch_name,fBranchValues[i]);
  }

  for(int i=0;i<ndim;i++){
    TString branch_name = Form("%s_%s_flag",
			       (fKeyContainer[index][i].first).Data(),
			       (fKeyContainer[index][i].second).Data());
    aOutput->ConstructTreeBranch(slope_tree_name[index],branch_name,fFlagValues[i]);
  }
}

void TaSuperCycle::FillSlopeTree(TaOutput *aOutput,Int_t index,
				 vector<Double_t> &fBranchValues,
				 vector<Double_t> &fFlagValues){

  Int_t ndim = fSlopeContainer[index].size();
  for(int i=0;i<ndim;i++){
    fBranchValues[i] = fSlopeContainer[index][i];
    fFlagValues[i] = fFlagContainer[index][i];
  }
  aOutput->FillTree(slope_tree_name[index]);

}
