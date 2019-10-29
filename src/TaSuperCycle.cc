#include "TaAccumulator.hh"
#include "TaSuperCycle.hh"
#include "TaPrinter.hh"

#include "TMath.h"
#include "TDecompLU.h"

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
      kScheme.push_back(fConfig->GetAnalysisParameter(ana_index,"scheme"));
      fcoil_list.push_back(fConfig->GetCoilList(ana_index));
      fmonitor_list.push_back(fConfig->GetMonitorList(ana_index));
      slope_tree_name.push_back(fConfig->GetAnalysisParameter(ana_index,"tree_name"));
    }
    iter_anatype++;
  }
}

void TaSuperCycle::CalcSlopes(){
  Int_t nMod = kScheme.size();
  for(int imod =0;imod<nMod;imod++){
    TMatrixD mrhs, mlhs;
    Bool_t kComplete1 = MakeMatrixFromList( fmonitor_list[imod],fcoil_list[imod], mrhs );
    Bool_t kComplete2 = MakeMatrixFromList( fDetectorList, fcoil_list[imod], mlhs );
    Bool_t isGood=kFALSE;
    TMatrixD sol(fDetectorList.size(), fmonitor_list[imod].size());
    if(kComplete2 && kComplete1)
      isGood = GetMatrixSolution( mlhs , mrhs, sol);
    isGoodSlopes.push_back(isGood);
    fSolutionArray.push_back(sol);
  }
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
      else
	return kFALSE;
    }
  }
  return isComplete;
}

Bool_t TaSuperCycle::GetMatrixSolution(TMatrixD lhs, TMatrixD rhs,
				       TMatrixD &sol){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__<< endl;
#endif
  Bool_t isGood;
  TMatrixD mX(rhs);
  TMatrixD mXT = rhs.T();
  TMatrixD mXsq = mX*mXT;
  TDecompLU lu(mXsq);
  if(!lu.Decompose()){
    cout << "LU Decomposition failed, rhs matrix may be singular" << endl;
    cout << "Fail to solve and Zero out solution" << endl;
#ifdef DEBUG  
    mXsq.Print();
#endif
    sol.Zero();
    isGood = kFALSE;
  }else{
    TMatrixD inv = (mXsq).Invert();
    sol = lhs*mXT*inv;
#ifdef NOISY  
    sol.Print();
#endif
    isGood = kTRUE;
  }
  return isGood;
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
  CalcSlopes();
  FillSlopes();
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

void TaSuperCycle::FillSlopes(){
  Int_t nMode = kScheme.size();
  for(int imod=0;imod<nMode;imod++){
    vector< pair<TString, TString> > fDMPairArray;
    vector< Double_t > fSlopeVector;
    vector<TString> monitor_list = fmonitor_list[imod];
    auto iter_mon = monitor_list.begin();
    while(iter_mon!=monitor_list.end()){
      Int_t imon  = iter_mon-monitor_list.begin();
      Int_t myIndex = fDVIndexMapByName[*iter_mon];
      TaDataElement* this_element = fDependentVarArray[ myIndex];
      auto iter_det = fDetectorList.begin();
      while(iter_det!=fDetectorList.end()){
	Int_t idet = iter_det - fDetectorList.begin();
	fDMPairArray.push_back(make_pair(*iter_det,*iter_mon));
	fSlopeVector.push_back(fSolutionArray[imod](idet,imon));
	iter_det++;
      }
      iter_mon++;
    }
    fSlopes.push_back(fSlopeVector);
    fDetMonPairArray.push_back(fDMPairArray);
  }
}

void TaSuperCycle::ConstructSlopeTreeBranch(TaOutput *aOutput,Int_t index,
					    vector<Double_t> &fBranchValues){
  Int_t ndim = fSlopes[index].size();
  fBranchValues.resize(ndim,0);
  
  for(int i=0;i<ndim;i++){
    TString branch_name = Form("%s_%s",
			       (fDetMonPairArray[index][i].first).Data(),
			       (fDetMonPairArray[index][i].second).Data());
    aOutput->ConstructTreeBranch(slope_tree_name[index],branch_name,fBranchValues[i]);
  }
}

void TaSuperCycle::FillSlopeTree(TaOutput *aOutput,Int_t index,
				 vector<Double_t> &fBranchValues){
  Int_t ndim = fSlopes[index].size();
  for(int i=0;i<ndim;i++){
    fBranchValues[i] = fSlopes[index][i];
  }
  aOutput->FillTree(slope_tree_name[index]);
}
