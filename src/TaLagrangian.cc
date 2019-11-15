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
  fRawDVlist = aConfig->GetRawDVlist();
  fChannelDefinition = aConfig->GetChannelDefintion();
  LoadConstraint(ana_index,aConfig);
}

TMatrixD TaLagrangian::GetDetMonCovMatrix(Int_t imini){
  Int_t nRawDV = fRawDVlist.size();
  Int_t nIV = fIndependentVar.size();

  TMatrixD retMatrix(nIV,nRawDV);
  for(int irow=0;irow<nIV;irow++)
    for(int icol=0;icol<nRawDV;icol++)
      retMatrix[irow][icol]=GetCovariance(fDVMaps[fRawDVlist[icol]],
					  fIndependentVar[irow],
					  imini);
#ifdef NOISY
  cout << __PRETTY_FUNCTION__<< endl;
  retMatrix.Print();
#endif
  return retMatrix;

}

void TaLagrangian::LoadConstraint(Int_t ana_index, TaConfig *aConfig){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  TString ext_filename = aConfig->GetAnalysisParameter(ana_index,"sens_input");
  TFile ext_file(ext_filename);
  if(!ext_file.IsOpen()){
    cout << ext_filename <<" doesn't exist!" << endl;
    return;
  }
  vector<TString> sCoilList = aConfig->GetCoilList(ana_index);
  
  vector<TString> *ext_dv_array =
    (vector<TString>*)ext_file.Get("dv_array");
  vector<TString> *ext_coil_array =
    (vector<TString>*)ext_file.Get("coil_array");

  TMatrixD* sens_matrix = (TMatrixD*) ext_file.Get("sens_matrix"); 
  nCoil = sCoilList.size();
  Int_t nDet = fRawDVlist.size();
  Int_t nMon = sIVlist.size();

  detConstraints.ResizeTo(nCoil,nDet);
  monConstraints.ResizeTo(nCoil,nMon);

  for(int icoil=0;icoil<nCoil;icoil++){
    auto iter_coil = find((*ext_coil_array).begin(),(*ext_coil_array).end(),
			  sCoilList[icoil]);
    if(iter_coil==(*ext_coil_array).end())
      cerr<< sCoilList[icoil] << " not found in matrix files." << endl;
    Int_t index_coil = iter_coil-(*ext_coil_array).begin();
    for(int idet=0;idet<nDet;idet++){
      auto iter_det = find((*ext_dv_array).begin(),(*ext_dv_array).end(),
			   GetBaseName(fRawDVlist[idet]));
      if(iter_det==(*ext_dv_array).end())
	cerr<< fRawDVlist[idet] << " not found in matrix files." << endl;
      Int_t index_det = iter_det-(*ext_dv_array).begin();
      detConstraints(icoil,idet)=(*sens_matrix)[index_det][index_coil];
    }

    for(int imon=0;imon<nMon;imon++){
      auto iter_mon = find((*ext_dv_array).begin(),(*ext_dv_array).end(),
			   GetBaseName(sIVlist[imon]));
      if(iter_mon==(*ext_dv_array).end())
	cerr<< sIVlist[imon] << " not found in matrix files." << endl;
      Int_t index_mon = iter_mon-(*ext_dv_array).begin();
      monConstraints(icoil,imon)=(*sens_matrix)[index_mon][index_coil];
    }
  }
#ifdef NOISY 
  monConstraints.Print();
  detConstraints.Print();
#endif  
  ext_file.Close();
}

vector<vector<Double_t> >  TaLagrangian::Solve(TMatrixD CovDM, TMatrixD CovMM){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  Int_t nMon = CovMM.GetNrows();
  Int_t nDet = CovDM.GetNcols();
  TMatrixD lhsM(nMon+nCoil,nMon+nCoil);    
  TMatrixD rhsM(nMon+nCoil,nDet);    
  TMatrixD monConstraintsTran(nMon,nCoil);
  monConstraintsTran.Transpose(monConstraints);
  lhsM.SetSub(0,nMon,monConstraintsTran);
  lhsM.SetSub(nMon,0,monConstraints);
  rhsM.SetSub(nMon,0,detConstraints);
  lhsM.SetSub(0,0,CovMM);
  rhsM.SetSub(0,0,CovDM);

  TMatrixD invlhsM = lhsM.Invert();
  TMatrixD solutionM = invlhsM*rhsM;
  TMatrixD slopeM(nMon,nDet);
  slopeM=solutionM.GetSub(0,nMon-1,0,nDet-1);
#ifdef DEBUG	
  cout << " -- Slopes Matrix (nMon x nDet) " << endl;
  slopeM.Print();
#endif
  vector<vector<Double_t> > fSlopesContainer;
  
  auto iter_dv = sDVlist.begin();
  while(iter_dv!=sDVlist.end()){
    auto iter_find = find(fRawDVlist.begin(),fRawDVlist.end(),
			  *iter_dv);
    if(iter_find!=fRawDVlist.end()){
      Int_t index = iter_find - fRawDVlist.begin();
      fSlopesContainer.push_back(GetColumnVector(slopeM,index));
    }else{
      auto iter_map = fChannelDefinition.find(*iter_dv);
      if(iter_map == fChannelDefinition.end()){
	cerr<< " -- Error: " <<  *iter_dv 
	    << "'s definition is not found " << endl;
	vector<Double_t> fZero(nMon,0.0);
	fSlopesContainer.push_back(fZero);
      }else{
	vector< pair<Double_t,TString> >  myDefArray = fChannelDefinition[*iter_dv];
	vector<Double_t> mySlope(nMon,0.0);
	auto iter_def = myDefArray.begin();
	while(iter_def!=myDefArray.end()){
	  Double_t weight =  (*iter_def).first;
	  TString element_name =  (*iter_def).second;
	  Int_t index = FindRawDVIndexFromList(element_name);
	  vector<Double_t> fElementSlope(nMon,0.0);
	  if(index!=-1)
	    fElementSlope = GetColumnVector(slopeM,index);
	  for(int imon=0;imon<nMon;imon++)
	    mySlope[imon]+=weight*fElementSlope[imon];
	  iter_def ++;
	} // end of loop over combiner element
	fSlopesContainer.push_back(mySlope);
      } // end of test if definition is found 
    } // end test if it is combined channel
    iter_dv++;
  } // end of dv loop
  return fSlopesContainer;
}

TString TaLagrangian::GetBaseName(TString in){
  in.ReplaceAll("diff_","");
  in.ReplaceAll("asym_","");
  return in;
}

Int_t TaLagrangian::FindRawDVIndexFromList(TString raw_name){
  Int_t index = -1;
  auto iter_find = find(fRawDVlist.begin(),fRawDVlist.end(),
			raw_name);
  if(iter_find!=fRawDVlist.end())
    index = iter_find - fRawDVlist.begin();

  return index;
}
