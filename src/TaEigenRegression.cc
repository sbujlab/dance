#include "TaEigenRegression.hh"
#include "TBranch.h"
#include "TLeaf.h"
#include "TEventList.h"
#include "TStopwatch.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TMath.h"
ClassImp(TaEigenRegression);

TaEigenRegression::TaEigenRegression(Int_t ana_index, TaConfig *fConfig){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  Init(ana_index,fConfig);
}

void TaEigenRegression::Process(TaOutput *fOutput){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif

  ConstructOutputs(fOutput);
  Int_t nMini =  minirun_range.size();
  const Int_t nDV = sDVlist.size();
  const Int_t nIV = sIVlist.size();

  vector<Double_t> fvector(nIV,0.0); // place-holder
  vector<Double_t> fvector_truncated(nIV,0.0); // place-holder
  Double_t coeff[nDV][nIV];
  TString branch_desc = Form("coeff[%d][%d]/D",nDV,nIV);
  fOutput->ConstructTreeBranch("mini_"+tree_name,"coeff",branch_desc,coeff);

  Double_t eigvec[nIV][nIV];
  branch_desc = Form("eigv[%d][%d]/D",nIV,nIV);
  fOutput->ConstructTreeBranch("mini_"+tree_name,"eigvec",branch_desc,eigvec);

  TString leaflist = "hw_sum/D:block0:block1:block2:block3";
  for(int ich=0;ich<nIV;ich++){
    TaChannel* aEigenVar = new TaChannel(tree_name,Form("diff_evMon%d",ich));
    fEigenVar.push_back(aEigenVar);
  }
  for(int ich=0;ich<nIV;ich++){
    if(!kOutputMiniOnly)
      fEigenVar[ich]->ConstructTreeBranch(fOutput,leaflist);
    fEigenVar[ich]->ConstructMiniTreeBranch(fOutput,"mini_"+tree_name);
    fEigenVar[ich]->ConstructSumTreeBranch(fOutput,"sum_"+tree_name);
  }
  for(int ich=0;ich<nIV;ich++){
    fEigenVar[ich]->ConnectChannels(fIndependentVar,fvector);
    fEigenVar[ich]->ConstructSlopeBranch(fOutput,"mini_"+tree_name); // A very nice trick
  }
  // >>> prepare for truncated corrections
  for(int ich=0;ich<nDV;ich++){
    TString myName  = sDVlist[ich];
    TaChannel *aOutputChannel = new TaChannel(tree_name+"_tr",branch_prefix+myName);
    TaChannel *aCorrection = new TaChannel(tree_name+"_tr","cor_"+myName);
    fOutputChannels_Truncated.push_back(aOutputChannel);
    fCorrection_Truncated.push_back(aCorrection);
  }
  for(int ich=0;ich<nDV;ich++)
    fOutputChannels_Truncated[ich]->DefineSubtraction(fDependentVar[ich],fCorrection_Truncated[ich]);

  for(int ich=0;ich<nDV;ich++){
    if(!kOutputMiniOnly){
      fOutputChannels_Truncated[ich]->ConstructTreeBranch(fOutput,leaflist);
      fCorrection_Truncated[ich]->ConstructTreeBranch(fOutput,leaflist);
    }
    fOutputChannels_Truncated[ich]->ConstructMiniTreeBranch(fOutput,"mini_"+tree_name+"_tr");
    fOutputChannels_Truncated[ich]->ConstructSumTreeBranch(fOutput,"sum_"+tree_name+"_tr");
  }

  // done with truncation setup <<< 

  for(int ich=0;ich<nDV;ich++){
    fCorrections[ich]->ConnectChannels(fEigenVar,fvector);
    fCorrections[ich]->ConstructSlopeBranch(fOutput,"mini_"+tree_name);
  }
  for(int ich=0;ich<nDV;ich++){
    fCorrection_Truncated[ich]->ConnectChannels(fEigenVar,fvector_truncated);
    fCorrection_Truncated[ich]->ConstructSlopeBranch(fOutput,"mini_"+tree_name+"_tr");
  }

  for(int imini=0; imini<nMini;imini++){
    TMatrixD CovDM = GetDetMonCovMatrix(imini);
    TMatrixD CovMM = GetMonMonCovMatrix(imini);
    TMatrixDSym symCovMM(nIV);
    for(int i=0;i<nIV;i++)
      for(int j=0;j<nIV;j++)
	symCovMM[i][j] = CovMM[i][j];
    TMatrixDSymEigen symCovMM_eig(symCovMM);
    TMatrixD eigen_vector = symCovMM_eig.GetEigenVectors();
    TVectorD eigen_values = symCovMM_eig.GetEigenValues();
#ifdef NOISY
    eigen_vector.Print();
#endif
    TMatrixD lambda(nIV,nIV);
    for(int i=0;i<nIV;i++)
      for(int j=0;j<nIV;j++)
	if(i==j)
	  lambda[i][j] = eigen_values[i];
	else
    	  lambda[i][j] = 0.0;
    TMatrixD CovDM_eigen(CovDM);
    TMatrixD eigen_vector_trans(eigen_vector);
    eigen_vector_trans.T();
    CovDM_eigen = eigen_vector_trans*CovDM;
    // vector<vector<Double_t> > fSlopes = Solve(CovDM_eigen,lambda);
    vector<vector<Double_t> > fSlopes_raw = Solve(CovDM,CovMM);
    vector<vector<Double_t> > fSlopes = RotateSlope(fSlopes_raw, eigen_vector_trans);
    
#ifdef NOISY
    lambda.Print();
#endif
    for(int idv=0;idv<nDV;idv++)
      for(int iiv=0;iiv<nIV;iiv++)
	coeff[idv][iiv] = fSlopes[idv][iiv];

    for(int ich=0;ich<nIV;ich++){
      vector<Double_t> fEigenVector(nIV);
      for(int icol=0;icol<nIV;icol++)
	fEigenVector[icol] = eigen_vector[icol][ich];
      fEigenVar[ich]->ConnectChannels(fIndependentVar,fEigenVector);
    }

    for(int ich=0;ich<nDV;ich++){
      vector<Double_t> fprefactors = fSlopes[ich];
      fCorrections[ich]->ConnectChannels(fEigenVar,fprefactors);
    }
    // >>>  Rank and Truncate Correction; 
    for(int ich=0;ich<nDV;ich++){
      vector<Double_t> fprefactors = fSlopes[ich];
      // sort and select the top five by correction size 
      vector<Double_t> fProduct(nIV); 
      vector<Int_t > fRank(nIV);
      for(int iiv=0;iiv<nIV;iiv++){
	fProduct[iiv] = TMath::Abs(fSlopes[ich][iiv]) * TMath::Sqrt(eigen_values[iiv]);
	fRank[iiv] = iiv;
      }
      for(int i=1;i<nIV;i++){
	for(int j=i-1;j>=0;j--){
	  if(fProduct[j]< fProduct[j+1] ){
	    Double_t swap = fProduct[j];
	    fProduct[j] = fProduct[j+1];
	    fProduct[j+1] = swap;
	    Int_t idx_swap = fRank[j];
	    fRank[j] = fRank[j+1];
	    fRank[j+1] = idx_swap;
	  }
	}
      }
      for(int i=5; i<nIV;i++) // only keep top five
	fprefactors[ fRank[i] ] = 0.0;
      fCorrection_Truncated[ich]->ConnectChannels(fEigenVar,fprefactors);
    }

    int istart = minirun_range[imini].first;
    int iend = minirun_range[imini].second;
    for(int ievt=istart;ievt<=iend;ievt++){
      GetEntry(ievt);
      for(int ich=0;ich<nIV;ich++)
	fEigenVar[ich]->CalcCombination();
      for(int ich=0;ich<nDV;ich++){
	fCorrection_Truncated[ich]->CalcCombination();
	fOutputChannels_Truncated[ich]->CalcCombination();
      }
      CalcCombination();
      if((fChannelCutFlag->fOutputValue.hw_sum)==1){
	for(int ich=0;ich<nIV;ich++){
	  fEigenVar[ich]->AccumulateMiniSum();
	  fEigenVar[ich]->AccumulateRunSum();
	}
	for(int ich=0;ich<nDV;ich++){
	  fCorrection_Truncated[ich]->AccumulateRunSum();
	  fOutputChannels_Truncated[ich]->AccumulateRunSum();
	  fCorrection_Truncated[ich]->AccumulateMiniSum();
	  fOutputChannels_Truncated[ich]->AccumulateMiniSum();
	}
      }

      AccumulateMiniSum();
      AccumulateRunSum();
      if(!kOutputMiniOnly){
	fOutput->FillTree(tree_name);
	fOutput->FillTree(tree_name+"_tr");
      }
    } // end of  event loop
    UpdateMiniStat();
    for(int ich=0;ich<nIV;ich++)
      fEigenVar[ich]->UpdateMiniStat();
    for(int ich=0;ich<nDV;ich++){
      fCorrection_Truncated[ich]->UpdateMiniStat();
      fOutputChannels_Truncated[ich]->UpdateMiniStat();
    }
    
    fOutput->FillTree("mini_"+tree_name);
    fOutput->FillTree("mini_"+tree_name+"_tr");
    ResetMiniAccumulator();
    for(int ich=0;ich<nIV;ich++)
      fEigenVar[ich]->ResetMiniAccumulator();
    for(int ich=0;ich<nDV;ich++){
      fCorrection_Truncated[ich]->ResetMiniAccumulator();
      fOutputChannels_Truncated[ich]->ResetMiniAccumulator();
    }
    
  } // end of mini run loop
  UpdateRunStat();
  for(int ich=0;ich<nIV;ich++)
    fEigenVar[ich]->UpdateRunStat();
  for(int ich=0;ich<nDV;ich++){
    fCorrection_Truncated[ich]->UpdateRunStat();
    fOutputChannels_Truncated[ich]->UpdateRunStat();
  }
    
  fOutput->FillTree("sum_"+tree_name);
  fOutput->FillTree("sum_"+tree_name+"_tr");
}


vector<vector<Double_t> > TaEigenRegression::Solve(TMatrixD CovDM, TMatrixD CovMM){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif

  TMatrixD lhsM(CovMM);
  TMatrixD rhsM(CovDM);
    
  TMatrixD invlhsM = lhsM.Invert();
  TMatrixD solutionM = invlhsM*rhsM;
#ifdef NOISY	
  cout << " -- Slopes Matrix (nMon x nDet) " << endl;
  solutionM.Print();
#endif
  vector<vector<Double_t> >  fSlopesContainer;
  Int_t ncol = solutionM.GetNcols();
  for(Int_t icol=0;icol<ncol;icol++)
    fSlopesContainer.push_back(GetColumnVector(solutionM,icol));
  return fSlopesContainer;
}


void TaEigenRegression::WriteSummary(){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif

}

Double_t TaEigenRegression::GetCovariance(TaChannel *c1, TaChannel *c2,Int_t imini){
  Double_t cov;
  Int_t start_evt = minirun_range[imini].first;
  Int_t end_evt = minirun_range[imini].second;
  TaAccumulator fCovAccumulator;
  
  for(int ievt=start_evt;ievt<=end_evt;ievt++){
    if(fChannelCutFlag->GetEntry(ievt))
      fCovAccumulator.Update(c1->GetEntry(ievt),c2->GetEntry(ievt));
  }
  cov = fCovAccumulator.GetM2()/fCovAccumulator.GetN();

  return cov;
}

TMatrixD TaEigenRegression::GetDetMonCovMatrix(Int_t imini){
  Int_t nDV = fDependentVar.size();
  Int_t nIV = fIndependentVar.size();

  TMatrixD retMatrix(nIV,nDV);
  for(int irow=0;irow<nIV;irow++)
    for(int icol=0;icol<nDV;icol++)
      retMatrix[irow][icol]=GetCovariance(fDependentVar[icol],
					  fIndependentVar[irow],
					  imini);
#ifdef NOISY
  cout << __PRETTY_FUNCTION__<< endl;
  retMatrix.Print();
#endif
  return retMatrix;
}

TMatrixD TaEigenRegression::GetMonMonCovMatrix(Int_t imini){
  Int_t nIV = fIndependentVar.size();
  TMatrixD retMatrix(nIV,nIV);

  for(int irow=0;irow<nIV;irow++)
    for(int icol=irow;icol<nIV;icol++){
      retMatrix[irow][icol]=GetCovariance(fIndependentVar[icol],
					  fIndependentVar[irow],
					  imini);
      if(icol!=irow)
      	retMatrix[icol][irow] = retMatrix[irow][icol];
    }
#ifdef NOISY
  cout << __PRETTY_FUNCTION__<< endl;
  retMatrix.Print();
#endif
  return retMatrix;
}

vector<Double_t> TaEigenRegression::GetColumnVector(TMatrixD matrix_in, Int_t icol ){
  Int_t nrow = matrix_in.GetNrows();
  vector<Double_t> ret_vector;
  for(int irow=0;irow<nrow;irow++)
    ret_vector.push_back(matrix_in(irow,icol));
  return ret_vector;
}

vector<Double_t> TaEigenRegression::GetRowVector(TMatrixD matrix_in, Int_t irow ){
  Int_t ncol = matrix_in.GetNcols();
  vector<Double_t> ret_vector;
  for(int icol=0;icol<ncol;icol++)
    ret_vector.push_back(matrix_in(irow,icol));
  return ret_vector;
}


vector<vector<Double_t> > TaEigenRegression::RotateSlope(vector<vector<Double_t> > input,
							 TMatrixD Q_T){
  vector<vector<Double_t> > fRet;
  Int_t nIV = sIVlist.size();
  Int_t nDV = sDVlist.size();
  TMatrixD slope_m(nIV,nDV);
  TMatrixD slope_rot(nIV,nDV);
  for(int iiv=0;iiv<nIV;iiv++)
    for(int idv=0;idv<nDV;idv++)
      slope_m[iiv][idv]  =  input[idv][iiv];

  slope_rot = Q_T*slope_m;

  for(int idv=0;idv<nDV;idv++){
    vector<Double_t> ftemp;
    for(int iiv=0;iiv<nIV;iiv++)
      ftemp.push_back(slope_rot[iiv][idv]);
    fRet.push_back(ftemp);
  }
  return fRet;
}
