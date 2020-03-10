#include "TSolver.hh"

ClassImp(TSolver);
TSolver::TSolver(){
}

Bool_t TSolver::SolveMatrix(){
  cout << " -- Using matrix solution " << endl;
  cout << " -- " << fCoilIndex.size() << " coils are loaded " << endl;
  if(!IsCoilSufficient()){
    cerr << " -- Array of Coil is not sufficient " << endl;
    kGoodSolution = kFALSE;
    return kFALSE ;
  }

  
  const Int_t nMon = fMonArray[0].size();
  cout << " -- " << nMon << " monitors are used " << endl;

  const Int_t nDet = 1;
  
  Int_t nCoil = fCoilIndex.size();
  TMatrixD det(nCoil,nDet);
  TMatrixD sol(nMon,nDet);
  TMatrixD mon(nCoil,nMon);
  TMatrixD monTrans(nMon,nCoil);

  for(int i=0;i<nCoil;i++){
    for(int j=0;j<nDet;j++)
      det[i][j]=fDetArray[i];

    for(int j=0;j<nMon;j++){
      mon[i][j]=fMonArray[i][j];
      monTrans[j][i]=fMonArray[i][j];
    }
  }
  TMatrixD invMTM = (monTrans*mon).Invert();
  sol = invMTM*monTrans*det;
  
  for(int i=0;i<nMon;i++)
    fSolution.push_back(sol[i][0]);

  kGoodSolution = kTRUE;
  return kTRUE;
}

void TSolver::LoadSensitivity(Int_t index, Double_t fDet, vector<Double_t> fMon){
  
  fCoilIndex.push_back(index);
  fDetArray.push_back(fDet);
  fMonArray.push_back(fMon);
}

Bool_t TSolver::IsCoilSufficient(){
  Bool_t kSufficient = kTRUE;
  auto iter_cond = fMinimumCondition.begin();
  while(iter_cond!=fMinimumCondition.end()){
    Int_t myIndex = *iter_cond;
    auto iter_find = find(fCoilIndex.begin(),fCoilIndex.end(),myIndex);
    if(iter_find==fCoilIndex.end()){
      kSufficient = kFALSE;
      break;
    }
    iter_cond++;
  }
  return kSufficient;
}
void TSolver::SetCondition(vector<Int_t> input){
  fMinimumCondition.clear();
  auto iter = input.begin();
  while(iter!=input.end()){
    fMinimumCondition.push_back(*iter);
    iter++;
  }
}

void TSolver::WriteSolution(vector<Double_t>  &fSlope){
  fSlope=fSolution;
}

// void TSolver::Print(){
//   cout << "\n ====== Results ====== " << endl;
//   Int_t nsize = fSolution.size();
//   cout << " --  slopes " << endl;
//   for(int i =0;i<nsize;i++){
//     cout << fSolution[i]*1e3<< endl;
//   }
//   Double_t residual=0.0;;
//   Int_t neq = fCoilIndex.size();
//   for(int icoil=0;icoil<neq;icoil++){
//     Double_t correction=0.0;
//     for(int i =0;i<nsize;i++){
//       correction+= fSolution[i]*fMonArray[icoil][i];
//     }

//     residual += pow(fDetArray[icoil]-correction,2);
//   }
//   cout << " -- residual : " << sqrt(residual)/neq << endl;
//   cout << " ====== END ======  \n" << endl;
// }
