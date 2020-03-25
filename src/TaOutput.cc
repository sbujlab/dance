#include "TaOutput.hh"
ClassImp(TaOutput)

TaOutput::TaOutput():nBranches(0){
  outputFile = new TFile("test.root","RECREATE");
  fPrinter = new TaPrinter("test.log");
  parity_scale.ppm=1e-6;
  parity_scale.ppb=1e-9;
  parity_scale.um=1e-3;
  parity_scale.nm=1e-6;
  run_number = 999999;
}

TaOutput::TaOutput(TaConfig* aConfig):nBranches(0){
  TString prefix = aConfig->GetConfigParameter("output_prefix");
  TString path = aConfig->GetConfigParameter("output_path");
  run_number = aConfig->GetRunNumber();
  int seg_number = aConfig->GetSegNumber();
  outputFile = new TFile(path+prefix+Form("%d.%03d.root",run_number,seg_number),"RECREATE");
  fPrinter = new TaPrinter(path+prefix+Form("%d_summary.txt",run_number));
  parity_scale.ppm=1e-6;
  parity_scale.ppb=1e-9;
  parity_scale.um=1e-3;
  parity_scale.nm=1e-6;

}

TaOutput::~TaOutput(){

}

void TaOutput::ConstructTreeBranch(TString treeName, 
				   TString branchName,
				   Double_t &value){
  outputFile->cd();
  if(fTreeArrayByName.find(treeName)==fTreeArrayByName.end()){
    TTree *newTree = new TTree(treeName,"");
    newTree->Branch("unit",&parity_scale,"ppm/D:ppb:um:nm");
    fTreeArrayByName[treeName]=newTree;
    fTreeArray.push_back(newTree);
  }
  
  Int_t myIndex = nBranches;
  nBranches++;
  TBranch *myBranch = fTreeArrayByName[treeName]->Branch(branchName,
  							 &value,
  							 "hw_sum/D");
  fBranchArray.push_back(myBranch);
  fBranchIndex[make_pair(treeName,branchName)]=myIndex;

}

void TaOutput::ConstructStatTreeBranch(TString treeName, 
				      TString branchName,
				      STAT &value){
  outputFile->cd();
  if(fTreeArrayByName.find(treeName)==fTreeArrayByName.end()){
    TTree *newTree = new TTree(treeName,"");
    newTree->Branch("unit",&parity_scale,"ppm/D:ppb:um:nm");
    fTreeArrayByName[treeName]=newTree;
    fTreeArray.push_back(newTree);
  }
  
  Int_t myIndex = nBranches;
  nBranches++;
  TString leaflist="mean/D:err:rms:m2:num_samples";
  TBranch *myBranch = fTreeArrayByName[treeName]->Branch(branchName,
  							 &value,
  							 leaflist);
  fBranchArray.push_back(myBranch);
  fBranchIndex[make_pair(treeName,branchName)]=myIndex;
}

void TaOutput::FillTree(TString treeName){
  outputFile->cd();
  if(fTreeArrayByName.find(treeName)==fTreeArrayByName.end()){
    cout << treeName << "not found " << endl;
    return;
  }
  fTreeArrayByName[treeName]->Fill();
}

void TaOutput::Write(){
  outputFile->cd();
  Int_t nTrees = fTreeArray.size();
  for(int i=0;i<nTrees;i++){
    fTreeArray[i]->Write();
  }

}
void TaOutput::Close(){
  outputFile->Close();
}
