/**********************************   
Dithering Sensitivity Calculator
	author: Tao Ye
	<tao.ye@stonybrook.edu>
	last update: Oct 2019
************************************/
#include "TaInput.hh"
#include "TaOutput.hh"
#include "TaConfig.hh"

#include "TaDitAna.hh"
#include "TaSuperCycle.hh"

#include "TStopwatch.h"
#include "TString.h"
#include "TROOT.h"
#include <iostream>

using namespace std;
int main(int argc, char** argv){
  TStopwatch tsw;

  tsw.Start();
  if(argc==1){
    cout << "\n*********************************************** " <<endl;
    cout << " Calculate It! (calcit) " << endl;
    cout <<  "\t - a Sensitivity Calculator "  << endl;
    cout<< "\t author: Tao Ye " << endl;
    cout<< "\t <tao.ye@stonybrook.edu>" << endl;
    cout<< "\t last update: Oct 2019" << endl;
    cout << "***********************************************\n " <<endl;
    return 0;
  }
  int opt;
  TString confFileName;
  TString japanFileName;
  Int_t run_number=0;
  while( (opt=getopt(argc,argv,":f:r:c:"))!=-1){
    switch(opt){
    case ':':
      cerr << argv[optind-1] << "requires value. " <<endl;
      return 1;
    case '?':
      cerr << "Unknown arguments " <<optopt << endl;
      return 1;
    case 'f':
      japanFileName=TString(optarg);
      break;
    case 'c':
      confFileName=TString(optarg);
      break;
    case 'r':
      run_number= atoi(optarg);
      break;
    }
  }

  TaConfig *fConfig= new TaConfig();
  fConfig->SetInputName(japanFileName);
  fConfig->ParseFile(confFileName);
  fConfig->SetRunNumber(run_number);

  TaInput *fInput = new TaInput(fConfig);
  fInput->LoadROOTFile();
  TaOutput *fOutput = new TaOutput(fConfig); 

  TaDitAna *fDithering = new TaDitAna(fConfig);
  fDithering->LoadModulationData(fInput);
  fDithering->WriteToTree(fOutput);
  // fDithering->ComputeSensitivities();
  // fDithering->PrintSummary();
  
  fOutput->Write();
  fOutput->Close();
  // vector<TString> det_array = fConfig->GetDetArray();
  // vector<TString> mon_array = fConfig->GetMonArray();
  // vector<TString> coil_array = fConfig->GetCoilArray();
  // output_rootfile->WriteObject(&det_array,"det_array");
  // output_rootfile->WriteObject(&mon_array,"mon_array");
  // output_rootfile->WriteObject(&coil_array,"coil_array");
  // output_rootfile->Close();

  cout <<" -- " ;
  tsw.Print();
  return 0;
}

