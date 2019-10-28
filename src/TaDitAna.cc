#include "TCut.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TEventList.h"

#include "TaDitAna.hh"

ClassImp(TaDitAna);

TaDitAna::TaDitAna(TaConfig* aConfig){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  tree_name = aConfig->GetConfigParameter("tree_name");
  TString bmodcut_str =  aConfig->GetConfigParameter("bmod_cut");
  bmod_cut = bmodcut_str.Data();
  RegisterRawDataElements(aConfig->GetRawElementArray());
  fDependentVarArray = BuildDataElementArray( aConfig->GetDependentVarArray() );

#ifdef DEBUG
  auto iter_dep = fDependentVarArray.begin();
  while(iter_dep!=fDependentVarArray.end()){
    cout << (*iter_dep)->GetName() << endl;;
    iter_dep++;
  }
#endif

  vector<TString> fCoilNameArray;
  for(int i=1;i<=7;i++)
    fCoilNameArray.push_back(Form("bmod_trim%d",i));
  RegisterRawDataElements(fCoilNameArray);
  for(int i=1;i<=7;i++)
    fCoilArray.push_back(fDataElementMap[Form("bmod_trim%d",i)]);

  ProcessDefinitions(aConfig->GetDataElementDefinitions());

  protoCycle.LoadDetectorList(aConfig->GetDetectorList());
  protoCycle.RegisterDependentVarArray(fDependentVarArray);
  protoCycle.RegisterCoilArray(fCoilArray);
  protoCycle.InitAccumulators();

  if(fDependentVarArray.size()==0){
    cerr << " Error: Empty Dependent Channel Array " << endl;
  }

}

TaDitAna::~TaDitAna(){}

Bool_t TaDitAna::LoadModulationData(TaInput *aInput){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  TTree *evt_tree = aInput->GetEvtTree();

  TEventList *elist = new TEventList("elist");
  cout << " -- beam mod event cuts"
       << bmod_cut.GetTitle() << endl;
  Int_t nGoodEvents = evt_tree->Draw(">>+elist",bmod_cut,"goff");
  if(nGoodEvents==0){
    cerr<< " Error: No Beam Modulation Data is found. " << endl;
    cerr<< " Error: Analysis Aborted!  " << endl;
    return kFALSE;
  }
  cout << " -- " << nGoodEvents
       << " good beam modulation events are found " << endl;
  
  Double_t cycle_id;
  Double_t CodaEventNumber;
  Double_t bmwobj;
  evt_tree->SetBranchAddress("bmwobj",&bmwobj);
  evt_tree->SetBranchAddress("bmwcycnum",&cycle_id);
  evt_tree->SetBranchAddress("CodaEventNumber",&CodaEventNumber);

  evt_tree->SetBranchStatus("*",0);
  evt_tree->SetBranchStatus("CodaEventNumber",1);
  evt_tree->SetBranchStatus("bmwcycnum",1);
  evt_tree->SetBranchStatus("bmwobj",1);
  //***
  RegisterBranchAddress(evt_tree);
  //***
  TaSuperCycle supercycle_buff = protoCycle;
  Double_t last_cycle_id=0;
  Double_t last_obj=0;
  Int_t cur_index=0;
  Double_t bmwobj_isLock=kFALSE;
  for(int i=0;i<nGoodEvents;i++){
    Int_t ievt=elist->GetEntry(i);
    evt_tree->GetEntry(ievt);
    if(cycle_id==0 || bmwobj==0)
      continue;
    if(cycle_id >last_cycle_id){
      if(last_cycle_id!=0){
      	fSuperCycleArray.push_back(supercycle_buff);
	supercycle_buff.CalcSensitivities();
      }

      supercycle_buff = protoCycle;
      cout <<" -- Found a new supercycle: ID = " << cycle_id << endl;
      cout <<" -- Starting at CodaEventNumber: " << CodaEventNumber << endl;
      last_cycle_id=cycle_id;
      supercycle_buff.SetCycleID(cycle_id);
    }

    if(bmwobj!=last_obj){
      bmwobj_isLock=kFALSE;
      for(int icoil=1;icoil<=7;icoil++){
	if(bmwobj==icoil){
	  cur_index = icoil-1;
	  last_obj=bmwobj;
	  bmwobj_isLock = kTRUE;
	  break;
	} 
      } // end of searching loop
    } // end of new bmwobj lock

    if(bmwobj_isLock)
      supercycle_buff.UpdateSamples(cur_index);

  } // end of Good Events loop
  supercycle_buff.CalcSensitivities();
  fSuperCycleArray.push_back(supercycle_buff);

  return kTRUE;
}

// Bool_t TaDitAna::CalcSensitivities(){
//   Int_t cycID;
//   Int_t run_number = fInput->GetRunNumber();
//   vector<Double_t> dummy_vec(nCoil);
//   Vec2D detsens_buff(nDet,dummy_vec);
//   Vec2D detsens_err_buff(nDet,dummy_vec);
//   Vec2D monsens_buff(nMon,dummy_vec);
//   Vec2D monsens_err_buff(nMon,dummy_vec);
//   sens_tree->Branch("cycID",&cycID);
//   sens_tree->Branch("run",&run_number);
//   for(int icoil=0;icoil<nCoil;icoil++){
//     for(int idet=0;idet<nDet;idet++){
//       sens_tree->Branch(fDetectorArray[idet]+Form("_coil%d",coil_index[icoil]),
// 		    &detsens_buff[idet][icoil]);
//       sens_tree->Branch(fDetectorArray[idet]+Form("_coil%d_err",coil_index[icoil]),
// 		    &detsens_err_buff[idet][icoil]);

//     }
//     for(int imon=0;imon<nMon;imon++){
//       sens_tree->Branch(fMonitorArray[imon]+Form("_coil%d",coil_index[icoil]),
// 		    &monsens_buff[imon][icoil]);
//       sens_tree->Branch(fMonitorArray[imon]+Form("_coil%d_err",coil_index[icoil]),
// 		    &monsens_err_buff[imon][icoil]);
//     }
//   }
//   Int_t nCycle = fSuperCycleArray.size();
//   DetSens.resize(nDet);
//   DetSens_err.resize(nDet);
//   MonSens.resize(nMon);
//   MonSens_err.resize(nMon);
//   for(int idet=0;idet<nDet;idet++){
//     DetSens[idet].resize(nCoil,0.0);
//     DetSens_err[idet].resize(nCoil,-1.0);
//   }
//   for(int imon=0;imon<nMon;imon++){
//     MonSens[imon].resize(nCoil,0.0);
//     MonSens_err[imon].resize(nCoil,-1.0);
//   }
//   for(int icyc=0;icyc<nCycle;icyc++){
//     fSuperCycleArray[icyc].ComputeSensitivities();
//     for(int imon=0;imon<nMon;imon++){
//       for(int icoil=0;icoil<nCoil;icoil++){
// 	double this_slope = fSuperCycleArray[icyc].GetMonSens(imon,icoil);
// 	double this_error = fSuperCycleArray[icyc].GetMonSens_err(imon,icoil);
// 	monsens_buff[imon][icoil] = this_slope;
// 	monsens_err_buff[imon][icoil] = this_error;
// 	if(this_error!=-1){
// 	  if(MonSens_err[imon][icoil]==-1){
// 	    MonSens[imon][icoil] = this_slope;
// 	    MonSens_err[imon][icoil] = this_error;
// 	  }
// 	  else{
// 	    double last_error = MonSens_err[imon][icoil];
// 	    double last_slope = MonSens[imon][icoil];
// 	    double weight = 1.0/pow(last_error,2)+1.0/pow(this_error,2);
// 	    MonSens[imon][icoil]=(1.0/weight)*(last_slope*1.0/pow(last_error,2)+
// 					       this_slope*1.0/pow(this_error,2));
// 	    MonSens_err[imon][icoil]=last_error*this_error/sqrt(pow(last_error,2)
// 								+pow(this_error,2));
// 	  }
// 	}
//       } // coil loop
//     }// monitor loop
//     for(int idet=0;idet<nDet;idet++){
//       for(int icoil=0;icoil<nCoil;icoil++){
// 	double this_slope = fSuperCycleArray[icyc].GetDetSens(idet,icoil);
// 	double this_error = fSuperCycleArray[icyc].GetDetSens_err(idet,icoil);
// 	detsens_buff[idet][icoil]=this_slope;
// 	detsens_err_buff[idet][icoil]=this_error;
// 	if(this_error!=-1){
// 	  if(DetSens_err[idet][icoil]==-1){
// 	    DetSens[idet][icoil] = this_slope;
// 	    DetSens_err[idet][icoil] = this_error;
// 	  }
// 	  else{
// 	    double last_error = DetSens_err[idet][icoil];
// 	    double last_slope = DetSens[idet][icoil];
// 	    double weight = 1.0/pow(last_error,2)+1.0/pow(this_error,2);
// 	    DetSens[idet][icoil]=(1.0/weight)*(last_slope*1.0/pow(last_error,2)+
// 					       this_slope*1.0/pow(this_error,2));
// 	    DetSens_err[idet][icoil]=last_error*this_error/sqrt(pow(last_error,2)
// 								+pow(this_error,2));
// 	  }
// 	}
//       }// coil loop
//     } // Detector loop
//     cycID = fSuperCycleArray[icyc].GetCycleID();
//     sens_tree->Fill();
//   } // Cycle loop
//   return kTRUE;
// }
// void TaDitAna::PrintSummary(){
// #ifdef NOISY
//   cout << __FUNCTION__ << endl;
// #endif

// }


void TaDitAna::ProcessDefinitions(vector<pair<TString,TString> > fDefinitions){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  vector<pair<TString,TString> >::iterator iter=fDefinitions.begin();
  while(iter!=fDefinitions.end()){
    TString myName=(*iter).first;
    TString definition=(*iter).second;
#ifdef DEBUG
    cout << "Combo Name: " << myName << endl;
#endif

    TaDataElement* tobedefined = new TaDataElement(myName);
    while( definition.Length()>0 ){
      Ssiz_t next_plus = definition.Last('+');
      Ssiz_t next_minus = definition.Last('-');
      Ssiz_t length=definition.Length();
      Ssiz_t head=0;
      if(next_minus>next_plus)
	head = next_minus;
      if(next_plus>next_minus)
	head = next_plus;

      TString extracted = definition(head,length-head);
#ifdef DEBUG
      cout <<  extracted << endl;
#endif
      if(extracted.Contains('*')){
	Ssiz_t aster_pos = extracted.First('*');
	Ssiz_t form_length = extracted.Length();
	Double_t factor = TString(extracted(0,aster_pos)).Atof();
	TString elementName = extracted(aster_pos+1,form_length-(aster_pos+1));
	TaDataElement* tobeconnected = fDataElementMap[elementName];
	tobedefined->AddElement(factor,tobeconnected);
#ifdef DEBUG
	cout << factor <<"\t" << elementName << endl;
#endif
      }
      else if(extracted(0,1)=="+"){
	TString bare_name = extracted(1,length-head-1);
  	TaDataElement* tobeconnected = fDataElementMap[bare_name];
	tobedefined->AddElement(1,tobeconnected);
#ifdef DEBUG
	cout << +1 <<"\t" << extracted << endl;
#endif
      } else if(extracted(0,1)=="-"){
	TString bare_name = extracted(1,length-head-1);
  	TaDataElement* tobeconnected = fDataElementMap[bare_name];
	tobedefined->AddElement(-1,tobeconnected);
#ifdef DEBUG
	cout << -1 <<"\t" << extracted << endl;
#endif
      }else{
  	TaDataElement* tobeconnected = fDataElementMap[extracted];
	tobedefined->AddElement(1,tobeconnected);
#ifdef DEBUG
	cout << 1 <<"\t" << extracted << endl;
#endif
      }
      definition.Remove(head,length-head);
    }

    iter++;
    if(fDataElementMap.find(myName)!=fDataElementMap.end())
      fDataElementMap[myName]=tobedefined;
  }
}

vector<TaDataElement*> TaDitAna::BuildDataElementArray( vector<TString> device_array ) {
  vector<TString>::iterator iter = device_array.begin();
  vector<TaDataElement*> fArray_ret;
  while(iter!=device_array.end()){
    if( fDataElementMap.find(*iter)!=fDataElementMap.end())
      fArray_ret.push_back(fDataElementMap[*iter]);
    iter++;
  }  
  return fArray_ret;
}


void TaDitAna::RegisterRawDataElements(vector<TString> device_array){
  vector<TString>::iterator iter = device_array.begin();
  while(iter!=device_array.end()){
    TaDataElement* new_ptr = new TaDataElement(*iter);
    if(fDataElementMap.find(*iter)==fDataElementMap.end()){
      fDataElementMap[*iter] = new_ptr; 
      fRawDataElementArray.push_back(new_ptr);
    }
    iter++;
  }
}

void TaDitAna::RegisterBranchAddress(TTree *fTree){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  vector<TaDataElement*>::iterator iter = fRawDataElementArray.begin();
  while(iter!=fRawDataElementArray.end()){
    TString myName=(*iter)->GetName();
    TBranch* fBranch = fTree->GetBranch(myName);
    if(fBranch==NULL){
      cout << "Branch " << myName  << " is not found" << endl;
      iter++;
      continue;
    }
    (*iter)->RegisterBranchAddress(fBranch);
    fTree->SetBranchStatus(myName,1);
    iter++;
  }
}
