#include "TCut.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TEventList.h"

#include "TaDithering.hh"

ClassImp(TaDithering);

TaDithering::TaDithering(TaInput* aInput){fInput = aInput;}
TaDithering::~TaDithering(){}

Bool_t TaDithering::LoadConfig(TaConfig *fConfig){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  TCut evt_cut(fConfig->GetEventCut());  
  cout << "-- Using Event Cut: " 
       << evt_cut.GetTitle() << endl;
  TCut bmw_cut ="bmwobj>0";
  bmod_cut = evt_cut && bmw_cut;
  TTree* evt_tree = fInput->GetEvtTree();
  det_array = fConfig->GetDetArray();
  mon_array = fConfig->GetMonArray();
  coil_array = fConfig->GetCoilArray();

  if(det_array.size()==0 ||
     mon_array.size()==0){
    cerr << " Error: Empty Detector or Monitor Array " << endl;
    cerr << " Error: Nothing I can do... Bye-bye! " << endl;
    return kFALSE;
  }
  if(coil_array.size()>7){
    cerr << " Error: Found " << coil_array.size() << " coils " << endl;
    cerr << " Error: Oversize Coil Array! " << endl;
    return kFALSE;
  }
  if(coil_array.size()>= mon_array.size()){
    cerr << " Error: Over-constraint condition ! " << endl;
    return kFALSE;
  }
  
  TObjArray* branch_array = evt_tree->GetListOfBranches();
  int ndet = det_array.size();
  for(int idet=0;idet<ndet;idet++){
    TObject *obj_ptr = branch_array->FindObject(det_array[idet]);
    if(obj_ptr==NULL){
      cerr<< " Error : Branch " 
	  << det_array[idet] 
	  << " is not found " << endl;
      return kFALSE;
    }
  }
  int nmon = mon_array.size();
  for(int imon=0;imon<nmon;imon++){
    TObject *obj_ptr = branch_array->FindObject(mon_array[imon]);
    if(obj_ptr==NULL){
      cerr<< " Error : Branch " 
	  << mon_array[imon] 
	  << " is not found " << endl;
      return kFALSE;
    }
  }
  int ncoil = coil_array.size();
  for(int icoil=0;icoil<ncoil;icoil++){
    TObject *obj_ptr = branch_array->FindObject(coil_array[icoil]);
    if(obj_ptr==NULL){
      cerr<< " Error : Branch " 
	  << coil_array[icoil] 
	  << " is not found " << endl;
      return kFALSE;
    }
  }

  nDet = det_array.size();
  nMon = mon_array.size();
  nCoil = coil_array.size();

  sens_tree = new TTree("sens","Sensitivities Tree");
  sens_tree->SetMarkerStyle(20);

  return kTRUE;
}

Bool_t TaDithering::LoadModulationData(){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  TTree *evt_tree = fInput->GetEvtTree();
  evt_tree->SetBranchStatus("*",0);
  evt_tree->SetBranchStatus("ErrorFlag",1);
  evt_tree->SetBranchStatus("CodaEventNumber",1);
  evt_tree->SetBranchStatus("bmw*",1);
  evt_tree->SetBranchStatus("bmod*",1);
  for(int idet=0; idet<nDet;idet++)
    evt_tree->SetBranchStatus(det_array[idet],1);
  for(int imon=0; imon<nMon;imon++)
    evt_tree->SetBranchStatus(mon_array[imon],1);

  TEventList *elist = new TEventList("elist");
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
  
  const int ndet = det_array.size();
  Double_t det_value[ndet];
  for(int idet=0; idet<ndet;idet++){
    TLeaf* this_leaf=evt_tree->GetBranch(det_array[idet])->GetLeaf("hw_sum");
    this_leaf->SetAddress(&det_value[idet]);
  }

  const int nmon = mon_array.size();
  Double_t mon_value[nmon];
  for(int imon=0; imon<nmon;imon++){
    TLeaf* this_leaf=evt_tree->GetBranch(mon_array[imon])->GetLeaf("hw_sum");
    this_leaf->SetAddress(&mon_value[imon]);
  }
  
  vector<TLeaf*> trimcard_leaf;
  int ncoil = coil_array.size();
  for(int icoil=0;icoil<ncoil;icoil++){
    Ssiz_t length=coil_array[icoil].Length();
    Int_t this_index= TString(coil_array[icoil][length-1]).Atoi();
    if(this_index<=0 || this_index>7){
      cout << " Warning: Coil ID make no sense !" << endl;
      cout << " Warning: will skip "
	   << coil_array[icoil] << endl;
      continue;
    }
    TLeaf* this_leaf=evt_tree->GetBranch(coil_array[icoil])->GetLeaf("value");
    trimcard_leaf.push_back(this_leaf);
    coil_index.push_back(this_index);
  }
  // FIXME: Should be ok if it is in sens mode
  ncoil = coil_index.size();
  if(ncoil>= nmon){
    cerr << __PRETTY_FUNCTION__ 
	 << " Error: Over-constraint condition ! " << endl;
    return kFALSE;
  }
  cout << " -- Using Coil ID: ";
  for(int icoil=0;icoil<ncoil;icoil++)
    cout << coil_index[icoil] << "\t";
  cout << endl;

  TaSuperCycle protoCycle;
  protoCycle.Resize(ndet,nmon,ncoil);
  TaSuperCycle sc_buff = protoCycle;
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
      if(last_cycle_id!=0)
	fSuperCycleArray.push_back(sc_buff);
      sc_buff = protoCycle;
      cout <<" -- Found a new supercycle: ID = " << cycle_id << endl;
      cout <<" -- Starting at CodaEventNumber: " << CodaEventNumber << endl;
      last_cycle_id=cycle_id;
      sc_buff.SetCycleID(cycle_id);
    }

    if(bmwobj!=last_obj){
      bmwobj_isLock=kFALSE;
      for(int icoil=0;icoil<ncoil;icoil++){
	if(bmwobj==coil_index[icoil]){
	  cur_index = icoil;
	  last_obj=bmwobj;
	  bmwobj_isLock = kTRUE;
	  break;
	} 
      } // end of searching loop
    } // end of new bmwobj lock

    if(bmwobj_isLock){
      Double_t trimcard_value=trimcard_leaf[cur_index]->GetValue();
      sc_buff.LoadCoilData(cur_index,trimcard_value/10.0);
      for(int idet=0;idet<ndet;idet++)
	sc_buff.LoadDetData(idet,cur_index,
			     det_value[idet],
			     trimcard_value/10.0);
      for(int imon=0;imon<nmon;imon++)
	sc_buff.LoadMonData(imon,cur_index,
			     mon_value[imon],
			     trimcard_value/10.0);
    } // end of Data Extraction
  } // end of Good Events loop
  fSuperCycleArray.push_back(sc_buff);
  return kTRUE;
}

Bool_t TaDithering::ComputeSensitivities(){
  Int_t cycID;
  Int_t run_number = fInput->GetRunNumber();
  vector<Double_t> dummy_vec(nCoil);
  Vec2D detsens_buff(nDet,dummy_vec);
  Vec2D detsens_err_buff(nDet,dummy_vec);
  Vec2D monsens_buff(nMon,dummy_vec);
  Vec2D monsens_err_buff(nMon,dummy_vec);
  sens_tree->Branch("cycID",&cycID);
  sens_tree->Branch("run",&run_number);
  for(int icoil=0;icoil<nCoil;icoil++){
    for(int idet=0;idet<nDet;idet++){
      sens_tree->Branch(det_array[idet]+Form("_coil%d",coil_index[icoil]),
		    &detsens_buff[idet][icoil]);
      sens_tree->Branch(det_array[idet]+Form("_coil%d_err",coil_index[icoil]),
		    &detsens_err_buff[idet][icoil]);

    }
    for(int imon=0;imon<nMon;imon++){
      sens_tree->Branch(mon_array[imon]+Form("_coil%d",coil_index[icoil]),
		    &monsens_buff[imon][icoil]);
      sens_tree->Branch(mon_array[imon]+Form("_coil%d_err",coil_index[icoil]),
		    &monsens_err_buff[imon][icoil]);
    }
  }
  Int_t nCycle = fSuperCycleArray.size();
  DetSens.resize(nDet);
  DetSens_err.resize(nDet);
  MonSens.resize(nMon);
  MonSens_err.resize(nMon);
  for(int idet=0;idet<nDet;idet++){
    DetSens[idet].resize(nCoil,0.0);
    DetSens_err[idet].resize(nCoil,-1.0);
  }
  for(int imon=0;imon<nMon;imon++){
    MonSens[imon].resize(nCoil,0.0);
    MonSens_err[imon].resize(nCoil,-1.0);
  }
  for(int icyc=0;icyc<nCycle;icyc++){
    fSuperCycleArray[icyc].ComputeSensitivities();
    for(int imon=0;imon<nMon;imon++){
      for(int icoil=0;icoil<nCoil;icoil++){
	double this_slope = fSuperCycleArray[icyc].GetMonSens(imon,icoil);
	double this_error = fSuperCycleArray[icyc].GetMonSens_err(imon,icoil);
	monsens_buff[imon][icoil] = this_slope;
	monsens_err_buff[imon][icoil] = this_error;
	if(this_error!=-1){
	  if(MonSens_err[imon][icoil]==-1){
	    MonSens[imon][icoil] = this_slope;
	    MonSens_err[imon][icoil] = this_error;
	  }
	  else{
	    double last_error = MonSens_err[imon][icoil];
	    double last_slope = MonSens[imon][icoil];
	    double weight = 1.0/pow(last_error,2)+1.0/pow(this_error,2);
	    MonSens[imon][icoil]=(1.0/weight)*(last_slope*1.0/pow(last_error,2)+
					       this_slope*1.0/pow(this_error,2));
	    MonSens_err[imon][icoil]=last_error*this_error/sqrt(pow(last_error,2)
								+pow(this_error,2));
	  }
	}
      } // coil loop
    }// monitor loop
    for(int idet=0;idet<nDet;idet++){
      for(int icoil=0;icoil<nCoil;icoil++){
	double this_slope = fSuperCycleArray[icyc].GetDetSens(idet,icoil);
	double this_error = fSuperCycleArray[icyc].GetDetSens_err(idet,icoil);
	detsens_buff[idet][icoil]=this_slope;
	detsens_err_buff[idet][icoil]=this_error;
	if(this_error!=-1){
	  if(DetSens_err[idet][icoil]==-1){
	    DetSens[idet][icoil] = this_slope;
	    DetSens_err[idet][icoil] = this_error;
	  }
	  else{
	    double last_error = DetSens_err[idet][icoil];
	    double last_slope = DetSens[idet][icoil];
	    double weight = 1.0/pow(last_error,2)+1.0/pow(this_error,2);
	    DetSens[idet][icoil]=(1.0/weight)*(last_slope*1.0/pow(last_error,2)+
					       this_slope*1.0/pow(this_error,2));
	    DetSens_err[idet][icoil]=last_error*this_error/sqrt(pow(last_error,2)
								+pow(this_error,2));
	  }
	}
      }// coil loop
    } // Detector loop
    cycID = fSuperCycleArray[icyc].GetCycleID();
    sens_tree->Fill();
  } // Cycle loop
  return kTRUE;
}
void TaDithering::PrintSummary(){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  Int_t ndet = det_array.size();
  Int_t nmon = mon_array.size();
  Int_t ncoil = coil_array.size();
  cout << "-- Coil # : " ;
  for(int ic=0;ic<nCoil;ic++)
    cout << coil_index[ic] << "\t";
  cout << endl;
  for(int idet=0;idet<ndet;idet++){
    cout << det_array[idet] << "\t";
    for(int ic=0;ic<ncoil;ic++){
      cout << DetSens[idet][ic]<<"+/-"
	   << DetSens_err[idet][ic]<<"\t";
    }
    cout << endl;
  }
  for(int imon=0;imon<nmon;imon++){
    cout << mon_array[imon] << "\t";
    for(int ic=0;ic<ncoil;ic++){
      cout << MonSens[imon][ic]<<"+/-"
	   << MonSens_err[imon][ic]<<"\t";
    }	
    cout << endl;
  }
  Int_t nCycle = fSuperCycleArray.size();
  for(int icyc=0;icyc<nCycle;icyc++){
    cout << "-- cycle ID : " 
  	 << fSuperCycleArray[icyc].GetCycleID() << endl;
    fSuperCycleArray[icyc].ComputeSensitivities();
    cout << "-- Coil # : " ;
    for(int ic=0;ic<nCoil;ic++)
      cout << coil_index[ic] << "\t";
    cout << endl;
    cout << "-- nSamples : " ;
    for(int ic=0;ic<nCoil;ic++)
      cout << fSuperCycleArray[icyc].GetNSamples(ic)<< "\t";
    cout << endl;
    for(int idet=0;idet<ndet;idet++){
      cout << det_array[idet] << "\t";
      for(int ic=0;ic<ncoil;ic++){
  	cout << fSuperCycleArray[icyc].GetDetSens(idet,ic)<<"+/-"
  	     << fSuperCycleArray[icyc].GetDetSens_err(idet,ic)<<"\t";
      }
      cout << endl;
    }
    for(int imon=0;imon<nmon;imon++){
      cout << mon_array[imon] << "\t";
      for(int ic=0;ic<ncoil;ic++){
  	cout << fSuperCycleArray[icyc].GetMonSens(imon,ic)<<"+/-"
  	     << fSuperCycleArray[icyc].GetMonSens_err(imon,ic)<<"\t";
      }	
      cout << endl;
    }
  }
}

TMatrixD TaDithering::GetDetSensMatrix(){
  TMatrixD detMatrix(nDet,nCoil);
  for(int irow=0;irow<nDet;irow++)
    for(int icol=0;icol<nCoil;icol++)
      detMatrix[irow][icol] = DetSens[irow][icol];
  return detMatrix;
}
TMatrixD TaDithering::GetMonSensMatrix(){
  TMatrixD monMatrix(nMon,nCoil);
  for(int irow=0;irow<nMon;irow++)
    for(int icol=0;icol<nCoil;icol++)
      monMatrix[irow][icol] =MonSens[irow][icol];
  return monMatrix;
}
