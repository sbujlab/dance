#include "TCut.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TEventList.h"

#include "TaDitAna.hh"
#include "TaPrinter.hh"

ClassImp(TaDitAna);

TaDitAna::TaDitAna(TaConfig* aConfig){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__ << endl;
#endif
  tree_name = aConfig->GetConfigParameter("tree_name");
  TString bmodcut_str =  aConfig->GetConfigParameter("bmod_cut");
  bmod_cut = bmodcut_str.Data();
  RegisterRawDataElements(aConfig->GetRawElementArray());

  vector<TString> fCoilNameArray;
  for(int i=1;i<=7;i++)
    fCoilNameArray.push_back(Form("bmod_trim%d",i));
  RegisterRawDataElements(fCoilNameArray);
  for(int i=1;i<=7;i++)
    fCoilArray.push_back(fDataElementMap[Form("bmod_trim%d",i)]);

  ProcessDefinitions(aConfig->GetDataElementDefinitions());
  fDependentVarArray = BuildDataElementArray( aConfig->GetDependentVarArray() );

  templateCycle.LoadDetectorList(aConfig->GetDetectorList());
  templateCycle.RegisterDependentVarArray(fDependentVarArray);
  templateCycle.RegisterCoilArray(fCoilArray);
  templateCycle.InitAccumulators();
  templateCycle.ConfigSlopesCalculation(aConfig);

  if(aConfig->GetConfigParameter("device_error_cut")=="on")
    templateCycle.EnableDeviceErrorCut();

  if(fDependentVarArray.size()==0){
    cerr << " Error: Empty Dependent Channel Array " << endl;
  }

#ifdef DEBUG
  auto iter_dep = fDependentVarArray.begin();
  while(iter_dep!=fDependentVarArray.end()){
    cout << (*iter_dep)->GetName() << endl;;
    iter_dep++;
  }
#endif

}

Bool_t TaDitAna::LoadModulationData(TaInput *aInput){
#ifdef NOISY
  cout << __FUNCTION__ << endl;
#endif
  TTree *evt_tree = aInput->GetEvtTree();

  TEventList *elist = new TEventList("elist");
  TEventList *elist_all = new TEventList("elist_all");
  cout << " -- beam mod event cuts"
       << bmod_cut.GetTitle() << endl;
  Int_t nGoodEvents = evt_tree->Draw(">>+elist",bmod_cut,"goff");
  Int_t nBeamModEvents = evt_tree->Draw(">>+elist_all","bmwcycnum>0","goff");
  if(nBeamModEvents==0){
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
  TaSuperCycle supercycle_buff = templateCycle;
  Double_t last_cycle_id=0;
  Double_t last_obj=0;
  Int_t cur_index=0;
  Double_t bmwobj_isLock=kFALSE;
  Bool_t kGoodEvent;
  for(int i=0;i<nBeamModEvents;i++){
    Int_t ievt=elist_all->GetEntry(i);
    evt_tree->GetEntry(ievt);
    if(elist->Contains(ievt))
      kGoodEvent = kTRUE;
    else
      kGoodEvent = kFALSE;

    if(cycle_id==0 || bmwobj==0)
      continue;
    if(cycle_id >last_cycle_id){
      if(last_cycle_id!=0){
      	fSuperCycleArray.push_back(supercycle_buff);
      }

      supercycle_buff = templateCycle;
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

    if(bmwobj_isLock && kGoodEvent)
      supercycle_buff.UpdateSamples(cur_index);

  } // end of Good Events loop
  fSuperCycleArray.push_back(supercycle_buff);
  return kTRUE;
}
void TaDitAna::Process(){
  auto iter = fSuperCycleArray.begin();
  while(iter!=fSuperCycleArray.end()){
    (*iter).CalcSensitivities();
    (*iter).CalcSlopes();
    iter++;
  }
}
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
    if(fDataElementMap.find(myName)==fDataElementMap.end())
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
    if(myName=="bmod_trim7")
      myName="bmod_ramp";
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

void TaDitAna::PrintSummary(TaOutput* aOutput){
  TaPrinter *fPrinter = aOutput->GetPrinter();
  auto iter=fSuperCycleArray.begin();
  while(iter!=fSuperCycleArray.end()){
    (*iter).WriteToPrinter(fPrinter);
    iter++;
  }
  fPrinter->PrintToFile();
  fPrinter->Print(std::cout);
  fPrinter->Close();
}
void TaDitAna::WriteToTree(TaOutput* aOutput){
#ifdef NOISY
  cout << __PRETTY_FUNCTION__<< endl;
#endif
  auto iter=fSuperCycleArray.begin();

  Int_t nDV = fDependentVarArray.size();
  Int_t nCoil =7;
  vector<Double_t> fSens(nDV*nCoil);
  vector<Double_t> fSensErr(nDV*nCoil);
  vector<Double_t> fNSamps(nDV*nCoil);
  for(int idv=0;idv<nDV;idv++){
    for(int ic=1;ic<=nCoil;ic++){
      TString dv_name = fDependentVarArray[idv]->GetName();
      TString coil_name = fCoilArray[ic-1]->GetName();
      Int_t myIndex = (*iter).GetIndex(make_pair(dv_name,coil_name));
      TString chName;
      chName = Form("%s_coil%d",dv_name.Data(),ic);
      aOutput->ConstructTreeBranch("sens",chName,fSens[myIndex]);
      chName = Form("%s_coil%d_err",dv_name.Data(),ic);
      aOutput->ConstructTreeBranch("sens",chName,fSensErr[myIndex]);
      chName = Form("%s_coil%d_nsamp",dv_name.Data(),ic);
      aOutput->ConstructTreeBranch("sens",chName,fNSamps[myIndex]);
    }
  }
  Double_t cycID;
  aOutput->ConstructTreeBranch("sens","cycID",cycID);
  Double_t run_number = aOutput->GetRunNumber();
  aOutput->ConstructTreeBranch("sens","run",run_number);

  while(iter!=fSuperCycleArray.end()){
    cycID =(Double_t) (*iter).GetCycleID();
    for(int i=0;i<nDV*nCoil;i++){
      fSens[i] = (*iter).GetSensitivity(i);
      fSensErr[i] = (*iter).GetErrorBar(i);
      fNSamps[i]=(*iter).GetNSamples(i);
    }
    aOutput->FillTree("sens");
    iter++;
  }

  // Write to slopes Tree
  Int_t nMod = templateCycle.GetNumberOfSlopeMode();
  for(int imod=0;imod<nMod;imod++){
    TString mod_name = templateCycle.GetTreeName(imod);
    Double_t cycID;
    aOutput->ConstructTreeBranch(mod_name,"cycID",cycID);
    Double_t run_number = aOutput->GetRunNumber();
    aOutput->ConstructTreeBranch(mod_name,"run",run_number);

    vector<Double_t> fBranchValues;
    vector<Double_t> fFlagValues;
    auto iter_cyc = fSuperCycleArray.begin();
    (*iter_cyc).ConstructSlopeTreeBranch(aOutput,imod,
					 fBranchValues,fFlagValues);
    while(iter_cyc!=fSuperCycleArray.end()){
      cycID = (*iter_cyc).GetCycleID();
      (*iter_cyc).FillSlopeTree(aOutput,imod,
				fBranchValues,fFlagValues);
      iter_cyc++;
    }
  }
}
