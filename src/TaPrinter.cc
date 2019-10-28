#include "TaPrinter.hh"


ClassImp(TaPrinter);

TaPrinter::TaPrinter(TString filename){
  logfile.open(filename.Data());
  nColumns=0;
  std::cout << __PRETTY_FUNCTION__ <<  ": opening " <<  filename << endl;
  column_sep ="|";
}
void TaPrinter::Close(){
  logfile.close();
}

void TaPrinter::AddHeader(vector<TString> input){
  Int_t ncol = input.size();
  NewLine();
  for(int icol=0;icol<ncol;icol++){
    AddStringEntry(input[icol]);
  }
  InsertHorizontalLine();
}

void TaPrinter::AddFloatEntry(Double_t finput){
  TString format_input = Form("%.2f",finput);
  current_row.AddEntry(format_input);
}

void TaPrinter::AddFloatEntryWithError(Double_t finput,Double_t error_bar){
  TString format_input = Form("%.2f +/- %.2f",finput,error_bar);
  current_row.AddEntry(format_input);
}

void TaPrinter::NewLine(){
  if(current_row.nEntries()!=0){
    fRows.push_back(current_row);
    if(current_row.nEntries()>nColumns)
      nColumns = current_row.nEntries();
    current_row.Clear();
  }
}
void TaPrinter::InsertHorizontalLine(){
  NewLine();
  AddStringEntry("hline");
  NewLine();
}

void TaPrinter::Print(ostream &out){
  NewLine();
  Int_t nrow = fRows.size();
  Int_t total_width=0;
  column_wid.resize(nColumns,0);
  
  for(int icol=0;icol<nColumns;icol++){
    for(int irow=0;irow<nrow;irow++){    
      Int_t this_colwid = fRows[irow].GetColumnWidth(icol);
      if(this_colwid>column_wid[icol])
	column_wid[icol] = this_colwid;
    }
    total_width +=(column_wid[icol]+column_sep.Length());
  }

  for(int irow=0;irow<nrow;irow++){
    if(fRows[irow].nEntries()==1 && fRows[irow].GetEntry(0)=="hline"){
      for(int i=0;i<total_width;i++)
	out<< "-";
      out << endl;
    }
    else{
      for(int icol=0;icol<nColumns;icol++){
	out<< setw(column_wid[icol]) << setiosflags(ios::left)
	   << Center(fRows[irow].GetEntry(icol), column_wid[icol]);
	out << column_sep;
      }
      out << endl;
    }
  }
}

void TaPrinter::PrintToFile(){
  Print(logfile);
}

TString TaPrinter::Center(TString input, Int_t this_colwid){

  Int_t string_length = input.Length();
  Int_t npad = floor( (this_colwid+1 -string_length)*0.5);
  TString a_single_space =" ";
  TString output;

  for(int i=0;i<npad;i++)
    output +=a_single_space;

  output+=input;

  return output;
}
