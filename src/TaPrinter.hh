#ifndef __TAPRINTER_HH__
#define __TAPRINTER_HH__

#include <ostream>
#include <iostream>
#include <fstream>
#include "TString.h"

using namespace std;
class TaRow{
public:
  TaRow(){};
  virtual ~TaRow(){};
  void AddEntry(TString entry){
    fEntries.push_back(entry);
    column_width.push_back(entry.Length());
  }
  
  void Clear(){
    fEntries.clear();
    column_width.clear();
  }

  Int_t nEntries(){ return fEntries.size();};

  TString GetEntry(int ie) {
    if(ie>=nEntries())
      return "n/a";
    else
      return fEntries[ie];};

  Int_t GetColumnWidth(int ie) {
    if(ie>=nEntries())
      return 3;
    else
      return column_width[ie];};

private:
  vector<TString> fEntries;
  vector<Int_t> column_width;
  ClassDef(TaRow,0);
};

using namespace std;
class TaPrinter{
public:
  TaPrinter(TString filename);
  virtual  ~TaPrinter(){};
  void AddHeader(vector<TString>);
  void AddStringEntry(TString input){ current_row.AddEntry(input);};
  void AddFloatEntry(Double_t finput);
  void AddFloatEntryWithError(Double_t finput, Double_t err_bar);

  void InsertHorizontalLine();
  void NewLine();
  void Close();
  TString Center(TString,Int_t);
  void Print(ostream &out);
  void PrintToFile();
private:
  ofstream logfile;
  TaRow current_row;
  vector<TaRow> fRows;
  vector<Int_t> column_wid;
  Int_t nColumns;
  TString column_sep ;
  ClassDef(TaPrinter,0);
};

#endif
