class TSolver{
public:
  TSolver();
  TSolver(Int_t id);
  virtual ~TSolver(){};
  void LoadSensitivity(Int_t index, Double_t fDet, vector<Double_t> fMon);
  Bool_t SolveMatrix();
  void WriteSolution(vector<Double_t> &fSlope);
  void SetCondition(vector<Int_t>);
  // void Print();
  // inline vector<Double_t> GetSolution() const {return fSolution;};

private:

  Bool_t kGoodSolution;
  vector<Int_t> fCoilIndex; //[icoil]
  vector<Double_t>  fDetArray;//[icoil]
  vector< vector<Double_t> > fMonArray; // [icoil][imon]

  vector<Int_t> fMinimumCondition; // a minimum requirement to trigger solver
  Bool_t IsCoilSufficient();
  vector<Double_t>  fSolution;
  ClassDef(TSolver,0);
};

