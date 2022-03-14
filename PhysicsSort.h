#ifndef PHYSICSSORT_hh_guard
#define PHYSICSSORT_hh_guard
#include <map>
#include <string>
#include <array>

// ROOT
#include <TFile.h>
#include <TTree.h>

#include "DetectorSort.h"


// Class to sort "physics level" data
// More high-level parameters than the detector data
class PhysicsSort {
public:
  static const UInt_t kNumSi = 4;
  
public:
  PhysicsSort(DetectorSort& detsort);
  ~PhysicsSort();

public:
  void Clear();
  void Fill()   { fTree->Fill();  }
  void Write()  { fTree->Write(); }
  void Calculate();
  void PrintEnergyTimeMismatches() const;
	struct Hit { double E, ES, T;	UInt_t Ch; };
	
private:
  template<class Type_t>
  void CreateBranch(const std::string& name, Type_t*& br)
  {  br = nullptr; fTree->Branch(name.c_str(), &br);  }

	bool ExtractHits(
		const std::string& detname_base,
		int chMin, int chMax, bool inclueShorts,
		std::vector<Hit>& vHits);

  void AddEnergyTimeMismatch(const std::string& detname);
  void CalculateSi();
  void CalculateSB();
  void CalculatePPAC();
  void CalculatePhoswich();
  void CalculateCoinc();
  double GetSiTheta(UInt_t ring);
  
private:
  // SI PARAMS
  double Si_E[kNumSi];
  double Si_T[kNumSi];
  UInt_t Si_Sector[kNumSi];
  UInt_t Si_Ring[kNumSi];
  double Si_E12;
  double Si_E123;
  double Si_Etot;
  double Si_ThetaLab;

  // SB PARAMS
  double SB_dE;
  double SB_E;

  // PPAC
  double PPAC_X[2];
  double PPAC_Y[2];
	double PPAC_T[2];
	double PPAC_E[2];

  // PHOSWICH
  double Phoswich_elong_L;
  double Phoswich_elong_R;
  double Phoswich_eshort_L;
  double Phoswich_eshort_R;
  double Phoswich_time_L;
  double Phoswich_time_R;
  double Phoswich_elong;
  double Phoswich_eshort;
  double Phoswich_time; 

  // COINCIDENCE ----
  double TOF_PPAC12;
  double TOF_PPAC_Phoswich;
  double TOF_Si_PPAC;
  double TOF_Si_Phoswich;
  
private:
  DetectorSort& fDetectorSort;
  TTree* fTree;
  std::map<std::string, Long64_t> fMisMatches;
};


#endif
