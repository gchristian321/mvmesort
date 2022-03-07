#ifndef MVMESORT_hh_guard
#define MVMESORT_hh_guard
#include <map>
#include <string>
#include <array>

struct ChannelInfo {
  std::string fDetectorName;
  UInt_t fDetectorChannel;
  std::array<double, 3> fEcal;
  std::array<double, 3> fTcal;
  std::array<double, 3> fEScal;
  double fThreshold;
};

class ChannelMap {
public:
  ChannelMap(const std::string& filename);
  ChannelInfo GetChannelInfo(const std::string& moduleName, UInt_t moduleCh);
  const ChannelInfo& GetChannelInfoByDetector(const std::string& detName, UInt_t detCh) const;
  void Print() const;
  void SetWarn(bool warn) { fWarn = warn; }
private:
  bool fWarn;
  std::map<std::string, std::map<UInt_t, ChannelInfo> > m_;
  std::map<std::pair<std::string, UInt_t>, ChannelInfo> channelMaps_;
};


typedef std::map<UInt_t, std::vector<double>*> BrMap;

class Detector {
public:
  BrMap fE;
  BrMap fT;
  BrMap fES;
  BrMap fThreshold;
  const std::vector<double>& GetEnergyHits(UInt_t channel) const
  { return GetHits(channel,0); }
  const std::vector<double>& GetTimeHits(UInt_t channel) const
  { return GetHits(channel,1); }
  const std::vector<double>& GetEnergyShortHits(UInt_t channel) const
  { return GetHits(channel,2); }
  const std::vector<double>& GetThreshold(UInt_t channel) const
  { return GetHits(channel,3); }
  void Clear();
private:
  const std::vector<double>& GetHits(UInt_t channel, UInt_t which) const;
};


// Class to unpack "detector channel" level data
// Stuff like Si1, Ring1, etc.
// Also applied calibnration to the data
class DetectorSort {
public:
  static const UInt_t kRings = 23;
  static const UInt_t kSectors = 16;
public:
  DetectorSort(const std::string& fileOutName,
	       const std::string& channelMapFile,
	       bool save = false);
  ~DetectorSort();

public:
  void SetWarnChannelMap(bool warn) { fChannelMap.SetWarn(warn); }
  void Clear();
  void Fill()   { if(fSave) fTree->Fill();  }
  void Write()  { if(fSave) fTree->Write(); }
  void CdFile() { fFile.cd(); }
  void AddData(UInt_t moduleCh, const std::string& module_name,
	       const std::string& storage_name, double value);
  const Detector& GetDetectorData(const std::string& detector_name) const;
  const ChannelMap& GetChannelMap() const { return fChannelMap; }
  
private:
  std::map<std::string, Detector> fDetectorData;
  
private:
  ChannelMap fChannelMap;
  TFile fFile;
  TTree* fTree;
  bool fSave;
};

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
  
private:
  template<class Type_t>
  void CreateBranch(const string& name, Type_t*& br)
  {  br = nullptr; fTree->Branch(name.c_str(), &br);  }

  void AddEnergyTimeMismatch(const string& detname, UInt_t channel);
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
  double PPAC_X1;
  double PPAC_X2;
  double PPAC_Y1;
  double PPAC_Y2;
  double PPAC_T1;
  double PPAC_T2;
  double PPAC_E1;
  double PPAC_E2;

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
  std::map<std::pair<std::string, UInt_t>, Long64_t> fMisMatches;
};


#endif
