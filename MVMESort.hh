#ifndef MVMESORT_hh_guard
#define MVMESORT_hh_guard
#include <map>
#include <string>
#include <array>

typedef std::map<UInt_t, std::vector<double>*> BrMap;

struct ChannelInfo {
  std::string fDetectorName;
  UInt_t fDetectorChannel;
  std::array<double, 3> fCalibration;  
  double fThreshold;
};

class ChannelMap {
public:
  ChannelMap(const std::string& filename);
  ChannelInfo GetChannelInfo(const std::string& moduleName, UInt_t moduleCh);
  void Print() const;
  void SetWarn(bool warn) { fWarn = warn; }
private:
  bool fWarn;
  std::map<std::string, std::map<UInt_t, ChannelInfo> > m_;
};

struct Detector {
  BrMap fE;
  BrMap fT;
  BrMap fES;
};


// Class to unpack "detector channel" level data
// Stuff like Si1, Ring1, etc.
// Also applied calibnration to the data
class MVMESort {
public:
  static const UInt_t kRings = 23;
  static const UInt_t kSectors = 16;
public:
  MVMESort(const std::string& fileOutName,
	   const std::string& channelMapFile,
	   bool save = false);
  ~MVMESort();

public:
  void SetWarnChannelMap(bool warn) { fChannelMap.SetWarn(warn); }
  void Clear();
  void Fill()  { if(fSave) fTree->Fill();  }
  void Write() { if(fSave) fTree->Write(); }
  void AddData(UInt_t moduleCh, const std::string& module_name,
	       const std::string& storage_name, double value);

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
  PhysicsSort();
  ~PhysicsSort();
};


#endif
