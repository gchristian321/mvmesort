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

class MVMESort {
public:
  static const UInt_t kRings = 23;
  static const UInt_t kSectors = 16;
public:
  MVMESort(const std::string& fileOutName,
	   const std::string& channelMapFile,
	   bool saveRaw = false);
  ~MVMESort();

public:
  void SetWarnChannelMap(bool warn) { fChannelMap.SetWarn(warn); }
  void Clear();
  void Fill()  { fTree->Fill();  }
  void Write() { fTree->Write(); }
  void AddData(UInt_t moduleCh, const std::string& module_name,
	       const std::string& storage_name, double value);
private:
  template<class T> void InitVector(const std::string& brname, std::vector<T>*& v)
  {
    v = nullptr;
    if(fSaveRaw) {
      fTree->Branch(brname.c_str(), &v);
    }
    else {
      v = new vector<T>();
    }
    if(!v) {
      throw runtime_error("Problem in MVMESort::InitVector");
    }
  }
  
private:
  std::map<std::string, Detector> fBrMap; // holds raw data
  
private:
  ChannelMap fChannelMap;
  TFile fFile;
  TTree* fTree;
  bool fSaveRaw;
};




#endif
