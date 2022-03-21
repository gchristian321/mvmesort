#ifndef DETECTORSORT_hh_guard
#define DETECTORSORT_hh_guard
#include <map>
#include <string>
#include <array>
#include <memory>

// ROOT
#include <TFile.h>
#include <TTree.h>

enum ADCMode {
	kPHA = 0,
	kPSD = 1
};

struct ChannelInfo {
  std::string fDetectorName;
  std::array<double, 3> fEcal;
  std::array<double, 3> fTcal;
  std::array<double, 3> fEScal;
  double fThreshold;
	ADCMode fMode;
};

class ChannelMap {
public:
  ChannelMap(const std::string& filename);
  ChannelInfo GetChannelInfo(const std::string& moduleName, UInt_t moduleCh);
  const ChannelInfo& GetChannelInfoByDetector(const std::string& detName) const;
	const std::vector<ChannelInfo>& GetAllChannelInfo() const
		{ return channelInfo_; }
  void Print() const;
  void SetWarn(bool warn) { fWarn = warn; }
private:
  bool fWarn;
	// m_ -->>  <module name, <moduleCh, ChInfo> >
  std::map<std::string, std::map<UInt_t, ChannelInfo> > m_;
  std::map<std::string, ChannelInfo> channelMaps_;
	std::vector<ChannelInfo> channelInfo_;
};


typedef std::map<UInt_t, std::vector<double>*> BrMap;

class Detector {
public:
	Detector(const std::string& name, ADCMode type);
	Detector(const std::string& name, ADCMode type, TTree*);
	~Detector();

public:
	void AddEnergyHit(double e)
		{ if(fE) fE->push_back(e); else throw std::runtime_error("(AddHit) fE is null"); }
	void AddTimeHit(double t)
		{ if(fT) fT->push_back(t); else throw std::runtime_error("(AddHit) fT is null"); }
	void AddEnergyShortHit(double es)
		{ if(fES) fES->push_back(es); else throw std::runtime_error("(AddHit) fES is null"); }
	
	const std::vector<double>& GetEnergyHits() const
		{ if(fE) return *fE; else throw std::runtime_error("(GetHit) fE is null"); }
	const std::vector<double>& GetTimeHits() const
		{ if(fT) return *fT; else throw std::runtime_error("(GetHit) fT is null"); }
	const std::vector<double>& GetEnergyShortHits() const
		{ if(fES) return *fES; throw std::runtime_error("(GetHit) fES is null");   }
	const std::string& GetName() const { return fName; }
	
	void Clear();
	void SetupBranches(TTree* tree);
	
private:
	std::string fName;
	ADCMode fADCMode;
	std::vector<double>* fE;
	std::vector<double>* fT;
	std::vector<double>* fES;
};


// Class to unpack "detector channel" level data
// Stuff like Si1, Ring1, etc.
// Also applied calibration to the data
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
  void Fill()   { if(fSave) fTree->Fill(); }
  void Write()  {	if(fSave) fTree->Write(); }
  void CdFile() { fFile.cd(); }
  void AddData(UInt_t moduleCh, const std::string& module_name,
	       const std::string& storage_name, double value);
  Detector& GetDetectorData(const std::string& detector_name);
  const ChannelMap& GetChannelMap() const { return fChannelMap; }
  
private:
	void AddDetectorData(const std::shared_ptr<Detector>& det);
	void AddDetectorData(const std::string& name, ADCMode type, TTree* t=0);	
  
private:
  std::map<std::string, std::shared_ptr<Detector> > fDetectorData;
  ChannelMap fChannelMap;
  TFile fFile;
  TTree* fTree;
  bool fSave;
};


#endif
