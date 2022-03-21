#include <iostream>
#include <sstream>
#include <fstream>
#include <TObjString.h>
#include "DetectorSort.h"
using namespace std;



ChannelMap::ChannelMap(const string& filename):
  fWarn(true)
{
  ifstream ifs(filename.c_str());
  if(!ifs.good()) {
    stringstream ss; ss << "Bad channel map file: " << filename;
    throw invalid_argument(ss.str());
  }

  string line;
  getline(ifs, line); // header
  while(getline(ifs,line)) {
    vector<string> ltok;
    {
      unique_ptr<TObjArray> tok(TString(line.c_str()).Tokenize(","));
      for(int i=0; i< tok->GetEntries(); ++i) {
				ltok.push_back (
					string(static_cast<TObjString*>(
									 tok->At(i))->GetString().Data()) );
      }
    }
    bool skip_line = false;
    bool read_line = ltok.size() == 15;
    if(read_line) {
      for(int i=0; i< ltok.size(); ++i) {
				if(ltok[i] == "") { read_line = false; }
      }
    }
    if(ltok.size() >= 3 && ltok.at(2) == "EMPTY") {
      skip_line = true;
    }

    if(read_line) {
      // module name, module ch, detector name, detector ch,
			// mode, ecal[0,1,2], thresh, tcal[0,1,2], escal[0,1,2]
      string moduleName = ltok.at(0);
      UInt_t moduleCh   = atoi(ltok.at(1).c_str());
      string detName    = ltok.at(2) + ltok.at(3);
			ADCMode mode      = static_cast<ADCMode>(atoi(ltok.at(4).c_str()));
			
      ChannelInfo chInfo;
      chInfo.fDetectorName = detName;
			chInfo.fMode = mode;
      for(int i=0; i< 3; ++i) {
				chInfo.fEcal[i]  = atof(ltok.at(i+5).c_str());
				chInfo.fTcal[i]  = atof(ltok.at(i+8).c_str());
				chInfo.fEScal[i] = atof(ltok.at(i+11).c_str());
      }
      chInfo.fThreshold = atof(ltok.at(8).c_str());
      
      auto emp = m_[moduleName].emplace(moduleCh, chInfo);
      if(!emp.second) {
				cerr << "WARNING: duplicate channel in mapping for (module, channel): "
						 <<	moduleName << ", " << moduleCh  << "\n";
				cerr << "  Previous ---> : " << emp.first->second.fDetectorName << endl;
				cerr << "  This -------> : " << detName << endl;
				cerr << "Keeping previous entry ONLY\n";
      } else {
				// also add to list for channel-level lookup
				channelMaps_.emplace(chInfo.fDetectorName, chInfo);
				channelInfo_.push_back(chInfo);
      }
    }
    else if(skip_line) {
      ;
    }
    else if(fWarn) {
      cerr << "WARNING: Skipping line in channel map file: \"" << line << "\"" << endl;
    }
  }
}

const ChannelInfo& ChannelMap::GetChannelInfoByDetector(
	const std::string& detName) const
{
  auto it = channelMaps_.find(detName);
  if(it != channelMaps_.end()) return it->second;
  else {
    throw invalid_argument
      (Form("Can't find Channel Map for detector \"%s\"",detName.c_str()));
  }
}

ChannelInfo ChannelMap::GetChannelInfo(const std::string& moduleName, UInt_t moduleCh)
{
  auto it = m_.find(moduleName);
  if (it != m_.end()) {
    auto it1 = it->second.find(moduleCh);
    if(it1 != it->second.end()) {
      return it1->second;
    }
    else if (fWarn) {
      cerr << "Can't find channel map for module " << moduleName <<
				" and channel " << moduleCh <<endl;
    }
  }
  else if (fWarn) {
    cerr << "Can't find channel map for module " << moduleName << endl;
  }
  ChannelInfo ch;
  ch.fDetectorName = "ERROR";
  return ch;
}

void ChannelMap::Print() const
{
  for (const auto& it : m_) {
    cout << "Module: " << it.first << "\n";
    for (const auto& it1 : it.second) {
      cout << "\tChannel >>detName: " << it1.first << " >>" <<
				it1.second.fDetectorName << "\n";
    }
  }
}

void Detector::Clear()
{
	if(fE && fE->size()) fE->clear();
	if(fT && fT->size()) fT->clear();
	if(fES && fES->size()) fES->clear();
	if(fChannel && fChannel->size()) fChannel->clear();
}

Detector::Detector(const std::string& name, ADCMode adcMode):
	fName(name), fADCMode(adcMode), fE(0), fT(0), fES(0), fChannel(0)
{   }

Detector::Detector(const std::string& name, ADCMode adcMode, TTree* tree):
	fName(name), fADCMode(adcMode), fE(0), fT(0), fES(0), fChannel(0)
{
	SetupBranches(tree);
}

Detector::~Detector()
{   }

void Detector::SetupBranches(TTree* tree)
{
	const Int_t bufsize = 32000;
	tree->Branch(Form("%s_E", GetName().c_str()), &fE, bufsize);
	tree->Branch(Form("%s_T", GetName().c_str()), &fT, bufsize);
	tree->Branch(Form("%s_Channel", GetName.c_str()), &fChannel, bufsize);
	if(fADCMode == kPSD) {
		tree->Branch(Form("%s_ES", GetName().c_str()), &fES, bufsize);
	}
}

DetectorSort::DetectorSort(const string& outFileName,
													 const string& channelMapFile,
													 bool save):
  fChannelMap(channelMapFile),
  fFile(outFileName.c_str(), "recreate"),
  fSave(save)
{  
  fFile.cd();
  fTree = new TTree("tdet", "Detector-level sorted events.");
	for(const auto& info : fChannelMap.GetAllChannelInfo()) {
		AddDetectorData(info.fDetectorName, info.fMode, fTree);
	}
}

DetectorSort::~DetectorSort()
{  }

void DetectorSort::AddDetectorData(
	const std::string& name, ADCMode mode, TTree* t)
{
	if(!t) {
		return AddDetectorData(make_shared<Detector>(name, mode));
	} else {
		return AddDetectorData(make_shared<Detector>(name, mode,t));
	}
}

void DetectorSort::AddDetectorData(const shared_ptr<Detector>& det)
{
	auto emp = fDetectorData.emplace(det->GetName(), det);
	if(!emp.second) {
		throw runtime_error(
			Form("Duplicate detector \"%s\"", det->GetName().c_str()));
	}			
}

Detector& DetectorSort::GetDetectorData(const string& detName)
{
  auto it = fDetectorData.find(detName);
  if(it != fDetectorData.end())  return *(it->second);
  else {
    throw invalid_argument(
			Form("DetectorSort::GetDetectorData(\"%s\");",detName.c_str()));
  }
}

void DetectorSort::Clear()
{
  for(auto& it : fDetectorData) {
    it.second->Clear();
  }
}

void DetectorSort::AddData(UInt_t moduleCh,
													 const string& module_name,
													 const string& storage_name,
													 double value)
{
  auto chInfo = fChannelMap.GetChannelInfo(module_name, moduleCh);
  if(chInfo.fDetectorName == "ERROR") { return; }

  const string det_name = chInfo.fDetectorName;
 
	Detector& detector = GetDetectorData(det_name);
  auto calibrate = 
    [&value](const std::array<double,3>& Pcal) {
			// calibration
			double cal_value =
				Pcal[2] * value * value +								\
				Pcal[1] * value  +											\
				Pcal[0];  
			return cal_value;
    };
	auto check_back_of_string =
    [](const string& s, const string& toFind) {
      if (string(s.substr(s.size() - toFind.size())) == toFind) return true;
      else return false;
    };
     
  if(check_back_of_string(storage_name, "amplitude") ||
		 check_back_of_string(storage_name, "integration_long"))
	{
		detector.AddEnergyHit(calibrate(chInfo.fEcal));
  }
  else if(check_back_of_string(storage_name, "channel_time"))
	{
		detector.AddTimeHit(calibrate(chInfo.fTcal));
  }
  else if(check_back_of_string(storage_name, "integration_short"))
	{
		detector.AddEnergyShortHit(calibrate(chInfo.fEScal));
  }
  else { // IGNORE module stuff
  }
}
