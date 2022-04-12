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
			// 0          , 1        , 2            , 3          , 4   , 5,6,7      , 8     , 9,10,11    , 12,13,14
      // module name, module ch, detector name, detector ch, mode, ecal[0,1,2], thresh, tcal[0,1,2], escal[0,1,2]
      string moduleName = ltok.at(0);
      UInt_t moduleCh   = atoi(ltok.at(1).c_str());
      string detName    = ltok.at(2) + ltok.at(3);
			ADCMode mode      = static_cast<ADCMode>(atoi(ltok.at(4).c_str()));
			
      ChannelInfo chInfo;
      chInfo.fDetectorName = detName;
			chInfo.fMode = mode;
      for(int i=0; i< 3; ++i) {
				chInfo.fEcal[i]  = atof(ltok.at(i+5).c_str()); //5,6,7
				chInfo.fTcal[i]  = atof(ltok.at(i+9).c_str()); //8,9,10
				chInfo.fEScal[i] = atof(ltok.at(i+12).c_str());//11,12,13
      }
      chInfo.fThreshold = atof(ltok.at(8).c_str());

			// Module Name,Module Channel,Detector Name,Detector Channel,"Mode (0 ADC, 1 PSD)",Energy Cal P0,Energy Cal P1,Energy Cal P2,Energy Threshold,Time Cal P0,Time Cal P1,Time Cal P2,Short Cal P0,Short Cal P1,Short Cal P2
			// mdpp32_1_scp_Sector_Si1And2,0,Si1_Sector,1,0,268.693,0.17483,0,0,0,1,0,0,1,0

			
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

	if(fEraw && fEraw->size())  fEraw->clear();
	if(fTraw && fTraw->size())  fTraw->clear();
	if(fESraw && fESraw->size()) fESraw->clear();
}

Detector::Detector(const std::string& name, ADCMode adcMode):
	fName(name), fADCMode(adcMode), fE(0), fT(0), fES(0)
{   }

Detector::Detector(const std::string& name, ADCMode adcMode, TTree* tree):
	fName(name), fADCMode(adcMode), fE(0), fT(0), fES(0)
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
	tree->Branch(Form("%s_Eraw", GetName().c_str()), &fEraw, bufsize);
	tree->Branch(Form("%s_Traw", GetName().c_str()), &fTraw, bufsize);
	if(fADCMode == kPSD) {
		tree->Branch(Form("%s_ES", GetName().c_str()), &fES, bufsize);
		tree->Branch(Form("%s_ESraw", GetName().c_str()), &fESraw, bufsize);
	}
}

namespace { inline void add_hit_(
	UInt_t rawval,
	vector<double>* ve,
	vector<UInt_t>* ver,
	const array<double,3>& Pcal,
	const string& what)
{
	if(!ve || !ver) throw runtime_error("(AddHit) "+what+"(raw) is null");
	const double cal_value =
		Pcal[2] * rawval * rawval +	\
		Pcal[1] * rawval  +				\
		Pcal[0];  
	ver->push_back(rawval);
	ve->push_back(cal_value);
} }

void Detector::AddEnergyHit(
	UInt_t eraw, const std::array<double,3>& Pcal)
{ add_hit_(eraw,fE,fEraw,Pcal,"fE"); }

void Detector::AddTimeHit(
	UInt_t traw, const std::array<double,3>& Pcal)
{ add_hit_(traw,fT,fTraw,Pcal,"fT"); }

void Detector::AddEnergyShortHit(
	UInt_t esraw, const std::array<double,3>& Pcal)
{ add_hit_(esraw,fES,fESraw,Pcal,"fES"); }


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

	auto check_back_of_string =
    [](const string& s, const string& toFind) {
      if (string(s.substr(s.size() - toFind.size())) == toFind) return true;
      else return false;
    };

	UInt_t valueInt = UInt_t(round(value));
	
  if(check_back_of_string(storage_name, "amplitude") ||
		 check_back_of_string(storage_name, "integration_long"))
	{ 
		detector.AddEnergyHit(valueInt, chInfo.fEcal);
  }
  else if(check_back_of_string(storage_name, "channel_time"))
	{ 
		detector.AddTimeHit(valueInt, chInfo.fTcal);
  }
  else if(check_back_of_string(storage_name, "integration_short"))
	{ 
		detector.AddEnergyShortHit(valueInt, chInfo.fEScal);
  }
  else { } // IGNORE module stuff
}
