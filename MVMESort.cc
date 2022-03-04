#include "MVMESort.hh"

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
	ltok.push_back ( string(static_cast<TObjString*>(tok->At(i))->GetString().Data()) );
      }
    }
    bool skip_line = false;
    bool read_line = ltok.size() == 8;
    if(read_line) {
      for(int i=0; i< ltok.size(); ++i) {
	if(ltok[i] == "") { read_line = false; }
      }
    }
    if(ltok.size() >= 3 && ltok.at(2) == "EMPTY") {
      skip_line = true;
    }

    if(read_line) {
      // module name, module ch, detector name, detector ch, cal0, cal1, cal2, thresh
      string moduleName = ltok.at(0);
      UInt_t moduleCh   = atoi(ltok.at(1).c_str());
      string detName    = ltok.at(2);
      UInt_t detCh      = atoi(ltok.at(3).c_str());
      
      ChannelInfo chInfo;
      chInfo.fDetectorName = detName;
      chInfo.fDetectorChannel = detCh;
      for(int i=0; i< 3; ++i) {
	chInfo.fCalibration[i] = atof(ltok.at(i+4).c_str());
      }
      chInfo.fThreshold = atof(ltok.at(7).c_str());
      
      auto emp = m_[moduleName].emplace(moduleCh, chInfo);
      if(!emp.second) {
	cerr << "WARNING: duplicate channel in mapping for (module, channel): " <<
	  moduleName << ", " << moduleCh  << "\n";
	cerr << "  Previous ---> (detector name, channel): " << emp.first->second.fDetectorName
	     << ", " <<  emp.first->second.fDetectorChannel << endl;
	cerr << "  This -------> (detector name, channel): " << detName << ", " << detCh << endl;
	cerr << "Keeping previous entry ONLY\n";
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
  ch.fDetectorChannel = 0xffffffff;
  return ch;
}

void ChannelMap::Print() const
{
  for (const auto& it : m_) {
    cout << "Module: " << it.first << "\n";
    for (const auto& it1 : it.second) {
      cout << "\tChannel >>detName, detCh: " << it1.first << " >>" <<
	it1.second.fDetectorName << ", " << it1.second.fDetectorChannel << "\n";
    }
  }
}


MVMESort::MVMESort(const string& outFileName,
		   const string& channelMapFile,
		   bool save):
  fChannelMap(channelMapFile),
  fFile(outFileName.c_str(), "recreate"),
  fSave(save)
{  
  fFile.cd();
  fTree = new TTree("tdet", "Detector-level sorted events.");

  auto CreateBranch =
    [&](const string& name, vector<double>*& v) {
      v = nullptr;
      fTree->Branch(name.c_str(), &v);
    };
  
  ///  CREATE SI BRANCHES ///
  for (int iSi = 1; iSi<= 4; ++iSi) {
    string sectorName = Form("Si%i_Sector",iSi);
    string ringName = Form("Si%i_Ring",iSi);

    fDetectorData.emplace(sectorName, Detector());
    for(UInt_t i=1; i<= kSectors; ++i) {
      vector<double> *vE = 0, *vT = 0;
      
      fDetectorData[sectorName].fE.emplace(i,vE);
      fDetectorData[sectorName].fT.emplace(i,vT);

      CreateBranch(Form("%s%i_E", sectorName.c_str(),i), fDetectorData[sectorName].fE[i]);
      CreateBranch(Form("%s%i_T", sectorName.c_str(),i), fDetectorData[sectorName].fT[i]);
    }

    if(iSi == 1 || iSi == 3) {
      fDetectorData.emplace(ringName, Detector());
      for(UInt_t i=1; i<= kRings; ++i) {
	vector<double> *vE = 0, *vT = 0;

	fDetectorData[ringName].fE.emplace(i,vE);
	fDetectorData[ringName].fT.emplace(i,vT);

	CreateBranch(Form("%s%i_E", ringName.c_str(),i), fDetectorData[ringName].fE[i]);
	CreateBranch(Form("%s%i_T", ringName.c_str(),i), fDetectorData[ringName].fT[i]);
      }
    }
  }

  // CREATE SB BRANCHES
  // 1--> DE; 2-->E
  for(int i=1; i<= 2; ++i) {
    vector<double> *vE=0, *vT=0;
    string brName = "SB";
    fDetectorData[brName].fE.emplace(i, vE);
    fDetectorData[brName].fT.emplace(i, vT);

    CreateBranch(Form("%s%i_E",brName.c_str(),i), fDetectorData[brName].fE[i]);
    CreateBranch(Form("%s%i_T",brName.c_str(),i), fDetectorData[brName].fT[i]);    
  }
  
  // CREATE PPAC BRANCHES
  for(int i=1; i<= 2; ++i) {
    for(auto& what : {"X","Y","E"}) {
	vector<double> *vE=0, *vT=0;
	string brName = Form("PPAC_%s",what);
	fDetectorData[brName].fE.emplace(i, vE);
	fDetectorData[brName].fT.emplace(i, vT);

	string what1 = what; if(what1 == "E") what1 = "Cathode";
	CreateBranch(Form("PPAC%i_%s_E",i,what1.c_str()), fDetectorData[brName].fE[i]);
	CreateBranch(Form("PPAC%i_%s_T",i,what1.c_str()), fDetectorData[brName].fT[i]);
    }
  }

  // CREATE PHOSWICH BRANCHES
  // 1--> Left; 2--> Right
  for(int iP=1; iP<= 2; ++iP) {
    vector<double> *vE=0, *vT=0, *vES=0;
    string brName = "Phoswich";
    fDetectorData[brName].fE.emplace(iP, vE);
    fDetectorData[brName].fT.emplace(iP, vT);
    fDetectorData[brName].fT.emplace(iP, vES);

    string side = iP == 1 ? "L" : "R";
    CreateBranch(Form("Phoswich_%s_E",  side.c_str()), fDetectorData[brName].fE[iP]);
    CreateBranch(Form("Phoswich_%s_T",  side.c_str()), fDetectorData[brName].fT[iP]);
    CreateBranch(Form("Phoswich_%s_ES", side.c_str()), fDetectorData[brName].fES[iP]);
  }
}

MVMESort::~MVMESort()
{  }

void MVMESort::Clear()
{
  auto do_clear =
    [](BrMap& brmap) {
      for(auto& it : brmap) {
	if(it.second) it.second->clear();
      }
    };
          
  for(auto& it : fDetectorData) {
    do_clear(it.second.fE);
    do_clear(it.second.fT);
    do_clear(it.second.fES);
  }
}

void MVMESort::AddData(UInt_t moduleCh,
		       const string& module_name,
		       const string& storage_name,
		       double value)
{
  auto chInfo = fChannelMap.GetChannelInfo(module_name, moduleCh);
  if(chInfo.fDetectorName == "ERROR") { return; }

  const string det_name = chInfo.fDetectorName;
  const UInt_t det_ch = chInfo.fDetectorChannel;
  auto it = fDetectorData.find(det_name);
  if(it == fDetectorData.end()) {
    cerr << "ERROR: coundn't find branch for detector \"" << det_name << "\"\n";
    return;
  }

  // calibration
  double cal_value =
    chInfo.fCalibration[2] * value * value + \
    chInfo.fCalibration[1] * value  + \
    chInfo.fCalibration[0];    
  
  Detector& detector = it->second;
  auto try_push_back =
    [&](BrMap& brmap) {
      auto it1 = brmap.find(det_ch);
      if(it1 == brmap.end()) {
	cerr << "ERROR: no detector channel number " << det_ch
	     << ", for detector " << det_name << "\n";
      }
      else {
	if(it1->second) {
	  it1->second->push_back(cal_value);
	} else {
	  cerr << "ERROR: NULL vector for detector channel number " << det_ch
	       << ", for detector " << det_name << "\n";
	}
      }
    };
  auto check_back =
    [](const string& s, const string& toFind) {
      if (string(s.substr(s.size() - toFind.size())) == toFind) return true;
      else return false;
    };
     
  if(check_back(storage_name, "amplitude")) {
    try_push_back(detector.fE);
  }
  else if(check_back(storage_name, "channel_time")) {
    try_push_back(detector.fT);
  }
  else {
    ;
  }
}
