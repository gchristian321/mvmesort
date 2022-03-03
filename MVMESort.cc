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
    bool isgood = ltok.size() == 7;
    if(isgood) {
      for(int i=0; i< 7; ++i) {
	if(ltok[i] == "") { isgood = false; }
      }
    }
    if(isgood) {
      // module name, module ch, detector name, detector ch
      UInt_t moduleCh = atoi(ltok[1].c_str()), detCh = atoi(ltok[3].c_str());
      string moduleName = ltok[0], detName = ltok[2];

      ChannelInfo chInfo;
      chInfo.fDetectorName = detName;
      chInfo.fDetectorChannel = detCh;
      for(int i=0; i< 3; ++i) {
	chInfo.fCalibration[i] = atof(ltok.at(i+4).c_str());
      }
      
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
		   bool saveRaw):
  fChannelMap(channelMapFile),
  fFile(outFileName.c_str(), "recreate"),
  fSaveRaw(saveRaw)
{  
  fFile.cd();
  fTree = new TTree("t", "Sorted MVME Tree");

  ///  CREATE SI BRANCHES ///
  for (int iSi = 1; iSi<= 4; ++iSi) {
    string sectorName = Form("Si%i_Sector",iSi);
    string ringName = Form("Si%i_Ring",iSi);

    fBrMap.emplace(sectorName, Detector());
    for(UInt_t i=1; i<= kSectors; ++i) {
      vector<double> *vE = 0, *vT = 0;
      
      fBrMap[sectorName].fE.emplace(i,vE);
      fBrMap[sectorName].fT.emplace(i,vT);

      InitVector(Form("%s%i_E", sectorName.c_str(),i), fBrMap[sectorName].fE[i]);
      InitVector(Form("%s%i_T", sectorName.c_str(),i), fBrMap[sectorName].fT[i]);
    }

    if(iSi == 1 || iSi == 3) {
      fBrMap.emplace(ringName, Detector());
      for(UInt_t i=1; i<= kRings; ++i) {
	vector<double> *vE = 0, *vT = 0;

	fBrMap[ringName].fE.emplace(i,vE);
	fBrMap[ringName].fT.emplace(i,vT);

	InitVector(Form("%s%i_E", ringName.c_str(),i), fBrMap[ringName].fE[i]);
	InitVector(Form("%s%i_T", ringName.c_str(),i), fBrMap[ringName].fT[i]);
      }
    }
  }

  // CREATE SB BRANCHES
  // 1--> DE; 2-->E
  for(int i=1; i<= 2; ++i) {
    vector<double> *vE=0, *vT=0;
    string brName = "SB";
    fBrMap[brName].fE.emplace(i, vE);
    fBrMap[brName].fT.emplace(i, vT);

    InitVector(Form("%s%i_E",brName.c_str(),i), fBrMap[brName].fE[i]);
    InitVector(Form("%s%i_T",brName.c_str(),i), fBrMap[brName].fT[i]);    
  }
  
  // CREATE PPAC BRANCHES
  for(int i=1; i<= 2; ++i) {
    for(auto& what : {"X","Y","E"}) {
	vector<double> *vE=0, *vT=0;
	string brName = Form("PPAC_%s",what);
	fBrMap[brName].fE.emplace(i, vE);
	fBrMap[brName].fT.emplace(i, vT);

	string what1 = what; if(what1 == "E") what1 = "Cathode";
	InitVector(Form("PPAC%i_%s_E",i,what1.c_str()), fBrMap[brName].fE[i]);
	InitVector(Form("PPAC%i_%s_T",i,what1.c_str()), fBrMap[brName].fT[i]);
    }
  }

  // CREATE PHOSWICH BRANCHES
  // 1--> Left; 2--> Right
  for(int iP=1; iP<= 2; ++iP) {
    vector<double> *vE=0, *vT=0, *vES=0;
    string brName = "Phoswich";
    fBrMap[brName].fE.emplace(iP, vE);
    fBrMap[brName].fT.emplace(iP, vT);
    fBrMap[brName].fT.emplace(iP, vES);

    string side = iP == 1 ? "L" : "R";
    InitVector(Form("Phoswich_%s_E",  side.c_str()), fBrMap[brName].fE[iP]);
    InitVector(Form("Phoswich_%s_T",  side.c_str()), fBrMap[brName].fT[iP]);
    InitVector(Form("Phoswich_%s_ES", side.c_str()), fBrMap[brName].fES[iP]);
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
          
  for(auto& it : fBrMap) {
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
  auto it = fBrMap.find(det_name);
  if(it == fBrMap.end()) {
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
