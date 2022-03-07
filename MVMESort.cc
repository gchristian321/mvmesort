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
    bool read_line = ltok.size() == 14;
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
	chInfo.fEcal[i]  = atof(ltok.at(i+4).c_str());
	chInfo.fTcal[i]  = atof(ltok.at(i+7).c_str());
	chInfo.fEScal[i] = atof(ltok.at(i+10).c_str());
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
      } else {
	// also add to list for channel-level lookup
	channelMaps_.emplace(make_pair(chInfo.fDetectorName, chInfo.fDetectorChannel), chInfo);
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

const ChannelInfo& ChannelMap::GetChannelInfoByDetector(const std::string& detName, UInt_t detCh) const
{
  auto it = channelMaps_.find(make_pair(detName,detCh));
  if(it != channelMaps_.end()) return it->second;
  else {
    throw invalid_argument
      (Form("Can't find Channel Map for detector \"%s\" and channel %i",detName.c_str(),detCh));
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



void Detector::Clear()
{
  auto do_clear =
    [](BrMap& brmap) {
      for(auto& it : brmap) {
	if(it.second) it.second->clear();
      }
    };
  do_clear(fE);
  do_clear(fT);
  do_clear(fES);
}

const std::vector<double>& Detector::GetHits(UInt_t channel, UInt_t which) const
{
  const BrMap* pM = 0;
  switch(which) {
  case 0:
    pM = &fE;
    break;
  case 1:
    pM = &fT;
    break;
  case 2:
    pM = &fES;
    break;
  case 3:
    pM = &fThreshold;
    break;
  default:
    throw invalid_argument(Form("BAD WHICH: %i",which));
    break;
  }
  auto it = pM->find(channel);
  if(it == pM->end()){
    throw runtime_error(Form("Detector::GetHits(): Bad Channel %i", channel));
  }
  return *(it->second);
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

DetectorSort::~DetectorSort()
{  }

const Detector& DetectorSort::GetDetectorData(const string& detName) const
{
  auto it = fDetectorData.find(detName);
  if(it != fDetectorData.end())  return it->second;
  else {
    stringstream ss; ss << "DetectorSort::GetDetectorData(\"" << detName << "\")";
    throw invalid_argument(ss.str());
  }
}

void DetectorSort::Clear()
{
  for(auto& it : fDetectorData) {
    it.second.Clear();
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
  const UInt_t det_ch = chInfo.fDetectorChannel;
  auto it = fDetectorData.find(det_name);
  if(it == fDetectorData.end()) {
    cerr << "ERROR: coundn't find branch for detector \"" << det_name << "\"\n";
    return;
  }
 
  Detector& detector = it->second;
  auto try_push_back =
    [&](BrMap& brmap, const std::array<double,3>& Pcal) {
      auto it1 = brmap.find(det_ch);
      if(it1 == brmap.end()) {
	cerr << "ERROR: no detector channel number " << det_ch
	     << ", for detector " << det_name << "\n";
      }
      else {
	if(it1->second) {
	  // calibration
	  double cal_value =
	    Pcal[2] * value * value +	\
	    Pcal[1] * value  + \
	    Pcal[0];  
	  it1->second->push_back(cal_value);
	} else {
	  cerr << "ERROR: NULL vector for detector channel number " << det_ch
	       << ", for detector " << det_name << "\n";
	}
      }
    };
  auto check_back_of_string =
    [](const string& s, const string& toFind) {
      if (string(s.substr(s.size() - toFind.size())) == toFind) return true;
      else return false;
    };
     
  if(check_back_of_string(storage_name, "amplitude")) {
    try_push_back(detector.fE, chInfo.fEcal);
  }
  else if(check_back_of_string(storage_name, "channel_time")) {
    try_push_back(detector.fT, chInfo.fTcal);
  }
  else if(check_back_of_string(storage_name, "integration_long")) {
    try_push_back(detector.fE, chInfo.fEcal);
  }
  else if(check_back_of_string(storage_name, "integration_short")) {
    try_push_back(detector.fES, chInfo.fEScal);
  }
  else { // IGNORE module stuff
  }
}


PhysicsSort::PhysicsSort(DetectorSort& detsort):
  fDetectorSort(detsort)
{
  fDetectorSort.CdFile();
  fTree = new TTree("tphys", "Physics level sorted data");

  for(int i=0; i< kNumSi; ++i) {
    int iSi = i+1;
    fTree->Branch(Form("Si%i_E",iSi), Si_E+i, Form("Si%i_E/D",iSi));
    fTree->Branch(Form("Si%i_T",iSi), Si_T+i, Form("Si%i_T/D",iSi));
    fTree->Branch(Form("Si%i_Sector",iSi), Si_Sector+i, Form("Si%i_Sector/i",iSi));
    if(iSi == 1 || iSi == 3) {
      fTree->Branch(Form("Si%i_Ring",iSi), Si_Ring+i, Form("Si%i_Ring/i",iSi));
    }
    fTree->Branch("Si_E12", &Si_E12, "Si_E12/D");
    fTree->Branch("Si_E123",&Si_E123,"Si_E123/D");
    fTree->Branch("Si_Etot",&Si_Etot,"Si_Etot/D");
    fTree->Branch("Si_ThetaLab",&Si_ThetaLab,"Si_ThetaLab/D");
  }

  // SB ---------
  fTree->Branch("SB_dE", &SB_dE, "SB_dE/D");
  fTree->Branch("SB_E", &SB_E, "SB_E/D");

  // PPAC -------
  fTree->Branch("PPAC_X1", &PPAC_X1,"PPAC_X1/D");
  fTree->Branch("PPAC_X2", &PPAC_X2,"PPAC_X2/D");
  fTree->Branch("PPAC_Y1", &PPAC_Y1,"PPAC_Y1/D");
  fTree->Branch("PPAC_Y2", &PPAC_Y2,"PPAC_Y2/D");
  fTree->Branch("PPAC_T1", &PPAC_T1,"PPAC_T1/D");
  fTree->Branch("PPAC_T2", &PPAC_T2,"PPAC_T2/D");
  fTree->Branch("PPAC_E1", &PPAC_E1,"PPAC_E1/D");
  fTree->Branch("PPAC_E2", &PPAC_E2,"PPAC_E2/D");

  // Phoswich
  fTree->Branch("Phoswich_elong_L", 	&Phoswich_elong_L,      "Phoswich_elong_L/D");      
  fTree->Branch("Phoswich_elong_R", 	&Phoswich_elong_R,      "Phoswich_elong_R/D");      
  fTree->Branch("Phoswich_eshort_L", 	&Phoswich_eshort_L,     "Phoswich_eshort_L/D");     
  fTree->Branch("Phoswich_eshort_R", 	&Phoswich_eshort_R,     "Phoswich_eshort_R/D");     
  fTree->Branch("Phoswich_elong", 	&Phoswich_elong,	       "Phoswich_elong/D");          
  fTree->Branch("Phoswich_eshort", 	&Phoswich_eshort,       "Phoswich_eshort/D");        
  fTree->Branch("Phoswich_time",&Phoswich_time,"Phoswich_time/D");

  // COINC
  fTree->Branch("TOF_PPAC12",&TOF_PPAC12,"TOF_PAPC12/D");
  fTree->Branch("TOF_PPAC_Phoswich",&TOF_PPAC_Phoswich,"TOF_PPAC_Phoswich/D");
  fTree->Branch("TOF_Si_PPAC",&TOF_Si_PPAC,"TOF_Si_PPAC/D");
  fTree->Branch("TOF_Si_Phoswich",&TOF_Si_Phoswich,"TOF_Si_Phoswich/D");
}

PhysicsSort::~PhysicsSort()
{ }

namespace { template<class T> void setnan(T& t)
{ t = -sqrt(-1); } }

void PhysicsSort::Clear()
{
  for(int i=0; i< kNumSi; ++i) {
    setnan(Si_E[i]);
    setnan(Si_T[i]);
    setnan(Si_Sector[i]);
    setnan(Si_Ring[i]);
  }
  setnan(Si_E12);
  setnan(Si_E123);
  setnan(Si_Etot);
  setnan(Si_ThetaLab);

  setnan(SB_dE);
  setnan(SB_E);

  setnan(PPAC_X1);
  setnan(PPAC_X2);
  setnan(PPAC_Y1);
  setnan(PPAC_Y2);
  setnan(PPAC_T1);
  setnan(PPAC_T2);
  setnan(PPAC_E1);
  setnan(PPAC_E2);

  setnan(Phoswich_elong_L);
  setnan(Phoswich_elong_R);
  setnan(Phoswich_eshort_L);
  setnan(Phoswich_eshort_R);
  setnan(Phoswich_elong);
  setnan(Phoswich_eshort);
  setnan(Phoswich_time);

  setnan(TOF_PPAC12);
  setnan(TOF_PPAC_Phoswich);
  setnan(TOF_Si_PPAC);
  setnan(TOF_Si_Phoswich);

}

void PhysicsSort::AddEnergyTimeMismatch(const string& detname, UInt_t channel)
{
  auto key = make_pair(detname, channel);
  auto it = fMisMatches.find(key);
  if( it != fMisMatches.end() ) {
    it->second += 1;
  }
  else {
    fMisMatches[key] = 1;
  }
}

void PhysicsSort::PrintEnergyTimeMismatches() const
{
  cout << "Summary of energy/time mismatches...\n";
  for(const auto& it : fMisMatches) {
    printf("%s, %i: %lli\n",
	   it.first.first.c_str(), it.first.second, it.second);
  }
}

void PhysicsSort::Calculate()
{
  CalculateSi();
  CalculateSB();
  CalculatePPAC();
  CalculatePhoswich();
  CalculateCoinc(); 
}

//// TODO Handle Multiple Hits Better ////
void PhysicsSort::CalculateSi()
{
  int MatchWindow = 20000; // channels --> TODO make this dynamic
  auto get_hits =
    [&] (const string& detname,
	 size_t iBegin, size_t iEnd,
	 vector<double> &vE, vector<double>& vT, vector<UInt_t>& vN)
    {
      const Detector& det = fDetectorSort.GetDetectorData(detname);
      
      for(size_t iCh = iBegin; iCh <= iEnd; ++iCh) {
	// energy
	//
	const vector<double>& energies = det.GetEnergyHits(iCh);
	const vector<double>& times = det.GetTimeHits(iCh);
	if(energies.size() != times.size()) {
	  AddEnergyTimeMismatch(detname, iCh);
	  // cerr << "WARNING: Silicon energy / time size mismatch: Detector, Channel: " << 
	  //   detname << ", " << iCh << "\n";
	}
	else {
 	  double thresh =
	    fDetectorSort.GetChannelMap().
	    GetChannelInfoByDetector(detname,iCh).fThreshold;
	  for(size_t iHit=0; iHit< energies.size(); ++iHit) {
	    if(energies.at(iHit) > thresh) {
	      vE.push_back(energies.at(iHit));
	      vT.push_back(times.at(iHit));
	      vN.push_back(iCh);
	    }
	  }
	}
      }
    };

  for(int iSi=1; iSi< kNumSi; ++iSi) {
    //
    // Sectors
    vector<double> vESector, vTSector;
    vector<UInt_t> vNSector;
    get_hits(Form("Si%i_Sector", iSi),
	     1, DetectorSort::kSectors,
	     vESector, vTSector, vNSector);
    
    //
    // Rings
    if(iSi == 1 || iSi == 3) {
      vector<double> vERing, vTRing;
      vector<UInt_t> vNRing;
      get_hits(Form("Si%i_Ring", iSi),
	       1, DetectorSort::kRings,
	       vERing, vTRing, vNRing);

      if(vTRing.size()) {
	// ring + sector matching
	auto iRing = min_element(vTRing.begin(), vTRing.end()) - vTRing.begin();
	const double ringE = vERing.at(iRing);
	const double ringT = vTRing.at(iRing);
	const UInt_t ringN = vNRing.at(iRing);

	vector<Int_t> ssort(vTSector.size());
	TMath::Sort(Int_t(vTSector.size()), &(vTSector[0]), &(ssort[0]), false);
	for(size_t iSector = 0; iSector< vESector.size(); ++iSector) {
	  auto iSectorSorted = ssort.at(iSector);
	  if(fabs(vESector.at(iSectorSorted) - ringE) < MatchWindow) {
	    Si_E[iSi-1] = ringE;
	    Si_T[iSi-1] = ringT;
	    Si_Ring[iSi-1] = ringN;
	    Si_Sector[iSi-1] = vNSector.at(iSectorSorted);
	    break;
	  }
	}
      }
    }
    else { // sectors only
      // take earliest hit
      if(vTSector.size()) {
	auto iMin = min_element(vTSector.begin(), vTSector.end()) - vTSector.begin();
	
	Si_E[iSi-1] = vESector.at(iMin);
	Si_T[iSi-1] = vTSector.at(iMin);
	Si_Sector[iSi-1] = vNSector.at(iMin);
      }
    }
  }
  
  if(!isnan(Si_E[0])) {
    Si_Etot = Si_E[0];
    Si_ThetaLab = GetSiTheta(Si_Ring[0]);
    
    if(!isnan(Si_E[1])) {
      Si_Etot += Si_E[1];
      Si_E12 = Si_Etot;

      if(!isnan(Si_E[2])) {
	Si_Etot += Si_E[2];
	Si_E123 = Si_Etot;

	if(!isnan(Si_E[3])) {
	  Si_Etot += Si_E[3];
	}
      }
    }
  }    
}

double PhysicsSort::GetSiTheta(UInt_t ring)
{
  // TODO --> Use actual values
  // these are just micron placeholders
  // (and maybe even wrong at that)
  const double kDist = 150; // 150 mm 
  const double rInner = 23.06/2;
  const double rOuter = 70./2;
  const double pitch = (rOuter-rInner)/48;
  double rMid;
  if(ring < 1 || ring > DetectorSort::kRings) {
    rMid = -sqrt(-1);
  }
  else if(ring == 1) { // middle ring 4
    rMid = rInner + pitch*3 + pitch/2;
  }
  else {
    rMid = rInner + pitch*4 + 2*pitch*(ring-1) + pitch/2;
  }
  return atan(rMid/kDist);
}

void PhysicsSort::CalculateSB()
{
  const Detector& det = fDetectorSort.GetDetectorData("SB");

  for(int i=1; i<= 2; ++i) {
    auto vE = det.GetEnergyHits(i);
    auto vT = det.GetTimeHits(i);
  
    if(vE.size() == vT.size()) {
      if(!vE.empty()) {
	auto iMin = min_element(vT.begin(), vT.end()) - vT.begin();
	if(i==1) SB_dE = vE.at(iMin);
	else     SB_E  = vE.at(iMin);
      }
    }
    else {
      AddEnergyTimeMismatch("SB",i);
    }
  }
}

void PhysicsSort::CalculatePPAC()
{
  auto get_data =
    [&](const string& detname, int iCh,
	double& E, double& T) 
    {
      auto det = fDetectorSort.GetDetectorData(detname);
      auto vE = det.GetEnergyHits(iCh);
      auto vT = det.GetTimeHits(iCh);
      bool returnval = false;
      if(vE.size() == vT.size()) {
	if(!vT.empty()) {
	  auto iMin = min_element(vT.begin(), vT.end()) - vT.begin();
	  E = vE.at(iMin);
	  T = vT.at(iMin);
	  returnval = true;
	}
      } else {
	AddEnergyTimeMismatch(detname,iCh);
      }
      return returnval;
    };

  for(int i=1;i<=2;++i) {
    double EX,TX,EY,TY,EC,TC;
    bool have_X = get_data("PPAC_X",i,EX,TX);
    bool have_Y = get_data("PPAC_Y",i,EY,TY);
    bool have_C = get_data("PPAC_E",i,EC,TC); // cathode

    if(have_C) {
      if(i==1) { PPAC_E1 = EC; PPAC_T1 = TC; }
      else     { PPAC_E2 = EC; PPAC_T2 = TC; }
	
      if(have_X) {
	if(i==1) PPAC_X1 = TX - TC;
	else     PPAC_X2 = TX - TC;
      }
      if(have_Y) {
	if(i==1) PPAC_Y1 = TY - TC;
	else     PPAC_Y2 = TY - TC;
      }
    }
  }
}

void PhysicsSort::CalculatePhoswich()
{
  auto get_data =
    [&](int iCh, double& E, double& T, double& ES) 
    {
      auto det = fDetectorSort.GetDetectorData("Phoswich");
      auto vE = det.GetEnergyHits(iCh);
      auto vT = det.GetTimeHits(iCh);
      auto vES= det.GetEnergyShortHits(iCh);
      bool returnval = false;
      if(vE.size() == vT.size() && vES.size() == vE.size()) {
	if(!vT.empty()) {
	  auto iMin = min_element(vT.begin(), vT.end()) - vT.begin();
	  E = vE.at(iMin);
	  T = vT.at(iMin);
	  ES= vES.at(iMin);
	  returnval = true;
	}
      } else {
	AddEnergyTimeMismatch("Phoswich",iCh);
      }
      return returnval;
    };

  get_data(1, Phoswich_elong_L, Phoswich_time_L, Phoswich_eshort_L);
  get_data(2, Phoswich_elong_R, Phoswich_time_R, Phoswich_eshort_R);

  if(!isnan(Phoswich_elong_L) && !isnan(Phoswich_elong_R)) {
    Phoswich_elong = sqrt(Phoswich_elong_L*Phoswich_elong_R);
  }
  if(!isnan(Phoswich_time_L)  && !isnan(Phoswich_time_R)) {
    Phoswich_time = (Phoswich_time_L + Phoswich_time_R)/2;
  }
  if(!isnan(Phoswich_eshort_L) && !isnan(Phoswich_eshort_R)) {
    Phoswich_eshort = sqrt(Phoswich_eshort_L*Phoswich_eshort_R);
  }
}

void PhysicsSort::CalculateCoinc()
{
  auto calcTOF =
    [](double t1,double t2,double& tof) {
      if(!isnan(t1) && !isnan(t2)) tof = t2-t2;
    };
  calcTOF(PPAC_T1,PPAC_T2,TOF_PPAC12);
  calcTOF(PPAC_T1,Phoswich_time,TOF_PPAC_Phoswich);
  calcTOF(Si_T[0],PPAC_T1,TOF_Si_PPAC);
  calcTOF(Si_T[0],Phoswich_time,TOF_Si_Phoswich);
}
