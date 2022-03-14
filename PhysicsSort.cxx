#include <iostream>
#include "PhysicsSort.h"
using namespace std;


PhysicsSort::PhysicsSort(DetectorSort& detsort):
  fDetectorSort(detsort)
{
  fDetectorSort.CdFile();
  fTree = new TTree("tphys", "Physics level sorted data");

	// Si ---------
  for(int i=0; i< kNumSi; ++i) {
		Si_E[i] = 0;
		Si_T[i] = 0;
		Si_Sector[i] = 0;
		Si_Ring[i] = 0;
		Si_RingSectorMatches[i] = 0;
		Si_ThetaLab[i] = 0;
		
    int iSi = i+1;
    fTree->Branch(Form("Si%i_E",iSi), &(Si_E[i]));
    fTree->Branch(Form("Si%i_T",iSi), &(Si_T[i]));
    fTree->Branch(Form("Si%i_Sector",iSi), &(Si_Sector[i]));
    if(iSi == 1 || iSi == 3) {
      fTree->Branch(Form("Si%i_Ring",iSi), &(Si_Ring[i]));
      fTree->Branch(Form("Si%i_RingSectorMatches",iSi),
										&(Si_RingSectorMatches[i]));
    } else {
			Si_Ring[i] = new vector<UInt_t>();
			Si_RingSectorMatches[i] = new vector<UInt_t>();
		}
		fTree->Branch("Si_ThetaLab",&(Si_ThetaLab[i]));
	}
	
	Si_E1 = 0;
	Si_E12 = 0;
	Si_E123 = 0;
	Si_Etot = 0; 
	fTree->Branch("Si_E1",  &Si_E1);
	fTree->Branch("Si_E12", &Si_E12);
	fTree->Branch("Si_E123",&Si_E123);
	fTree->Branch("Si_Etot",&Si_Etot);

  // SB ---------
  fTree->Branch("SB_dE", &SB_dE, "SB_dE/D");
  fTree->Branch("SB_E", &SB_E, "SB_E/D");

  // PPAC -------
  fTree->Branch("PPAC_X1", PPAC_X+0,"PPAC_X1/D");
  fTree->Branch("PPAC_X2", PPAC_X+1,"PPAC_X2/D");
  fTree->Branch("PPAC_Y1", PPAC_Y+0,"PPAC_Y1/D");
  fTree->Branch("PPAC_Y2", PPAC_Y+1,"PPAC_Y2/D");
  fTree->Branch("PPAC_T1", PPAC_T+0,"PPAC_T1/D");
  fTree->Branch("PPAC_T2", PPAC_T+1,"PPAC_T2/D");
  fTree->Branch("PPAC_E1", PPAC_E+0,"PPAC_E1/D");
  fTree->Branch("PPAC_E2", PPAC_E+1,"PPAC_E2/D");

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
    Si_E[i]->clear();
    Si_T[i]->clear();
    Si_Sector[i]->clear();
    Si_Ring[i]->clear();
    Si_RingSectorMatches[i]->clear();
		Si_ThetaLab[i]->clear();
	}
	Si_E1->clear();
	Si_E12->clear();
	Si_E123->clear();
	Si_Etot->clear();

  setnan(SB_dE);
  setnan(SB_E);

	for(int i=0; i< 2; ++i) {
		setnan(PPAC_X[i]);
		setnan(PPAC_Y[i]);
		setnan(PPAC_T[i]);
		setnan(PPAC_E[i]);
	}
	
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

void PhysicsSort::AddEnergyTimeMismatch(const string& detname)
{
  auto it = fMisMatches.find(detname);
  if( it != fMisMatches.end() ) {
    it->second += 1;
  }
  else {
    fMisMatches[detname] = 1;
  }
}

void PhysicsSort::PrintEnergyTimeMismatches() const
{
  cout << "Summary of energy/time mismatches...\n";
  for(const auto& it : fMisMatches) {
    printf("%s: %lli\n",
					 it.first.c_str(), it.second);
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

bool PhysicsSort::ExtractHits(
	const string& detname_base,
	int chMin, int chMax, bool includeShort,
	vector<PhysicsSort::Hit>& vHits)
{
	vHits.clear();
	
	for(int iCh = chMin; iCh <= chMax; ++iCh) {
		string detname = Form("%s%i",detname_base.c_str(),iCh);
		const Detector& det = fDetectorSort.GetDetectorData(detname);

		const vector<double>& energies = det.GetEnergyHits();
		const vector<double>& times    = det.GetTimeHits();
		if(energies.size() != times.size()) {
			AddEnergyTimeMismatch(detname);
			return false;
		}
		if(includeShort && det.GetEnergyShortHits().size() != energies.size()) {
			AddEnergyTimeMismatch(detname + " [SHORT/LONG MISMATCH]");
			return false;
		}
		// equal hit length
		//
		const double thresh = fDetectorSort.GetChannelMap().
			GetChannelInfoByDetector(detname).fThreshold;
		
		for(size_t iHit=0; iHit< energies.size(); ++iHit) {
			if(energies.at(iHit) > thresh) {

				vHits.push_back(Hit());
				Hit& h = vHits.back();
				h.Ch = iCh;
				h.E  = energies.at(iHit);
				h.T  = times.at(iHit);
				h.ES = includeShort ? det.GetEnergyShortHits().at(iHit) : -sqrt(-1);
			}
		}
	}
	return true;
}


namespace {
bool minTime(const PhysicsSort::Hit& l, const PhysicsSort::Hit& r) {
	return l.T < r.T;
}; }

//// TODO Better Hit Matching ////
void PhysicsSort::CalculateSi()
{
  int MatchWindow = 20000; // channels --> TODO make this dynamic
	
  for(int iSi=1; iSi< kNumSi; ++iSi) {
    //
    // Sectors
		vector<Hit> hitSector;
		ExtractHits(
			Form("Si%i_Sector", iSi),
			1, DetectorSort::kSectors, false, hitSector);
    
    //
    // Rings
    if(iSi == 1 || iSi == 3) {
			vector<Hit> hitRing;
      ExtractHits(
				Form("Si%i_Ring", iSi),
				1, DetectorSort::kRings, false, hitRing);

      if(hitRing.size()) {
				// ring + sector matching
				sort(hitRing.begin(), hitRing.end(), minTime);
				sort(hitSector.begin(), hitSector.end(), minTime);

				// loop rings
				for(const auto& hR : hitRing) {				
					Si_E[iSi-1]->push_back(hR.E);
					Si_T[iSi-1]->push_back(hR.T);
					Si_Ring[iSi-1]->push_back(hR.Ch);
					Si_ThetaLab[iSi-1]->push_back(GetSiTheta(hR.Ch));

					// search for sector match
					Si_RingSectorMatches[iSi-1]->push_back(0);
					Si_Sector[iSi-1]->push_back(255);
					UInt_t& num_matches = Si_RingSectorMatches[iSi-1]->back();
					UInt_t& sector_no = Si_Sector[iSi-1]->back();
					
					for(const auto& hS : hitSector) {
						if(fabs(hS.E - hR.E) < MatchWindow) {
							++num_matches;
							sector_no = num_matches == 1 ? hS.Ch : 255;
						}
					}
				}
      }
    }
    else { // sectors only
      // take earliest hit
      if(hitSector.size()) {
				sort(hitSector.begin(), hitSector.end(), minTime);
				for(const auto& h : hitSector) {
					Si_E[iSi-1]->push_back(h.E);
					Si_T[iSi-1]->push_back(h.T);
					Si_Sector[iSi-1]->push_back(h.Ch);
				}
      }
    }
  }

	for(size_t i1 = 0; i1< Si_E[0]->size(); ++i1) {
		Si_E1->push_back(Si_E[0]->at(i1));
		Si_Etot->push_back(Si_E[0]->at(i1));
		
		if(i1 < Si_E[1]->size()) {
			Si_Etot->back() += Si_E[1]->at(i1);
			Si_E12->push_back(Si_E[0]->at(i1) + Si_E[1]->at(i1));

			if(i1 < Si_E[2]->size()) {
				Si_Etot->back() += Si_E[2]->at(i1);
				Si_E123->push_back(
					Si_E[0]->at(i1) + Si_E[1]->at(i1) + Si_E[2]->at(i1));

				if(i1 < Si_E[3]->size()) {
					Si_Etot->back() += Si_E[3]->at(i1);
				}
			}
		}
	}

}

void PhysicsSort::CalculateSB()
{
  for(int i=1; i<= 2; ++i) {
		vector<Hit> hits;
		ExtractHits("SB",i,i,false,hits);

		if(hits.size()) {
			auto iMin = min_element(hits.begin(), hits.end(), minTime);
			if(i==1) SB_dE = iMin->E;
			else     SB_E  = iMin->E;
		}
    else {
      AddEnergyTimeMismatch("SB");
    }
  }
}

void PhysicsSort::CalculatePPAC()
{
  for(int i=1;i<=2;++i) {
		vector<Hit> hitX, hitY, hitCathode;
		ExtractHits("PPAC_X", i, i, false, hitX);
		ExtractHits("PPAC_Y", i, i, false, hitY);
		ExtractHits("PPAC_Cathode", i, i, false, hitCathode);
		
		if(hitCathode.size()) {
			auto iCathode = min_element(
				hitCathode.begin(), hitCathode.end(), minTime);

			PPAC_E[i-1] = iCathode->E;
			PPAC_T[i-1] = iCathode->T;
			
			if(hitX.size()) {
				auto iX = min_element(
					hitX.begin(), hitX.end(), minTime);

				PPAC_X[i-1] = iX->T - iCathode->T;
			}
			if(hitY.size()) {
				auto iY = min_element(
					hitY.begin(), hitY.end(), minTime);

				PPAC_Y[i-1] = iY->T - iCathode->T;
			}
		}
	}
}

void PhysicsSort::CalculatePhoswich()
{
	vector<Hit> hitL, hitR;
	ExtractHits("Phoswich",1,1,true,hitL);
	ExtractHits("Phoswich",2,2,true,hitR);

	if(hitL.size() && hitR.size()) {
		auto iL = min_element(hitL.begin(), hitL.end(), minTime);
		auto iR = min_element(hitR.begin(), hitR.end(), minTime);

		Phoswich_time_L   = iL->T;
		Phoswich_time_R   = iR->T;
		Phoswich_elong_L  = iL->E;
		Phoswich_elong_R  = iR->E;
		Phoswich_eshort_L = iL->ES;
		Phoswich_eshort_R = iR->ES;

		Phoswich_elong = sqrt(iL->E * iR->E);
		Phoswich_eshort = sqrt(iL->ES * iR->ES);
		Phoswich_time = 0.5 * (iL->T + iR->T);
	}
}

void PhysicsSort::CalculateCoinc()
{
  auto calcTOF =
    [](double t1,double t2,double& tof) {
      if(!isnan(t1) && !isnan(t2)) tof = t2-t2;
    };
  calcTOF(PPAC_T[0],PPAC_T[1],TOF_PPAC12);
  calcTOF(PPAC_T[0],Phoswich_time,TOF_PPAC_Phoswich);

	if(Si_T[0]->size()) {
		calcTOF(Si_T[0]->at(0),PPAC_T[0],TOF_Si_PPAC);
		calcTOF(Si_T[0]->at(0),Phoswich_time,TOF_Si_Phoswich);
	}
}
