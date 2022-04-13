#include <iterator>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <TLorentzVector.h>
#include "catima/structures.h"
#include "catima/catima.h"
#include "PhysicsSort.h"
using namespace std;


// some experiment-specific constants
namespace {
// M(21Ne), M(p), M(t), M(19Ne) [MeV/c^2]
const array<double,4> kReactionMasses = {19550.534,938.27203,2808.9210,17695.029};
const double kBeamEnergy = 840;
const double kAMU_p = 1.0072765;
const double kAMU_d = 2.0135532;
const double kAMU_t = 3.0155007;
}

PhysicsSort::PhysicsSort(DetectorSort& detsort):
  fDetectorSort(detsort)
{
	fMatchWindow = numeric_limits<double>().max();
																								 
  fDetectorSort.CdFile();
  fTree = new TTree("tphys", "Physics level sorted data");

	// Si ---------
  for(int i=0; i< kNumSi; ++i) {
		Si_E[i] = 0;
		Si_T[i] = 0;
		Si_Sector[i] = 0;
		Si_Ring[i] = 0;
		Si_Mult[i] = 0;
		
    int iSi = i+1;
    fTree->Branch(Form("Si%i_E",iSi), &(Si_E[i]));
    fTree->Branch(Form("Si%i_T",iSi), &(Si_T[i]));
    fTree->Branch(Form("Si%i_Sector",iSi), &(Si_Sector[i]));
    if(iSi == 1 || iSi == 3) {
      fTree->Branch(Form("Si%i_Ring",iSi), &(Si_Ring[i]));
    } else {
			Si_Ring[i] = new vector<UInt_t>();
		}
		fTree->Branch(Form("Si%i_Mult",iSi), &(Si_Mult[i]));
	}
	
	Si_E1 = 0;
	Si_E12 = 0;
	Si_E123 = 0;
	Si_Etot = 0;
	Si_Etot_p = 0;
	Si_Etot_d = 0;
	Si_Etot_t = 0;
	Si_Ecm_pt = 0;
	Si_ThetaLab = 0;
	fTree->Branch("Si_E1",  &Si_E1);
	fTree->Branch("Si_E12", &Si_E12);
	fTree->Branch("Si_E123",&Si_E123);
	fTree->Branch("Si_Etot",&Si_Etot);
	fTree->Branch("Si_Etot_p",&Si_Etot_p);
	fTree->Branch("Si_Etot_d",&Si_Etot_d);
	fTree->Branch("Si_Etot_t",&Si_Etot_t);
	fTree->Branch("Si_ThetaLab",&Si_ThetaLab);
	fTree->Branch("Si_Ecm_pt",&Si_Ecm_pt);

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
  fTree->Branch("TOF_PPAC12",&TOF_PPAC12,"TOF_PPAC12/D");
  fTree->Branch("TOF_PPAC_Phoswich",&TOF_PPAC_Phoswich,"TOF_PPAC_Phoswich/D");
  fTree->Branch("TOF_Si_PPAC",&TOF_Si_PPAC,"TOF_Si_PPAC/D");
  fTree->Branch("TOF_Si_Phoswich",&TOF_Si_Phoswich,"TOF_Si_Phoswich/D");
}

PhysicsSort::~PhysicsSort()
{ }

namespace {

template<class T> void setnan(T& t)
{ t = -sqrt(-1); }

template<class T> void doclear(T* t)
{ if(t->size()) t->clear(); }

}

void PhysicsSort::Clear()
{
  for(int i=0; i< kNumSi; ++i) {
    doclear(Si_E[i]);
    doclear(Si_T[i]);
    doclear(Si_Sector[i]);
    doclear(Si_Ring[i]);
		Si_Mult[i] = 0;
	}
	doclear(Si_E1);
	doclear(Si_E12);
	doclear(Si_E123);
	doclear(Si_Etot);
	doclear(Si_Etot_p);
	doclear(Si_Etot_d);
	doclear(Si_Etot_t);
	doclear(Si_Ecm_pt);
	doclear(Si_ThetaLab);

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

vector<PhysicsSort::Hit> PhysicsSort::ExtractHits(
	const string& detname_base,
	int chMin, int chMax, bool includeShort)
{
	vector<Hit> vHits;
	vHits.reserve(chMax - chMin + 1);
	
	for(int iCh = chMin; iCh <= chMax; ++iCh) {
		string detname = Form("%s%i",detname_base.c_str(),iCh);
		const Detector& det = fDetectorSort.GetDetectorData(detname);

		const vector<UInt_t>& raw_energies = det.GetRawEnergyHits();
		const vector<double>& energies = det.GetEnergyHits();
		const vector<double>& times    = det.GetTimeHits();
		if(energies.size() != times.size()) {
			AddEnergyTimeMismatch(detname);
		}
		else if(includeShort && det.GetEnergyShortHits().size() != energies.size()) {
			AddEnergyTimeMismatch(detname + " [SHORT/LONG MISMATCH]");
		}
		else {
			// equal hit length
			//
			const double thresh = fDetectorSort.GetChannelMap().
				GetChannelInfoByDetector(detname).fThreshold;
		
			for(size_t iHit=0; iHit< energies.size(); ++iHit) {
				if(raw_energies.at(iHit) > 0x1 && 
					 raw_energies.at(iHit) < 0xffff &&
					 energies.at(iHit) > thresh)
				{						 
					Hit h;
					h.Ch = iCh;
					h.E  = energies.at(iHit);
					h.T  = times.at(iHit);
					h.ES = includeShort ? det.GetEnergyShortHits().at(iHit) : -sqrt(-1);
					vHits.emplace_back(h);
				}
			}
		}
	}
	return vHits;
}


namespace { bool minTime(
	const PhysicsSort::Hit& l,
	const PhysicsSort::Hit& r) {
	return l.T < r.T;
}; }

UInt_t PhysicsSort::MatchRingSector(
	const PhysicsSort::Hit& hR,
	vector<PhysicsSort::Hit>& hitSector)
{
	UInt_t sectorMatch = 255;
	auto itLastMatch = partition(
		hitSector.begin(), hitSector.end(),
		[&](const Hit& h) {
			return fabs(h.E - hR.E) < fMatchWindow;
		});
	
	if(itLastMatch != hitSector.begin()) {
		auto itMinTime = min_element(
			hitSector.begin(), itLastMatch,
			[&](const Hit& l, const Hit& r) {
				return fabs(l.T - hR.T) < fabs(r.T - hR.T);
			});

		sectorMatch = itMinTime->Ch;
		hitSector.erase(itMinTime);
	}
	
	return sectorMatch;
}
				
		
//// TODO Better Hit Matching ////
void PhysicsSort::CalculateSi()
{
  for(int iSi=1; iSi<= kNumSi; ++iSi) {
    //
    // Sectors
		auto nameS = Form("Si%i_Sector", iSi);
		auto hitSector =
			ExtractHits(nameS,	1, DetectorSort::kSectors, false);
    
    //
    // Rings
    if(iSi == 1 || iSi == 3) {
			auto nameR = Form("Si%i_Ring", iSi);
			vector<Hit> hitRing =
				ExtractHits(nameR, 1, DetectorSort::kRings, false);

			// loop rings (sorted by time)
			sort(hitRing.begin(), hitRing.end(), minTime);
			for(const auto& hR : hitRing) {
				Si_E[iSi-1]->push_back(hR.E);
				Si_T[iSi-1]->push_back(hR.T);
				Si_Ring[iSi-1]->push_back(hR.Ch);						
				Si_Sector[iSi-1]->push_back(
					MatchRingSector(hR, hitSector));
			}
		}
		else { // sectors only
			// take earliest hit
			sort(hitSector.begin(), hitSector.end(), minTime);
			for(const auto& hS : hitSector) {
				Si_E[iSi-1]->push_back(hS.E);
				Si_T[iSi-1]->push_back(hS.T);
				Si_Sector[iSi-1]->push_back(hS.Ch);
			}
		}

		// multiplicity
		// (store as vector for TTree::Draw matching)
		Si_Mult[iSi-1] = Si_E[iSi-1]->size();		
	}

	for( const auto& ringNo : *(Si_Ring[0]) ) {
		Si_ThetaLab->push_back(GetSiTheta(ringNo));
	}	
	
	for(size_t iHit = 0; iHit< Si_E[0]->size(); ++iHit) {
		int depth = 0;

		if(iHit < Si_E[0]->size()) {
			depth = 1;
			Si_E1->push_back(Si_E[0]->at(iHit));
			Si_Etot->push_back(Si_E[0]->at(iHit));
		
			if(iHit < Si_E[1]->size()) {
				depth = 2;
				Si_Etot->back() += Si_E[1]->at(iHit);
				Si_E12->push_back(Si_E[0]->at(iHit) + Si_E[1]->at(iHit));

				if(iHit < Si_E[2]->size()) {
					depth = 3;
					Si_Etot->back() += Si_E[2]->at(iHit);
					Si_E123->push_back(
						Si_E[0]->at(iHit) + Si_E[1]->at(iHit) + Si_E[2]->at(iHit));
				
					if(iHit < Si_E[3]->size()) {
						depth = 4;
						Si_Etot->back() += Si_E[3]->at(iHit);
					}
				}
			}
		}

		Si_Etot_p->push_back( AddBackDeadLayers(kAMU_p,1,iHit,depth) );
		Si_Etot_d->push_back( AddBackDeadLayers(kAMU_d,1,iHit,depth) );
		Si_Etot_t->push_back( AddBackDeadLayers(kAMU_t,1,iHit,depth) );
		Si_Ecm_pt->push_back(
			CalcEcm(kReactionMasses,
							kBeamEnergy,
							Si_Etot_t->back(),
							Si_ThetaLab->at(iHit)
				)
			);
	}
}

double PhysicsSort::CalcEcm(
	const array<double,4>& M, double ebeam, double elab, double thetalab)
{
	auto ekin2p =
		[](double ekin, double M)
			{ return sqrt((ekin+M)*(ekin+M) - M*M); };

	TLorentzVector lv1(0,0,ekin2p(ebeam,M[0]),ebeam+M[0]);
	TLorentzVector lv2(0,0,0,M[1]);
	TLorentzVector lv0 = lv1+lv2;
	
	const double p3 = ekin2p(elab,M[2]);
	TLorentzVector lv3(p3*sin(thetalab),0,p3*cos(thetalab),elab+M[2]);
	TLorentzVector lv4 = lv0 - lv3;
	return lv4.M() - M[3];
}

double PhysicsSort::AddBackDeadLayers(
	double A, int Z, size_t iHit, int depth)
{
	const array<double,4> thick = {477,1537,1496,1544};
	const double dead_Al = 0.3;
	const double dead_Si = 0.5;
	const double target_thick = 12;
	
	catima::Projectile p(A,Z);
	catima::Material Si  = catima::get_material(14);
	catima::Material Al  = catima::get_material(13);
	catima::Material CH2 = catima::get_material(catima::material::CH2);
	
	auto get_energy =
		[&](int idet) {
			try { return Si_E.at(idet-1)->at(iHit); }
			catch (exception& e) {
				cout << "Problem in AddBackDeadLayers::get_energy\n";
				throw (e);
			}
		};
	double theta;
	try { theta = Si_ThetaLab->at(iHit); }
	catch (exception& e) {
		cerr << "AddBackDeadLayers:: Problem in accessing ThetaLab["<<iHit<<"]\n";
		throw e;
	}

	auto set_thickness =
		[&theta](catima::Material& m, double thick_um) {
			const double thick_cm = thick_um/1e4;
			m.thickness_cm(thick_cm/cos(theta));
		};
	
	set_thickness(Si, dead_Si);
	set_thickness(Al, dead_Al*2);	
	p.T = 0;
	for(int idet = depth; idet >0; --idet) {
		p.T += get_energy(idet)/A;
		if(idet > 1) {
			p.T = catima::energy_in(p,Si);
			p.T = catima::energy_in(p,Al);
			p.T = catima::energy_in(p,Si);
		}
		else {
			set_thickness(Al, dead_Al);
			p.T = catima::energy_in(p,Si);
			p.T = catima::energy_in(p,Al);
			break;
		}
	}

	// half target thickness
	set_thickness(CH2, 0.5*target_thick);
	p.T = catima::energy_in(p,CH2);
	
	return p.T * A;
}

void PhysicsSort::CalculateSB()
{
  for(int i=1; i<= 2; ++i) {
		auto hits = ExtractHits("SB",i,i,false);

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
		auto hitX = ExtractHits("PPAC_X", i, i, false);
		auto hitY = ExtractHits("PPAC_Y", i, i, false);
		auto hitCathode = ExtractHits("PPAC_Cathode", i, i, false);
		
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
	auto hitL = ExtractHits("Phoswich",1,1,true);
	auto hitR = ExtractHits("Phoswich",2,2,true);

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
      if(!isnan(t1) && !isnan(t2)) tof = t2-t1;
    };
  calcTOF(PPAC_T[0],PPAC_T[1],TOF_PPAC12);
  calcTOF(PPAC_T[0],Phoswich_time,TOF_PPAC_Phoswich);

	if(Si_T[0]->size()) {
		calcTOF(Si_T[0]->at(0),PPAC_T[0],TOF_Si_PPAC);
		calcTOF(Si_T[0]->at(0),Phoswich_time,TOF_Si_Phoswich);
	}
}
