#include <algorithm>
#include <iostream>
#include <numeric>
#include <fstream>
#include <string>
#include <memory>
#include <array>
#include <map>

#include <TSpline.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;

//\todo DEAD LAYERS

namespace {

const array<double,4> Thick = {
	500, 1500, 1500, 1500 // Si thicknesses, micron
};
const array<double,5> AMU = {
	1.00728, 2.01355, 3.0155, 3.01493, 4.00151
};
const array<string,5> particles = {
	"p","d","t","3He","4He"
};

// string get_particleName(const string& filename) {
// 	auto pos_ = filename.rfind("_")+1;
// 	string particleName =
// 		filename.substr(pos_, filename.rfind(".") - pos_);
// 	return particleName;
// }

double get_mass(const string& particleName) {
	auto it = find(particles.begin(),particles.end(),particleName);
	if(it == particles.end()) {
		throw invalid_argument(
			"Bad particle name: " + particleName);
	}
	return AMU.at(it - particles.begin());
}

unique_ptr<TSpline3> read_file(const string& particleName)
{
	const string filename =
		string(EXECPATH) + "/Si_dedx/" +  "dedx_si_" + particleName + ".txt";
	
	auto bail =
		[&filename](const TString& what) {
			throw invalid_argument(
				"Bad LISE++ energy loss file: \""+filename+"\": " + what);
		};

	ifstream ifs(filename);
	if(!ifs.good()) bail("ifs not good");
	double M = get_mass(particleName);

	vector<double> e,dedx;
	array<double,12> v;
	string header; getline(ifs,header);
	while(ifs>>v[0]>>v[1]>>v[2]>>v[3]>>v[4]>>v[5]>>v[6]>>v[7]>>v[8]>>v[9]>>v[10]>>v[11]) {
		e.push_back(v[4] * M); // -->> ATIMA
		dedx.push_back(v[5]);
	}
	if(e.empty() || dedx.empty())
		bail("Read nothing from file!");

	return make_unique<TSpline3>(
		"edespline",e.data(),dedx.data(),int(e.size()));
}

void calc_ede(const string& particleName,
							const unique_ptr<TSpline3>& spl,
							double elow,
							double ehigh,
							const unique_ptr<array<double,4> >& thresh,
							const unique_ptr<array<double,4> >& saturate,
							int color = 1) {
	auto t = make_unique<TTree>(
		Form("ede_%s",particleName.c_str()),
		Form("SiStack energy losses for particle: %s",particleName.c_str()) );

	double eion,etot;
	array<double,4> de;
	array<int,4> sat;
	t->Branch("eion",&eion);
	t->Branch("etot",&etot);
	for(int i=0; i< 4; ++i){
		t->Branch(Form("de%i",i+1),de.data()+i);
		t->Branch(Form("sat%i",i+1),sat.data()+i);
	}

	double enow;
	eion = elow;
	while(eion <= ehigh) {
		fill(de.begin(),de.end(),0);
		fill(sat.begin(),sat.end(),0);

		enow = eion;
		for(int i=0; i< 4; ++i) {
			de[i] = spl->Eval(enow) * Thick[i];
			if(de[i] > enow) {
				de[i] = enow;
				break;
			}
			enow -= de[i];
		}
		for(int i=0; i< 4; ++i) {
			if(thresh.get() && de[i] < thresh->at(i)) {
				de[i] = 0;
			}
			if(saturate.get() && de[i] > saturate->at(i)) {
				de[i] = 0;
				sat[i] = 1;
			}
		}
		etot = accumulate(de.begin(),de.end(),0.);
		t->Fill();
		eion += 0.1;
	}
	{ // to make saturation show up
		eion=etot=0;
		fill(de.begin(),de.end(),0);
		fill(sat.begin(),sat.end(),1);
		t->Fill();
	}

	t->SetLineWidth(3);
	t->SetMarkerStyle(6);
	t->SetLineColor(color);
	t->SetMarkerColor(color);
	t->Write();
}

struct args_t {
	string fout;
	unique_ptr<array<double,4> > thresh;
	unique_ptr<array<double,4> > saturate;
	array<double,2> erange;
	
	args_t(): fout("ede-sistack.root"),
						thresh(nullptr),
						saturate(nullptr),
						erange({1.,200.})
		{}
};

args_t handle_args(int argc, char** argv, bool& help) {
	help = false;
	
	args_t out;
	if(argc == 1) return out;

	for(int i=1; i< argc; ++i) {
		if(0);
		else if (string(argv[i]) == "-o" && i+1 < argc)
			out.fout = argv[++i];
		else if (string(argv[i]) == "-thresh" && i+4 < argc) {
			out.thresh.reset(new array<double,4>);
			for(int j=0;j<4;++j) out.thresh->at(j) = atof(argv[++i]);
		}
		else if (string(argv[i]) == "-sat" && i+4 < argc) {
			out.saturate.reset(new array<double,4>);
			for(int j=0;j<4;++j) out.saturate->at(j) = atof(argv[++i]);
		}
		else if (string(argv[i]) == "-erange" && i+2 < argc) {
			for(int j=0;j<2;++j) out.erange.at(j) = atof(argv[++i]);
		}
		else if (string(argv[i]) == "-h" || string(argv[i]) == "--help") {
			help = true;
		}
		else {
			cout << "Skipping invalid argument " << argv[i] << endl;
			help = true;
		}
	}
	return out;
}

int give_help()
{
	cout << "usage: ede [-o <output file name>] [-erange <low> <high>] [-thresh <e1> <e2> <e3> <e4>] [-sat <e1> <e2> <e3> <e4>] [-h --help]\n";
	return 1;
}

} // namespace

int main(int argc, char** argv)
{
	bool help;
	auto args = handle_args(argc,argv,help);
	if(help) return give_help();
	
	auto f = make_unique<TFile>(args.fout.c_str(),"recreate");

	int i=0;
	array<int,5> colors = {
		1,2,3,4,6
	};
	for(auto p : particles) {
		auto spl = read_file(p);
		calc_ede(
			p,spl,
			args.erange[0],
			args.erange[1],
			args.thresh,
			args.saturate,
			colors.at(i++)
			);
	}
}
