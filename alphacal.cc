namespace alphacal {

double peakFunction(double*x, double*p)
{
	double H = p[0];
	double c = p[1];
	double FWHM = p[2];
	
	double R = p[3];
	double beta = p[4] * c;

	double side = p[5];
	
	double GausPart =
		H*(1-R) * TMath::Gaus(*x,c,(FWHM/2.355));

	double dx = side <= 0 ? *x-c : c-*x;
	double SkewPart;
	if(R == 0) { SkewPart = 0; }
	//else if(dx > 0) { SkewPart = 0; }
	else {
		SkewPart = H*R * TMath::Exp( dx/beta ) * \
			TMath::Erfc( dx/(sqrt(2)*(FWHM/2.355)) + (FWHM/2.355)/(sqrt(2)*beta) );
	}

	return GausPart + SkewPart;
}

double peakFunction3(double*x, double*p)
{
	return peakFunction(x,p) + peakFunction(x,p+6) + peakFunction(x,p+12);
}

TF1* create_fpeak(int npeaks, string name = "")
{
	if(name == "") name = Form("fpeak%i",npeaks);
	auto f = dynamic_cast<TF1*>(gROOT->GetListOfFunctions()->FindObject(name.c_str()));
	if(f) return f;

	TF1* fpeak = 0;
	if(npeaks == 3) fpeak = new TF1(name.c_str(), peakFunction3, 0, 65535, 18);
	else if(npeaks == 1) fpeak = new TF1(name.c_str(), peakFunction, 0, 65535, 6);
	else return fpeak;
	
	fpeak->SetNpx(10000);
	for(int i=0; i< npeaks; ++i) {
		fpeak->FixParameter(5+6*i, -1);
	}
	
	return fpeak;
}

array<double, 3> find_peaks(
	TTree* t, int det, int channel, const string& ringOrSector, TH1** hout=0) {
	array<double,3> out = {0,0,0};
	string param = Form("Si%i_%s%i_E",det,ringOrSector.c_str(),channel);
	TH1D *halpha = new TH1D("halpha","", 512, 0, 65535);
	t->Project("halpha",param.c_str());
	static TSpectrum s;
	s.Search(halpha,2,"nodraw");
	auto npk = s.GetNPeaks();
	if(npk >= 3) {
		vector<int> srt(npk);
		TMath::Sort(npk,s.GetPositionY(),&srt[0]);

		for(int i=0; i< 3; ++i) {
			out[i] = s.GetPositionX()[srt.at(i)];
		}
	}
	if(hout) {
		vector<double> bc(halpha->GetNbinsX());
		for(size_t i=0; i< bc.size(); ++i) {
			bc.at(i) = halpha->GetBinContent(i+1);
		}
		int nb = halpha->GetNbinsX();
		double lo = halpha->GetBinLowEdge(1);
		double hi = halpha->GetBinLowEdge(nb+1);
		halpha->Delete();
		*hout = new TH1D("halpha","",nb,lo,hi);
		for(size_t i=0; i< bc.size(); ++i) {
			(*hout)->SetBinContent(i+1,bc.at(i));
		}
	}
	return out;
}

struct FitResults {
	array<double,3> H;
	array<double,3> C;
	array<double,3> FWHM;
	array<double,3> R;
	array<double,3> BETA;
	bool success;
};

FitResults fit_peaks(
	TTree* t, int det, int channel, const string& ringOrSector, bool batch = false) {
	TH1* halpha = 0;
	auto pk = find_peaks(t,det,channel,ringOrSector,&halpha);
	if(pk[0] == 0 && pk[1] == 0 && pk[2] == 0) {
		FitResults nogood;
		nogood.success = false;
		return nogood;
	}

	TF1* fpk = create_fpeak(3);

	double R = 0.01;
	double beta = 0.001;
	for(int i=0; i< 3; ++i) {
		double c = pk[i];
		double H = halpha->GetBinContent(halpha->FindBin(c));
		double FWHM = (40/5e3) * c;

		fpk->SetParameter(0 + i*6, H);
		fpk->SetParameter(1 + i*6, c);
		fpk->SetParameter(2 + i*6, FWHM);
		fpk->FixParameter(3 + i*6, 0);
		fpk->FixParameter(4 + i*6, 0);
	}

	auto fr = halpha->Fit(fpk,"qnsw"); // gaus only
	for(int i=0; i< 3; ++i) {
		fpk->FixParameter(0 + i*6, fr->GetParams()[0+i*6]);
		fpk->FixParameter(1 + i*6, fr->GetParams()[1+i*6]);
		fpk->FixParameter(2 + i*6, fr->GetParams()[2+i*6]);
		fpk->ReleaseParameter(3 + i*6);
		fpk->ReleaseParameter(4 + i*6);
		fpk->SetParameter(3 + i*6, R);
		fpk->SetParameter(4 + i*6, beta);
	}
	fr = halpha->Fit(fpk,"qns"); // fixed gaus, fits skew
	for(int i=0; i< 18; ++i) {
		if(i == 5 || i == 5+6*1 || i == 5+6*2) {
			continue;
		}
		else {
			fpk->ReleaseParameter(i);
		}
	}
	array<double,3> centers = {
		fr->GetParams()[1],fr->GetParams()[1+6],fr->GetParams()[1+12]
	};
	auto ilo = min_element(centers.begin(),centers.end())-centers.begin();
	auto ihi = max_element(centers.begin(),centers.end())-centers.begin();
	double lo = fr->GetParams()[1 + ilo*6] - fr->GetParams()[2 + ilo*6]*3;
	double hi = fr->GetParams()[1 + ihi*6] + fr->GetParams()[2 + ihi*6]*3;
	string opt = batch ? "qns" : "s";
	fr = halpha->Fit(fpk,opt.c_str(),"",lo,hi); // fits everything
	if(!batch) {
		halpha->GetXaxis()->SetRangeUser(lo,hi);
		halpha->Draw("E");
	} else {
		halpha->Delete();
	}

	double cfit[3];
	int ifit[3];
	for(int i=0; i< 3; ++i) {
		cfit[i] = fr->GetParams()[1+i*6];
	}
	TMath::Sort(3,cfit,ifit,false);
		
	FitResults out;
	out.success = (fr->Status() == 0);
	for(int i=0; i< 3; ++i) {
		out.H[i]    = fr->GetParams()[0+ifit[i]*6];
		out.C[i]    = fr->GetParams()[1+ifit[i]*6];
		out.FWHM[i] = fr->GetParams()[2+ifit[i]*6];
		out.R[i]    = fr->GetParams()[3+ifit[i]*6];
		out.BETA[i] = fr->GetParams()[4+ifit[i]*6];
	}
	return out;
};

FitResults fit_peaks1(
	TTree* t, int det, int channel, const string& ringOrSector, bool batch = false) {
	TH1* halpha = 0;
	auto pk = find_peaks(t,det,channel,ringOrSector,&halpha);
	if(pk[0] == 0 && pk[1] == 0 && pk[2] == 0) {
		FitResults nogood;
		nogood.success = false;
		return nogood;
	}

	array<TFitResultPtr,3> fr;
	array<TF1*,3> fpk;
	for(int i=0; i< 3; ++i) {
		fpk[i] = create_fpeak(1,Form("fpeak_%i",i));
	}

	double R = 0.01;
	double beta = 0.001;
	for(int i=0; i< 3; ++i) {
		double c = pk[i];
		double H = halpha->GetBinContent(halpha->FindBin(c));
		double FWHM = (40/5e3) * c;

		fpk[i]->SetParameter(0, H);
		fpk[i]->SetParameter(1, c);
		fpk[i]->SetParameter(2, FWHM);
		fpk[i]->FixParameter(3, 0);
		fpk[i]->FixParameter(4, 0);

		fr[i] = halpha->Fit(fpk[i],"qnsw"); // gaus only
	}

	array<double,3> lo,hi;
	for(int i=0; i< 3; ++i) {
		fpk[i]->FixParameter(0, fr[i]->GetParams()[0]);
		fpk[i]->FixParameter(1, fr[i]->GetParams()[1]);
		fpk[i]->FixParameter(2, fr[i]->GetParams()[2]);
		fpk[i]->ReleaseParameter(3);
		fpk[i]->ReleaseParameter(4);
		fpk[i]->SetParameter(3, R);
		fpk[i]->SetParameter(4, beta);

		lo[i] = fr[i]->GetParams()[1] - fr[i]->GetParams()[2]*2.5;
		hi[i] = fr[i]->GetParams()[1] + fr[i]->GetParams()[2]*2.5;
		
		fr[i] = halpha->Fit(fpk[i],"qns","",lo[i],hi[i]); // fixed gaus, fits skew
	}

	array<double,3> centers;
	for(int i=0;i<3;++i) {
		for(int j=0; j< 5; ++j) {
			fpk[i]->ReleaseParameter(j);
		}
		string opt = batch ? "qns" : "s";
		fr[i] = halpha->Fit(fpk[i],opt.c_str(),"",lo[i],hi[i]);
		centers[i] = fr[i]->GetParams()[1];
		fpk[i]->SetRange(lo[i],hi[i]);
	}


	int icsort[3];
	TMath::Sort(3,&centers[0],icsort,false);
		
	FitResults out;
	out.success = true;
	for(int i=0; i< 3; ++i) {
		out.H[i]    = fr[icsort[i]]->GetParams()[0];
		out.C[i]    = fr[icsort[i]]->GetParams()[1];
		out.FWHM[i] = fr[icsort[i]]->GetParams()[2];
		out.R[i]    = fr[icsort[i]]->GetParams()[3];
		out.BETA[i] = fr[icsort[i]]->GetParams()[4];
		if(fr[icsort[i]]->Status() != 0) out.success = false;
	}

	if(!batch) {
		halpha->GetXaxis()->SetRangeUser(lo[icsort[0]],hi[icsort[2]]);
		halpha->SetMarkerStyle(20);
		halpha->SetLineWidth(2);
		halpha->SetMarkerColor(kBlack);
		halpha->SetLineColor(kBlack);
		halpha->Draw("E");
		for(int i=0; i< 3; ++i) fpk[i]->Draw("SAME");
	} else {
		halpha->Delete();
	}

	return out;
};


pair<double,double> fit_calibration(const array<double, 3>& channels) {
	array<double, 3> energies = { 5.157, 5.486, 5.805 };
	TGraph gr(3,&channels[0],&energies[0]);
	auto fr = gr.Fit("pol1","qns");
	return make_pair(fr->GetParams()[0], fr->GetParams()[1]);
}

map<int, pair<double,double> > fit_all_peaks(
	TTree* t, int det, string ringOrSector, bool fit) {
	map<int, pair<double, double> > m;

	int iMax = ringOrSector == "Ring" ? 23 : 16;
	for(int i=1; i<= iMax; ++i) {
		if(fit) {auto pk = fit_peaks1(t,det,i,ringOrSector,true);
			if(pk.success) {
				m.emplace(i, fit_calibration(pk.C));
			}
			else {
				m.emplace(i, make_pair(0.,0.));
			}
		}
		else {
			auto pk = find_peaks(t,det,i,ringOrSector);
			m.emplace(i, fit_calibration(pk));
		}
	}
	return m;
}

void save_all_peaks(
	ostream& out, TTree* t,
	int det, string ringOrSector, bool fit) {
	auto pk = fit_all_peaks(t,det,ringOrSector,fit);
	for(const auto& p : pk) {
		out << p.first << "," <<
			p.second.first << "," <<
			p.second.second << endl;
	}
}

}
