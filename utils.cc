namespace mvmeutils {

const int kRings = 23;
const int kSectors = 16;


TH1* StripHist(
	TTree* t, int detector, const string& ringOrSector,
	const char* gate, const char* name, const char* title,
	int bins, double low, double high)
{
	TString rs = TString(ringOrSector.c_str()); rs.ToLower();
	const int nstrips = rs == "ring" ? kRings : kSectors;
	TH2D* hst = new TH2D(name,title,bins,low,high,nstrips,1,nstrips+1);

	for(int i=1; i<= nstrips; ++i) {
		string rs1 = rs == "ring" ? "Ring" : "Sector";
		auto nevt = t->Draw(Form("Si%i_%s%i_E",detector,rs1.c_str(),i),gate,"goff");
		for(int j=0; j< nevt; ++j) {
			hst->Fill(t->GetV1()[j], i);
		}
	}
	hst->Draw("colz");
	return hst;
}

}
