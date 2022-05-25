/*
	Notes about sizes.
	From original file 340 Mb, ~24 s to read
	
	a. For each "storage", create 2x branches:
	  data (UShort_t)
    channel (UShort_t)
		--> 199 M

	b. For each "storage", create 1x branches:
	  data&storage (UInt_t)
		--> 133 M

  c. One branch for data, channel, storage ID, each UShort_t
    map<int, string for storage ID>
		--> 108 M

	d. One branch for data, channel, storage ID, bitpacked into single UInt_t
    map<int, string for storage ID>
		--> 67 M

*/
#ifndef WRITE_SPARSE_hh_guard
#define WRITE_SPARSE_hh_guard
#include <map>
#include <string>
#include <array>
#include <memory>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>

class MVMEExperiment;
class Storage;

struct MeasurementType {
	std::vector<UInt_t>* fPackedData;
};


class WriteSparse {
public:
  WriteSparse(MVMEExperiment* experiment, const std::string& fileOutName);
  ~WriteSparse();

  void Clear();
  void Fill()   {	for(auto& t : fTrees) t->Fill(); }
  void Write();
  void CdFile() { fFile.cd(); }
  void AddData(UInt_t event, const Storage*, UInt_t channel, double paramValue);
  
private:
	// std::map<const Storage*, std::shared_ptr<MeasurementType> > fStorageMap;
	std::vector<TObjArray> fStorageId;
	std::vector< std::map<std::string, UShort_t> > fStorageIdByName;
  TFile fFile;
	std::vector<TTree*> fTrees;
	std::vector<std::shared_ptr<MeasurementType> > fPackedData;
};


#endif
