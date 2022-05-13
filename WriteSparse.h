#ifndef WRITE_SPARSE_hh_guard
#define WRITE_SPARSE_hh_guard
#include <map>
#include <string>
#include <array>
#include <memory>

// ROOT
#include <TFile.h>
#include <TTree.h>


class MVMEExperiment;
class Storage;

struct MeasurementType {
	std::vector<UInt_t>* fData;
	std::vector<UInt_t>* fChannel;
};


class WriteSparse {
public:
  WriteSparse(MVMEExperiment* experiment, const std::string& fileOutName);
  ~WriteSparse();

  void Clear();
  void Fill()   {	for(auto& t : fTrees) t->Fill(); }
  void Write()  {	for(auto& t : fTrees) t->Write(); }
  void CdFile() { fFile.cd(); }
  void AddData(const Storage*, UInt_t channel, double paramValue);
  
private:
	std::map<const Storage*, std::shared_ptr<MeasurementType> > fStorageMap;
  TFile fFile;
	std::vector<TTree*> fTrees;
};


#endif