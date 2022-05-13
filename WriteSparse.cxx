#include <iostream>

#include "WriteSparse.h"
#include "Experiment_mvme.h"
using namespace std;



WriteSparse::WriteSparse(MVMEExperiment* experiment, const string& outFileName):
  fFile(outFileName.c_str(), "recreate")
{
	fFile.cd();
	for(const auto& event : experiment->GetEvents()) {
		fTrees.emplace_back(
			new TTree(event->GetName(), 
								Form("Sparse tree for \"%s\"", event->GetName())
				)
			);

		for(auto& storage : event->GetDataSourceStorages()) {
			auto measurement = make_shared<MeasurementType>();
			measurement->fData = measurement->fChannel = nullptr;
			fTrees.back()->Branch(storage.name.c_str(), &(measurement->fData));
			fTrees.back()->Branch(Form("%s_channel",storage.name.c_str()),
														&(measurement->fChannel));
			fStorageMap.emplace(&storage, measurement);
		}	
	}
}

WriteSparse::~WriteSparse()
{  }

void WriteSparse::Clear()
{
  for(auto& it : fStorageMap) {
    it.second->fData->clear();
		it.second->fChannel->clear();
  }
}

void WriteSparse::AddData(const Storage* storage,
													UInt_t channel,
													double paramValue)
{
	if(std::isnan(paramValue)) return;
	auto it = fStorageMap.find(storage);
	if(it == fStorageMap.end()) {
		throw invalid_argument(
			Form("WriteSparse::AddData: Couldn't find storage at %p",storage)
			);
	}
	it->second->fChannel->push_back(channel);
	it->second->fData->push_back(UInt_t(paramValue));
}
