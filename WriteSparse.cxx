#include <iostream>
#include "TObjString.h"
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
		fStorageId.emplace_back(TObjArray());
		fStorageIdByName.emplace_back(map<std::string,UShort_t>());

		for(auto& storage : event->GetDataSourceStorages()) {
			fStorageIdByName.back().emplace(storage.name, fStorageId.back().GetEntries());
			fStorageId.back().Add(new TObjString(storage.name.c_str()));
		}
		if(fStorageId.back().GetEntries() > 1023) {
			throw runtime_error(
				Form("Too many storages [%i], max is 1023!", fStorageId.back().GetEntries()));
		}
		auto measurement = make_shared<MeasurementType>();
		measurement->fPackedData = nullptr;
		fTrees.back()->Branch("packed_data", &(measurement->fPackedData));
		fPackedData.emplace_back(measurement);
	}
}

WriteSparse::~WriteSparse()
{  }

void WriteSparse::Clear()
{
	for(auto& it : fPackedData) {
		it->fPackedData->clear();
	}
}

void WriteSparse::AddData(UInt_t event,
													const Storage* storage,
													UInt_t channel,
													double paramValue)
{
	if(std::isnan(paramValue)) return;
	
	UInt_t val = 0;
	auto it = fStorageIdByName.at(event).find(storage->name);
	if(it == fStorageIdByName.at(event).end()) {
		throw invalid_argument(Form("Don't recognize storage: %s", storage->name.c_str()));
	}
	const UInt_t storageID = it->second;
	const UInt_t param = UInt_t(paramValue);
	val = ((channel >> 16) & 0x3f) | val;
	val = ((storageID >> 22) & 0x3ff) | val;
	val = param | val;
	
	fPackedData.at(event)->fPackedData->push_back(val);
}

void WriteSparse::Write()
{
	for(size_t i=0; i< fTrees.size(); ++i) {
		fTrees.at(i)->Write();
		fStorageId.at(i).Write(
		 	Form("%s_StorageID", fTrees.at(i)->GetName())
			);
	}
}
