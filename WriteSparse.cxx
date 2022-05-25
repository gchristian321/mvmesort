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
		fStorageIdByName.emplace_back(map<std::string,UShort_t>());

		UInt_t sid = 0;
		for(auto& storage : event->GetDataSourceStorages()) {
			fStorageIdByName.back().emplace(storage.name, sid++);
		}
		if(sid > 0x7ff) {
			throw runtime_error(
				Form("Too many storages [%i], max is 0x7ff!", sid));
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
	if(storage->name.find("timestamp") < storage->name.size()) return;
	// ^^ timestamps are larger than 32 bits and need special treatment
	//    we will not write them, for now
	
	UInt_t val = 0;
	auto it = fStorageIdByName.at(event).find(storage->name);
	if(it == fStorageIdByName.at(event).end()) {
		throw invalid_argument(Form("Don't recognize storage: %s", storage->name.c_str()));
	}
	const UInt_t storageID = it->second;
	const UInt_t param = UInt_t(paramValue);
	val = (channel << 16) | val;
	val = (storageID << 21) | val;
	val = param | val;
	
	fPackedData.at(event)->fPackedData->push_back(val);
}

void WriteSparse::Write()
{
	for(size_t i=0; i< fTrees.size(); ++i) {
		fTrees.at(i)->Write();
		TObjArray storage(fStorageIdByName.at(i).size());
		for(const auto& it : fStorageIdByName.at(i)) {
			storage.AddAt(new TObjString(it.first.c_str()), it.second);
			//cout << "Adding " << it.second << " : " << it.first << endl;
		}
		storage.Write(
		 	Form("%s_StorageID", fTrees.at(i)->GetName()),
			TObject::kSingleKey );
	}
}
