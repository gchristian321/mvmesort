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
		fStorageModuleIdByName.emplace_back(StorageMap_t());

		UInt_t sid = 0;
		for(auto& module : event->GetModules()) {
			for(auto& storage : module->GetDataStorages()) {
				auto names = make_pair<string,string>(
					module->GetName(), string(storage.name));
				fStorageModuleIdByName.back().emplace(names, sid++);

				if(sid > 0x7ff) {
					throw runtime_error(
						Form("Too many storages [%i], max is 0x7ff!", sid));
				}
			}
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
													const string& module_name,
													const Storage* storage,
													UInt_t channel,
													double paramValue)
{
	if(std::isnan(paramValue)) return;
	if(storage->name.find("timestamp") < storage->name.size()) return;
	// ^^ timestamps are larger than 32 bits and need special treatment
	//    we will not write them, for now
	
	UInt_t val = 0;
	auto it = fStorageModuleIdByName.at(event).find(
		make_pair(module_name, storage->name));
	if(it == fStorageModuleIdByName.at(event).end()) {
		throw invalid_argument(
			Form("Don't recognize module, storage: %s, %s",
					 module_name.c_str(), storage->name.c_str()));
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

		TObjArray storage(fStorageModuleIdByName.at(i).size());
		for(const auto& it : fStorageModuleIdByName.at(i)) {
			storage.AddAt(new TObjString(it.first.second.c_str()), it.second);
		}
		storage.Write(
		 	Form("%s_StorageID", fTrees.at(i)->GetName()),
			TObject::kSingleKey );
		
		TObjArray module(fStorageModuleIdByName.at(i).size());
		for(const auto& it : fStorageModuleIdByName.at(i)) {
			module.AddAt(new TObjString(it.first.first.c_str()), it.second);
		}
		module.Write(
		 	Form("%s_ModuleID", fTrees.at(i)->GetName()),
			TObject::kSingleKey );
	}
}
