// Modified from code in
// ~/Packages/mvme-1.4.8-3-Linux-x64/share/mvme_root_client/mvme_root_client.cc
#include <iostream>
#include <fstream>
#include <string>

#include <dlfcn.h> // dlopen, dlsym, dlclose
#include <getopt.h>
#include <signal.h>
#include <sys/stat.h>

// ROOT
#include <TFile.h>
#include <TH1D.h>
#include <TRandom.h> // gRandom
#include <TROOT.h>   // gROOT
#include <TSystem.h> // gSystem

// mvme
#include "mvme_root_event_objects.h" // base classes for generated experiment ROOT objects
#include "DetectorSort.h"
#include "PhysicsSort.h"

using namespace std;
using MVMERunInfo = std::map<std::string, std::string>;


//
// replay from ROOT file
//
int run_mvmesort(const std::string& inputFilename,
								 const std::string& outputFilename,
								 const std::string& channelMapFile,
								 double matchWindow,
								 bool saveRaw,
								 bool channelMapWarn)
{
  // Replay from ROOT file.
  TFile f(inputFilename.c_str(), "read");

  cout << ">>> Reading MVMERunInfo from " << inputFilename << endl;

  auto runInfoPtr = reinterpret_cast<MVMERunInfo *>(f.Get("MVMERunInfo"));

  if (!runInfoPtr)
	{
		cout << "Error: input file " << inputFilename << " does not contain an MVMERunInfo object"
				 << endl;
		return 1;
	}

  auto runInfo = *runInfoPtr;


  std::string expName = runInfo["ExperimentName"];
  std::string expStructName = expName;
  std::string runId = runInfo["RunID"];
  std::string libName = "lib" + expName + "_mvme.so";

  cout << ">>> Run info: ExperimentName = " << expName << ", RunID=" << runId << endl;

	std::string cmd = "new " + expStructName + "();";
	auto exp = std::unique_ptr<MVMEExperiment>(
		reinterpret_cast<MVMEExperiment *>(
			gROOT->ProcessLineSync(cmd.c_str())));


  if (!exp)
	{
		cout << "Error creating an instance of the experiment class '"
				 << expStructName << "'" << endl;
		return 1;
	}

  
  
  
  // Setup tree branch addresses to point to the experiment subobjects.
  auto eventTrees = exp->InitTrees(&f);

  if (eventTrees.size() != exp->GetNumberOfEvents())
	{
		cout << "Error: could not read experiment eventTrees from input file "
				 << inputFilename << endl;
		return 1;
	}

  // Replay tree data
  if(eventTrees.size() != 1) {
    cerr << "ERROR: eventTrees.size() == " << eventTrees.size() <<
      ", don't know how to handle this!\n";
    return 1;
  }
  
  DetectorSort detSort(outputFilename.c_str(), channelMapFile, saveRaw);
  PhysicsSort physSort(detSort);
	if(matchWindow > 0) {
		physSort.SetRingSectorMatchWindow(matchWindow);
	}
	
  detSort.SetWarnChannelMap(channelMapWarn);

  for(size_t eventIndex = 0; eventIndex < eventTrees.size(); eventIndex++) {
    auto tree = eventTrees[eventIndex];
    auto event = exp->GetEvent(eventIndex);
    //      auto analyzeFunc = analysis.eventFunctions[eventIndex];

    cout << "Replaying data from tree '" << tree->GetName() << "'..." << std::flush;

    const auto entryCount = tree->GetEntries();

    cout << "\nPercent complete: 0... ";
    flush(cout);
    for(int64_t entryIndex = 0; entryIndex < entryCount; entryIndex++) {
      detSort.Clear();
      physSort.Clear();

      // Fills the event and its submodules with data read from the ROOT tree.
      tree->GetEntry(entryIndex);

      // Read values from the generated array members of the events
      // module classes and fill the raw histograms.
      for(auto& m : event->GetModules()) {
				for(const auto& storage : m->GetDataStorages()) {
					for(size_t moduleCh = 0; moduleCh < storage.size; moduleCh++) {
						double paramValue = storage.ptr[moduleCh];
						if (!std::isnan(paramValue)) {
							detSort.AddData(moduleCh, m->GetName(), storage.name, paramValue);
						}
					} // module Ch
				} // storages
      } // modules
      detSort.Fill();

      physSort.Calculate();
      physSort.Fill();
      
      // Progress counter //
      static int64_t last = 0;
      if(int64_t(10.*entryIndex/entryCount) > last) {
				last = int64_t(10.*entryIndex/entryCount);
				cout << 10*last << "... ";
				flush(cout);
      }
    } // entryIndex
    cout << "100\nDone! Processed " << entryCount << " events\n";
  } // treeIndex

  detSort.Write();
  physSort.Write();

  physSort.PrintEnergyTimeMismatches();
  
  return 0;
}


int main(int argc, char** argv)
{
	auto usage =
		[&]() {
			cerr << "usage: mvmesort <input file> <output file> <channel map> [--no-save-raw] [--no-channel-map-warn]" << endl;
			cerr << "OPTIONAL ARGUMENTS\n";
			cerr << "  --match-window=<window> -->> set Si ring + sector energy matching window in MeV (default: DBL_MAX [very large])\n";
			cerr << "  --no--save-raw -->> disable saving of \"raw\" (detector-level) data to output file (default IS to save)\n";
			cerr << "  --no-channel-map-warn  -->> turn off warnings about problems with the channel map file (default IS to warn)\n";
			return 1;
		};

	if(argc < 4) return usage();

	auto check_begin =
		[](const char* arg, const string& key) {
			return string(arg).substr(0, key.size()) == key;
		};
	auto extract_end =
		[](const char* arg, const string& key) {
			return string(string(arg).substr(key.size()));
		};
	
	bool saveRaw=true, warnChannelMap=true;
	double matchWindow = -1;
	for(int i=4; i< argc; ++i) {
		if(string(argv[i]) == "--no-save-raw") saveRaw = false;
		else if(string(argv[i]) == "--no-channel-map-warn") warnChannelMap = false;
		else if(check_begin(argv[i],"--si-match-window=")) {
			matchWindow = atof(extract_end(argv[i], "--si-match-window=").c_str());
			cout << "Set Si ring+sector energy match window to " << matchWindow << "MeV\n";
		}
		else return usage();
	}

	int err = run_mvmesort(
		argv[1], argv[2], argv[3],
		matchWindow, saveRaw, warnChannelMap);
	if(!err) { // set output permissions to have group write
		gSystem->Exec(Form("chmod g+w %s",argv[2]));
	}
	return err;
}
