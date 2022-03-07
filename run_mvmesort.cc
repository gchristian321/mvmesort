// Modified from code in
// ~/Packages/mvme-1.4.8-3-Linux-x64/share/mvme_root_client/mvme_root_client.cc
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

// // mvme
// #include <mvme/Mustache/mustache.hpp> // mustache template engine
// #include <mvme/event_server/common/event_server_lib.h> // mvme event_server protocol parsing and socket handling
#include "mvme_root_event_objects.h" // base classes for generated experiment ROOT objects

#include "MVMESort.cc"

using namespace std;
using MVMERunInfo = std::map<std::string, std::string>;

std::unique_ptr<MVMEExperiment> compile_and_load_experiment_library(
    const std::string &expName)
{
    const std::string expStructName = expName;
    const std::string libName = "lib" + expName + "_mvme.so";

    // Run the ROOT pre-make macro
    {
        cout << ">>> Executing ROOT pre-make macro file mvme_root_premake_hook.C ..." << endl;
        auto res = gROOT->ProcessLineSync(".x mvme_root_premake_hook.C");
        cout << endl << "<<< End of output from ROOT premake hook" << endl << endl;

        if (res)
            return {};
    }

    // Run make
    {
        cout << ">>> Running make ..." << endl;
        int res = gSystem->Exec("make");
        cout << endl << "<<< End of make output" << endl;

        if (res)
            return {};
    }

    // Load experiment library
    {
        cout << ">>> Loading experiment library " << libName << endl;
        int res = gSystem->Load(libName.c_str());

        if (res != 0 && res != 1)
        {
            cout << "Error loading experiment library " << libName << endl;
            return {};
        }
    }

    // Create and return an instance of the generated experiment class cast to
    // the base class type. (Only the base class type is known to this code
    // here, not the concrete subclass from the generated code).
    std::string cmd = "new " + expStructName + "();";

    return std::unique_ptr<MVMEExperiment>(
        reinterpret_cast<MVMEExperiment *>(
            gROOT->ProcessLineSync(cmd.c_str())));
}


//
// replay from ROOT file
//
int run_mvmesort(const std::string& inputFilename,
		 const std::string& outputFilename,
		 const std::string& channelMapFile,
		 bool saveRaw = true,
		 bool channelMapWarn = false)
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

  auto exp = compile_and_load_experiment_library(expName);

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
