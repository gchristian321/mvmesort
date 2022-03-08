

#include "Experiment_mvme.h"
#include <iostream>

extern "C"
{

// Called once after loading this library. Additional command line parameters
// given to mvme_root_client (appearing after '--') are passed in the args
// variable.
bool init_analysis(const std::vector<std::string> &args)
{
    std::cerr << __FUNCTION__ << "(), args: ";
    for (const auto &arg: args)
        std::cerr << arg << " ";
    std::cerr << std::endl;

    return true;
}

// Called before the mvme_root_client program exits.
bool shutdown_analysis()
{
    std::cerr << __FUNCTION__ << std::endl;
    return true;
}

// Called each time a new run stats.
bool begin_run(const std::string &inputSource, const std::string &runId, bool isReplay)
{
    std::cerr << __FUNCTION__ << std::endl;
    return true;
}

// Called each time a run ends.
bool end_run()
{
    std::cerr << __FUNCTION__ << std::endl;
    return true;
}

bool analyze_event0(Event_event0 *event)
{
    // ===== Custom analysis code for Event_event0 goes here =====

    return true;
}

} // end extern "C"

