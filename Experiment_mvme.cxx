
#include "Experiment_mvme.h"

namespace event0_modules
{

Module_mdpp32_1_scp_Sector_Si1And2::Module_mdpp32_1_scp_Sector_Si1And2()
    : MVMEModule("mdpp32_1_scp_Sector_Si1And2", "Module mdpp32_1_scp_Sector_Si1And2")

{
    RegisterDataStorage(mdpp32_1_scp_Sector_Si1And2_amplitude, 32, 16, "mdpp32_1_scp_Sector_Si1And2_amplitude",
        {  });
    RegisterDataStorage(mdpp32_1_scp_Sector_Si1And2_channel_time, 32, 16, "mdpp32_1_scp_Sector_Si1And2_channel_time",
        {  });
    RegisterDataStorage(mdpp32_1_scp_Sector_Si1And2_trigger_time, 2, 16, "mdpp32_1_scp_Sector_Si1And2_trigger_time",
        {  });
    RegisterDataStorage(mdpp32_1_scp_Sector_Si1And2_module_timestamp, 1, 30, "mdpp32_1_scp_Sector_Si1And2_module_timestamp",
        {  });
    RegisterDataStorage(mdpp32_1_scp_Sector_Si1And2_extended_timestamp, 1, 16, "mdpp32_1_scp_Sector_Si1And2_extended_timestamp",
        {  });
}

Module_mdpp32_1_scp_Sector_Si1And2::~Module_mdpp32_1_scp_Sector_Si1And2()
{}


Module_mdpp32_2_scp_Sector_Si3And4::Module_mdpp32_2_scp_Sector_Si3And4()
    : MVMEModule("mdpp32_2_scp_Sector_Si3And4", "Module mdpp32_2_scp_Sector_Si3And4")

{
    RegisterDataStorage(mdpp32_2_scp_Sector_Si3And4_amplitude, 32, 16, "mdpp32_2_scp_Sector_Si3And4_amplitude",
        {  });
    RegisterDataStorage(mdpp32_2_scp_Sector_Si3And4_channel_time, 32, 16, "mdpp32_2_scp_Sector_Si3And4_channel_time",
        {  });
    RegisterDataStorage(mdpp32_2_scp_Sector_Si3And4_trigger_time, 2, 16, "mdpp32_2_scp_Sector_Si3And4_trigger_time",
        {  });
    RegisterDataStorage(mdpp32_2_scp_Sector_Si3And4_module_timestamp, 1, 30, "mdpp32_2_scp_Sector_Si3And4_module_timestamp",
        {  });
    RegisterDataStorage(mdpp32_2_scp_Sector_Si3And4_extended_timestamp, 1, 16, "mdpp32_2_scp_Sector_Si3And4_extended_timestamp",
        {  });
}

Module_mdpp32_2_scp_Sector_Si3And4::~Module_mdpp32_2_scp_Sector_Si3And4()
{}


Module_mdpp32_3_scp_Ring_Si1::Module_mdpp32_3_scp_Ring_Si1()
    : MVMEModule("mdpp32_3_scp_Ring_Si1", "Module mdpp32_3_scp_Ring_Si1")

{
    RegisterDataStorage(mdpp32_3_scp_Ring_Si1_amplitude, 32, 16, "mdpp32_3_scp_Ring_Si1_amplitude",
        {  });
    RegisterDataStorage(mdpp32_3_scp_Ring_Si1_channel_time, 32, 16, "mdpp32_3_scp_Ring_Si1_channel_time",
        {  });
    RegisterDataStorage(mdpp32_3_scp_Ring_Si1_trigger_time, 2, 16, "mdpp32_3_scp_Ring_Si1_trigger_time",
        {  });
    RegisterDataStorage(mdpp32_3_scp_Ring_Si1_module_timestamp, 1, 30, "mdpp32_3_scp_Ring_Si1_module_timestamp",
        {  });
    RegisterDataStorage(mdpp32_3_scp_Ring_Si1_extended_timestamp, 1, 16, "mdpp32_3_scp_Ring_Si1_extended_timestamp",
        {  });
}

Module_mdpp32_3_scp_Ring_Si1::~Module_mdpp32_3_scp_Ring_Si1()
{}


Module_mdpp32_4_scp_Ring_Si3::Module_mdpp32_4_scp_Ring_Si3()
    : MVMEModule("mdpp32_4_scp_Ring_Si3", "Module mdpp32_4_scp_Ring_Si3")

{
    RegisterDataStorage(mdpp32_4_scp_Ring_Si3_amplitude, 32, 16, "mdpp32_4_scp_Ring_Si3_amplitude",
        {  });
    RegisterDataStorage(mdpp32_4_scp_Ring_Si3_channel_time, 32, 16, "mdpp32_4_scp_Ring_Si3_channel_time",
        {  });
    RegisterDataStorage(mdpp32_4_scp_Ring_Si3_trigger_time, 2, 16, "mdpp32_4_scp_Ring_Si3_trigger_time",
        {  });
    RegisterDataStorage(mdpp32_4_scp_Ring_Si3_module_timestamp, 1, 30, "mdpp32_4_scp_Ring_Si3_module_timestamp",
        {  });
    RegisterDataStorage(mdpp32_4_scp_Ring_Si3_extended_timestamp, 1, 16, "mdpp32_4_scp_Ring_Si3_extended_timestamp",
        {  });
}

Module_mdpp32_4_scp_Ring_Si3::~Module_mdpp32_4_scp_Ring_Si3()
{}


Module__mdpp16_scp_PPAC_And_SB::Module__mdpp16_scp_PPAC_And_SB()
    : MVMEModule("_mdpp16_scp_PPAC_And_SB", "Module _mdpp16_scp_PPAC_And_SB")

{
    RegisterDataStorage(_mdpp16_scp_PPAC_And_SB_amplitude, 16, 16, "_mdpp16_scp_PPAC_And_SB_amplitude",
        {  });
    RegisterDataStorage(_mdpp16_scp_PPAC_And_SB_channel_time, 16, 16, "_mdpp16_scp_PPAC_And_SB_channel_time",
        {  });
    RegisterDataStorage(_mdpp16_scp_PPAC_And_SB_trigger_time, 2, 16, "_mdpp16_scp_PPAC_And_SB_trigger_time",
        {  });
    RegisterDataStorage(_mdpp16_scp_PPAC_And_SB_module_timestamp, 1, 30, "_mdpp16_scp_PPAC_And_SB_module_timestamp",
        {  });
    RegisterDataStorage(_mdpp16_scp_PPAC_And_SB_extended_timestamp, 1, 16, "_mdpp16_scp_PPAC_And_SB_extended_timestamp",
        {  });
}

Module__mdpp16_scp_PPAC_And_SB::~Module__mdpp16_scp_PPAC_And_SB()
{}


Module__mdpp16_qdc_Phoswich::Module__mdpp16_qdc_Phoswich()
    : MVMEModule("_mdpp16_qdc_Phoswich", "Module _mdpp16_qdc_Phoswich")

{
    RegisterDataStorage(_mdpp16_qdc_Phoswich_channel_time, 16, 16, "_mdpp16_qdc_Phoswich_channel_time",
        {  });
    RegisterDataStorage(_mdpp16_qdc_Phoswich_integration_long, 16, 16, "_mdpp16_qdc_Phoswich_integration_long",
        {  });
    RegisterDataStorage(_mdpp16_qdc_Phoswich_integration_short, 16, 16, "_mdpp16_qdc_Phoswich_integration_short",
        {  });
    RegisterDataStorage(_mdpp16_qdc_Phoswich_trigger_time, 2, 16, "_mdpp16_qdc_Phoswich_trigger_time",
        {  });
    RegisterDataStorage(_mdpp16_qdc_Phoswich_module_timestamp, 1, 30, "_mdpp16_qdc_Phoswich_module_timestamp",
        {  });
    RegisterDataStorage(_mdpp16_qdc_Phoswich_extended_timestamp, 1, 16, "_mdpp16_qdc_Phoswich_extended_timestamp",
        {  });
}

Module__mdpp16_qdc_Phoswich::~Module__mdpp16_qdc_Phoswich()
{}

} // end namespace event0_modules

Event_event0::Event_event0()
    : MVMEEvent("event0", "Storage for event 'event0'")
{
    AddModule(&mdpp32_1_scp_Sector_Si1And2);
    AddModule(&mdpp32_2_scp_Sector_Si3And4);
    AddModule(&mdpp32_3_scp_Ring_Si1);
    AddModule(&mdpp32_4_scp_Ring_Si3);
    AddModule(&_mdpp16_scp_PPAC_And_SB);
    AddModule(&_mdpp16_qdc_Phoswich);
}

Event_event0::~Event_event0()
{}


Experiment::Experiment()
    : MVMEExperiment("Experiment", "")
{
    AddEvent(&event0);
}

Experiment::~Experiment()
{}

