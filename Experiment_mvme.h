#ifndef __MVME_EVENT_EXPORT_GUARD_Experiment__
#define __MVME_EVENT_EXPORT_GUARD_Experiment__

// Contains declarations of the MVMEExperiment, MVMEEvent and MVMEModule base
// classes
#include "mvme_root_event_objects.h"

namespace event0_modules
{
    struct Module_mdpp32_1_scp_Sector_Si1And2: public MVMEModule
    {
        Module_mdpp32_1_scp_Sector_Si1And2();
        virtual ~Module_mdpp32_1_scp_Sector_Si1And2();



        double mdpp32_1_scp_Sector_Si1And2_amplitude[32];
        double mdpp32_1_scp_Sector_Si1And2_channel_time[32];
        double mdpp32_1_scp_Sector_Si1And2_trigger_time[2];
        double mdpp32_1_scp_Sector_Si1And2_module_timestamp[1];
        double mdpp32_1_scp_Sector_Si1And2_extended_timestamp[1];

        ClassDef(Module_mdpp32_1_scp_Sector_Si1And2, 1);
    };

    struct Module_mdpp32_2_scp_Sector_Si3And4: public MVMEModule
    {
        Module_mdpp32_2_scp_Sector_Si3And4();
        virtual ~Module_mdpp32_2_scp_Sector_Si3And4();



        double mdpp32_2_scp_Sector_Si3And4_amplitude[32];
        double mdpp32_2_scp_Sector_Si3And4_channel_time[32];
        double mdpp32_2_scp_Sector_Si3And4_trigger_time[2];
        double mdpp32_2_scp_Sector_Si3And4_module_timestamp[1];
        double mdpp32_2_scp_Sector_Si3And4_extended_timestamp[1];

        ClassDef(Module_mdpp32_2_scp_Sector_Si3And4, 1);
    };

    struct Module_mdpp32_3_scp_Ring_Si1: public MVMEModule
    {
        Module_mdpp32_3_scp_Ring_Si1();
        virtual ~Module_mdpp32_3_scp_Ring_Si1();



        double mdpp32_3_scp_Ring_Si1_amplitude[32];
        double mdpp32_3_scp_Ring_Si1_channel_time[32];
        double mdpp32_3_scp_Ring_Si1_trigger_time[2];
        double mdpp32_3_scp_Ring_Si1_module_timestamp[1];
        double mdpp32_3_scp_Ring_Si1_extended_timestamp[1];

        ClassDef(Module_mdpp32_3_scp_Ring_Si1, 1);
    };

    struct Module_mdpp32_4_scp_Ring_Si3: public MVMEModule
    {
        Module_mdpp32_4_scp_Ring_Si3();
        virtual ~Module_mdpp32_4_scp_Ring_Si3();



        double mdpp32_4_scp_Ring_Si3_amplitude[32];
        double mdpp32_4_scp_Ring_Si3_channel_time[32];
        double mdpp32_4_scp_Ring_Si3_trigger_time[2];
        double mdpp32_4_scp_Ring_Si3_module_timestamp[1];
        double mdpp32_4_scp_Ring_Si3_extended_timestamp[1];

        ClassDef(Module_mdpp32_4_scp_Ring_Si3, 1);
    };

    struct Module__mdpp16_scp_PPAC_And_SB: public MVMEModule
    {
        Module__mdpp16_scp_PPAC_And_SB();
        virtual ~Module__mdpp16_scp_PPAC_And_SB();



        double _mdpp16_scp_PPAC_And_SB_amplitude[16];
        double _mdpp16_scp_PPAC_And_SB_channel_time[16];
        double _mdpp16_scp_PPAC_And_SB_trigger_time[2];
        double _mdpp16_scp_PPAC_And_SB_module_timestamp[1];
        double _mdpp16_scp_PPAC_And_SB_extended_timestamp[1];

        ClassDef(Module__mdpp16_scp_PPAC_And_SB, 1);
    };

    struct Module__mdpp16_qdc_Phoswich: public MVMEModule
    {
        Module__mdpp16_qdc_Phoswich();
        virtual ~Module__mdpp16_qdc_Phoswich();



        double _mdpp16_qdc_Phoswich_channel_time[16];
        double _mdpp16_qdc_Phoswich_integration_long[16];
        double _mdpp16_qdc_Phoswich_integration_short[16];
        double _mdpp16_qdc_Phoswich_trigger_time[2];
        double _mdpp16_qdc_Phoswich_module_timestamp[1];
        double _mdpp16_qdc_Phoswich_extended_timestamp[1];

        ClassDef(Module__mdpp16_qdc_Phoswich, 1);
    };

} // end namespace event0_modules

struct Event_event0: public MVMEEvent
{
    Event_event0();
    virtual ~Event_event0();

    event0_modules::Module_mdpp32_1_scp_Sector_Si1And2 mdpp32_1_scp_Sector_Si1And2;
    event0_modules::Module_mdpp32_2_scp_Sector_Si3And4 mdpp32_2_scp_Sector_Si3And4;
    event0_modules::Module_mdpp32_3_scp_Ring_Si1 mdpp32_3_scp_Ring_Si1;
    event0_modules::Module_mdpp32_4_scp_Ring_Si3 mdpp32_4_scp_Ring_Si3;
    event0_modules::Module__mdpp16_scp_PPAC_And_SB _mdpp16_scp_PPAC_And_SB;
    event0_modules::Module__mdpp16_qdc_Phoswich _mdpp16_qdc_Phoswich;

    ClassDef(Event_event0, 1);
};

struct Experiment: public MVMEExperiment
{
    Experiment();
    ~Experiment();

    Event_event0 event0;

    ClassDef(Experiment, 1);
};

#endif /* __MVME_EVENT_EXPORT_GUARD_Experiment__ */
