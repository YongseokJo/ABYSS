//
// Created by Giuliano Iorio on 2020-02-21.
//

#ifndef SEVN_STATIC_MAIN_H
#define SEVN_STATIC_MAIN_H

#include <property.h>
#include <BinaryProperty.h>
#include <Orbit.h>
//Processes
#include <Processes.h>
#include <Collision.h>
#include <Hardening.h>
#include <Circularisation.h>
#include <MTstability.h>
#include <Windaccretion.h>
//SN
//SN models
#include <supernova.h>
#include <rapid.h>
#include <delayed.h>
#include <compactness.h>
#include <directcollapse.h>
#include <deathmatrix.h>
//Kick models
#include <kicks.h>
#include <zeros.h>
#include <specialkicks.h>
#include <hobbs.h>
#include <unified.h>
//Pair instability models
#include <pairinstability.h>
#include <pimapelli20.h>
#include <piiorio22.h>
#include <pifarmer19.h>
//Neutron star models
#include <neutrinomassloss.h>
#include <nmllattimer89.h>

//TODO At the moment the order in which the processes and properties are evolved depend on the order we define them here.
//TODO It is better to put some check in the derived properties to be sure that the properties that should have been evolved before have been already evolved



/*******************************************
 ******** PROPERTY STATIC INIT *************
 *******************************************/

    /** Instruction for adding new property  PROP
     * 1- Add the size_t PROP::ID
     * 2- Initialise the static (fake) instance PROP PROP::_prop.
     * NB: take care of the order of initilisation. The properies will be evolved exactly in the same order
     * in which they are initialised here
     */

    //TODO In all the other cases _size is called just size. Change _size to size?
    /** Define:
     * @size: total number of instances of Property and derived classes (both fake and real).
     * @all: vector containing all the (pointers) to the Property.
     * @PrintMap: map containing the pair (Property_name, Property_id) for output purpose.
     * Note all and PrintMap are filled during the instantiation of the fake Processes (see below)
     * */
    size_t Property::_size = 0;  //setting the initial size at 0.
    vector<Property*> Property::all;
    Property::_PrintMap Property::PrintMap;

    /** Define all the Properties ID
    * */
    size_t Mass::ID;
    size_t Radius::ID;
    size_t Luminosity::ID;
    size_t Inertia::ID;
    size_t NextOutput::ID;
    size_t dRdt::ID;
    size_t Phase::ID;
    size_t RemnantType::ID;
    size_t Localtime::ID;
    size_t Timestep::ID;
    size_t Temperature::ID;
    size_t MCO::ID;
    size_t MHE::ID;
    size_t RCO::ID;
    size_t RHE::ID;
    size_t Hsup::ID;
    size_t HEsup::ID;
    size_t Csup::ID;
    size_t Nsup::ID;
    size_t Osup::ID;
    size_t Qconv::ID;
    size_t Tconv::ID;
    size_t Depthconv::ID;
    size_t Spin::ID;
    size_t Xspin::ID;
    size_t dMdt::ID;
    size_t dMHEdt::ID;
    size_t dMCOdt::ID;
    size_t dMcumul_binary::ID;
    size_t dMcumul_RLO::ID;
    size_t dMcumulacc_wind::ID;
    size_t Bmag::ID;
    size_t OmegaRem::ID;
    size_t Rs::ID;
    size_t Worldtime::ID;
    size_t Lambda::ID;
    size_t Ebind::ID;
    size_t PhaseBSE::ID;
    size_t Zams::ID;
    size_t Event::ID;
    size_t AngMomSpin::ID;
    size_t OmegaSpin::ID;
    size_t Zmet::ID;
    size_t dMRLOdt::ID;
    size_t dMaccwinddt::ID;
    size_t Plife::ID;


    //TODO to be removed
    size_t NSsalpha::ID;


    /// To be added
    //size_t Mconv::ID;
    //size_t Rconv::ID;
    //size_t Tconv::ID;

    /**
     * Initialisation of the static (fake) instance
     * NB: take care of the order of initialisation. The properties will be evolved exactly in the same order
     * in which they are initialised here
     */
    Localtime Localtime::_localtime; //special evolve function (outside the main evolve() loop)
    Phase Phase::_phase;//special evolve function (outside the main evolve() loop)
    NextOutput NextOutput::_nextoutput; //special evolve function (outside the main evolve() loop)
    Worldtime Worldtime::_worldtime; //time evolution of the simulation (it always starts from zero, independently of the localtime of the star)

    //fundamental stellar parameters (They evolve from lookup tables)
    //WARNING: Mass and Radius needs to be the first two!
    Mass Mass::_mass;
    Radius Radius::_radius;
    MHE MHE::_masshe;
    MCO MCO::_massco;
    Luminosity Luminosity::_luminosity;
    ///////////////////////////////////
    //Optional  (they can be table or derived from above properties)
    RHE RHE::_rhe;
    RCO RCO::_rco;
    Inertia Inertia::_inertia;
    Hsup Hsup::_hsup;
    HEsup HEsup::_hesup;
    Csup Csup::_csup;
    Nsup Nsup::_nsup;
    Osup Osup::_osup;
    Qconv Qconv::_qconv;
    Tconv Tconv::_tconv;
    Depthconv Depthconv::_depthconv;

    //Non tab property
    RemnantType RemnantType::_remnanttype;
    Bmag Bmag::_bmag;
    OmegaRem OmegaRem::_omegarem; //Notice OmegaRem evolves after Bmag since it depends on the updated value of Bmag


    //derived stellar parameters. NB: They have to evolve after fundamental stellar parameters.
    Xspin Xspin::_xspin;
    dMdt dMdt::_dmdt;
    dMHEdt dMHEdt::_dmhedt;
    dMCOdt dMCOdt::_dmcodt;
    dRdt dRdt::_drdt;
    Temperature Temperature::_temperature;
    Rs Rs::_rs;
    AngMomSpin AngMomSpin::_angmomspin;
    OmegaSpin OmegaSpin::_omegaspin;
    Spin Spin::_spin;

    //Properties that are updated only in binary
    dMcumul_binary dMcumul_binary::_dMcumul_binary;
    dMcumul_RLO    dMcumul_RLO::_dMcumul_RLO;
    dMcumulacc_wind dMcumulacc_wind::_dMcumulacc_wind;

    //JIT property
    Lambda Lambda::_lambda;
    Ebind Ebind::_ebind;
    PhaseBSE PhaseBSE::_phasebse;
    Zams Zams::_zams;
    Event Event::_event;
    Zmet Zmet::_zmet;
    dMRLOdt dMRLOdt::_dmorlodt;
    dMaccwinddt dMaccwinddt::_dmaccwinddt;
    Plife Plife::_plife;

    //TODO to be removed;
    NSsalpha NSsalpha::_nssalpha;

    //predict the next time step based on the calculated variations of stellar parameters
    //WARNING:Timestep needs to be the last one!
    Timestep Timestep::_timestep; //this needs be the very last element to evolve
/***************************************************************************************/

/**************************************************
 *********** SN PRESCRIPTIONS   *******************
 **************************************************/
    //SN explosions
    delayed delayed::_delayed;
    delayed_gauNS _delayed_gauns;
    rapid rapid::_rapid;
    rapid_gauNS _rapid_gauns;
    compactness compactness::_compactness;
    directcollapse directcollapse::_directcollapse;
    DeathMatrix DeathMatrix::_deathmatrix;
    //PISN prescriptions
    PIDisabled _ppidisabled;
    PIMapelli20 _pimapelli20;
    PIIorio22 PIIorio22::_piiorio22;
    PIIorio22Limited PIIorio22Limited::_piiorio22limited;
    PIFarmer19 PIFarmer19::_pifarmer19;
    // SN kick prescriptions
    Hobbs Hobbs::_hobbs;
    HobbsPure HobbsPure::_hobbspure;
    Unified Unified::_unified;
    EC15CC265 EC15CC265::_ec15cc265;
    Zeros Zeros::_zeros;
    CC15 CC15::_cc15;
    ECUS30 ECUS30::_ecus30;
    //SN presciptions
    NMLDisabled NMLDisabled::_nmldisabled;
    NMLLattimer89 NMLLattimer89::_nmllattimer89;
/***************************************************************************************/

/*************************************************************
 ****** MASS TRANSFER STABILITY PRESCRIPTIONS ****************
 *************************************************************/
    /***                    Always stable              ***/
    MT_Stable                          MT_Stable::_mt_stable;
    /***                    Always unstable              ***/
    MT_UnStable                        MT_UnStable::_mt_unstable;
    /***                   Qcrit style stability              ***/
    Qcrit_Hurley                      Qcrit_Hurley::_qcrit_hurley;
    Qcrit_Hurley_Webbink              Qcrit_Hurley_Webbink::_qcrit_hurley_webbink;
    Qcrit_Hurley_Webbink_Shao         Qcrit_Hurley_Webbink_Shao::_qcrit_hurley_webbink_shao;
    Qcrit_COSMIC_Neijssel             Qcrit_COSMIC_Neijssel::_qcrit_cosmic_neijssel;
    Qcrit_COSMIC_Claeys               Qcrit_COSMIC_Claeys::_qcrit_cosmic_claeys;
    Qcrit_StarTrack                   Qcrit_StarTrack::_qcrit_startrack;
    Qcrit_Radiative_Stable            Qcrit_Radiative_Stable::_qcrit_radiative_stable;
    Qcrit_HRadiative_Stable           Qcrit_HRadiative_Stable::_qcrit_hradiative_stable;


/**************************************************
 ******** BINARY PROPERTY STATIC INIT *************
 **************************************************/


    /** Instruction for adding new property  BPROP
     * 1- Add the size_t BPROP::ID
     * 2- Initialise the static (fake) instance BPROP BPROP::_bprop.
     * NB: take care of the order of initilisation. The properties will be evolved exactly in the same order
     * in which they are initialised here
     */

    /** Define:
     * @size: total number of instances of BinaryProperty and derived classes (both fake and real).
     * @all: vector containing all the (pointers) to the BinaryProperty.
     * @PrintMap: map containing the pair (BinaryProperty_name, BinaryProperty_id) for output purpose.
     * Note all and PrintMap are filled during the instantiation of the fake Processes (see below)
     * */
     size_t BinaryProperty::size = 0;



//setting the initial size at 0.
     std::vector<BinaryProperty*> BinaryProperty::all;
     BinaryProperty::_PrintMap BinaryProperty::PrintMap;

    /** Define all the BinaryProperty ID
     */
    size_t Eccentricity::ID;
    size_t Semimajor::ID;
    size_t BTimestep::ID;
    size_t BWorldtime::ID;
    size_t dadt::ID;
    size_t dedt::ID;
    size_t AngMom::ID;
    size_t Period::ID;
    size_t GWtime::ID;
    size_t RL0::ID;
    size_t RL1::ID;
    size_t BEvent::ID;

    /**
    * Initialisation of the static (fake) instance
    * NB: take care of the order of initialisation. The binary properties will be evolved exactly in the same order
    * in which they are initialised here
    */

    //Stuff that does not evolve but just use special_evolve
    BWorldtime BWorldtime::_bworldtime;

    //fundamental binary parameters (they need to be evolved before other properties)
    Eccentricity Eccentricity::_eccentricity;
    Semimajor Semimajor::_semimajor;

    //derived binary parameters from fundamental
    dadt dadt::_dadt;
    dedt dedt::_dedt;
    AngMom AngMom::_angmom;
    Period Period::_period;
    GWtime GWtime::_gwtime;
    RL0 RL0::_rl0;
    RL1 RL1::_rl1;
    BEvent BEvent::_bevent;


    //predict the next time step based on the calculated variations of stellar and binary parameters
    BTimestep BTimestep::_btimestep; //this needs be the very last element to evolve
/***************************************************************************************/


/*******************************************
******** PROCESS STATIC INIT *************
*******************************************/

    /** Instruction for adding new property  PROC
     * 1- Add the size_t PROC::ID
     * 2- Initialise the static (fake) instance PROC PROC::_proc.
     * NB: take care of the order of initialisation. The processes will be evolved exactly in the same order
     * in which they are initialised here
     */

    /** Define:
     * @size: total number of instances of Process and derived classes (both fake and real).
     * @all: vector containing all the (pointers) to the processes.
     * @PrintMap: map containing the pair (process_name, process_id) for output purpose.
     * Note all and PrintMap are filled during the instantiation of the fake Processes (see below)
     * */
    size_t Process::size = 0;
    std::vector<Process*> Process::all;
    Process::_PrintMap Process::PrintMap;

    /** Define all the Process ID
     */
    size_t CommonEnvelope::ID;
    size_t RocheLobe::ID;
    size_t SNKicks::ID;
    size_t GWrad::ID;
    size_t Windaccretion::ID;
	size_t Tides::ID;
	size_t Mix::ID;
    size_t Kollision::ID;
    size_t Hardening::ID;
    size_t Circularisation::ID;


    /**
    * Initialisation of the static (fake) instance
    * NB: take care of the order of initialisation. The processes will be evolved exactly in the same order
    * in which they are initialised here
    */
    Windaccretion Windaccretion::_windaccretion;
    SNKicks SNKicks::_snkicks;
    RocheLobe RocheLobe::_rochelobe;
    GWrad GWrad::_gwrad;
	Tides Tides::_tides;
	Mix Mix::_mix;
    Kollision Kollision::_kollision;
    Hardening Hardening::_hardening;
    Circularisation Circularisation::_circularisation;
    CommonEnvelope CommonEnvelope::_commonenvelope;

    //Processes that have been already refactored
    WindaccretionHurley WindaccretionHurley::_windaccretionhurley;
    WindaccretionDisabled WindaccretionDisabled::_windaccretiondisabled;
    WindaccretionHurleyTT24mod WindaccretionHurleyTT24mod::_windaccretionhurleytt24mod;

    KollisionDisabled KollisionDisabled::_kollisiondisabled;
    KollisionHurley KollisionHurley::_kollisionhurley;

    HardeningDisabled HardeningDisabled::_hardeningdisabled;
    HardeningFastCluster HardeningFastCluster::_hardeningfastcluster;

    CircularisationDisabled CircularisationDisabled::_circularisationdisabled;
    CircularisationPeriastron CircularisationPeriastron::_circularisationperiastron;
    CircularisationPeriastronFull CircularisationPeriastronFull::_circularisationperiastronfull;
    CircularisationAngMom CircularisationAngMom::_circularisationangmom;
    CircularisationSemimajor CircularisationSemimajor::_circularisationsemimajor;

/*******************************************
 ****** ORBIT CHANGE STATIC INIT ***********
 *******************************************/

/** Instruction for adding new orbit changes ORBC
 * Note this is a bit different with respect to Properties, BinaryProperties and Processes.
 * 1- Just call ORBC ORBC::_orbc
 * 2- Add the name of the new orb_change in the related map name (follow the windsmap_name example)
 */

///Tides
//map_name defined in lookup_and_phases.h
const Lookup::TIDESMAP_NAME  Lookup::tidesmap_name{
        {TidesMode::_Tsimple, "simple"},
        {TidesMode::_Tsimple_notab, "simple_notab"},
        {TidesMode::_Tdisabled, "disabled_tides"}
};

disabled_tides disabled_tides::_disabled_tides;
Tides_simple Tides_simple::_Tides_simple;
Tides_simple_notab Tides_simple_notab::_Tides_simple_notab;

///GW
//map_name defined in lookup_and_phases.h
const Lookup::GW_NAME  Lookup::gwmap_name{
        {GWMode::_GWdisabled, "disabled_GW"},
        {GWMode::_GWPeters, "Peters_GW"}
};

disabled_gw disabled_gw::_disabled_gw;
Peters_gw Peters_gw::_peters_gw;

//RL
//map_name defined in lookup_and_phases.h
const Lookup::RL_NAME Lookup::rlmap_name{
        {RLMode::_RLdisabled,"disabled_RL"},
        {RLMode::_RLHurley,"Hurley_RL"},
        {RLMode::_RLHurleyBSE,"Hurley_RL_BSE"},
        {RLMode::_RLHurleymod,"Hurley_mod_RL"},
        {RLMode::_RLHurleyprop,"Hurley_RL_prop"},
};

disabled_rl disabled_rl::_disabled_rl;
Hurley_rl Hurley_rl::_Hurley_rl;
Hurley_rl_bse Hurley_rl_bse::_Hurley_rl_bse;
Hurley_mod_rl Hurley_mod_rl::_Hurley_mod_rl;
Hurley_rl_propeller Hurley_rl_propeller::_RLHurleyprop;

//Mix
const Lookup::MIX_NAME Lookup::mixmap_name{
        {MixMode::_Mixdisabled,"disabled_mix"},
        {MixMode::_Mixsimple,"simple_mix"},
};

disabled_mix disabled_mix::_disabled_mix;
simple_mix simple_mix::_simple_mix;

//SN
const Lookup::SNK_NAME Lookup::snkmap_name{
        {SNKickMode::_SNKickdisabled,"disabled_snk"},
        {SNKickMode::_SNKickHurley,"hurley_snk"},
};

disabled_SNKicks disabled_SNKicks::_disabled_SNKicks;
Hurley_SNKicks Hurley_SNKicks::_Hurley_SNKicks;

//CE
const Lookup::CE_NAME Lookup::cemap_name{
        {CEMode::_CEdisabled,"disabled_ce"},
        {CEMode::_CEEnergy,"energy_ce"},
};

disabled_CE  disabled_CE::_disabled_CE;
energy_CE   energy_CE::_energy_CE;




#endif //SEVN_STATIC_MAIN_H
