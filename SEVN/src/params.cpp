//
// Created by iorio on 16/03/20.
//

#include <params.h>
#include <fstream>
#define PAR(value, documentation) std::make_pair(value,documentation)

//STATIC
std::size_t SEVNpar::instance_counter=0;

//SETTABLE VALUES
std::set<std::string> SEVNpar::public_settable{};
std::set<std::string> SEVNpar::private_settable{"star_lambda","star_lambda_fth","star_lambda_pureHe", "star_tshold_WR_envelope", //Star properties
                                                "jtrack_tshold_dm_rel","jtrack_h_err_rel_max","jtrack_max_dm_factor",
                                                "jtrack_min_dm_factor","jtrack_max_iteration","jtrack_dm_step","turn_WR_to_pureHe", //Change track
                                                "gw_tshold", "gw_onlybco", //GW
                                                "eddington_factor", //general accretion
                                                "rlo_f_mass_accreted","rlo_eps_nova", //RLO
                                                "rlo_gamma_angmom","rlo_QHE","rlo_stability","rlo_max_nuclearmt","rlo_enable_collision","rlo_mtstable_ms",  //RLO
                                                "sn_co_lower_sn","sn_co_lower_ecsn","sn_co_lower_sn_pureHe", "sn_co_lower_ecsn_pureHe",
                                                "sn_max_ns_mass","sn_Mchandra","sn_min_vkick", "sn_kick_velocity_stdev",
                                                "sn_kicks","sn_pairinstability", "sn_neutrinomloss", //SN additional options
                                                "sn_compact_csi25_tshold", "sn_compact_fallback",
                                                "sn_Mremnant_average_NS", "sn_Mremnant_std_NS", //NS mass remnant
                                                "sn_Mejected_average", "sn_Mremnant_average",//SN
                                                "w_alpha","w_beta", //WInds
                                                "ts_maximum_variation","ts_min_points_per_phase","ts_min_dt","ts_max_dt", //Timesteps
                                                "ts_check_spin", "ts_check_spin_bse", "ts_check_NSspin", //Timesteps
                                                "io_literal_phases","io_logfile", //Options
                                                "ce_alpha", "ce_knce", "ce_kce", //CE
                                                "hard_rhoc","hard_sigma","hard_xi","hard_kappa","hard_mass_average", //Hardening options
                                                "wmode", "rlmode", "tmode", "gwmode", "mixmode", "collmode", "kmode", "hardmode", "cemode", "snmode","inertiamode","circmode", //Process options
                                                "Z","spin", //overwrite Z and spin
                                                "xspinmode","xspin_sigma_maxwell","xspin_bavera",//Xspin BBH
                                                "tabuse_rhe", "tabuse_rco", "tabuse_inertia", "tabuse_envconv","tabuse_Xsup", //Use extra tables options
                                                "dtout","tini","tf", //Global evolution parameters
                                                "ev_setwM", "ev_setwM_log","ev_setwM_tphase","ev_R_linear_interp", //Evolution
                                                "ev_max_repetitions", "ev_Nchunk", "ev_naked_tshold","ev_set_maxCO","ev_set_minHE",
                                                "nthreads", "name_prefix","check_stalling","check_stalling_time", //Systems
                                                "omode",
                                                "initerror_stop",
                                                "log_level", //logging
                                                "rseed", //random seed
                                                "ns_detailed_accretion","ns_magnetic_tscale","ns_magnetic_mscale", //neutron stars
                                                "list", "ibmode", "tables", "tables_HE", "o", "scol", "bcol",
                                                "use_thg_hurley", "optimistic_scenario_hg" //Others

};


///Double params
int SEVNpar::default_value_num() {

    ///STAR
    params_num["star_lambda"]              = PAR(-1,"if >0 Constant Lambda in binding energy (Eq. 69 in Hurley02). If -1 use Lambda from Claeys et al. 2014 (Appendix A).");
    params_num["star_lambda_pureHe"]       = PAR(0.5,"Constant lambda to use for the pureHe stars. Notice: some Lambda model have their own estimate of Lamdba_he (for example option -4), in this case this value is not considered");
    params_num["star_lambda_fth"]          = PAR(1,"Fraction of internal energy that goes to the binding energy. Used only if star_lambda<0. Notice that some Lambda models do not include an option for the fraction of internal energu (e.g. option -5)");
    params_num["star_tshold_WR_envelope"]  = PAR(0.021, "Relative difference threshold between envelope (Mass-MHE) and total mass to define a star as Wolf Rayet");


    //Neutron stars
    params_num["ns_magnetic_tscale"] = PAR(1000,"Magnetic field decay timescale in Myr");
    params_num["ns_magnetic_mscale"] = PAR(0.15,"Magnetic field decay mass-scale in Msun");

    //Black holes
    params_num["xspin_sigma_maxwell"]    = PAR(0.1, "options for standard deviation of Maxwellian distribution for Xspin");


    ///Tracks check
    //determined when reading the look-up tables
    params_num["max_zams"]             = PAR(2e10, "Max Zams Mass in the loaded tables");  /*!< Max Zams Mass in the loaded tables */
    params_num["min_zams"]             = PAR(1e10, "Min Zams Mass in the loaded tables");   /*!< Min Zams Mass in the loaded tables */
    params_num["max_z"]                = PAR(2e10, "Max Z Mass in the loaded tables");   /*!< Max Z  in the loaded tables */
    params_num["min_z"]                = PAR(1e10, "Min Z Mass in the loaded tables");   /*!< Min Z  in the loaded tables */
    params_num["max_zams_he"]          = PAR(2e10, "Max Zams Mass in the loaded pureHE tables");  /*!< Max Zams Mass in the loaded tables */
    params_num["min_zams_he"]          = PAR(1e10, "Min Zams Mass in the loaded pureHE tables");
    params_num["max_z_he"]             = PAR(2e10, "Max Z Mass in the loaded  pureHE tables");
    params_num["min_z_he"]             = PAR(1e10, "Min Z Mass in the loaded pureHE tables");

    ///SN
    params_num["sn_co_lower_sn"]                  = PAR(1.44, "Minimum value of the CO core  Mass  to explode as SN");
    params_num["sn_co_lower_ecsn"]                = PAR(1.38, "Minimum value of the CO core  Mass  to explode as electron capture SN");
    params_num["sn_co_lower_sn_pureHe"]           = PAR(-1, "Minimum value of the CO core  Mass  to explode as SN for pureHe star, if -1 use the same of H star");
    params_num["sn_co_lower_ecsn_pureHe"]         = PAR(-1, "Minimum value of the CO core  Mass  to explode as electron capture SN for pureHe star, if -1 use the same of H star");
    params_num["sn_max_ns_mass"]                  = PAR(3.00, "Maximum mass allowed for a NS");
    params_num["sn_Mchandra"]                     = PAR(1.44, "Chandrasekar mass limit for WD");
    params_num["sn_compact_csi25_tshold"]         = PAR(0.35,  "csi25 parameter threshold for explosion/implosion decision, if -1 use a stochastic threshold based on the results of  Patton&Sukhbold20");
    params_num["sn_compact_fallback"]             = PAR(0.9,  "Fallback fraction for implosions in the compact SN option");
    params_num["sn_Mremnant_average_NS"]          = PAR(1.33, "Mean value of the Gaussian sampling the remnant mass for NS. Notice it is not used in all SNmodel. If SN model is compact this value overwrite the parameter is used also as sn_Mremnant_average.");
    params_num["sn_Mremnant_std_NS"]              = PAR(0.09, "Std value of the Gaussian sampling the remnant mass for NS. Notice it is not used in all SNmodel");
    params_num["sn_Mejected_average"]             = PAR(-1, "Average remnant mass to use in the Unified kick model.  -1 use the default SN model value (if any)");
    params_num["sn_Mremnant_average"]             = PAR(-1,  "Average ejected mass to use in the Unified kick model. -1 use the default SN model value (if any)");
    params_num["sn_min_vkick"]                    = PAR(0,  "Minimum Vkick after a SN explotsion in km/s");
    params_num["sn_kick_velocity_stdev"]          = PAR(265.0, "Standard deviation  of the Maxwellian distribution of kick velocity (Used in the Hobbs and Unified SN kick model)");

    ///Jump track parameters
    params_num["jtrack_tshold_dm_rel"] = PAR(0.01,"relative accumulated mass threshold to change tracks");   /*!< relative accumulated mass threshold to change tracks */
    params_num["jtrack_h_err_rel_max"] = PAR(5.0e-3, "relative difference in total mass to consider the convergence reached "); /*!< relative difference in total mass to consider the convergence reached */
    params_num["jtrack_max_dm_factor"] = PAR(1.2, "The maximum new zams that will be tested is Mzams +  jtrack_DM_factor_max*dM_accumul");   /*!< The maximum new zams that will be tested is Mzams +  jtrack_DM_factor_max*dM_accumul */
    params_num["jtrack_min_dm_factor"] = PAR(0.0, "The minimum new zams that will be tested is Mzams +  jtrack_DM_factor_min*dM_accumul");    /*!< The minimum new zams that will be tested is Mzams +  jtrack_DM_factor_min*dM_accumul */
    params_num["jtrack_max_iteration"] = PAR(10.0, "Maximum, number of iteration to found the convergence");   /*!< Maximum, number of iteration to found the convergence */
    params_num["jtrack_dm_step"]       = PAR(0.1, "Mass increment when trying to matche the core ans in case the binding energy"); /*!< Mass increment when trying to matche the core ans in case the binding energy */

    ///GWcheck parameters
    params_num["gw_tshold"]            = PAR(1.0, "Time in units of Hubble time. If GWtime<GW_TSHOLD*tHubble, enable GW process");    /*!< Time in units of Hubble time. If GWtime<GW_TSHOLD*tHubble, enable GW process*/

    ///Winds parameters
    params_num["w_alpha"]              = PAR(1.5, "factor to tune mass accretion through winds (Eq.6 Hurley+02)");
    params_num["w_beta"]               = PAR(0.125, "factor to tune stellar wind velocity (Eq.9 in Hurley+02)");

    ///Timestep
    params_num["ts_maximum_variation"]    = PAR(0.05,"Relative maximum variation of stellar and binary properties used in the adaptive time step");
    params_num["ts_min_points_per_phase"] = PAR(10,"Set the maximum time step so that we have at least N points evaluated in each phase");
    params_num["ts_min_dt"]               = PAR(-1,"Force the adaptive timestep to be larger than this value, it will it has the priority on any other option, -1 means that the option is disabled");
    params_num["ts_max_dt"]               = PAR(-1,"Force the adaptive timestep to be smaller than this value, it has the priority on any other option, -1 means that the option is disabled");

    ///Roche Lobe
    params_num["rlo_f_mass_accreted"]   = PAR(0.5, "Fraction of mass lost trough the  RLO that is accreted on the other star");
    params_num["rlo_eps_nova"]          = PAR(0.001, "Fraction of accreted matter retained in nova eruption");
    params_num["rlo_gamma_angmom"]      = PAR(-2, "Parameter to manage the ang mom loss in hurley_rlo, if -1 "
                                                  "angmom is lost from the primary (Jeans mode), -2 angmom is lost from the secondary (re-isotropic emission),"
                                                  "otherwise if>0 gamma is the fraction of ang mom lost from the system.");
    params_num["rlo_max_nuclearmt"]     = PAR(5,"Max value of the mass to use in the normalisation of the nuclear mass transfer (Eq. 59 Hurley+02)");

    //General accretion
    params_num["eddington_factor"]      = PAR(1.0, "Eddington factor to limit accretion on a compact object. 1 means accretion exactly at the Eddington limit, >1 is super Eddington accretion");


    ///Kicks parameters
    //params_num["sn_Mejected_average"] = PAR(10.45,"Average ejected mass for the delayed model");
    //params_num["sn_Mremnant_average"] = PAR(1.36,"Average remnant mass for the delayed model");
    //params_num["sn_Mejected_average"]   = PAR(10.90,"Average ejected mass for the rapid model");
    //params_num["sn_Mremnant_average"]   = PAR(1.27,"Average remnant mass for the rapid model");


    ///CE
    params_num["ce_alpha"]  = PAR(3,"alpha in binding energy (Eq. 73 in Hurley02)");
    params_num["ce_knce"]    = PAR(1, "Fraction of non core mass not participating to the CE (e.g. a MS star) retained after the CE coalescence."
                                      "If -1, use the eq. 77 in Hurley 2002 (ce_kce is ignored unless it is -1, in this case the SEVN1 binding energy method is used)");
    params_num["ce_kce"]   = PAR(1, "Fraction of non core mass  participating to the CE (e.g. envelope of giants) retained after the CE coalescence."
                                    "If -1, use a rescaled version of eq. 77 In Hurley just for che CE mass");


    ///Hardening
    params_num["hard_rhoc"]              =  PAR(3.99e4,"central density of the cluster in Msun/pc^3");
    params_num["hard_sigma"]             =  PAR(5.0, "3D velocity dispersion of the cluster in km/s");
    params_num["hard_xi"]                =  PAR(3.0, "xi parameter for fastcluster hardening option (Eq. 11 Mapelli+21)");
    params_num["hard_kappa"]             =  PAR(0.1, "kappa parameter for fastcluster hardening option (Eq. 13 Mapelli+21)");
    params_num["hard_mass_average"]      =  PAR(1, "Average mass of the perturber stars in the environment producing the hardening in Msun");

    ///Evolution
    params_num["ev_max_repetitions"]  = PAR(100, "Maximum number of repetitions allowed in the sse and bse. If we reach this number an error is raised");
    params_num["ev_naked_tshold"]     = PAR(1E-4,"Mass difference threshold (Msun) between envelope and core to set a star as nakedHe or nakedCO.");
    params_num["ev_Nchunk"]           = PAR(1000, "Evolve Nchunk at time");




    //Systems
    params_num["nthreads"]              = PAR(1, "Number of threads to be used");
    params_num["check_stalling_time"]   = PAR(5, "Time in seconds (units granularity) after which the system is considered stalled and the evolution is stopped and an error is raised)");



    return EXIT_SUCCESS;
}

int SEVNpar::default_value_str() {

    //SN KICKS
    params_str["sn_kicks"]               = PAR("unified", "SN kick model");
    params_str["sn_pairinstability"]     = PAR("mapelli20", "SN pair instability model");
    params_str["sn_neutrinomloss"]       = PAR("lattimer89", "SN neutrino mass loss model");


    ///INPUT OPTIONS
    params_str["wmode"]                  = PAR("hurleyTT24mod", "Option for Wind mass transfer Process");
    params_str["rlmode"]                 = PAR("hurley_rl", "Option for Roche Lobe mass transfer Process");
    params_str["tmode"]                  = PAR("tides_simple", "Option for Tides  Process");
    params_str["gwmode"]                 = PAR("peters", "Option for GW decay  Process");
    params_str["mixmode"]                = PAR("simple", "Option for stellar mix Process");
    params_str["collmode"]               = PAR("hurley", "Option for stellar collision at periastron");
    params_str["hardmode"]               = PAR("disabled", "Option for hardening in clusters");
    params_str["circmode"]               = PAR("angmom", "Option for orbit circularisation at the onset of RLO");
    params_str["kmode"]                  = PAR("hurley_kick", "Option for SN kick Process");
    params_str["cemode"]                 = PAR("energy", "Option for Common Envelope Process");
    params_str["snmode"]                 = PAR("list", "if list use the snmode in input, otherwise ust this option for all the stars.");
    params_str["inertiamode"]            = PAR("Hurley","option for inertia estimate when tabuse_inertia is false");
    params_str["xspinmode"]              = PAR("disabled","options for black hole spin models");

    ///Roche Lobe
    params_str["rlo_stability"]          = PAR("qcrit_Hradiative_stable", "Option for RLO mass transfer stability");

    ///Overwrite Z
    params_str["Z"]                      = PAR("list", "if list use the Z in the input list, otherwise ust this option for all the stars.");
    //Overwrite Spin
    params_str["spin"]                   = PAR("list", "if list use the spin in the input list, otherwise ust this option for all the stars.");

    ///Evolution
    params_str["dtout"]                  = PAR("list", "If list use the dtout reported in the input list, otherwise use this value for all the stars and binaries.");
    params_str["tini"]                   = PAR("list", "If list use the tini reported in the input list, otherwise use this value for all the stars.");
    params_str["tf"]                     = PAR("list", "If list use the tf reported in the input list, otherwise use this value for all the stars and binaries.");

    params_str["ev_setwM"]               = PAR("linear", "Option to set weights for mass interpolation");
    params_str["ev_setwM_log"]           = PAR("log", "Option to set weights for mass interpolation for log properties (i.e. Radius,Luminosity,Inertia)");
    params_str["ev_setwM_tphase"]        = PAR("rational", "Option to set weights for estimating the stellar phases times");

    ///System
    params_str["name_prefix"]            = PAR("","prefix to add to the name of the systems");

    ///Output
    params_str["omode"]                  = PAR("csv", "Define the results output format (ascii or csv)");
    params_str["o"]                      = PAR("sevn_output", "Complete path to the output folder. It will be created by the code");
    params_str["scol"]                   = PAR("Worldtime:Mass:MHE:MCO:Radius:Phase:RemnantType", "Additional columns to print in the output file for single evolution runs. Default is empty, but any property of single stars can be added (check names in the Property class)");
    params_str["bcol"]                   = PAR("Semimajor:Eccentricity:BEvent", "Additional columns to print in the output file for binary evolution runs. Default is empty, but any property of binary stars can be added (check names in the Property class)");

    ///Input
    params_str["ibmode"]                 = PAR("new", "Input file format for binaries [new*] [legacy] [sevn1]");


    ///Folders
    params_str["tables"]                 = PAR(utilities::get_absolute_SEVN_path()+"/tables/SEVNtracks_parsec_ov05_AGB", "Complete path to look-up tables");
    params_str["tables_HE"]              = PAR(utilities::get_absolute_SEVN_path()+"/tables/SEVNtracks_parsec_pureHe36", "Complete path to look-up tables of pure-HE stars");
    params_str["list"]                   = PAR("", "Complete path to input file (list of binaries or single stars)");
    params_str["myself"]                 = PAR(utilities::get_absolute_SEVN_path(), "Complete path to the SEVN folder");



    ///LOG
    #ifdef DEBUG
    params_str["log_level"]               = PAR("debug", "Log output level: debug, info, warning, error");
    #else
    params_str["log_level"]               = PAR("error", "Log output level: debug, info, warning, error");
    #endif

    return EXIT_SUCCESS;
}

int SEVNpar::default_value_bool(){

    //Timesteps
    params_bool["ts_check_spin"]          = PAR(false,"If true take into account the variation (SSE only) of OmegaSpin in the adaptive timestep");
    params_bool["ts_check_spin_bse"]      = PAR(false,"If true take into account the variation (BSE only) of OmegaSpin in the adaptive timestep");
    params_bool["ts_check_NSspin"]        = PAR(false,"If true take into account the variation of OmegaRem for NS in the adaptive timestep. It should be set to true if interested on pulsars");



    //run options
    params_bool["initerror_stop"]         = PAR(false, "If true terminate the run when a error on a system initialisation is thrown");

    ///IO options
    params_bool["io_literal_phases"]      = PAR(false, "If true print the phase as literal otherwise as integer");
    params_bool["io_logfile"]             = PAR(true,"If true output the logfile");

    //Interpolation
    params_bool["ev_R_linear_interp"]        = PAR(false, "if true interpolate the Radius using a linear scaling (default uses log scale)");
    params_bool["turn_WR_to_pureHe"]         = PAR(false, "if true transform the star to pureHe as soon as it is flagged as WR (see parameter star_tshold_WR_envelope)");


    ///Evolution
    params_bool["ev_set_maxCO"]           = PAR(false,"If true the first time a star develops a CO core, we set the maximum CO core Mass for SSE as the last value of the interpolating tracks");
    params_bool["ev_set_minHE"]           = PAR(false,"If true the first time a star develops a CO core, we set the minimum HE core Mass for SSE as the last value of the interpolating tracks");

    ///GWcheck parameters
    params_bool["gw_onlybco"]              = PAR(false, "If true activate the GW orbital decay on for binary compact objects");


    ///RLO options
    params_bool["rlo_QHE"]                = PAR(false, "If true enable the Quasi Homogeneous Evolution  after a RLO mass transfer following Elrdige&Stanway11");
    params_bool["rlo_enable_collision"]   = PAR(false, "If true allow collision at periastron during RLO");
    params_bool["rlo_mtstable_ms"]        = PAR(false,  "If true mass transfer from radiative MS and pureHE MS are always stable");


    ///Extra tables options
    //Option to use or not extra tables of stellar properties
    params_bool["tabuse_rhe"]              = PAR(true, "If true interpolate RHE from tables (if present)");
    params_bool["tabuse_rco"]              = PAR(true, "If true interpolate RCO from tables (if present)");
    params_bool["tabuse_inertia"]          = PAR(false, "If true use the inertia from tables (if present)");
    params_bool["tabuse_envconv"]          = PAR(true, "If true estimate the properties of the convective envelope using the tables (xxxconv.dat, optional tables)");
    params_bool["tabuse_Xsup"]             = PAR(false, "If true use the surface abundance tables (xxxsup.dat, optional tables)");

    //Use HG time from Hurley
    params_bool["use_thg_hurley"]          = PAR(false, "If true estimate the HG time from the Hurley+00 functional forms instead of using the convective envelope");
    ///Optimistic scenario
    params_bool["optimistic_scenario_hg"]  =  PAR(false, "If true enable optimistic scenario for HG, i.e. allow to start a CE after a unstable mass transfer");

    //NS accretion details
    params_bool["ns_detailed_accretion"]   = PAR(false, "If true activate a detailed handling of the mass accretion on Neutron star accounting for the variation of spin and magnetic field");

    //SEED options
    params_bool["rseed"]                   = PAR(false, "If true the random seed is given in input in the source list (last column)");

    //BH options
    params_bool["xspin_bavera"]                  = PAR(false, "If true the Bavera correction for the black-hole spin is applied");


    params_bool["check_stalling"]                = PAR(true, "If true check stalling stars. If the elapsed evolution time is larger than 20s an error is thrown");


    return 0;

}


///load
int SEVNpar::load(int n, char **val, bool initialise){

    if (initialise) init(); //Start from init() to be sure that the values that are not defined got the default value

    //Here make a check, if everything it is ok we should have a pair of value for each
    // runtime option (option-value) + the name of the executable, so if everything is ok it has to be odd.
    // Therefore if it is even raise an error.
    if (n%2==0){
        svlog.critical("Oops it seems there is an error in the list of runtime parameter: the total number"
                       "of option-value arguments is an odd number (maybe a value or a space separator missing?), check it.",
                       __FILE__,__LINE__,sevnstd::sevnio_error(""));
    }


    for (int i = 1; i<n; i+=2) {
        std::string argv_name = std::string(&val[i][1]); //starts from 1 to remove the - sign
        std::string svalue = std::string(&val[i+1][0]);
        set_from_string(argv_name,svalue);
    }
    if (!check()) //Check error in the inserted values
    {
        svlog.critical("Parameters check failed",__FILE__,__LINE__,sevnstd::params_error());
    }


    return EXIT_SUCCESS;

}