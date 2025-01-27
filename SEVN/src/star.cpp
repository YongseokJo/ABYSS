//
// Created by Giuliano Iorio on 2020-03-11.
//

//
// Created by Giuliano Iorio on 2020-03-11.
//

#include <star.h>
#include <string>
#include <params.h>
#include <supernova.h>
#include <remnant.h>
#include <BSEintegrator.h>

//0-JUMP_CONVERGE 1-JUMP 2-NO_JUMP
utilities::jump_convergence Star::find_new_track(bool is_merger){


    utilities::jump_convergence outcome=utilities::NO_JUMP;

    //If remnant or empty just exit
    if (amiremnant() or amiempty())
        return  outcome;


    ///1- CHECK IF THE STAR HAS BECOME A pureHE
    //1 - pureHE because it became WR (only if turn_WR_to_pureHe is enabled)
    bool has_become_WR = get_svpar_bool("turn_WR_to_pureHe") and amiWR() and !aminakedhelium() and !aminakedco() and !isempty; //Check if the star has become just now a naked helium star
    //1a - pureHE because it became a Naked Helium
    bool has_become_He_naked = getp(MHE::ID)+get_svpar_num("ev_naked_tshold")>=getp(Mass::ID) and !aminakedhelium() and !aminakedco() and !isempty; //Check if the star has become just now a naked helium star
    //1b - pureHE due to the QHE evolution
    bool evolve_after_MS_QHE = amifollowingQHE() and getp(Phase::ID)==TerminalMainSequence and getp_0(Phase::ID)==MainSequence;
    ///2- CHECK IF THE STAR HAS BECOME a NAKED CO
    bool has_become_CO_naked = false;
    //Enter here only if the star has not already become a HE naked, otherwise it could happens that the star becomes immediately a naked CO
    //even if it should pass before through the pureHE phase.
    //Notice that it possible to enter in this check even if abs(MCO-MHE)<get_svpar_num("ev_naked_tshold"),
    //in this case we can force the star to become direcly nakedCO withouth passing through the pureHE phase
    if (!has_become_He_naked or getp(MCO::ID)+get_svpar_num("ev_naked_tshold")>=getp(MHE::ID)){
        has_become_CO_naked = getp(MCO::ID)+get_svpar_num("ev_naked_tshold")>=getp(Mass::ID) and !aminakedco() and !isempty;
    }



    std::string log_mess;


    //1- Check if Naked CO (First check do not jump
    if (has_become_CO_naked){
        log_mess = log_mess_COnaked();
        set_COnaked();
        outcome=utilities::JUMP_CONVERGE;
    }
    //2- Check if Naked He
    else if (has_become_He_naked or has_become_WR){
        outcome=jump_to_pureHE_tracks();
        log_mess = log_mess_HEnaked(outcome);
    }
    //3- Check if it becames a Naked Helium after QHE evolution
    else if (evolve_after_MS_QHE){
        //Before to jump make all the Mass as MHE
        properties[MHE::ID]->copy_V_from(properties[Mass::ID]);
        properties[RHE::ID]->copy_V_from(properties[Radius::ID]);
        reset_QHE();
        //utilities::hardwait("PLPL",amifollowingQHE(),getp(MHE::ID),getp(RHE::ID),__FILE__,__LINE__);
        outcome=jump_to_pureHE_tracks();
        log_mess = log_mess_HEnaked(outcome);
    }
    //4- If we are changing not to become a NakedHe or Naked CO, check the times, if we are very close to the end do not jump
    else if (std::abs(plife()-1)<utilities::TINY and !force_jump){
        svlog.warning("Stars did not change tracks because it is very close to  a change of phase (plife approx 1).",__FILE__,__LINE__);
    }
    //5- If a core is present (and it is not a pureHE) enter  here only if the jump has been forced.
    else if ( !aminakedco() and (getp(MCO::ID)>1e-3 or (getp(MHE::ID)>1e-3 and !aminakedhelium()) )  ){
        if (force_jump){
            double best_rel_err=1e30;
            outcome=match_core(best_rel_err, is_merger);
            log_mess = log_mess_jumped(outcome);
        }
    }
    //6- The normal case: No naked Helium, no naked CO, no Core, no close to the end,  no remnant just a change during the MS
    else if (!amiremnant() and !aminakedco()){


        double dMcumul_binary_frac=getp(dMcumul_binary::ID)/getp(Mass::ID);

        //Check If the accumulated mass (positive or negative) is larger than a certain percent of the total Mass
        //Enter also if the jump has been forced
        if (std::abs(dMcumul_binary_frac)>=get_svpar_num("jtrack_tshold_dm_rel") or force_jump){

            double old_zams=get_zams();
            double old_lt = getp(Localtime::ID);

            double best_rel_err=1e30;

            //auto t1 = std::chrono::high_resolution_clock::now();
            outcome=match_M(best_rel_err, is_merger);
            //auto t2 = std::chrono::high_resolution_clock::now();
            //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
            //utilities::wait("Duration",duration);
            log_mess = log_mess_jumped(outcome);

            std::string outcome_label= outcome==0 ? "Converge=True" : "Converge=False";

            svlog.pdebug("Stars ID:", get_ID(), "Name:", get_name(),
                        "changed tracks at Worldtime",getp(Worldtime::ID),
                        "\n Accumulated mass was",dMcumul_binary_frac*getp(Mass::ID),"(",100*dMcumul_binary_frac,"% Mtot)",
                        "\n New Mzams=", get_zams(), "(old zams=",old_zams,")",
                        "\n New LT=", getp(Localtime::ID), "(old LT=",old_lt,")",
                        "\n dM_frac=", (get_zams()-old_zams)/(dMcumul_binary_frac*getp(Mass::ID)),
                        "\n Relative error=",100*best_rel_err,"%", outcome_label);
        }

    }


    force_jump = false; // Reset force_jump
    //In case of jump
    if (outcome!=utilities::NO_JUMP){
        properties[dMcumul_binary::ID]->reset(); //If the star jumps to a new track reset the cumulative dM counter
        //set_rzams(); //If the star jumps to a new track estimate the new rzams, but do not change rzams0
    }
    if (!log_mess.empty()){
        this->print_to_log(log_mess);
    }

    return outcome;

}

utilities::jump_convergence Star::jump_to_pureHE_tracks(){


    /**************** Some checks ******************/
    if (getp(MHE::ID)<=0)
        svlog.critical("Called jump_to_pureHE_tracks for a star with no HE core at time "+
        utilities::n2s(getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::jtrack_error());
    /*
    else if (std::abs(getp(Mass::ID)-getp(MHE::ID))>+get_svpar_num("ev_naked_tshold"))
        svlog.critical("Called jump_to_pureHE_tracks for a star with a total Mass ( "+
                       utilities::n2s(getp(Mass::ID),__FILE__,__LINE__) +") that is larger with respect to the Helium Mass more than the threshold (" +
                       utilities::n2s(getp(MHE::ID),__FILE__,__LINE__)+"). Worldtime "+
                       utilities::n2s(getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::jtrack_error());
    */

    //else if (getp(Phase::ID)>4)
    //    svlog.critical("Called jump_to_pureHE_tracks for a star with a Phase more evolved that CHEB "+
    //                   utilities::n2s(getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::jtrack_error());
    /**********************************************/





    //Handle plife and phase
    //Stars in SEVN have a He core already in Terminal Main Sequence and Hshellburning phase, but the pureHE track start from CoreHeburning.
    //If we have to jump to pureHE for these stars, we just assume that it has to jumpt at the beginning of phase 4 (plife0).
    double _Mratio_target = getp(Mass::ID)/getp(MCO::ID);
    double _plife=plife();
    int _phase=getp(Phase::ID);
    if (getp(Phase::ID)<4){
        _plife=utilities::TINY;
        _phase=4;
    }

    //Set Mass limits
    double mlow, mhigh;
    double min_mass = get_svpar_num("min_zams_he"), max_mass = get_svpar_num("max_zams_he");


    //GI: 26/04/2023
    //TODO: Is this the best implementation we can have? It is true that a stripped low mass star never ignite Helium?
    //Following Hurley, a HeWD is formed when a low-mass star (a star with MZAMS<MHEf, mass helium flash) loses the envelope when it is
    //has not yet started the core Helium burning (see setting kw = 10 in hrdiag.f).
    //We are replicating the Hurley implementation forcing the new generated pureHe stars to be turned into a remnant with a MCO=0
    //see supernova::WDformation
    double Z_over_Zsun = std::log10(get_Z()/0.02);
    double MHE_flash = 1.995 + 0.25*Z_over_Zsun + 0.087*Z_over_Zsun*Z_over_Zsun;  //Eq. 2 Hurley https://academic.oup.com/mnras/article/315/3/543/972062
    if (getp(Phase::ID)<Lookup::Phases::CoreHeBurning and get_zams()<=MHE_flash){
        //Notice, the order here is important: We have to 1) Advance the Localtime 2)Let the Phase change 3)turno into remnant
        //Notice, the Phase is updated only in Phase special evolve the remnant type inside turno_into_remnant
        properties[Localtime::ID]->update_from_binary(this, 1.001*(get_tphase()[Remnant] - getp(Localtime::ID)));
        properties[Phase::ID]->special_evolve(this);
        turn_into_remnant();
        return utilities::JUMP;
    }
    //In the following case the mass is below the minimum HE ZAMS in the table, but
    //the total mass is less than 0.4 Msun, so we assume that the stars cannot ignite He burning
    //and they end their life as a HeWD
    //TODO: this optiin should be enabled by a runtime parameter (maybe the value of the mininum mass, if it is 0 the option is disabled)
    else if(getp(Mass::ID)<min_mass and getp(Mass::ID)<0.4){
        properties[Localtime::ID]->update_from_binary(this, 1.001*(get_tphase()[Remnant] - getp(Localtime::ID)));
        properties[Phase::ID]->special_evolve(this);
        turn_into_remnant();
        return utilities::JUMP;
    }
    else if(getp(Mass::ID)<min_mass){
        svlog.critical("Trying to jump to a pureHE, but the Mass ("+utilities::n2s(getp(Mass::ID),__FILE__,__LINE__)
                       +") is lower than the minimum allowed by   the pureHE tracks ("+utilities::n2s(min_mass,__FILE__,__LINE__)+")",
                       __FILE__,__LINE__,sevnstd::jtrack_error());

        /* An example on how to turn such stars in Zombie, it is quite verbose, we should create turn_into_zombie as we have already done
         * for turn_into_remnant (we cannote use such funciton because it automatically calls SN->main())
        properties[Localtime::ID]->update_from_binary(this, 1.001*(get_tphase()[Remnant] - getp(Localtime::ID)));
        properties[Phase::ID]->special_evolve(this);
        set_staremnant(new Zombierem(this,getp(Mass::ID)));
        set_remnant();
        //TODO Set_remnant and set_empty are different because set_empty calls also the set_empy of properteis, we should make a decision and uniformise them
        //Evolve the properties to store the right values (if set_empty already happened)
        if (!isempty){
            remnant();
        }
         */
        return utilities::JUMP;
    }
    else if (getp(Mass::ID)>max_mass){
        svlog.critical("Trying to jump to a pureHE, but the Mass ("+utilities::n2s(getp(Mass::ID),__FILE__,__LINE__)
        +") is larger than the maximum allowed by   the pureHE tracks ("+utilities::n2s(max_mass,__FILE__,__LINE__)+")",
        __FILE__,__LINE__,sevnstd::jtrack_error());
    }




    double best_rel_err=1E30;
    utilities::jump_convergence convergence_M=utilities::NO_JUMP, convergence_core=utilities::NO_JUMP,_convergence;
    double best_zams_M,best_zams_core=utilities::NULL_DOUBLE;
    double Mratio_M, Mratio_core;


    //Find Best zams for Phase 4 or Phase>4 if MC0 is very small.
    // It can happens sometime that at the start of the TerminalHeCore Burning MC0 is small or even 0.
    // Here we use a linear interpolation
    if (_phase==4 or getp(MCO::ID)<1e-3){
        mlow  = std::max(getp(Mass::ID),min_mass);
        mhigh = std::min(2*getp(Mass::ID),max_mass);
        best_zams_M=find_mass_linear(mlow,mhigh,max_mass,min_mass,_plife,_phase,true,best_rel_err,convergence_M,Mass::ID);
    }
    // Find Best zams matching CO (when it is significant( for Phase 5 and 6. Here we use a bisection.
    // GI 20/11, sometimes matching only the CO can cause problems, so we decide to match in any case the total mass.
    //TODO This is an experimental feature, it it works we can use this criterium of Mratio_M also in match_core.
    else{
        mlow  = min_mass;
        mhigh = max_mass;
        best_zams_M=find_mass_bisection(mlow,mhigh,max_mass,min_mass,_plife,_phase,true,best_rel_err,convergence_M,Mass::ID);
        best_zams_core=find_mass_bisection(mlow,mhigh,max_mass,min_mass,_plife,_phase,true,best_rel_err,convergence_core,MCO::ID);

    }

    if (convergence_M==utilities::NO_JUMP and convergence_core==utilities::NO_JUMP)
        svlog.critical("Not jumped to a new track in jump_to_pureHE_tracks",__FILE__,__LINE__,sevnstd::jtrack_error());

    //Jump
    std::string tini_ss = utilities::make_pfile_str(_plife,_phase);
    //best_zams=0.5;
    Star s_best_M(this, 0, best_zams_M, get_Z(), tini_ss,true,get_rseed());
    Mratio_M = s_best_M.getp(Mass::ID)/s_best_M.getp(MCO::ID);
    _convergence = convergence_M;



    //If best_zams_core is not null check which one gives the best match to the Mass-Mcore ration
    //Do not enter also if was not possible to find a tracks (i.e. convergnece !=utilities NO JUMP)
    if (best_zams_core!=utilities::NULL_DOUBLE and convergence_core!=utilities::NO_JUMP){
        Star s_best_core(this, 0, best_zams_core, get_Z(), tini_ss,true,get_rseed());

        Mratio_core = s_best_core.getp(Mass::ID)/s_best_core.getp(MCO::ID);
        if (std::abs(Mratio_core-_Mratio_target) < std::abs(Mratio_M-_Mratio_target) ){
            jump_tracks(&s_best_core);
            _convergence=convergence_core;
        }
        else
            jump_tracks(&s_best_M);
    }
    else
        jump_tracks(&s_best_M);

    //Now we have successfully jumped (notice we cannot be here without a jump given the condition in row 201)
    //Therefore just be sure that MHE and RHE are the same of Mass and Radius.
    //For example, we can be here because Mass and MHE are within a given threshold of tollerance to become pureHE
    //so after the jumping we want to be sure to set Mass=MHE and Radius=RHE, this is not directly done
    //in jump tracks because it checks if the stars is a naked helium but the stars become a naked helium
    //just after this function call.
    copy_property_from(MHE::ID,Mass::ID);
    copy_property_from(RHE::ID,Radius::ID);


    //utilities::wait("Best zams",plife(),s_best.get_zams(),s_best.getp(MHE::ID),s_best.getp(MCO::ID),getp(MHE::ID),getp(MCO::ID),best_rel_err,convergence,__FILE__,__LINE__);
    //if (_convergence==utilities::NO_JUMP)
    //    svlog.critical("Not jumped to a new track in jump_to_pureHE_tracks",__FILE__,__LINE__);



    return _convergence;
}

utilities::jump_convergence Star::jump_to_normal_tracks(){

    if (!ami_on_pureHE_track)
        svlog.critical("Trying to jump to normal tracks for a non pureHE star",__FILE__,__LINE__,sevnstd::jtrack_error());
    if (getp(Mass::ID)<=getp(MHE::ID))
        svlog.critical("Trying to jump to normal tracks for a star without Hydrogen",__FILE__,__LINE__,sevnstd::jtrack_error());

    ami_on_pureHE_track=false;
    force_jump=true;

    utilities::jump_convergence outcome=find_new_track(true); //Use is_a_merger true, so we remove boundaries from MIN MAX ZAMS


    if (outcome==utilities::NO_JUMP)
        svlog.critical("Not jumped to a new track in jump_to_normal_tracks",__FILE__,__LINE__);

    return outcome;
}

utilities::jump_convergence Star::find_track_after_CE_Ebinding(double Ebind, double Min_Mass, double Max_mass, bool pureHE){

    int outcome=utilities::NO_JUMP;
    double best_rel_err=1e30;
    outcome = match_HE_and_binding(Ebind,best_rel_err, Min_Mass, Max_mass,pureHE);

    //We are on a new track Reset the cumulative dM counter
    properties[dMcumul_binary::ID]->reset();

    return outcome;

}

double  Star::find_mass_linear(double mlow, double mhigh, const double MAX_MASS, const double MIN_MASS, const double plife, const size_t Phase, const bool pureHE, double & best_rel_err, utilities::jump_convergence& convergence, const size_t property_id, const double z_look){


    //Check the Property ID
    // Only Mass properties can be used in this function
    if (property_id!=Mass::ID and property_id!=MHE::ID and property_id!=MCO::ID)
        svlog.critical("Property ID "+utilities::n2s(property_id,__FILE__,__LINE__)+" cannot be used as input in function find_mass_linear.",
                __FILE__,__LINE__,sevnstd::notimplemented_error());

    //Check the maximum and minimum mass
    double MAX_ZAMS = pureHE ? get_svpar_num("max_zams_he") : get_svpar_num("max_zams");
    double MIN_ZAMS = pureHE ? get_svpar_num("min_zams_he") : get_svpar_num("min_zams");

    if (MAX_MASS<=MIN_MASS)
        svlog.critical("The MAX_MASS in input ("+utilities::n2s(MAX_MASS,__FILE__,__LINE__)+
        ") is equal or smaller than MIN_MASS in input ("+utilities::n2s(MAX_MASS,__FILE__,__LINE__)+")");
    if (MAX_MASS>MAX_ZAMS)
        svlog.critical("The MAX_MASS in input ("+utilities::n2s(MAX_MASS,__FILE__,__LINE__)+
        ") is larger than the MAX ZAMS ("+utilities::n2s(MAX_ZAMS,__FILE__,__LINE__)+").",__FILE__,__LINE__,sevnstd::sevnerr());
    else if(MIN_MASS<MIN_ZAMS)
        svlog.critical("The MIN_MASS in input ("+utilities::n2s(MIN_MASS,__FILE__,__LINE__)+
                       ") is smaller than the MIN ZAMS ("+utilities::n2s(MIN_ZAMS,__FILE__,__LINE__)+").",__FILE__,__LINE__,sevnstd::sevnerr());
    else if (mhigh<=mlow)
        svlog.critical("The initial  mlow guess ("+utilities::n2s(mlow,__FILE__,__LINE__)+
                       ") is equal or smaller than the initial mhigh guess ("+utilities::n2s(mhigh,__FILE__,__LINE__)+")");
    else if (mlow<MIN_MASS)
        svlog.critical("The initial  mlow guess ("+utilities::n2s(mlow,__FILE__,__LINE__)+
                       ") is smaller than MIN_MASS ("+utilities::n2s(MIN_MASS,__FILE__,__LINE__)+").",__FILE__,__LINE__,sevnstd::sevnerr());
    else if (mhigh>MAX_MASS)
        svlog.critical("The initial mhigh guess ("+utilities::n2s(mhigh,__FILE__,__LINE__)+
                       ") is larger than the MAX MASS ("+utilities::n2s(MAX_MASS,__FILE__,__LINE__)+").",__FILE__,__LINE__,sevnstd::sevnerr());


    double err; //variable to temporarly store the relative error
    double const REL_ERR_TSHOLD = get_svpar_num("jtrack_h_err_rel_max"); //relative error threshold for convergence
    double best_zams=std::min(std::max(get_zams(),MIN_ZAMS),MAX_ZAMS);



    //TODO, Write a constructor that takes directly the percentange of life at a given phase instead to use the string initialiser
    std::string tini_ss = utilities::make_pfile_str(plife,Phase);


    /****** Check the first two guesses *********/
    //Check low
    Star s_mlow(this, 0, mlow, z_look, tini_ss, pureHE,get_rseed());
    err = utilities::rel_difference(getp(property_id), s_mlow.getp(property_id));
    if (err <= REL_ERR_TSHOLD) {
        best_rel_err = err;
        convergence=utilities::JUMP_CONVERGE;
        return mlow;
    } else if (err < best_rel_err) {
        best_zams = s_mlow.get_zams();
        best_rel_err = err;
    }



    //Check high
    Star s_mhigh(this, 1, mhigh, z_look, tini_ss, pureHE,get_rseed());
    err = utilities::rel_difference(getp(property_id), s_mhigh.getp(property_id));

    if (err <= REL_ERR_TSHOLD) {
        best_rel_err = err;
        convergence=utilities::JUMP_CONVERGE;
        return mhigh;
    } else if (err < best_rel_err) {
        best_zams = s_mhigh.get_zams();
        best_rel_err = err;
    }
    /******************************************************/



    /****** Iterative search *********/
    //variable for the iteraive section
    double slope, intercept, new_mass, old_mass = 1e30;
    const int  MAX_ITERATION=get_svpar_num("jtrack_max_iteration"); //maximum number of iterations
    int Niter=0;
    while(Niter<MAX_ITERATION && err>REL_ERR_TSHOLD){
        //for (size_t i = 0; i < MAX_ITERATION; i++) {

        //Find slope and intercept of the line passing through (x1,y1)=(dM_factor_low, s_mlow.Mass) and (x2,y2)=(dM_factor_high, s_mhigh.Mass)
        utilities::find_line(mlow, mhigh, s_mlow.getp(property_id), s_mhigh.getp(property_id), slope,intercept);

        //Estimate new dM_factor and mass
        new_mass = (getp(property_id) - intercept) / slope;
        //If outside the limits use one of the limits.
        new_mass = std::min(std::max(new_mass, MIN_MASS), MAX_MASS);

        //If the value is equal to the last one exit (e.g. this can happen when the values are outside the limits)
        if (new_mass == old_mass) break;
        else old_mass = new_mass;



        //set new dM_factor_low or dM_factor_high
        if (new_mass <= mlow) {
            s_mlow.reset(new_mass, z_look, tini_ss);
            err = utilities::rel_difference(getp(property_id), s_mlow.getp(property_id));
            mlow=new_mass;
            if (err < best_rel_err) {
                best_zams = s_mlow.get_zams();
                best_rel_err = err;
            }
        } else {
            s_mhigh.reset(new_mass, z_look, tini_ss);
            err = utilities::rel_difference(getp(property_id), s_mhigh.getp(property_id));
            mhigh=new_mass;
            if (err < best_rel_err) {
                best_zams = s_mhigh.get_zams();
                best_rel_err = err;
            }
        }



        //Jump tracks
        //Create a star from the best zams
        //Star s_best(this, 0, best_zams, utilities::NULL_DOUBLE, tini_ss.str());
        //jump_tracks(&s_best);
        //Check the convergence
        //if (best_rel_err <= REL_ERR_TSHOLD)
        //    return utilities::JUMP_CONVERGE;
        Niter++;
    }
    /******************************************************/

    /// If the convergence has not been reached jump to the tracks with the best relative error.
    convergence = best_rel_err <= REL_ERR_TSHOLD ? utilities::JUMP_CONVERGE : utilities::JUMP;

    return best_zams;

}

double  Star::find_mass_bisection(double mlow, double mhigh, const double MAX_MASS, const double MIN_MASS, const double plife, const size_t Phase, const bool pureHE, double & best_rel_err, utilities::jump_convergence& convergence, const size_t property_id, const double z_look){

    //Check the Property ID
    // Only Mass properties can be used in this function
    if (property_id!=Mass::ID and property_id!=MHE::ID and property_id!=MCO::ID)
        svlog.critical("Property ID "+utilities::n2s(property_id,__FILE__,__LINE__)+" cannot be used as input in function find_mass_linear.",
                       __FILE__,__LINE__,sevnstd::notimplemented_error());

    //Check the maximum and minimum mass
    const double MAX_ZAMS = pureHE ? get_svpar_num("max_zams_he") : get_svpar_num("max_zams");
    const double MIN_ZAMS = pureHE ? get_svpar_num("min_zams_he") : get_svpar_num("min_zams");


    if (MAX_MASS<=MIN_MASS)
        svlog.critical("The MAX_MASS in input ("+utilities::n2s(MAX_MASS,__FILE__,__LINE__)+
                       ") is equal or smaller than MIN_MASS in input ("+utilities::n2s(MAX_MASS,__FILE__,__LINE__)+")");
    if (MAX_MASS>MAX_ZAMS)
        svlog.critical("The MAX_MASS in input ("+utilities::n2s(MAX_MASS,__FILE__,__LINE__)+
                       ") is larger than the MAX ZAMS ("+utilities::n2s(MAX_ZAMS,__FILE__,__LINE__)+").",__FILE__,__LINE__,sevnstd::sevnerr());
    else if(MIN_MASS<MIN_ZAMS)
        svlog.critical("The MIN_MASS in input ("+utilities::n2s(MIN_MASS,__FILE__,__LINE__)+
                       ") is smaller than the MIN ZAMS ("+utilities::n2s(MIN_ZAMS,__FILE__,__LINE__)+").",__FILE__,__LINE__,sevnstd::sevnerr());
    else if (mhigh<=mlow)
        svlog.critical("The initial  mlow guess ("+utilities::n2s(mlow,__FILE__,__LINE__)+
                       ") is equal or smaller than the initial mhigh guess ("+utilities::n2s(mhigh,__FILE__,__LINE__)+")");
    else if (mlow<MIN_MASS)
        svlog.critical("The initial  mlow guess ("+utilities::n2s(mlow,__FILE__,__LINE__)+
                       ") is smaller than MIN_MASS ("+utilities::n2s(MIN_MASS,__FILE__,__LINE__)+").",__FILE__,__LINE__,sevnstd::sevnerr());
    else if (mhigh>MAX_MASS)
        svlog.critical("The initial mhigh guess ("+utilities::n2s(mhigh,__FILE__,__LINE__)+
                       ") is larger than the MAX MASS ("+utilities::n2s(MAX_MASS,__FILE__,__LINE__)+").",__FILE__,__LINE__,sevnstd::sevnerr());



    //Bisection alg variables
    //Starting time at plife and given phase
    std::string tini_ss = utilities::make_pfile_str(plife,Phase);
    double best_zams=get_zams();
    const int  MAX_ITERATION=get_svpar_num("jtrack_max_iteration"); //maximum number of iterations
    double const REL_ERR_TSHOLD = get_svpar_num("jtrack_h_err_rel_max"); //relative error threshold for convergence
    double Mcore_target=getp(property_id);
    double Mcore_err = 1e30;
    double Mcore_test;
    double test_zmass = std::max(Mcore_target,MIN_ZAMS);



    Star s_test(this, get_ID(), 0.5*(MIN_ZAMS+MAX_ZAMS), z_look, tini_ss, pureHE,get_rseed()); //Aux star used to change track
    //Notice the initial mass  0.5*(MIN_ZAMS+MAX_ZAMS) has not meaning, s_test here is just initiliased, we want just to be sure
    //to choice a mass allowd by the tables


    /*******************************************************************/
    ///Start bisection guess
    //Here we estimate the initial zams interval where to apply the bisection method.
    //Our very first initial guest are between zams_low=max(Mcore,1.001*MIN_ZAMS) and zams_up=0.999*max_zams.
    //Sometimes it is not guarantee that a given plife  a larger zams has a larger core (e.g. MHE can even decrease in massive stars due to drege-up).
    //Therefore we check if both expected core are larger or lower than the target one. If both are larger we increase the value of the smaller zams, otherwise we reduce
    //the value of the larger one

    double Mcore_test_low, Mcore_err_low, Mcore_test_high, Mcore_err_high;
    double mass_step=5;
    unsigned int Niter_iguess=0;
    for (;;){


        //1-First guesses for the interval
        best_rel_err=1e30;
        //Check m_low
        s_test.reset(mlow,z_look,tini_ss);
        Mcore_test_low=s_test.getp(property_id);
        Mcore_err_low = fabs(Mcore_test_low - Mcore_target)/Mcore_target;
        if (Mcore_err_low<best_rel_err){
            best_rel_err = Mcore_err_low;
            best_zams = mlow;
        }



        //Check m_up
        s_test.reset(mhigh,z_look,tini_ss);
        Mcore_test_high = s_test.getp(property_id);
        Mcore_err_high = fabs(Mcore_test_high - Mcore_target)/Mcore_target;
        if (Mcore_err_high<best_rel_err){
            best_rel_err = Mcore_err_low;
            best_zams = mhigh;
        }



        if (Niter_iguess==(unsigned)MAX_ITERATION){
            svlog.warning("Max number of iteration searching for bisection initial interval. Star will not jump",__FILE__,__LINE__);
            convergence = utilities::NO_JUMP;
            return best_zams; //Return best zams anyway, handling of the results outside has to depends also in the convergene value
        }
        if (Mcore_test_high>Mcore_target and Mcore_test_low>Mcore_target){
            //utilities::wait("ML",m_low,__FILE__,__LINE__);
            mlow=std::min(mlow+mass_step,0.999*MAX_ZAMS);
        }
        else if (Mcore_test_high<Mcore_target and Mcore_test_low<Mcore_target){
            //utilities::wait("MH",m_high,__FILE__,__LINE__);
            mhigh=std::max(mhigh-mass_step,1.001*MIN_ZAMS);
        }
        else
            break;

        Niter_iguess++;

    }
    /*******************************************************/


    Mcore_err=std::min(Mcore_err_low,Mcore_err_high);


    //Bisection algorithm
    unsigned int  Niter=0;
    while(Niter<(unsigned)MAX_ITERATION and Mcore_err>REL_ERR_TSHOLD){

        test_zmass=0.5*(mhigh+mlow);
        s_test.reset(test_zmass,z_look,tini_ss);
        Mcore_test = s_test.getp(property_id);
        Mcore_err = fabs(Mcore_test - Mcore_target)/Mcore_target;

        if (Mcore_err<best_rel_err){
            best_rel_err = Mcore_err;
            best_zams = test_zmass;
        }

        //Update extreme points
        if (Mcore_test>Mcore_target)
            mhigh  = test_zmass;
        else
            mlow = test_zmass;

        Niter++;
    }

    if (best_rel_err<=REL_ERR_TSHOLD)
        convergence=utilities::JUMP_CONVERGE;
    else
        convergence=utilities::JUMP;



    return best_zams;
}

utilities::jump_convergence Star::match_M(double &best_rel_err,bool is_merger){

    /**** First of all check that we are not trying to change track mathich the total Mass when we have a core ****/
    if(!aminakedhelium() and getp(MHE::ID)>1e-3)
        svlog.critical("A star (ID: " + utilities::n2s(get_ID(),__FILE__,__LINE__) +
                       ")  is trying to change tracks (LT:" +
                       utilities::n2s(getp(Localtime::ID),__FILE__,__LINE__) + ", WT:" +
                       utilities::n2s(getp(Worldtime::ID),__FILE__,__LINE__) +
                       ") matching the total mass but it has a non negligible core (MHE= "+getps(MHE::ID)+
                       ")",__FILE__,__LINE__,sevnstd::jtrack_error());

    if(aminakedhelium() and getp(MCO::ID)>1e-3)
        svlog.critical("A pureHE star (ID: " + utilities::n2s(get_ID(),__FILE__,__LINE__) +
                       ")  is trying to change tracks (LT:" +
                       utilities::n2s(getp(Localtime::ID),__FILE__,__LINE__) + ", WT:" +
                       utilities::n2s(getp(Worldtime::ID),__FILE__,__LINE__) +
                       ") matching the total mass but it has a non negligible core (MCO= "+getps(MCO::ID)+
                       ")",__FILE__,__LINE__,sevnstd::jtrack_error());

    /**** Now make a consistency check to be sure that we have a negligble core but we are at very evolved phase ****/
    if (!aminakedhelium() and getp(Phase::ID)>3){
        svlog.critical("A star (ID: " + utilities::n2s(get_ID(),__FILE__,__LINE__) +
                      ") that is trying to change tracks (LT:" +
                      utilities::n2s(getp(Localtime::ID),__FILE__,__LINE__) + ", WT:" +
                      utilities::n2s(getp(Worldtime::ID),__FILE__,__LINE__) +
                      ") has  not a core (MHE=MCO=0), but it is in a more evolved"
                      " phase than the TMS (current phase: "+get_phase_str()+").\nNotice that this could be an hint that something is wrong with the "
                                                                            "loaded stellar tables.", __FILE__,__LINE__,sevnstd::jtrack_error());
        return utilities::NO_JUMP;
    }
    else if (aminakedhelium() and getp(Phase::ID)>5){
        svlog.critical("A pureHE star (ID: " + utilities::n2s(get_ID(),__FILE__,__LINE__) +
                      ") that is trying to change tracks (LT:" +
                      utilities::n2s(getp(Localtime::ID),__FILE__,__LINE__) + ", WT:" +
                      utilities::n2s(getp(Worldtime::ID),__FILE__,__LINE__) +
                      ") is in a more evolved phase than CHEB"
                      " (current phase: "+get_phase_str()+").", __FILE__,__LINE__,sevnstd::jtrack_error());
        return utilities::NO_JUMP;
    }
    /*********************************************************/



    /**** Define some variables ****/
    double zlook = get_Z();
    double dM = getp(dMcumul_binary::ID); //dM accumulated for binary processes (can be positive or negative)
    double zams = get_zams(); //Value of the zams before the tracks change
    _UNUSED double err; //variable to temporarly store the relative err

    //Define the Max dM_factor allowed given the loaded tables (dM factor are always positive)
    const double MAX_ZAMS = get_max_zams();
    const double MIN_ZAMS = get_min_zams();
    //FInd the difference with the table boundary of interest: if dM>0 the boundary is MAX_ZAMS otherwise is MIN_ZAMS
    double table_boundary_diff = dM>=0 ? (MAX_ZAMS-zams) : (zams - MIN_ZAMS);
    //If the difference is less than 0.01 just does not jump


    if (table_boundary_diff<1e-2)
        return utilities::NO_JUMP;


    /****** Define and check the first two guesses and the mass boundary*********/
    double minzams,maxzams;
    double mhigh_guess, mlow_guess;
    double dm_jactor_max_mass = dM<0 ? get_svpar_num("jtrack_min_dm_factor") : get_svpar_num("jtrack_max_dm_factor");
    double dm_jactor_min_mass = dM<0 ? get_svpar_num("jtrack_max_dm_factor") : get_svpar_num("jtrack_min_dm_factor");





    //If it is  marger relax boundaries
    minzams       =  is_merger ? MIN_ZAMS : std::max(MIN_ZAMS,get_zams() + dm_jactor_min_mass * dM);
    maxzams       =  is_merger ? MAX_ZAMS : std::min(MAX_ZAMS,get_zams() + dm_jactor_max_mass * dM);


    //Handle limits
    if (dM==0){
        mlow_guess = MIN_ZAMS;
        mhigh_guess = MAX_ZAMS;
    }
    else{
        mlow_guess    =  dM<0 ? std::min(maxzams,std::max(minzams, get_zams() + dM))     : std::min(maxzams,std::max(minzams, get_zams() + 0.2*dM));
        mhigh_guess   =  dM<0 ? std::max(minzams,std::min(maxzams, get_zams() + 0.2*dM)) : std::max(minzams,std::min(maxzams, get_zams() + dM));
    }



    //TODO, in these cases maybe just do not jump?
    //This means mlow_guess=mhigh_guess=maxzams
    if (mlow_guess==maxzams) mlow_guess=minzams;
    //This means mlow_guess=mhigh_guess=minzams
    else if(mhigh_guess==minzams) mhigh_guess=maxzams;




    /********* Find the best zams ************/
    utilities::jump_convergence outcome;
    double best_zams=find_mass_linear(mlow_guess,mhigh_guess,maxzams,minzams,plife(),(int)getp(Phase::ID),aminakedhelium(),best_rel_err,
                             outcome,Mass::ID,zlook);



    ///jump to the tracks with the best relative error.
    //TODO, Write a constructor that takes directly the percentange of life at a given phase instead to use the string initialiser
    std::string tini_ss = utilities::make_pfile_str(plife(),(int)getp(Phase::ID));
    Star s_best(this, 0, best_zams, zlook, tini_ss,aminakedhelium(),get_rseed());
    jump_tracks(&s_best);

    return outcome;
}

utilities::jump_convergence Star::match_HE_and_binding(double Ebind, double & best_rel_err, double Min_Mass, double Max_Mass, bool pureHE){

    //TODO We are not using plife here, so the final star can have very different plife even if the convergence in Ebind has not been reached
    //TODO Since seems difficult to find the conv in Ebind, should be better to find the Mass with a bisection (not in time) and then check Ebind?


    //TODO put parameters in the parameter class
    //PARAMTERS
    double z_look = get_Z();
    double mass_step = 0.1;
    const double MAX_ZAMS = pureHE ? get_svpar_num("max_zams_he") : get_svpar_num("max_zams");
    _UNUSED const double MIN_ZAMS = pureHE ? get_svpar_num("max_zams_he") : get_svpar_num("max_zams");

    double test_zmass = getp(MHE::ID)-mass_step;

    int MAX_iterations = 50;
    double  MHE_target = getp(MHE::ID);
    double HE_err_max = 5.0e-3;
    double Ebind_err_max = 5.0e-2;
    double err_HE = 1e30;
    double err_Ebind = 1e30;
    double Ebind_test;

    double tlow,tup,t_test;
    double MHE_test;
    double M_test;

    double best_zams = 0.0;
    double best_time = 0.0;


    //std::string filename = "CE" + utilities::n2s(z_look,__FILE__,__LINE__)+ "_EbindMethod.txt";
    //std::string filename = "CE" + utilities::n2s(z_look,__FILE__,__LINE__)+ ".txt";


    //CEfile.open(filename);

    //CEfile<<"Z,Mzams,Mt,MHEt,Ebindt,Mmin,Mmax,MHE,Ebind"<<std::endl;
    best_rel_err = 1e30;
    int NmatchHe=0;

    while( (test_zmass<(MAX_ZAMS-mass_step)) and (err_Ebind > Ebind_err_max)){


        test_zmass += mass_step;

        Star s_test(this, get_ID(), test_zmass, z_look, 0.01,pureHE);
        if (getp(MCO::ID)>0){
            tlow = s_test.costart();
            tup  = 0.999*s_test.tphase[Remnant];
        }
        else{
            tlow = s_test.hestart();
            tup  = s_test.tphase[TerminalCoreHeBurning]; //Start of the CO
        }



        for (int i=0; i<MAX_iterations; i++){

            t_test = 0.5*(tlow+tup);
            s_test.reset(utilities::NULL_DOUBLE, z_look, t_test);

            MHE_test = s_test.getp(MHE::ID);
            M_test = s_test.getp(Mass::ID);

            err_HE = fabs(MHE_test - MHE_target)/MHE_target;

            //TODO, We match always the He core, should we consider the CO when it is present?
            if (err_HE< HE_err_max && M_test>=Min_Mass && M_test<=Max_Mass){
                Ebind_test = s_test.getp(Ebind::ID);
                err_Ebind = fabs(Ebind_test-Ebind)/fabs(Ebind);
                NmatchHe++;
                if (err_Ebind<best_rel_err){
                    best_rel_err = err_Ebind;
                    best_zams = test_zmass;
                    best_time = t_test;
                }

                //CEfile<<z_look<<","<<test_zmass<<","<<M_test<<","<<MHE_test<<","<<Ebind_test<<","<<Min_Mass<<","<<Max_Mass<<","
                //      <<MHE_target<<","<<Ebind<<","<<err_Ebind<<","<<Ebind_err_max<<std::endl;
                break;
            }



            if (MHE_test>MHE_target)
                tup  = t_test;
            else
                tlow = t_test;

        }

    }



    //If not match with core He has been found return an error
    if (NmatchHe==0)
        return utilities::NO_JUMP;
        //svlog.critical("No match with He core found for star with: "
        //              "\n MHE=" + utilities::n2s(getp(MHE::ID),__FILE__,__LINE__) +
        //             "\n MCO=" + utilities::n2s(getp(MCO::ID),__FILE__,__LINE__),
        //             __FILE__,__LINE__,sevnstd::jtrack_error());

    /// If the convergence has not been reached jump to the tracks with the best relative error.
    Star s_best(this, 0, best_zams, z_look, best_time, pureHE);
    //utilities::hardwait("AAA",plife(),getp(MHE::ID),getp(MCO::ID),getp(Ebind::ID),Ebind,get_zams(),__FILE__,__LINE__);
    //utilities::hardwait("AAA",s_best.plife(),s_best.getp(MHE::ID),s_best.getp(MCO::ID),s_best.getp(Ebind::ID),s_best.get_zams(),best_rel_err,__FILE__,__LINE__);
    jump_tracks(&s_best);
    //Jump tracks does not change the masses, including the total mass, here we want to change it
    update_from_binary(Mass::ID,s_best.getp(Mass::ID)-getp(Mass::ID));
    //If The star is or becomes a pureHE synch MHE and RHE with Mass and Radius
    if(aminakedhelium()){
        properties[MHE::ID]->copy_V_from(properties[Mass::ID]);
        properties[RHE::ID]->copy_V_from(properties[Radius::ID]);
    }



    //TODO what happens to other properties as for example the Inertia?
    //TODO In SEVN1, total masses larger than M1+M2 were allowed, but then they force an internal change_track wid DM such that the final mass is equal to M1 + M2
    //CEfile<<best_zams<<" "<<best_time<<"  "<<best_rel_err<<std::endl;
    //CEfile.close();

    if (best_rel_err<=Ebind_err_max){
        return utilities::JUMP_CONVERGE;
    }

    return utilities::JUMP;
}

utilities::jump_convergence Star::match_core(double &best_rel_err, _UNUSED bool is_merger){


    //TODO put parameters in the parameter class
    //PARAMETERS
    _UNUSED double z_look = get_Z();
    _UNUSED const int  MAX_ITERATION=get_svpar_num("jtrack_max_iteration");
    const double min_zams = aminakedhelium() ? get_svpar_num("min_zams_he") : get_svpar_num("min_zams");
    const double max_zams = aminakedhelium() ? get_svpar_num("max_zams_he") : get_svpar_num("max_zams");
    double const REL_ERR_TSHOLD = get_svpar_num("jtrack_h_err_rel_max"); //relative error threshold for convergence
    _UNUSED double Mcore_err = 1e30;
    _UNUSED double Mcore_test;

    const double MIN_MASS_CORE = 1e-3; //TODO it should be a parameter
    //MInimum Tshold to try to match the core, very small core means that are very close to the beginning of the phase
    //but here we can found some numer

    /*******************************************************************/
    ///Check if the innercore is too small do not jump
    unsigned int idCore;
    if ( (getp(MCO::ID)==0 and aminakedhelium()) or (getp(MHE::ID)==0 and getp(MCO::ID)==0))
        svlog.critical("You are trying to change track considering a core, but this star has not a core yet",__FILE__,__LINE__,sevnstd::jtrack_error());
    else if (getp(MCO::ID)>MIN_MASS_CORE)
        idCore=MCO::ID;
    else if (getp(MHE::ID)>MIN_MASS_CORE and !aminakedhelium())
        idCore=MHE::ID;
    else{
        svlog.warning("You are trying to change track considering a core, but this star has very tiny core (MCO=" + utilities::n2s(getp(MCO::ID),__FILE__,__LINE__) +
        ", MHE=" +  utilities::n2s(getp(MHE::ID),__FILE__,__LINE__)  + ", MIN_MASS_CORE= " + utilities::n2s(MIN_MASS_CORE,__FILE__,__LINE__) +")  so it does not jump",
        __FILE__,__LINE__);
        return utilities::NO_JUMP;
    }
    /*******************************************************/


    /********* Find the best zams ************/
    //TODO We should use min_zams and max_zams without multiplication
    //Start with the innermost core
    double m_low=1.001*min_zams;
    double m_high=0.999*max_zams;
    utilities::jump_convergence outcome;
    double best_zams;


    best_zams = find_mass_bisection(m_low,m_high,max_zams,min_zams,plife(),getp(Phase::ID),aminakedhelium(),best_rel_err,outcome,idCore,z_look);



    //If not jumped and we tried with MCO, try with MHE if possible
    if (outcome==utilities::NO_JUMP and idCore==MCO::ID and !aminakedhelium() and getp(MHE::ID)>MIN_MASS_CORE){
        best_zams = find_mass_bisection(m_low,m_high,max_zams,min_zams,plife(),getp(Phase::ID),aminakedhelium(),best_rel_err,outcome,MHE::ID, z_look);
    }
    //Very last change only if is_merger is true
    if (outcome==utilities::NO_JUMP and is_merger){
        //best_zams = find_mass_bisection(m_low,m_high,max_zams,min_zams,plife(),getp(Phase::ID),aminakedhelium(),best_rel_err,outcome,Mass::ID);
        //Notice for the total mass we use find_mass_linear since it will always found at list a track even if the convergene is not reached
        //while the bisection can get stuck to the initial choice of the interval that bound the mass.
        best_zams = find_mass_linear(m_low,m_high,max_zams,min_zams,plife(),getp(Phase::ID),aminakedhelium(),best_rel_err,outcome,Mass::ID, z_look);
    }


    //If stars does not jump just exit
    if (outcome==utilities::NO_JUMP)
        return utilities::NO_JUMP;




    ///jump to the tracks with the best relative error.
    //TODO, Write a constructor that takes directly the percentange of life at a given phase instead to use the string initialiser
    std::string tini_ss = utilities::make_pfile_str(plife(),(int)getp(Phase::ID));
    Star s_best(this, 0, best_zams, z_look, tini_ss,aminakedhelium(),get_rseed());
    jump_tracks(&s_best);




    if (best_rel_err<=REL_ERR_TSHOLD){
        return utilities::JUMP_CONVERGE;
    }

    return utilities::JUMP;

}

utilities::jump_convergence Star::jump_tracks(Star *s){

    ///Set the new parameters
    ami_on_pureHE_track=s->ami_on_pureHE_track;

    ///First of all change the properties (this is the first thing to do, before s is changed by the jump)
    //Handle the change of the properties. E.g. Mass is the same, Luminosity and Radius are equal to the new stars
    //Inertia is the weighted mean.
    for (auto & prop : properties)
        prop->changed_track(this,s);


    ///Init stars
    set_mzams(s->get_zams(),__FILE__,__LINE__);
    set_Z(s->get_Z(),__FILE__,__LINE__);
    init_on_lookup(); //Handle the pointer to the right table
    tracktimes(); //Find the right position in the tables
    needsinit=true; //Signal that the next lookup should reinitialised, it will be reset to false inside lookup
    update_jtrack_counter();
    //Set new rzams
    rzams=s->get_rzams(); //Set new rzams

    return EXIT_SUCCESS;
}

void Star::default_initialiser(bool isareset, bool pureHE) {



    /*** Initialise some class members ***/
    initialise_as_remnant=Lookup::Remnants::NotARemnant;
    needsinit = true; ///needsinit=true is needed in lookup position to reinisialise the binary search of the position in the lookup tables
    break_at_remnant = false;
    print_all_steps = false;
    //TODO I should eliminate print per phase in favour or "print_every_n".. like in gnuplot
    print_per_phase = false;
    print_only_end = false;
    print_events = false;
    repeatstep = false;
    iskicked = false;
    isempty = isconaked = once_conaked = isremnant = changedphase = false;
    initialised_on_zams = false;
    force_jump = false;
    ami_on_pureHE_track = pureHE;
    ami_following_QHE = false;
    evolution_step_completed = false;
    for (auto& v : vkick){
        v=std::nan("");
    }


    /*** Init mzamz, Z, sntype, Spin (init_1); tini, tf, dtout, LocalTime, Timestep (init_2); interpolating tracks (init_on_lookup) ***/
    init(init_params, isareset);

    /*** SN istance ***/
    if (!isareset) {
        SN = supernova::Instance(sntype, this);
        //utilities::hardwait("Called SN instance", get_name(), sntype, __FILE__, __LINE__);
        if (SN == nullptr)
            svlog.critical("Unknown SN explosion model for a star: [" + sntype + "]", __FILE__, __LINE__);
    }


    /*** Evolve the star to the initial time (without any print) or create a remnant depending on the
     * starting condition. Notice that if the code is starting directly with a remnant, mzams is the initial mass of the remnant
     * not of the progenitor.
     * ***/
    if (initialise_as_remnant==Lookup::Remnants::NotARemnant){
        evolve_to_tini();
    }
    else{
        evolve_to_tini_as_remnant();
    }


    //Set Z0 and mzams0
    Z0=get_Z();
    mzams0=get_zams();


    needsinit = false;

}
void Star::set_rzams() {
    if (initialise_as_remnant!=NotARemnant or getp(Phase::ID)!=NotARemnant){
        rzams = getp(Radius::ID);
    }
    else if ( initialised_on_zams and getp(Worldtime::ID)==0.){
        rzams = getp(Radius::ID);
    }
    else if (aminakedhelium()){
        Star aux_star(this,0,get_zams(),get_Z(),"cheb",aminakedhelium(),get_rseed());
        rzams=aux_star.getp(Radius::ID);
    }
    else{
        Star aux_star(this,0,get_zams(),get_Z(),"zams",aminakedhelium(),get_rseed());
        rzams=aux_star.getp(Radius::ID);
    }

}

void Star::reset_staremnant(){
    if (staremnant!= nullptr){
        delete staremnant;
        staremnant= nullptr;
    }
}

void Star::evolve_to_tini_as_remnant(){

    svlog.debug("Initialising star as compact object",__FILE__,__LINE__);

    ///Just propose a very large timestep
    properties[Timestep::ID]->resynch(1e6,false);

    //Evolve Localtime and Phase, but not Worldtime since this is just a fake pre-evolution phase to bring
    //the stars to the correct initial position as compact
    properties[Localtime::ID]->special_evolve(this);
    properties[Phase::ID]->special_evolve(this);
    //Initialise remnant properties
    SN->initialise_remnant(this,get_zams(),initialise_as_remnant);

    //initialize also the first value for the next output, very important otherwise the check isoutputime, will be always
    //false and the nextoutpu special evolve will be never called
    properties[NextOutput::ID]->init(get_dtout());

    //Upate all the properties if not empty
    if (initialise_as_remnant!=Lookup::Remnants::Empty)
        remnant();

}

utilities::evolution Star::evolve() {

    //First of all check if the star has to become a pureHe star because it satisfies the WR condition
    //This will be possible only if the option turn_WR_to_pureHe is enabled
    //Notice the aminakedco and amirmenant are needed to avoid awkward cases where a remnant or a naked CO could try to jump to pureHE
    //ATM this should not be possibile for remnants since MHE=0, but it could be for nakedCO where MHE=MASS=MCO
    if (get_svpar_bool("turn_WR_to_pureHe") and !aminakedhelium() and !aminakedco() and !amiremnant()  and amiWR()){
        //TODO Here we have two options: use the find_new_track that can automatically handle also the WR transformation or
        //be super specific and handle here the jump_to_pureHe, we are using the second for now
        //find_new_track();
        auto outcome = jump_to_pureHE_tracks();
        auto log_mess = log_mess_HEnaked(outcome);
        this->print_to_log(log_mess);
    }


    //Reset signal
    evolution_step_completed=false;


    int _count=0;
    for(;;){

        //TODO Transform this to a Warning and let the evolution continue (stopping the repetitions).
        if (_count>get_svpar_num("ev_max_repetitions"))
            svlog.critical("The number of repetitions in sse (Star ID:" +
                                   utilities::n2s(int(get_ID()),__FILE__,__LINE__) +") reaches the maximum allowed ("+
                    utilities::n2s(get_svpar_num("ev_max_repetitions"),__FILE__,__LINE__)+"). You can bypass this error"
                                                                                          "increasing the parameter ev_max_repetitions. However"
                                                                                          " a large number of repetitions is usually an hint that something is broken.",
                    __FILE__,__LINE__,sevnstd::sse_error());


        svlog.debug("Inside star "+utilities::n2s(ID,__FILE__,__LINE__));
        double ffff = properties[Localtime::ID]->get();
        svlog.debug("Localtime "+utilities::n2s(ffff,__FILE__,__LINE__));
        svlog.debug("Timestep "+utilities::n2s(properties[Timestep::ID]->get(),__FILE__,__LINE__));
        svlog.debug("Timestep 0 "+utilities::n2s(properties[Timestep::ID]->get_0(),__FILE__,__LINE__));
        svlog.debug("Nextout "+utilities::n2s(properties[NextOutput::ID]->get(),__FILE__,__LINE__));
        //utilities::wait();


        properties[Localtime::ID]->special_evolve(this);
        properties[Phase::ID]->special_evolve(this);
        properties[Worldtime::ID]->special_evolve(this);

        //if (get_ID()==0 and getp(Worldtime::ID)!=0 and getp(Worldtime::ID)<=40)
        //    std::cout<<getp(Worldtime::ID)<<" "<<getp_0(Localtime::ID)<<" "<<getp(Localtime::ID)<<" "<<getp(Timestep::ID)<<" "<<tphase[Nphases-1]<<" "<<getp(Phase::ID)<<__FILE__<<" "<<__LINE__<<std::endl;


        svlog.debug("Inside star "+utilities::n2s(ID,__FILE__,__LINE__));
        ffff = properties[Localtime::ID]->get();
        svlog.debug("Localtime after ev"+utilities::n2s(ffff,__FILE__,__LINE__));
        svlog.debug("Timestep "+utilities::n2s(properties[Timestep::ID]->get(),__FILE__,__LINE__));
        //utilities::wait();

        ///If empty, force V0=V and return
        if (isempty){
            //TODO Really needed to set everytime?
            for (auto &prop : properties) {
                prop->evolve_empty(this); //V0=V=nan
            }
            //Signal
            evolution_step_done();
            return utilities::SINGLE_STEP_EVOLUTION;
        }


        ///If time to explode and the star is not already a remnant turno into remnant and return.
        if (properties[Localtime::ID]->get() >= tphase[Remnant] and !isremnant) {
            turn_into_remnant();
            //Signal
            evolution_step_done();
            return utilities::SINGLE_STEP_EVOLUTION;
        }

        //if (get_ID()==1)
        //    utilities::hardwait("Check",isremnant,__FILE__,__LINE__);


        ///If it is remnant of conaked just evolve accordingly the properties
        if (isremnant){
            ///Here the star is a remnant
            /// so just evolve the time and put V0=V
            for (auto &prop : properties) {
                prop->evolve_remnant(this); //V0=V with some exception
            }
            //Signal
            evolution_step_done();
            return utilities::SINGLE_STEP_EVOLUTION;
        }
        else if (isconaked){
            ///Here the star is a remnant
            /// so just evolve the time and put V0=V
            for (auto &prop : properties) {
                prop->evolve_nakedco(this); //V0=V, with some exceptions
            }
            //Signal
            evolution_step_done();
            return utilities::SINGLE_STEP_EVOLUTION;
        }



        tracktimes();
        lookuppositions();



        for (auto &prop : properties) {

            prop->evolve(this);

        }

#ifdef DEBUG
        io->print_timelog(this);
            svlog.debug("COUNT "+std::to_string(_count));
#endif

        _count++;

        if (repeatstep)
            restore();
        else
            break;

    }

    //Signal
    evolution_step_done();

    if (_count>1) return utilities::REPEATED_EVOLUTION;
    else return utilities::SINGLE_STEP_EVOLUTION;
}

void Star::turn_into_remnant() {
    set_remnant();
    //Reset conaked and henaked
    //if (aminakedco()){
    //    isconaked = false; //Now this star is a remnant, reset nakedco
    //    once_conaked = true; //If the star was a nakedCO before to become a remnant, take memory of this (useful in case of repetition)
    //}
    //isremnant = true;
    //Let SN explodes (main will set isremnant to true)
    SN->main(this);
    //TODO Set_remnant and set_empty are different because set_empty calls also the set_empy of properteis, we should make a decision and uniformise them
    //Evolve the properties to store the right values (if set_empty already happened)
    if (!isempty){
        remnant();
    }
}

void Star::explode_as_SNI() {
    tobe_SNIA=false;
    auto wlog = log_message_SNIa();
    SN->explosion_SNI(this); //Star is set to empty
    print_to_log(wlog);
}

void Star::tobe_exploded_as_SNI(){
    tobe_SNIA=true;
}

bool Star::is_exploding_as_SNI() const{
    return tobe_SNIA;
}



void Star::explode_as_SNI(Binstar *b) {
    SN->explosion_SNI(this,b); //Star is set to empty
}


void Star::crem_transition_to_NS() {
    reset_staremnant(); //Reset staremant pointer
    set_staremnant(new NSCCrem(this,getp(Mass::ID))); //Set new remnant
    remnant(); //Call set_remnant on all the properties so that the the V are set to the correct new value
}

void Star::crem_transition_to_BH() {
    reset_staremnant(); //Reset staremant pointer
    set_staremnant(new BHrem(this,getp(Mass::ID))); //Set new remnant
    remnant(); //Call set_remnant on all the properties so that the the V are set to the correct new value
}

Material Star::whatamidonating() const  {


    //Not a remnant
    if (getp(Phase::ID)!=7){
        if(aminakedco())
            return Material::CO;
        else if(aminakedhelium())
            return Material::He;
        else
            return Material::H;
    } //Remnant
    else {

        double remtype = getp(RemnantType::ID);
        if (remtype==Remnants::HeWD)
            return Material::He;
        else if (remtype==Remnants::COWD)
            return Material::CO;
        else if (remtype==Remnants::ONeWD)
            return Material::ONe;
        else if (remtype==Remnants::NS_CCSN || remtype==Remnants::NS_ECSN)
            return Material::Neutron;
        else if (remtype==Remnants::BH)
            return Material::BHmaterial;
        else if (remtype==Remnants::Zombie)
            return Material::Unknownmaterial;
    }

    return Material::Unknownmaterial; //TODO this needs to be changed... in this case the function should return an ``error material''


}

Material Star::whatamidonating_0() const  {


    //Not a remnant
    if (getp_0(Phase::ID)!=7){
        if(aminakedco())
            return Material::CO;
        else if(aminakedhelium())
            return Material::He;
        else
            return Material::H;
    } //Remnant
    else {
        double remtype = getp_0(RemnantType::ID);
        if (remtype==Remnants::HeWD)
            return Material::He;
        else if (remtype==Remnants::COWD)
            return Material::CO;
        else if (remtype==Remnants::ONeWD)
            return Material::ONe;
        else if (remtype==Remnants::NS_CCSN || remtype==Remnants::NS_ECSN)
            return Material::Neutron;
        else if (remtype==Remnants::BH)
            return Material::BHmaterial;
        else if (remtype==Remnants::Zombie)
            return Material::Unknownmaterial;
    }

    return Material::Unknownmaterial; //TODO this needs to ba changed to Material::unknown

}

int Star::_main_get_bse_phase(bool old) const {

    const auto& getP = [this,&old](size_t _ID){
        if (old)
            return this->getp_0(_ID);
        else
            return this->getp(_ID);
    };


    int sevn_phase =  int(getP(Phase::ID));

    //First check remnant
    if (sevn_phase == Phases::Remnant){

        if(getP(RemnantType::ID) == Remnants::HeWD) return 10; // He White dwarf
        else if (getP(RemnantType::ID) == Remnants::COWD) return 11; //CO White dwarf
        else if (getP(RemnantType::ID) == Remnants::ONeWD) return 12; //Oxygen-Neon White dwarf
        else if (getP(RemnantType::ID) == Remnants::NS_CCSN or getP(RemnantType::ID) == Remnants::NS_ECSN) return 13; //Neutron star
        else if (getP(RemnantType::ID) == Remnants::BH) return 14; //Black hole
        else if (getP(RemnantType::ID) == Remnants::Empty) return 15; //Massless remnant
        else if (getP(RemnantType::ID) == Remnants::Zombie) return -2; //Zombie stars are special stars that SEVN cannot handle
    }
    //If nakedco treat it as a WD (CO or ONe) or still a NakedHelium
    else if ( (aminakedco() or once_conaked) and getp(Mass::ID)<get_svpar_num("sn_Mchandra")){
        //Mzams 1.6 limit from Section 6 Hurley+00
        if (get_zams()<1.6) return  11; //HeWD
        else return  12; //ONeWD
    }
    else if ( (aminakedco() or once_conaked)){
        return 7; //Consider it a Naked Helium without envelope //TODO (it this right?)
    }
    else if (sevn_phase == Phases::PreMainSequence)
        return -1;
        //GI 22/11/2020 In order to avoid to have phase 2 or 3 star with a very small He core we consider still MS if the HE content is small
        //TODO 1e-3 as HE threshold is ok?
    else if(sevn_phase == Phases::MainSequence or (sevn_phase == Phases::TerminalMainSequence and getp(MHE::ID)<1e-5)){
        ///Separate from Deep of fully convective low mass Ms and high mass MS
        //Use Convective table if present
        if(get_svpar_bool("tabuse_envconv")){
            //TODO The threshold used to define a convective envelope here is hardcoded (better use the radius?), we should test it
            return getP(Qconv::ID)>0.8 ? 0 : 1; //Deep and alrge convective envelope return 0 otherwise 1
        }
            //Use Mass threshold in Hurley+02
        else{
            return get_zams()<0.7 ? 0 : 1; //Low Mass Ms return 0 otherwise 1 (Hurley+02)
        }

    }
    else if(sevn_phase == Phases::TerminalMainSequence || sevn_phase == Phases::ShellHBurning){ //helium core

        if (get_svpar_bool("use_thg_hurley")){
            Tbgb tbgb(get_Z());
            Tms tms(get_Z());
            double Thg = tbgb(get_zams()) - tms(get_zams());
            if (getp(Localtime::ID) - *(get_tphase()+2) <=Thg){
                return 2;
            } else {
                return 3;
            }
        }

        double M_convective_envelope = getP(Qconv::ID)* getP(Mass::ID);
        double M_H_envelope = getP(Mass::ID) - getP(MHE::ID);

        //In Hurley00 (Sec. 5.1) the base of the GB is defined as the point where
        //Mcnv=0.4Menvelope for M<=Mheliumf_flash
        //Mcnv=0.33333333Menvelope for M>Mheliumf_flash
        //Mhelium_flash=1.995 + 0:25 * log(Z/0.02) +  0.087 log(Z/0.02)^2
        //TODO the following should be made consistent with the tracks of the look-up tables (we KNOW which are the stars that do not ignite He!!)
        //double _Z = std::log10(get_Z()/0.02);
        //double MHE_flash = 1.995 + 0.25*_Z + 0.087*_Z*_Z;  //Eq. 2 Hurley https://academic.oup.com/mnras/article/315/3/543/972062
        //double conv_fraction= get_zams()<MHE_flash ? 0.4 : 0.33333333333; //Hurley00 SEc. 5.1
        //We should use the code above to be consistent, but since most of the star in SEVN have M>Mhe_flash
        //we avoid to spend computational resources in estimating MHE_flash and use 0.333333333 as conv_fraction
        double conv_fraction=0.3333333333333;

        if(M_convective_envelope < conv_fraction*M_H_envelope) return 2; //Hertzsprung gap: small convective envelopes //formulas for M_CE in Sec.5.1 of https://arxiv.org/pdf/astro-ph/0001295.pdf
        else return 3; //First Giant Branch: large convective envelopes
    }
        //GI 22/11/2020 In order to avoid to have phase 5 or 6 or 8 star with a very small CO core we consider still CHEB if the CO content is small
        // TODO 1e-3 as HE threshold is ok?
    else if (sevn_phase == Phases::CoreHeBurning or (sevn_phase == Phases::TerminalCoreHeBurning and getP(MCO::ID) < 1e-5)) {
        if(!aminakedhelium())
            return 4; //Core Helium burning
        else
            return 7; // Naked Helium (MS)
    }
    else if (!aminakedhelium() and (sevn_phase == Phases::TerminalCoreHeBurning or sevn_phase == Phases::ShellHeBurning)) {

        double M_convective_envelope = getP(Qconv::ID)* getP(Mass::ID);
        double M_H_envelope = getP(Mass::ID) - getP(MHE::ID);
        double conv_fraction=0.4;

        if(M_convective_envelope < conv_fraction*M_H_envelope){
                return 4; //Core Helium burning
        }
        else {
            return 5; //Early AGB
        }
    }
    else if(aminakedhelium() and sevn_phase == Phases::TerminalCoreHeBurning){
        return 7;
    }
    else if(aminakedhelium() and sevn_phase == Phases::ShellHeBurning){
        //We assume all the pureHE are radiative.
        //We check with parsec tracks that this is totally consistent:
        //all the pureHE are completely radiative or have a very light convective envelope (Qconv<1e-3) toward the end
        //of their life if they are metal poor. Therefore, we use the bse phase 8 (pureHE HG assuming that GB is for stars with a significant convective layer)
        return 8;
    }
    else if (std::isnan(getP(Mass::ID))){
        return 15; //Massless remnant
    }
    else
        svlog.critical("Cannot set BSE phase",__FILE__,__LINE__,sevnstd::sse_error());

    return -1.0; //return negative number if cannot set the BSE phase

}

int Star::get_bse_phase() const {
    return _main_get_bse_phase();
}

int Star::get_bse_phase_0() const {
    return _main_get_bse_phase(true);
}

void Star::set_staremnant(Staremnant* remnant){

    //if (get_ID()==0) std::cout<<"SR "<<getp(Phase::ID)<<" "<<getp(Worldtime::ID)<<std::endl;
    if (staremnant != nullptr){
        svlog.error("Staremnant is already set "+utilities::n2s(getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__,false);
        //Discard old remnant
        delete staremnant;
        staremnant = nullptr;
    }

    staremnant = remnant;
    remnant = nullptr;
}

double Star::get_staremnant(size_t ID)  {
    return staremnant->get(this,ID);
}


std::string Star::log_message_SNIa() {
    std::string w = utilities::log_print("SNIa",this,
                                         getp(Mass::ID),getp(Radius::ID),getp(Luminosity::ID),
                                         getp(RemnantType::ID));
    return w;
}




double Star::inspect_param_dt(std::string p){


    ///Check if the global variable is set
    if (get_svpar_str("dtout")!="list"){
        p=get_svpar_str("dtout");
    }
    else if(get_svpar_str("dtout")=="list" and p==utilities::PLACEHOLDER){
        svlog.critical("Using placeholder option for dtout but the parameter dtout is"
                       " set to list",__FILE__,__LINE__,sevnstd::params_error());
    }


    if(utilities::string_is_number<double>(p)){
        auto dt = utilities::s2n<double>(p, __FILE__, __LINE__);
        dtout_generator= utilities::ListGenerator::make_unique(dt);
        dtout_type = Lookup::DTout_type::_DTVAL;
        return dt;
    }
        //TODO USe regex instead
    else if (p.front()=='{' and p.back()=='}'){
        std::vector<std::string> svalues = utilities::split(p.substr(1,p.size()-1),',');
        std::vector<double> dvalues(svalues.size());
        std::transform(svalues.begin(),svalues.end(),dvalues.begin(),[](std::string& val){return utilities::s2n<double>(val,__FILE__,__LINE__);});
        dtout_generator= utilities::ListGenerator::make_unique(dvalues);
        dtout_type = Lookup::DTout_type::_DTLIST;
    }
    else if (p.find(':')!=std::string::npos){
        std::vector<std::string> svalues = utilities::split(p,':');
        if (svalues.size()!=3){
            svlog.critical("The x:x:x dtout option need exactly 3 parameters.",__FILE__,__LINE__,sevnstd::params_error());
        }

        auto tmin=utilities::s2n<double>(svalues[0],__FILE__,__LINE__);
        auto tmax=utilities::s2n<double>(svalues[1],__FILE__,__LINE__);
        auto tstep=utilities::s2n<double>(svalues[2],__FILE__,__LINE__);
        dtout_generator=utilities::ListGenerator::make_unique(tstep,tmax,tmin);
        dtout_type = Lookup::DTout_type::_DTGENINTERVAL;
    }
    else if (p=="all"){
        print_all_steps = true;
        dtout_type = Lookup::DTout_type::_DTALL;
    }
    else if(p=="end"){
        print_only_end = true;
        dtout_type = Lookup::DTout_type::_DTEND;
    }
    else if(p=="events" or p=="eventsrlo"){
        print_events = true;
        dtout_type = p=="eventsrlo" ? Lookup::DTout_type::_DTEVNTRLO : Lookup::DTout_type::_DTEVNT;
    }
    else
        svlog.critical("Option "+p+" not implemented for dtout",__FILE__,__LINE__,sevnstd::notimplemented_error());

    return 1E30;
}


void Star::default_destructor(){
    ///Here we have to delete all the pointer to stuff allocated on the heap

    //Delete properties (properties are allocated on the heap, see property.h)
    if (!properties.empty()){
        for (auto& prop : properties){
            delete prop;
            prop = nullptr;
        }
    }

    //Delete SN
    if (SN!= nullptr){
        //Delete SN
        delete SN;
        SN=nullptr;
    }

    //Delete staremnant
    if (staremnant!= nullptr){
        delete staremnant;
        staremnant = nullptr;
    }
}

Star::~Star(){
    ///Here we have to delete all the pointer to stuff allocated on the heap
    default_destructor();
}






