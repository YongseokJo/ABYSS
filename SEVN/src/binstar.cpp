//
// Created by mario on 04/12/18.
//

#include <binstar.h>
#include <supernova.h>
#include <Collision.h>
#include <Hardening.h>
#include <Circularisation.h>
#include <Windaccretion.h>
#include <remnant.h>

#define FOR_BPROC for(size_t _i = 0; _i < process.size(); _i++)
#define FOR_BPROP for(size_t _i = 0; _i < property.size(); _i++)


/**********************************
 *  Evolve functions
 **********************************/

utilities::bse_evolution Binstar::evolve_binary() {

    ///1-Reset some flags
    force_tiny_dt=false; //Reset, so that a new force_tiny_dt can be called if needed.
    repeatstep=false; //Reset repeastep, in case it will be set to true again when the BTimestep is evolved.

    ///2-Special evolve
    property[BWorldtime::ID]->special_evolve(this);


    ///3-Check if the bynary still esist
    //It this is not a binary we do not have to evolve (except the time special evolve)
    if (broken or onesurvived or empty){
        property[BTimestep::ID]->special_evolve(this); //handle Timestep step T_0=T, T=1e30;
        return utilities::BIN_EV_NOT_DONE;
    }


    ///4-Check if some stellar process has destroied the systems
    //If some stellar process, e.g. a PISN SN, destroy a star set broken and return
    if (star[0]->isempty and star[1]->isempty){
        set_broken(); set_empty();
        return utilities::BIN_EV_SETBROKEN;
    }
    else if (star[0]->isempty or star[1]->isempty){
        set_broken(); set_onesurvived();
        return utilities::BIN_EV_SETBROKEN;
    }

    ///5-Reset variations   due to processes (all to zero... VB and VS)
    for (auto &proc : process) proc->restore();

    ///6-Outside the parallel update, we have here  a sequential update due to the special evolve function of the processes
    //Update 19/07/21 GI: Before this sequential checks directly returned only if the system was broken or destroyed.
    //This is in general not a problem for comenv, mix and is_swallowed since the timestep of the run is tiny and the binary process
    //cannot really do anything significant.  However, there could be some critical cases. For example, a star becomes a WD, so it enters
    //inside  the just_exploded check (the name is maybe misleading since it is true also for a WD formation), but at the same time during the
    // binary evolution a collision with a mix is triggered. Therefore, the evolution is repeated but even with a small timestep the star becomes
    //again a WD. So in the same step we have to handle a remnant transformazion and a mix. Waiting for a better idea, it is safer to
    //jut assume that after one of these special evolution step the binary evolution is skipped.
    if (star[0]->just_exploded() or star[1]->just_exploded() ){
        //Before to apply the kick to the binary, check if we have to apply the bavera Xspin correction
        //It uses the value of the Period (that has not been changed yet and the mass of the remnant progenitor)
        check_and_set_bavera_xspin();
        process[SNKicks::ID]->special_evolve(this);
        //Kick is an impulsive change of properties, so for safety just reset all the flags
        reset_evolution_flags();
        reset_all_processes_alarms();

        if (check_and_set_broken())  return utilities::BIN_EV_SETBROKEN; //THis is not a binary anymore, return.
    }
    //TODO When the circularisation is triggered, a special timestep with a timestep of 1e-15 is done to take into account only the circularisation
    //so we have a timestep in the ouput (actually the first one in which we set RLO start) in which the mass transfer through RLO is 0 because we consider only the circularisation
    else if (processes_alarm(Circularisation::ID)){
        process[Circularisation::ID]->special_evolve(this); //Circularise the orbit
        process[Circularisation::ID]->special_evolution_alarm_switch_off(); //Reset alarm
        if (check_and_set_broken()) svlog.critical("System broken after the Circularisation special evolve, this is not expected",
                                                   __FILE__,__LINE__, sevnstd::sanity_error(""));
    }
    else if (comenv){
        process[CommonEnvelope::ID]->special_evolve(this); //CE, comenv is set to false inside
        if (check_and_set_broken()) return utilities::BIN_EV_SETBROKEN; //THis is not a binary anymore, return.
    }
    else if (mix){
        process[Mix::ID]->special_evolve(this); //mix
        if (check_and_set_broken()) return utilities::BIN_EV_SETBROKEN; //THis is not a binary anymore, return.
    }
    else if (is_swallowed[0] or is_swallowed[1]){
        process[RocheLobe::ID]->special_evolve(this); //dynamic_swallowing
        if (check_and_set_broken()) return utilities::BIN_EV_SETBROKEN; //THis is not a binary anymore, return.
    }
    ///8-Main evolve
    else{

        ///8a-reset special flags (mix, comenv, is_swallowed)
        reset_evolution_flags();

        //TODO we are switching to handling the special evolution through processes alarm
        //instead of using evolution flags, at certain point we have to include also mix, comenv and is swallowed in the new framework
        reset_all_processes_alarms();

        //TODO when all the property and processes are in place use just a for cycle
        //TODO The evolve should return a flag if a special condition has been raised (e.g. a merger in RLO) so that we can check it and avoid to consider all the process in these specific cases
        ///8b-Evolve all the Processes
        process[GWrad::ID]->evolve(this);
        process[RocheLobe::ID]->evolve(this);
        process[Circularisation::ID]->evolve(this); //Notice: Circularisation after RocheLobe
        process[Windaccretion::ID]->evolve(this);
        process[Tides::ID]->evolve(this);
        process[Hardening::ID]->evolve(this);


        //The Collision at periastron is checked and taken into account only if the star is not in a RLO, or
        // if the runtime option rlo_enable_collision is set to true
        if (!process[RocheLobe::ID]->is_process_ongoing() or get_svpar_bool("rlo_enable_collision")){
            process[Kollision::ID]->evolve(this);
        }



        ///8c-Check if a special state has been flagged or is the system ahs been broken
        if (mix or comenv or is_swallowed[0] or is_swallowed[1] or processes_alarm(Circularisation::ID)){
            //If the flags are set the BTimestep evolve will call a repeat with a small timestep
            //TODO Here we have to evolve all the properties, otherwisw the BTimestep step call repeat and
            //the Binary property will return to the tstep-2 because we have not evolved them yet.
            for (auto& prop : property)
                prop->evolve(this);
            //property[BTimestep::ID]->evolve(this);
            return utilities::BIN_EV_DONE;
        }
        //Check if the system is not broken, some process can broken it, e.g. HeWd-HeWD RLO
        if (check_and_set_broken()) return utilities::BIN_EV_SETBROKEN; //THis is not a binary anymore, return.



        ///8d-Update the Single stellar properties due to the processes. This is a parallel update
        for (size_t k = 0;k < Property::all.size(); k++) { //the k-th property has varied because the i-th process has made that vary

            for (auto s : star){
                if (s->isempty){
                    utilities::wait("Is empty",__FILE__,__LINE__);
                    continue;
                }
                if (k!=dMcumul_RLO::ID){
                    s->update_from_binary(k, process[Windaccretion::ID]->get_var(s->get_ID(), k),this);
                    s->update_from_binary(k, process[Tides::ID]->get_var(s->get_ID(), k),this);
                    s->update_from_binary(k, process[GWrad::ID]->get_var(s->get_ID(), k),this);
                    s->update_from_binary(k, process[Hardening::ID]->get_var(s->get_ID(),k),this);
                }
                //dMcumul_RLO should uptade only the RL DM
                s->update_from_binary(k, process[RocheLobe::ID]->get_var(s->get_ID(), k),this);

            }

        }




        ///8e-Check if the binary evolution creates a nakedHe or nakedCO star
        ///WARNING: This function needs to be exactly there, after the update of the single star properties and before the update of the binary property
        //TODO Myabe we have to call this before to update sse and check if the new total mass will be unphysical
        check_nakedHe_or_nakedCO_after_binary_evolution();




        ///8f-Check outcome of accretion on compact objects
        //Check after possible correction of the accreted mass
        check_accretion_on_compact();
        if (check_and_set_broken()) return utilities::BIN_EV_SETBROKEN; //THis is not a binary anymore, return.


        ///8h-Update derived properties (basic properties could have been changes during binary processes)
        //GI: 17/03/22, the call to update_derived_properties_star is needed because when Inertia is
        //estimated not using the table we have to update it before to call check AnMomSPin to limit the spin to critical.
        //In principle we could call just the evolve of Inertia, but it could be dangerous (to check!)
        update_derived_properties_star(getstar(0));
        update_derived_properties_star(getstar(1));



        ///8h-Check AngMomSpin (reduce it if it has reached the maximum spin)
        //TODO At a certain point we have to analyse the effect of reaching critical velocities (maybe we have to remove mass?)
        //TODO When we pass to AngMom binary as fundamental parameter and a as derived, here we have to add to the AngMom
        //all the angular momentum in excess with respect to the one related to the critical velocity
        check_AngMomSpin_after_binary_evolution();


    }



    ///9-If system is not broken update the binary properties
    if (!check_and_set_broken()) {
        // Update the Binary  properties due to the processes
        for (auto &prop: property)
            prop->evolve(this);


        ///9b-Check if the Eccentricity is larger than 1
        //If eccentricity (just updated) is larger than one and a repeat has not been triggered
        // trigger  a mix at the end of the evolution and return.
        //Notice this is an exception wrt to the fact that the special evolution
        //are handled at the end. We are assuming that in the time step the processes
        //smoothly increase the eccentricity until a head-on merger is triggered at the end of the timestep
        //This is necessary to avoid to start the next step with Eccentricity>1 and make a mess with
        //other processes (that expect Eccentricity<1).
        if(getp(Eccentricity::ID)>=1 and !repeatstep){
            //TODO Maybe we have to create a proper special_evolve in Hardening
            //Merge the stars
            mix=true; //Enable mixing
            process[Mix::ID]->special_evolve(this); //mix
            if (check_and_set_broken()) return utilities::BIN_EV_SETBROKEN; //THis is not a binary anymore, return.
        }


        // Note. The evolution of BTimestep can set repeatstep to true.
        // The re-evolution is handled in  binstar->evolve.
        return utilities::BIN_EV_DONE;
    }
    else{
        return utilities::BIN_EV_NOT_DONE;
    }



}

void Binstar::synchronise_dt_star(){

    //TODO Add other case as isempty, is naked.
    if ( (star[0]->getp(Worldtime::ID)!=star[1]->getp(Worldtime::ID))||(star[0]->getp(Worldtime::ID)!=getp(BWorldtime::ID)))
        svlog.critical("Stellar and/or Binary Worldtimes are not sync S0_T0=" +std::to_string(star[0]->getp_0(Worldtime::ID)) +
        " S0_T=" +std::to_string(star[0]->getp(Worldtime::ID)) + " S1_T0="+std::to_string(star[1]->getp_0(Worldtime::ID)) +
                               " S1_T="+std::to_string(star[1]->getp(Worldtime::ID)) +
                               " B_T0=" + std::to_string(getp_0(BWorldtime::ID)) + " B_T=" +std::to_string(getp(BWorldtime::ID)),__FILE__,__LINE__);

    //Now the remnant tmax is handled directly in Timestep
    //if (star[0]->isremnant) star[0]->sync_with(tf-star[0]->getp(Worldtime::ID));
    //if (star[1]->isremnant) star[1]->sync_with(tf-star[1]->getp(Worldtime::ID));
    //We can leave this check because if broken the Binary properties are not evolved.
    if (broken) property[BTimestep::ID]->resynch(tf-star[0]->getp(Worldtime::ID));
    //property[BTimestep::ID]->resynch(1e30, false);  //TODO I have to evolve Btimestep
    double mindt;
    mindt = std::min({star[0]->getp(Timestep::ID), star[1]->getp(Timestep::ID),property[BTimestep::ID]->get()});
    star[0]->sync_with(mindt);
    star[1]->sync_with(mindt);
    property[BTimestep::ID]->resynch(mindt, false);


}

void Binstar::check_and_sync_sse() {

    //GI changed (17/10) the old if with a while  since the initial assumpation that using a smaller timestep guarantee
    //that the stars will not repeat the evolution seems not to hold in some cases (e.g. due to the oscillating evolution of the Radius).
    //The while guarantee that exiting from this function the two stars will be perfectly syncronised.
    int _count=0;

    while(star[0]->getp_0(Timestep::ID)!=star[1]->getp_0(Timestep::ID)){
        /// Re-evolve the star that did the longest timestep
        Star * to_resynch= star[0]->getp_0(Timestep::ID)>star[1]->getp_0(Timestep::ID)? star[0] : star[1];
        Star * other= star[1]->getp_0(Timestep::ID)>star[0]->getp_0(Timestep::ID)? star[0]  : star[1];
        //Restore data and resynch time step (handled by resynch(dt))
        to_resynch->resynch(other->getp_0(Timestep::ID));
        to_resynch->evolve();

        //Write  a warning if we have to repeat more the one. This means that a star signaled a repeat
        //after using a smaller timestep with respect to a larger one that was accepted.
        //This can be due to the oscillating evolution of some properties (e.g. the Radius).
        if (_count>0)
            svlog.pwarning("The star with ID",to_resynch->get_ID(),"repeated the evolution",
                    "after an original repetition made with a dt smaller than the last accepted one. This can be an hint"
                    "of an oscillating evolution of some properties (e.g. Radius) that is not followed well enough.");

        _count++;
    }

    ///Now update the BinaryTimestep to be equal to the actual stellar evolution timestep
    property[BTimestep::ID]->resynch(star[0]->getp_0(Timestep::ID), false);


}

void Binstar::check_and_sync_bse() {

    //We changed the old if(repeatstep) to while(repeatstep) because it is not sure that we reached convergence after the first repetition.
    unsigned _count=0;

    while (repeatstep){

        if (_count>get_svpar_num("ev_max_repetitions"))
            svlog.critical("The number of repetitions in bse reaches the maximum allowed ("+
                           utilities::n2s(get_svpar_num("ev_max_repetitions"),__FILE__,__LINE__)+"). You can bypass this error"
                                                                                                 " increasing the parameter ev_max_repetitions. However"
                                                                                                 " a large number of repetitions is usually an hint that something is broken.",
                           __FILE__,__LINE__,sevnstd::bse_error());


        //utilities::wait("BSE repeat step",__FILE__,__LINE__, getp(BTimestep::ID));

        ///Safety check on timesteps: the new proposed timestep should be smaller than the stellar one.
        if ( (getp(BTimestep::ID)>star[0]->getp_0(Timestep::ID)) || (getp(BTimestep::ID)>star[1]->getp_0(Timestep::ID))  ){
            utilities::wait("Check",broken,onesurvived,__FILE__,__LINE__);
            svlog.critical("The new proposed Binary timestep [" + utilities::n2s(getp(BTimestep::ID),__FILE__,__LINE__) +"] is larger than stellar timesteps [" +
                           utilities::n2s(star[0]->getp_0(Timestep::ID),__FILE__,__LINE__) +", " + utilities::n2s(star[1]->getp_0(Timestep::ID),__FILE__,__LINE__)+"]",__FILE__,__LINE__);
        }

        ///We have to evolve again the stars and the binary
        ///Re-evolve the stars with the new timestep propsed by the last call of BTimestep.evolve
        for (auto& s : star){
            //We use getp and not getp_0, as in check_and_sync_sse,
            // because if repeatstep is true,the new propesed timestep is the one to use in the repetition
            svlog.debug("Inside bse " + std::to_string(s->getp_0(Timestep::ID)) + " " + std::to_string(s->getp(Timestep::ID)));
            //if (getp(BWorldtime::ID)>4.5668373678e+01)
            //    utilities::wait("Before re-evolve",getp(BTimestep::ID),__FILE__,__LINE__);

            //Restore data and resynch time step (handled by resynch(dt))
            s->resynch(getp(BTimestep::ID));

            //CHeck if the star repeat the evolutin, this is not expected since the new BTimestep is smaller than the last accepted one.
            if(s->evolve()==utilities::REPEATED_EVOLUTION)
                svlog.warning("Star ID "+utilities::n2s(s->get_ID(),__FILE__,__LINE__)+ " in binary ID " +
                              utilities::n2s(get_ID(),__FILE__,__LINE__)+"("+get_name()+"):\n repeated evolution "
                               "after an original repetition made with a dt smaller than the last accepted one. This can be an hint "
                               "of an oscillating evolution of some properties (e.g. Radius) that is not followed well enough.", __FILE__,__LINE__);

            svlog.debug("Inside bse " + std::to_string(s->getp_0(Timestep::ID)) + " " + std::to_string(s->getp(Timestep::ID)));
            //if (getp(BWorldtime::ID)>4.5668373678e+01)
            //    utilities::wait("After re-evolve",getp(BTimestep::ID),__FILE__,__LINE__);
        }

        //TODO Since we add this, now it seems we do not need anymore chec_and_sync bse
        //We can instead use a single cycle in evolve of the type:
        // for(;;)
        //      evolve stars
        //      check_and_sync_sse
        //      eolve binary
        //      if (repeat step)
        //          resync stars, resync binaries
        //      else
        //            break

        //GI Added this 17/10. This is needed because some time the adaptive time step skip
        //some part where the star is varying and if binary want to use a smaller time step
        //we can end in this varying part and the star ask for another smaller time step.
        check_and_sync_sse();


        //if (star[0]->getp_0(Timestep::ID)!=star[1]->getp_0(Timestep::ID))
        //    utilities::wait("WEEE",star[0]->getp_0(Timestep::ID),star[1]->getp_0(Timestep::ID),__FILE__,__LINE__);

        /// Restore binary properties and revolve (processes are restored in evolve_binary)
        //We use getp and not getp_0, as in check_and_sync_sse,
        // because if repeatstep is true, the new proposed timestep is the one to use in the repetition
        resynch(getp(BTimestep::ID));
        evolve_binary();

        _count++;
    }

}

void Binstar::check_and_set_bavera_xspin(){
    if (get_svpar_bool("xspin_bavera")){
        for (const unsigned i : {0,1}){
            if (star[i] -> just_exploded() and star[i] -> amiBH() and star[i] -> amiWR_0() and star[1-i] -> amiBH()){
                //Apply Bavera correction
                //1- we are sure that the remnant is a BH so cast staremnant to the BH type
                auto brem = dynamic_cast<BHrem*>(star[i]->staremnant);
                //2-Apply the bavera correction
                //Notice the Period is the Period of the binary before of the SN kick
                //The mass is the mass of the remnant progenitor
                brem->apply_Bavera_correction_to_Xspin(getp(Period::ID), star[i] -> getp_0(Mass::ID)); //Before it was getp(Mass::ID) but this is the mass of the remnant ath this point
                //3- Deal with current modification of Xspin (in the nex steps this is not necessary since xspin has been changed internally in the staremant class)
                auto dv = brem->get(star[i],Xspin::ID) - star[i] -> getp(Xspin::ID);
                star[i] -> update_from_binary(Xspin::ID, dv);
            }
        }
    }
}


/*************************************/

/**********************************
 *  Special Evolve functions
 **********************************/

bool Binstar::check_and_set_QHE(Star *accretor){


    ///Check the preliminary condition for QHE:

    ///Very first check
    //The parameter rlo_QHE has to be turnend on
    if (!accretor->get_svpar_bool("rlo_QHE"))
        return false;

    ///Second preliminary check
    //1- The star should be in the MainSequence and should not have changed phase in the last evolution step
    //2- The star should not be already flagged as amifollowingQHE
    if(accretor->getp_0(Phase::ID)!=accretor->getp(Phase::ID) or accretor->getp(Phase::ID)!=MainSequence
        or accretor->amifollowingQHE())
        return false;

    ///Third Check, is RLO happening and is the Star the accretor?
    double DM_cumulative_RLO = accretor->getp(dMcumul_RLO::ID);
    //If the stars is not the accretor or RLO is not happening or rlo_QHE is disabled just return
    if (DM_cumulative_RLO<=0.)
        return false;


    ///Prelimary Checks passed, now do the main check based on Eldridge&Stanway11
    //TODO should be parameters?
    //Prameters from Sec. 2.2 of Eldridge&Stanway11 (https://arxiv.org/abs/1109.0288)
    double Zlim = 0.004;  //Condition on Metallicity
    double Mtot_threshold  = 10; //Condition on total Mass after accretion
    double fMass_threshold = 0.05;  //Condition on fraction of accreted Mass

    //The metallicity has to be lower or equatl to Zlim;
    bool trigger = accretor->get_Z()<=Zlim;
    //The Total Mass after evolution has to be larger than Mtot_threshold
    trigger = trigger and accretor->getp(Mass::ID)>Mtot_threshold;
    //The accreted mass has to be larger than fMass_threshold of the initial Mass.
    trigger = trigger and DM_cumulative_RLO>fMass_threshold*accretor->getp_0(Mass::ID);

    ///If trigger is true, set the QHE
    if (trigger)
        accretor->set_QHE();

    //utilities::hardwait("INSIDE",accretor->get_ID(),DM_cumulative_RLO,trigger,accretor->amifollowingQHE(),__FILE__,__LINE__);


    return trigger;

}

void Binstar::update_derived_properties_star(Star *s){
    //Update derived properties (basic properties could have been changes during binary processes)
    if (!s->amiempty()){
        for (auto& prop: s->properties){
            if (prop->amiderived()) prop->update_derived(s);
        }
    }
    //DO nothing for empty star (No binary processes in this case
}


utilities::sn_explosion Binstar::check_accretion_on_compact(size_t donorID, size_t accretorID, double DMaccreted){

    Star *donor    = getstar(donorID);
    Star *accretor = getstar(accretorID);

    double Maccretor = accretor->getp(Mass::ID); //Use current mass after evolution of stellar properties
    double remnant = accretor->getp(RemnantType::ID);//accretor->get_supernova()->get_remnant_type();
    Material donor_material = donor->whatamidonating_0(); //Check the material before the evolution of stellar properties

    //If it is not  a remnant return
    if (remnant==Remnants::NotARemnant)
        return utilities::SN_NOT_EXPLODE;

    //HeWD can only accrete He-rich matieral up to a mass of 0.7 or is it destroyed as SN Ia.
    if (remnant==Remnants::HeWD and donor_material==Material::He and Maccretor>=0.7 and DMaccreted>0){
        //trigger Type 1a SN explosion
        ///trigger Type 1a SN explosion
        //accretor->explode_as_SNI(this); //Set the star to empty
        accretor->explode_as_SNI(); //Set the star to empty
        //The star is now empty, the binary system becomes broken or empty
        set_broken();
        set_onesurvived(); //If onesurvived is already true, the system becomes empty

        utilities::wait("SN1a triggered by tripleHa",__FILE__,__LINE__);

        return utilities::SNIA_EXPLODE;
    }

    /* CO and ONeWDs accrete helium-rich material until the accumulated
    * material exceeds a mass of 0.15 when it ignites. For a COWD with
    * mass less than 0.95 the system will be destroyed as an ELD in a
    * possible Type 1a SN. COWDs with mass greater than 0.95 and ONeWDs
    * will survive with all the material converted to ONe (JH 30/09/99).
    * Now changed to an ELD for all COWDs when 0.15 accreted (JH 11/01/00).
    */
    if (remnant==Remnants::COWD and donor_material==Material::He and DMaccreted>=0.15){
        ///trigger Type 1a SN explosion
        //accretor->explode_as_SNI(this); //Set the star to empty
        accretor->explode_as_SNI(); //Set the star to empty

        //The star is now empty, the binary system becomes broken or empty
        set_broken();
        set_onesurvived(); //If onesurvived is already true, the system becomes empty


        utilities::wait("SN1a triggered by accretion on a COWD",__FILE__,__LINE__);
        return utilities::SNIA_EXPLODE;
    }

    //Possible transition due to increase of mass over the chandrasekhar or vto mass limits are handled elsewhere

    return utilities::SN_NOT_EXPLODE;

}

/*************************************/


/**********************************
 *  Init functions
 **********************************/

void Binstar::init(const std::vector<std::string> &params){

    //Check if dimension is correct. The check is made here since init is publicly exposed (init_stars and init_binary are privates).
    check_init_param(params);



    /***********************************
    * Init binary properties
    ************************************/
    //Init binary property (a,e)
    init_binary_properties(params);


    //Random seed
    /*** Random number generator and random name ***/
    if(io->rseed_provided()){
        std::string rseed_provided = params[inputmapbin.at(io->binput_mode).second]; //very last element of the params vector
        set_rseed(utilities::s2n<unsigned long>(rseed_provided, __FILE__, __LINE__), __FILE__,__LINE__); //7th column is the random seed, if provided
    }
    else if(rseed==0){
        set_rseed(utilities::gen_rseed(), __FILE__, __LINE__);
    }

    //initialize the random number generator with the seed (to generate IDs for binary systems)
    utilities::mtrand.seed(rseed);

    //Generate name
    name = get_svpar_str("name_prefix")+utilities::random_keygen(&utilities::mtrand);
    /*********************************************/






    /***********************************
    * Init single star properties
    ************************************/
    //GI Note: no need to use Lookup::OutputOption..ecc because we have defined using namespace Lookup in IO.h
    int inchoice=inputmapbin.at(io->binput_mode).first;
    //Init stars based on input type
    switch(inchoice){
        case InputBinaryOption::_new:
            init_stars(params);
            //init_binary(params);
            break;
        case InputBinaryOption::_legacy:
            init_stars_legacy(params);
            //init_binary_legacy(params);
            break;
        default:
            svlog.critical("Error on input binary choice ", __FILE__, __LINE__,sevnstd::sevnio_error());
    }
    /*********************************************/


    /***********************************
    * Init the derived properties. These depends both on the binary and stellar properties,
    * so we init them at the end.
    ************************************/
    for (auto& prop : property)
        prop->init_derived(this);
    /*********************************************/


    /***********************************
    * Init other properties
    ************************************/
    //Notice this need to be the last thing because both binary and stellar properties have to be already initiliased
    init_other_params();
    /*********************************************/


}

void Binstar::call_stars_constructor(std::vector<std::string> &params_star1, std::vector<std::string> &params_star2){

    size_t id1 = 0;
    size_t id2 = 1;

    if (io->rseed_provided()){
        params_star1.push_back(utilities::n2s(get_rseed(),__FILE__,__LINE__));
        params_star2.push_back(utilities::n2s(get_rseed(),__FILE__,__LINE__));
    }



    star[0] = new Star(io, params_star1, id1, false,get_rseed());
    std::ostringstream ss_star;
    ss_star<<*star[0];
    svlog.debug("Stars 1 initialised: IDbin: " + utilities::n2s(id1, __FILE__,__LINE__) + " " + ss_star.str());

    star[1] = new Star(io, params_star2, id2, false,get_rseed());
    ss_star.str("");
    ss_star<<*star[1];
    svlog.debug("Stars 2 initialised: IDbin: " + utilities::n2s(id1, __FILE__,__LINE__) + " " + ss_star.str());
}

void Binstar::init_stars(const std::vector<std::string> &params){

    //Common variable
    std::string  tf{params[12]}, dtout{params[13]};

    set_break_at_broken(tf);


    //Initialise vectors and assign fist elements (m,z,spin,sn, tini)
    std::vector<std::string> param_star1(params.begin(), params.begin()+5);
    std::vector<std::string> param_star2(params.begin()+5, params.begin()+10);
    //Insert the rest of the star parameters
    param_star1.insert(param_star1.end(), {tf, dtout});
    param_star2.insert(param_star2.end(), {tf, dtout});


    call_stars_constructor(param_star1, param_star2);

}

void Binstar::init_stars_legacy(const std::vector<std::string> &params){

    //Common variable
    std::string tf{params[8]}, t_ini{params[9]}, dtout{params[13]};

    set_break_at_broken(tf);

    //Mass
    std::string m1{params[0]}, m2{params[1]};
    //Z
    std::string z1{params[2]}, z2{params[3]};
    //Spin
    std::string sp1{params[4]}, sp2{params[5]};
    //sn_type
    std::string sn1{params[11]}, sn2{params[12]};

    //Initialise vectors
    std::vector<std::string> param_star1{m1, z1, sp1, sn1, t_ini, tf, dtout};
    std::vector<std::string> param_star2{m2, z2, sp2, sn2, t_ini, tf, dtout};

    call_stars_constructor(param_star1, param_star2);
}

void  Binstar::init_other_params(){


    //In this function we set the parameters that are already set in each star when we call init_stars or init_stars_legacy.
    //These parameters are: break_at_remnant, print_all_steps and print_per_phase, tf, dtout.
    //It is not useful to make another specialised function to assign them also in this class.
    // This function just check that the values are the same for the stars and then use the values from one of them to set also
    // the corrispondent attribute in this class.

    //TODO the print stuff should be totally remove the from the class star and binarstar, it should be handlded by an external class
    //the star and binary should just evolve not deal with storing data, this should be a separate class
    //dtout


    if (star[0]->get_dtout()==star[1]->get_dtout()) set_dtout(star[0]->get_dtout(), __FILE__, __LINE__);
    else svlog.critical("dtout attribute is not the same for the two stars",__FILE__,__LINE__,sevnstd::sevnio_error());


    //break_at_remnant
    if (star[0]->break_at_remnant and star[1]->break_at_remnant){
        break_at_remnant=star[0]->break_at_remnant;
        set_tf(utilities::LARGE,__FILE__,__LINE__);
        star[0]->set_tf(utilities::LARGE,__FILE__,__LINE__);
        star[1]->set_tf(utilities::LARGE,__FILE__,__LINE__);
    } else if (star[0]->get_tf()==star[1]->get_tf()) {
        set_tf(star[0]->get_tf(), __FILE__, __LINE__);
    } else if (star[0]->get_tf()!=star[1]->get_tf()){
        svlog.critical("tf attribute is not the same for the two stars",__FILE__,__LINE__,sevnstd::sevnio_error());
    } else{
        svlog.critical("break_at_remnant attribute is not the same for the two stars",__FILE__,__LINE__,sevnstd::sevnio_error());
    }


    //print_all_steps
    if (star[0]->print_all_steps==star[1]->print_all_steps) print_all_steps=star[0]->print_all_steps;
    else svlog.critical("print_all_steps attribute is not the same for the two stars",__FILE__,__LINE__,sevnstd::sevnio_error());

    //print_only_end
    if (star[0]->print_only_end==star[1]->print_only_end) print_only_end=star[0]->print_only_end;
    else svlog.critical("print_only_end attribute is not the same for the two stars",__FILE__,__LINE__,sevnstd::sevnio_error());

    //print_per_phase
    if (star[0]->print_per_phase==star[1]->print_per_phase) print_per_phase=star[0]->print_per_phase;
    else svlog.critical("print_per_phase attribute is not the same for the two stars",__FILE__,__LINE__,sevnstd::sevnio_error());
    //TODO We have to take into account that the two stars can have different phases times. For now we just disable the print_per_phase option
    if (print_per_phase) svlog.critical("print_per_phase is currently disabled in binary stellar evolution",__FILE__,__LINE__);

    //print_events
    if (star[0]->printevents()==star[1]->printevents()){
        print_events=star[0]->printevents();
        //Check if dtout option is eventsrlo
        //const auto dtout_current=utilities::s2n<double>(star[0]->get_initparams()[Star::Init_params::_Dtout],__FILE__,__LINE__);
        if (star[0]->get_dtout_type()==Lookup::DTout_type::_DTEVNTRLO){
               print_rlo=true;
        }
    }
    else svlog.critical("print_events attribute is not the same for the two stars",__FILE__,__LINE__,sevnstd::sevnio_error());


}

void Binstar::limit_and_correct_mass_transfer_for_donor_from_binary(Star *donor, Star *accretor,
                                                                    std::size_t mass_limiting_property_id){
    //TODO THis function it bit a mess, we should refactor it a t some point


    if (donor->getp(mass_limiting_property_id)<donor->getp(Mass::ID)){
        svlog.critical("The function limit_and_correct_mass_transfer_for_donor_from_binary can "
                       "work only when the limiting mass is lower than the total mass",__FILE__,__LINE__,sevnstd::bse_error());
    }

    ///Notice:
    //This function is called inside check_nakedHe_or_nakedCO_after_binary_evolution() that is
    //called after the update of the single stellar property and before the binary property.
    //Therefore, the single stellar property are corrected directly using update_from_binary, while
    //the binary property are corrected modifying the values in the VB tables of the processes using
    //the functions modify_SemimajorDV_by_a_factor and  modify_EccentricityDV_by_a_factor.
    ///

    double DM_RL_donor=0, DM_RL_accretor=0;


    svlog.pdebug("Inside limit and correct mass transfer","Mdonor",donor->getp(Mass::ID),"Maccretor",accretor->getp(Mass::ID),__FILE__,__LINE__);

    ///Estimate the DM from the processes
    //Consider all the process, if we enter in this part of the code for sure Mtot<MHe, therefore for sure DM for the donor is negative/
    //and it is positive for the accretor. Cycle over all the processes to get the total DM variation
    for (auto &proc : process){
        DM_RL_donor    += proc->get_var(donor->get_ID(), Mass::ID);
        DM_RL_accretor += proc->get_var(accretor->get_ID(), Mass::ID);
    }

    //We are not here due to BSE evolution
    if (utilities::areEqual(DM_RL_donor,0.0))
        return;

    ///Estimate the actual fraction of mass accred on the accretor (considering all the processes)
    //This value will be used to rescale the accreted mass once the donated mass has been accreted.
    double f_accreted = std::abs(DM_RL_accretor/DM_RL_donor);



    ///Estimate the correcting term
    //The effective mass donated is the entire envelope (before the processes evolution) is Mtot (before evolution) = DM + Mtot(now)   - MHE (now)
    double DM_RL_donor_effective     =    std::abs(DM_RL_donor) + donor->getp(Mass::ID) - donor->getp(mass_limiting_property_id); //Mass that has to be
    ///Notice
    //DM_RL_donor_effective is smaller than 0 means that while the star was losing mass, the core has growth so that now is larger
    //than the total mass even accounting for the mass loss. Here we make a very simple decision, we assume no mass transfer has happened.
    //TODO A better stuff to do is to estimate what is the crossing time between the mass loss and the core growth
    //ad use this value to estimate the DM_RL_donor_effective (as done in MHE::correct_interpolation_error, for example), however in this case
    //we have to manually change also the properties of the MHE core, a better option will be to raise in some way a repetion flag
    //proposing as new timestep something a bit smaller than the Mass-MHE crossing time.
    DM_RL_donor_effective=std::max(0.,DM_RL_donor_effective);



    //The effective mass accreted is f_accreted times the effective DM
    double DM_RL_accretor_effective  =    f_accreted*DM_RL_donor_effective;
    //In order to correct the current values we add the old DM and substract the real mass lost
    double DM_RL_donor_correction    =    std::abs(DM_RL_donor)-DM_RL_donor_effective; //Add the mass lost in the process evolution and subtract the effective mass loss (all are absolue values)
    //In order to correct the current values we subctrad the old accred mass and we add the real accreted mass.
    double DM_RL_accretor_correction =    DM_RL_accretor_effective - DM_RL_accretor; //Subtract the mass accreted in the process evolution and add the effective accreted mass.

    //GI 20/12/22: GI, before was without the abs, but this shold be the correct form
    double correction_fraction       =  std::abs(DM_RL_accretor_effective/DM_RL_donor);


    ///Correct mass
    //Donor
    //std::cout<<" Before donor "<< donor->getp(Mass::ID)<<" "<<donor->getp(dMcumul_binary::ID)<<" "<<donor->getp(dMcumul_RLO::ID)<<std::endl;
    donor->update_from_binary(Mass::ID,DM_RL_donor_correction);
    donor->update_from_binary(dMcumul_binary::ID,DM_RL_donor_correction);
    const auto correction_dMdonor=std::abs(DM_RL_donor_effective/DM_RL_donor);
    const auto DVRLO_donor= process[RocheLobe::ID]->get_var(donor->get_ID(),dMcumul_RLO::ID);
    donor->update_from_binary(dMcumul_RLO::ID,DVRLO_donor*correction_dMdonor - DVRLO_donor,this);
    //std::cout<<" After donor "<< donor->getp(Mass::ID)<<" "<<donor->getp(dMcumul_binary::ID)<<" "<<donor->getp(dMcumul_RLO::ID)<<std::endl;


    //Accretor
    //std::cout<<" Before accretor "<< accretor->getp(Mass::ID)<<" "<<accretor->getp(dMcumul_binary::ID)<<" "<<accretor->getp(dMcumul_RLO::ID)<<" "<<accretor->getp(dMcumulacc_wind::ID)<<" "<<DM_RL_accretor_correction<<std::endl;
    accretor->update_from_binary(Mass::ID, DM_RL_accretor_correction);
    accretor->update_from_binary(dMcumul_binary::ID,DM_RL_accretor_correction);
    if (DM_RL_accretor > 0){
        const auto correction_dMccretor=std::abs(DM_RL_accretor_effective/DM_RL_accretor) ;
        const auto DVRLO_accretor= process[RocheLobe::ID]->get_var(accretor->get_ID(),dMcumul_RLO::ID);
        const auto DVwind_accretor= process[Windaccretion::ID]->get_var(accretor->get_ID(),dMcumulacc_wind::ID);
        accretor->update_from_binary(dMcumul_RLO::ID,DVRLO_accretor*correction_dMccretor - DVRLO_accretor,this);
        accretor->update_from_binary(dMcumulacc_wind::ID,DVwind_accretor*correction_dMccretor - DVwind_accretor);
    }
    //std::cout<<" After accretor "<< accretor->getp(Mass::ID)<<" "<<accretor->getp(dMcumul_binary::ID)<<" "<<accretor->getp(dMcumul_RLO::ID)<<" "<<accretor->getp(dMcumulacc_wind::ID)<<std::endl;






    //Correct nakedHelium
    if (donor->aminakedhelium()) donor->copy_property_from(MHE::ID,Mass::ID);
    if (accretor->aminakedhelium()) accretor->copy_property_from(MHE::ID,Mass::ID);

    //Estimate new values for NS
    //In this part we correct the values evolved in NS (Bmag and OmegaRem) due to the mass accretion
    if (accretor->amiNS() and DM_RL_accretor_effective!=0){
        //1-Step: estimate the new cumlative DBmag and DOmegaRem  based on the amount of mass accrete.
        //Since in principle each process can have its own way to estimate these values, we
        //cycle over all processes, checking if they belog to the MaccretionProcess family, if this is the case
        //we estimate new values of DB and Domega using the corrected effective DM assuming that the old proprotionaly in the mass
        //accreted between the processes remains the same (E.g. if before the RLO accounts for the 80% of the accreted mass and the wind for the 20%, this
        // fraction remains the same when considering the new corrected mass)
        double DBmag_effective=0;
        double DOmegaRem_effective=0;
        double DBmag_old=0;
        double DOmegaRem_old=0;
        for (auto &proc : process){

            //Check the amount of variation in Bmag and OmegaRem, it is 0 for both, skip to the next proc
            double DBmag_old_proc = proc->get_var(accretor->get_ID(), Bmag::ID);
            double DOmegaRem_old_proc = proc->get_var(accretor->get_ID(), OmegaRem::ID);
            if (DBmag_old_proc==0 and DOmegaRem_old_proc==0){
                continue;
            }


            MaccretionProcess* Mproc = dynamic_cast<MaccretionProcess*>(proc);
            if (Mproc!= nullptr){
                DBmag_old += DBmag_old_proc;
                DOmegaRem_old += DOmegaRem_old_proc;
                DBmag_effective += Mproc->NS_DBmag_accretion(this,accretor,DM_RL_accretor_effective);
                DOmegaRem_effective += Mproc->NS_DOmegaRem_accretion(this,accretor,DM_RL_accretor_effective);
            }
        }
        //Update the values removing the old DV.
        accretor->update_from_binary(Bmag::ID,DBmag_effective-DBmag_old);
        accretor->update_from_binary(OmegaRem::ID,DOmegaRem_effective-DOmegaRem_old);
    }

    ///Correct orbital properties
    //The DV of the orbital properties set by the process are estimated using the old (no correct) mass transfer
    //In this part we correct the variation of Semimajor and Eccentricity (if circularises is True the DV eccentricity is not changed)
    //Notice here we correct the DV instead of the values because at this point we have still to call the evolve of binary properties.
    //TODO The correction here assumes that DA and DE change linearly with DM but this is not really the case.
    for (auto &proc : process){
        //Check if the current process is causing a mass transfer, if yes correct the estimated DV by the correction fraction
        if ( proc->is_mass_transfer_happening()){
            proc->modify_SemimajorDV_by_a_factor(correction_fraction); //DA_new=DA_old*correction_fraction
            proc->modify_EccentricityDV_by_a_factor(correction_fraction); //DE_new=DE_old*correction_fraction
        }
    }
}

void Binstar::check_nakedHe_or_nakedCO_after_binary_evolution(){

    double RLs[2] = {getp(RL0::ID),getp(RL1::ID)};

    //Set star pointer (swapped below), donor and accretor are just placeholder
    Star *donor = star[0];
    Star *accretor = star[1];

    for (size_t i=0; i<2; i++) {
        /** Now start to consider different cases **/
        /* if the star is a remnant, just skip, nothing to check */
        if (!donor->amiremnant()){
            //Flag to check if the star has just become a naked Hw or a naked CO (considering the ev_naked_tshold flag)
            bool is_become_He_naked = donor->getp(MHE::ID)>0 and donor->getp(MHE::ID)+get_svpar_num("ev_naked_tshold")>=donor->getp(Mass::ID) and !donor->aminakedhelium(); //Check if the star has become just now a naked helium star
            bool is_become_CO_naked = donor->getp(MCO::ID)>0 and donor->getp(MCO::ID)+get_svpar_num("ev_naked_tshold")>=donor->getp(Mass::ID) and !donor->aminakedco(); //Check if the star has become just now a naked helium star

            // If at least one of the two flag is true, it means that binary the combination
            // of stellar and binary processes removed more mass than what necessary to produce a naked  (He or CO) star
            // In this case we have to take care of the naked star producetion, otherwise we jsut skipe.
            if (is_become_He_naked or is_become_CO_naked){
                int donor_id = donor->get_ID();
                // Define other two flag to check if the Roche Lobe radius is within the core.
                // This is the only situation in which we allow to actually remove mass from the core.
                bool is_RL_inside_COcore = donor->getp(RCO::ID) > RLs[donor_id];
                bool is_RL_inside_core = donor->getp(RHE::ID) > RLs[donor_id];

                /** Check now different cases **/
                //1-Case, MTOT<MCO+tshold or MTOT<MHE+tshold, but RL>RHe and RL>RCO, so the star loses just the envelope.
                //so the outcome is a Naked Helium star or a Naked CO star. In this case we remove the excess of mass lost
                //from the core.
                if (!is_RL_inside_core and !is_RL_inside_COcore){
                    //Estimate the real Mass loss (and accretion) correcting for the Mass exceding the envelope
                    /** Notice
                     * we can enter here both because the total Mass is lower than the MHE (or MCO) Mass or because
                     * the difference is within the tshol.
                     * - In the first case we have to correct the masses and the orbit evolution.
                     * - In the second case we can go ahead and just create the nakedHe or nakedCO
                     */

                    //TODO The function limit_and_correct_mass_transfer_for_donor_from_binary works only if the donor is
                    // losing more mass than MHE, we should make this function general to be used also when we end here
                    // because  the difference is within the tshold
                    if (is_become_CO_naked){
                        // if MCO is currently larger than Mass correct the mass excess
                        if (donor->getp(MCO::ID)>donor->getp(Mass::ID))
                            limit_and_correct_mass_transfer_for_donor_from_binary(donor, accretor, MCO::ID);
                        // Set the nakedCO properties, the proper creation of the nakedCO is triggered and the end of the evolution step
                        donor->properties[Radius::ID]->copy_V_from(donor->properties[RCO::ID]);
                        donor->properties[Mass::ID]->copy_V_from(donor->properties[MCO::ID]);
                        donor->properties[RHE::ID]->copy_V_from(donor->properties[RCO::ID]);
                        donor->properties[MHE::ID]->copy_V_from(donor->properties[MCO::ID]);
                    } else if (is_become_He_naked){
                        // if MHE is currently larger than Mass correct the mass excess
                        if (donor->getp(MHE::ID)>donor->getp(Mass::ID))
                            limit_and_correct_mass_transfer_for_donor_from_binary(donor, accretor, MHE::ID);
                        // Set the nakedHE properties, the proper creation of the nakedCO is triggered and the end of the evolution step
                        donor->properties[Radius::ID]->copy_V_from(donor->properties[RHE::ID]);
                        donor->properties[Mass::ID]->copy_V_from(donor->properties[MHE::ID]);
                    }
                }
                //2-Case, MCO>MTOT and RCO>RL, we are actually removing mass from the core (both He and CO core)
                //The results ia naked CO star
                else if(is_RL_inside_COcore and is_become_CO_naked){
                    donor->properties[Radius::ID]->copy_V_from(donor->properties[RCO::ID]);
                    donor->properties[Mass::ID]->copy_V_from(donor->properties[MCO::ID]);
                    donor->properties[RHE::ID]->copy_V_from(donor->properties[RCO::ID]);
                    donor->properties[MHE::ID]->copy_V_from(donor->properties[MCO::ID]);
                    //TODO THis is a very interesting situation to investigate at a certain point (SEVN MEETING 25/05/2020)
                    //so we raise a non-blocking error just to signal it (myabe not so interesting anymore).
                    svlog.error(
                            "HEYY: We are in situation where the Radius of the CO Core is filling the ROche Lobe and we are transferring mass "
                            "both from the envelope and from the core (Binary name= "+
                            get_name()+").", __FILE__, __LINE__,false);
                }
                //3-Case, MCO>MTOT and RHE>RL>RCO>RL, we limit the mass loss to the current the CO MASS
                //The result is a naked CO star
                else if (is_RL_inside_core and is_become_CO_naked) {
                    // if MCO is currently larger than Mass correct the mass excess
                    if (donor->getp(MCO::ID)>donor->getp(Mass::ID))
                        limit_and_correct_mass_transfer_for_donor_from_binary(donor, accretor, MCO::ID);
                    // Set the nakedCO properties, the proper creation of the nakedCO is triggered and the end of the evolution step
                    donor->properties[Radius::ID]->copy_V_from(donor->properties[RCO::ID]);
                    donor->properties[Mass::ID]->copy_V_from(donor->properties[MCO::ID]);
                    donor->properties[RHE::ID]->copy_V_from(donor->properties[RCO::ID]);
                    donor->properties[MHE::ID]->copy_V_from(donor->properties[MCO::ID]);

                }
                //4-Case, MHE>MTOT and RHE>RLRL<RHE, the star loses mass from envelope and HE core
                else if (is_RL_inside_core and is_become_He_naked) {
                    donor->properties[Radius::ID]->copy_V_from(donor->properties[RHE::ID]);
                    donor->properties[MHE::ID]->copy_V_from(donor->properties[Mass::ID]); //MHE>=Mtot
                    //TODO THis is a  interesting situation to investigate at a certain point (SEVN MEETING 25/05/2020)
                    //so we raise a non-blocking error just to signal it (myabe not so interesting anymore).
                    svlog.warning(
                            "HEYY: We are in situation where the Radius of the Core is filling the ROche Lobe and we are trasnferring mass "
                            "both from the envelope and from the core.", __FILE__, __LINE__);
                }
                //Other cases? The cases above should be enough to take into account all the possible cases
                //Here we put an error to check if some case is outside our selection
                else {
                    svlog.critical(
                            "Where is your mind? Here we have a case that has not been considered in the if/else if cases",
                            __FILE__, __LINE__, sevnstd::bse_error());
                }

            }
        }
        // Since we do not knwo a priori which star is the donor one, swap the stars and check again
        utilities::swap_stars(donor,accretor);
    }

}

void Binstar::check_AngMomSpin_after_binary_evolution() {

    for (auto s : {getstar(0), getstar(1)}){
        double Lcrit = s->getp(Inertia::ID)*utilities::omega_crit(s->getp(Mass::ID),s->getp(Radius::ID));
        if (s->getp(AngMomSpin::ID)>Lcrit){
            s->update_from_binary(AngMomSpin::ID,Lcrit-s->getp(AngMomSpin::ID));
            update_derived_properties_star(s);
        } else if(s->getp(AngMomSpin::ID)<0.){
            s->update_from_binary(AngMomSpin::ID,fabs(s->getp(AngMomSpin::ID))+1e-10);
            update_derived_properties_star(s);
        }
    }
}
/*************************************/


/**********************************
 *  Other auxiliary functions
 **********************************/
double Binstar::Radx(size_t starID)  {

    if (starID!=0 and starID!=1){
        svlog.critical("starID can be only 0 or 1, you are using "+
        utilities::n2s(starID,__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::sanity_error(""));
    }

    const Star* s = getstar(starID);
    double Rc = s->aminakedhelium() ? s->getp_0(RCO::ID) : s->getp_0(RHE::ID); //beginning of the timestep
    double RL = starID==0? getp(RL0::ID) :  getp(RL1::ID); //binary is not getp_0 because it is only updated at the end

    return std::max(Rc,RL);
}

/**
 * Check all the processes alarm
 * @return True if at least one of the processes has a special_evolution_alarm switched on, False otherwise
 */
bool Binstar::process_alarm_any() const {

    for (const auto& proc: process){
        if (proc->is_special_evolution_alarm_on()) return  true;
    }

    return false;
}

bool Binstar::isoutputtime() {

        if (print_rlo and process[RocheLobe::ID]->is_process_ongoing() and !broken){
            return true;
        }
        else if (print_events){
            return property[BEvent::ID]->get()!=Lookup::EventsList::NoEvent;
        }

        return star[0]->isoutputtime();
}


Binstar::~Binstar(){
    default_destructor();
}

void Binstar::default_destructor(){

    ///Here we have to delete all the pointer to stuff allocated on the heap

    //Processes
    if (!process.empty()){
        for (auto& proc : process){
            delete proc;
            proc= nullptr;
        }
    }

    //Properties
    if (!property.empty()){
        for (auto& prop : property){
            delete prop;
            prop = nullptr;
        }
    }

    //Stars
    for (auto& _star : star){
        if (_star!= nullptr){
            delete _star;
            _star = nullptr;
        }
    }
}
