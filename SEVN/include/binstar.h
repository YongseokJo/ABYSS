//
// Created by mario on 04/12/18.
//

#ifndef SEVN_BINSTAR_H
#define SEVN_BINSTAR_H

#include <BinaryProperty.h>
#include <Processes.h>
#include <star.h>
#include <cmath>
#include <utilities.h>
#include <string>
#include <sstream>
#include <params.h>

#include <sevnlog.h>
using sevnstd::SevnLogging;

#include <lookup_and_phases.h>
#include "Orbit.h"

using namespace Lookup;


class Binstar {

public:
    ///MEMBERS///
    bool broken; /*!< If true the binary is broken*/
    bool onesurvived;  /*!< true if just one star survived else either both stars are present or no stars are present*/
    bool empty;
    bool break_at_remnant;  /*!< if true stop the evolution when both stars are remnant*/
    bool break_at_broken; /*!< if true stop the evolution when the binary breaks*/
    bool repeatstep;  /*!< Control if it is needed to repeat a step of the binary evolution. It is mostly used in BTimestep::evolve*/
    bool mix; /*!< If true the binary is mixing*/
    bool comenv; /*!< If true the binary is doing a  ce*/
    bool is_swallowed[2]; /*!< Array of bool corresponding to the ID of stars in binary. If true, the star has been swallowed*/
    double last_Timestep; /*!< it stores the last value of the Timestep used in the evolution */

    //TODO the print options should be removed from the classes stars and binaries and handled by another class only devoted to that
    /** Print options **/
    bool print_all_steps;  /*!<  If true print in output all the steps.*/
    bool print_per_phase;  /*!< If true print in output in each phase change*/
    bool print_only_end;  /*!< If true do not print intermediate steps */
    bool print_events; /*!< If true print all the events */
    bool print_rlo; /*!< If true print all the steps during a rlo */

    //Check disable, some process can force the e or a to jump to some values
    //in that case we have to disable the check for adaptive time step.
    bool disable_e_check; /*!< If true, the e_check in adpative timestep is disabled */
    bool disable_a_check; /*!< If true, the a_check in adpative timestep is disabled */
    bool disable_DM_check; /*!< If true, the DM_check (from binary evolution) in adpative timestep is disabled */
    bool disable_OmegaRem_NS_check; /*!< If true, the OmegaRem_NS_check (from binary evolution) in adpative timestep is disabled */
    bool disable_stellar_rotation_check; /*!< If true, the OmegaSpin (from binary evolution) in adpative timestep is disabled */
    bool force_tiny_dt; /*!< If true, it forces the next Binary time step to be equal to utilities::TINY */

    //Other variables
    double CE_E_tomatch = utilities::NULL_DOUBLE; /*!< Envelope binding energy of the system after coalesce during a CE  */

    //Guard
    bool evolution_step_completed;

    ////////////////

    ///CONSTRUCTORS///
    ///Constructors for fake objects
    Binstar(){}

    //Proper constructors
    Binstar(IO *_io, std::vector<std::string> &params, size_t &_ID, unsigned long _rseed=0){


        //GI 03/02/2023: Based on a issue regarding data leak after failed initialisation, I implement
        //this solution to robustly clean the heap allocated variables after a failed object initilisation.
        //Indeed, in such cases, the object is never build, so the desctructor is not available.
        //So, I have implemented a private method default_destructor that is called from the destructor to clean the heap
        //but can be called also here. So we use a try-catch-rethrow schema to check if the initilisation
        //was successful, otherwise we clean the heap and then rethrow the error.
        //We use this schema also for the other two constructors.

        //Try the initialisation
        try {
            /***********************************
            * Some variable assignment
            ************************************/
            io = _io;
            ID = _ID;
            break_at_remnant = false;
            broken = onesurvived = empty = false;
            break_at_remnant = break_at_broken = false;
            init_params = params;

            /** Print options **/
            print_all_steps = false;
            print_per_phase = false;
            print_events = false;
            print_rlo = false;

            last_Timestep = 0.0;
            repeatstep = false;
            mix = false;
            comenv = false;
            evolution_step_completed = false;
            is_swallowed[0] = is_swallowed[1] = false;
            disable_e_check = disable_a_check = disable_DM_check = disable_OmegaRem_NS_check = disable_stellar_rotation_check = force_tiny_dt = false;
            rseed = _rseed;
            /************************************/


            state.resize(io->printcolumns_binary.size());

            for (size_t i = 0; i < BinaryProperty::all.size(); i++) {
                property.push_back(BinaryProperty::all[i]->Instance()); //create instances for all properties
            }

            for (size_t i = 0; i < Process::all.size(); i++) {
                //std:cout<<"Procss "<< Process::all[i]->name() << " " << i <<std::endl;
                process.push_back(Process::all[i]->Instance(_io)); //create instances for all properties
            }



            /***********************************
            * Init Stars and binary params
            ************************************/
            init(init_params);
            /************************************/
        }
        //Catch any kind of error, dangerous in general, but we use it just to clean the heap and then
        //the error is rethrown and handled somewhere else.
        catch (...){ //Catch all errors
            default_destructor();  //Clean the heap associated to class members
            throw; //Rethrown the exact same  error
        }


    }
    ////////////////

    ///METHODS///
    inline double getp(const size_t &id) const {return property[id]->get(this);}
    inline double getp_0(const size_t &id) const {return property[id]->get_0(this);}
    inline  std::string get_name() const {return name;}
    inline size_t get_ID() const {return ID;}
    inline Star *getstar(const size_t &id)  {return star[id];}
    inline Process *getprocess(const size_t &id) {return process[id];}
    inline const std::vector<Process*>& getprocesses() const {return process;}
    inline const std::vector<double> & getstate() {return state;}
    /**
     * Method that return the Timestep of the last successful evolution step
     * @return the Timestep of the last successful evolution step in Myr
     */
    inline double get_last_Timestep() const {return last_Timestep;}

    /**
     * Return the id  of the  more evoveld star
     * @return Cases:
     *              - If both are empty: 0
     *              - If one star is empty: the id of other star
     *              - If stars have the same Phase: the id with the star with larger plife
     *              - In the other cases: the id of the star with larger Phase
     */
    inline unsigned  int get_id_more_evolved_star() const{
        //If one of the star is empty return the other one
        if (star[0]->isempty and !star[1]->isempty)
            return 1;
        else if (star[1]->isempty)
            return 0;
        //If both are remnant the older is the more compact following the order HeWD,COWD,ONEWD,NS,BH
        else if (star[0]->isremnant and star[1]->isremnant and star[0]->getp(RemnantType::ID)!=star[1]->getp(RemnantType::ID))
            return star[0]->getp(RemnantType::ID)>star[1]->getp(RemnantType::ID) ? 0 : 1;
        //If both are remnant of the same type, the more evolved is the one with the largest mass
        else if (star[0]->isremnant and star[1]->isremnant)
            return star[0]->getp(Mass::ID)>star[1]->getp(Mass::ID) ? 0 : 1;
        //In the stars have the same phase (not remnant handled before), return the more evolved is the one with the larger plife
        else if (star[0]->getp(Phase::ID)==star[1]->getp(Phase::ID))
            return star[0]->plife()>=star[1]->plife() ? 0 : 1;
        //Finally, if the stars are not empty, not remnant and don't have the same phase
        //the more evolved is the one with the largest phase.
        else
            return star[0]->getp(Phase::ID)>=star[1]->getp(Phase::ID) ? 0 : 1;
    }

    /**
     * Return t the  more evoveld star
     * @return Cases:
     *              - If both are empty: star[0]
     *              - If one star is empty: the other star
     *              - If stars have the same Phase: the star with larger plife
     *              - In the other cases: star with larger Phase
     */
    inline Star *get_more_evolved_star() {

        unsigned int _id=get_id_more_evolved_star();

        return getstar(_id);
    }

    /**
     * Estimate the effective radius during a RLO, i.e. Max(Rc,RL) for a given star.
     * @param starID  Id of the star (0 or 1)
     * @return Radx=max(Rc,RL)
     * Notice: the stellar properties are the one at the beginning of the timestep, the Rc is the CO radius
     * for nakedhelium stars.
     */
    double Radx(size_t starID);



    ///Get the value of the parameters store in the svpar class in io
    inline double get_svpar_num(std::string name) { return io->svpar.get_num(name);};
    inline std::string get_svpar_str(std::string name) { return io->svpar.get_str(name);};
    inline bool get_svpar_bool(std::string name) { return io->svpar.get_bool(name);};
    inline unsigned long get_rseed(){
        if (rseed==0)
            svlog.critical("rseed has not been initialised",__FILE__,__LINE__,sevnstd::bse_error());
        return rseed;
    }

    //TODO Move to the src file
    virtual void evolve() {

        evolution_step_completed=false;

        svlog.debug(std::to_string(property[BTimestep::ID]->get()));
        svlog.debug(std::to_string(property[BTimestep::ID]->get_0()));
        //utilities::wait(__PRETTY_FUNCTION__);


        //double twait=4.5668373678e+01;
        double twait=1E30;
        if (getp(BWorldtime::ID)>twait){
            svlog.debug("S0 T= " + utilities::n2s(star[0]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[0]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[0]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("S1 T= " + utilities::n2s(star[1]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[1]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[1]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("B T= " + utilities::n2s(getp(BTimestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(getp_0(BTimestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(getp(BWorldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            utilities::wait("Before synch_dt");
        }



        /*************************************
        * Synchronise proposed timesteps
        *************************************/
        synchronise_dt_star();



        if (getp(BWorldtime::ID)>twait){
            svlog.debug("S0 T= " + utilities::n2s(star[0]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[0]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[0]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("S1 T= " + utilities::n2s(star[1]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[1]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[1]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("B T= " + utilities::n2s(getp(BTimestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(getp_0(BTimestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(getp(BWorldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            utilities::wait("After synch_dt");
        }

        /*************************************
        * Single stellar evolution
        *************************************/
        //At this point timesteps are in sync
        star[0]->evolve();
        //timestep of star 0 may have changed at this point
        star[1]->evolve();
        //timestep of star 1 may have changed at this point



        //std::cout<<" 0: "<<star[0]->getp_0(Timestep::ID)<<" "<<star[0]->getp(Timestep::ID)<<std::endl;
        //std::cout<<" 1: "<<star[1]->getp_0(Timestep::ID)<<" "<<star[1]->getp(Timestep::ID)<<std::endl;

        if (getp(BWorldtime::ID)>twait){
            svlog.debug("S0 T= " + utilities::n2s(star[0]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[0]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[0]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("S1 T= " + utilities::n2s(star[1]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[1]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[1]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("B T= " + utilities::n2s(getp(BTimestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(getp_0(BTimestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(getp(BWorldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.pdebug("S0 type",star[0]->getp(Phase::ID),__FILE__,__LINE__);
            utilities::wait("After evolve star");
        }

        //Check if the the two stars have been evolved with the same timestep.
        //If not repeat the evolution of the star with the largest used timestep (V0).
        //Plus update the new timestep of the binary.
        check_and_sync_sse();



        if (getp(BWorldtime::ID)>twait){
            svlog.debug("S0 T= " + utilities::n2s(star[0]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[0]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[0]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("S1 T= " + utilities::n2s(star[1]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[1]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[1]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("B T= " + utilities::n2s(getp(BTimestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(getp_0(BTimestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(getp(BWorldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            utilities::wait("After check and synch");
        }


        /*************************************
        * Binary stellar evolution
        *************************************/
        //Here the timestep of the binary is equal to the one just done by the single stars (s->getp_0(Timestep::ID)).
        //So it is equal to one defined in synchronise_dt_star if the stars have not repeated the evoluion in evolve.
        //In that case it is equal to the actual timestep of the single stellar evolution
        evolve_binary();
        //NOT NEEDED ANYMORE
        //binary_evolution_status = evolve_binary();
        //Now if the the Binary properties have not been uptated, update in any case the BTimestep to be syncrhonised with stars Timestepgit
        //if (binary_evolution_status==utilities::BIN_EV_NOT_DONE){
         //   std::cout<<" WE "<<std::endl;
            //property[BTimestep::ID]->evolve(this);
        //    std::cout<<" WE2 "<<std::endl;
        //}

        if (getp(BWorldtime::ID)>twait){
            svlog.debug("S0 T= " + utilities::n2s(star[0]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[0]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[0]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("S1 T= " + utilities::n2s(star[1]->getp(Timestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(star[1]->getp_0(Timestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(star[1]->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            svlog.debug("B T= " + utilities::n2s(getp(BTimestep::ID),__FILE__,__LINE__) + " T0=" + utilities::n2s(getp_0(BTimestep::ID),__FILE__,__LINE__) + " WT=" + utilities::n2s(getp(BWorldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
            utilities::wait("After evolv binary");
        }



        //Check if the binary evolution has re-evolved the system with a smaller timestep
        //In that case just re-evolve the stars with this timestep and again the binary.
        //NB: This check assumes that the stars are already synchronised (check_and_sync_sse) and that
        //the binary has been evolved with the actual time step of the stellar evolution. Therefore,
        // the new proposed binary dt can be only smaller than the stellar dt. If it is larger
        check_and_sync_bse();
        //property[BTimestep::ID]->resynch(1e30);
        //At this point stars and binary processes have been evolved with the same timestep

        evolution_step_done();


        /****************************************************************
        * FINAL STUFF. All the updates that need to be done after a complete evolution step (all the repetitions already done) go here
        * Last checks, change tracks and check mass limits for remnants
        * Plus,  evolve the derived property so that they are updated considering processes and change of tracks
        *****************************************************************/
        //TODO The points below can be completely skipped if the binary is broken
        ///Check for new_tracks for stars and check mass_limits
        Star * this_star = star[0];
        Star * other_star = star[1];

        for (int i=0; i<2; i++){
            //First of all check if it has to explode
            //TODO this is a temporary solution, we have to think about how to handle SNIa
            if (this_star->is_exploding_as_SNI()){
                    //TODO Can we remove the difference between explode_as_SNI(binstar *b) and explode_as_SNI(), i.e. between Star::set_empty and Star::set_empty_in_bse
                    this_star->explode_as_SNI(this); //Use explode_as_SNI(Binstar* b) so that timestep_0 is reset to the current value
                    set_broken();
                    set_onesurvived();
                    //TODO, Here we are using a arbitrary process just to set the event, we should have a proper way to set the event without calling the preoceese
                    process[RocheLobe::ID]->set_event((double)Lookup::EventsList::SNIa);
            }
            check_and_set_QHE(this_star); //Check if we have to flag QHE evolution
            this_star->find_new_track(); //Check if we have to change tracks
            this_star->check_remnant_mass_limits(); //Check if we have to performe a Crem transition (e.g. from NS to BH)
            //If the Mass is extremely small set it to empty.
            if (this_star->getp(Mass::ID)<1e-3){this_star->set_empty();}

            //Update derived properties (basic properties could have been changes during binary processes)
            //TODO Now this check is made also inside binary_evolution. I have to check  if we can remove this call from here
            update_derived_properties_star(this_star);

            //update Maximum CO and minimum HE (only after the CO core has been formed)
            if (get_svpar_bool("ev_set_maxCO") and std::isnan(this_star->get_MCO_max()) and this_star->getp(Phase::ID) >= TerminalCoreHeBurning and !this_star->amiremnant()){
                this_star->set_MCO_max();
            }
            if (get_svpar_bool("ev_set_minHE") and std::isnan(this_star->get_MHE_min()) and this_star->getp(Phase::ID) >= TerminalCoreHeBurning and !this_star->amiremnant()){
                this_star->set_MHE_min();
            }
            utilities::swap_stars(this_star,other_star);
        }

        //Very last thing Binary event special evolve
        property[BEvent::ID]->special_evolve(this);

    }

    //TODO tO save space maybe we can flush the allstates data to the output file after a given number of stored elements.
    //TODO use a binary property NextOutput
    inline void recordstate_w_timeupdate() {

        ///Check that the evolution has been done
        evolution_guard(__FILE__,__LINE__);

        //TODO check that they are in sync
        star[0]->recordstate_w_timeupdate();
        star[1]->recordstate_w_timeupdate();

        recordstate_binary();

    }
    inline void recordstate() {

        //TODO check that they are in sync
        star[0]->recordstate();
        star[1]->recordstate();

        recordstate_binary();
    }
    inline void recordstate_binary(){

        for(size_t i = 0; i < io->printcolumns_binary.size(); i++){
            //Concerning the Timestep, put the used one not the proposed one
            if ((size_t)io->printIDs_binary[i]==BTimestep::ID)
                state[i] = getp_0((size_t)io->printIDs_binary[i]);
                //state[i] = star[0]->getp_0(Timestep::ID);
            else
                state[i] = getp((size_t)io->printIDs_binary[i]);

        }


        combinedstate = star[0]->getstate();
        combinedstate.insert(combinedstate.end(), star[1]->getstate().begin(), star[1]->getstate().end());
        combinedstate.insert(combinedstate.end(), state.begin(), state.end());

        allstates.push_back(combinedstate);

        //TODO The single recordstate of the two stars are now stored in combinedstate. Can we reset and free the memory of star.allstates?
        for (auto & cc: combinedstate)
            svlog.debug("Combined state: "+utilities::n2s(cc,__FILE__,__LINE__));
    }

    inline void set_broken(){
        broken = true;
        check_and_set_broken();
        //std::string w = utilities::log_print("BROKEN", this);
        //print_to_log(w);
    }
    inline void set_empty(){onesurvived=false; empty=true;}
    inline void set_onesurvived(){
        //set_onesurvived is triggered when a star in the system becomes empty:
        //if oneservived is already true, the binasry system is now empty
        //if empty the binary is empty and onesurvived has to be false.
        if (empty)
            onesurvived = false;
        //if true and not empty, not the binary system becomes empty and onesurvived is reset to false
        if (onesurvived){
            onesurvived = false;
            empty = true;
        }
        //If binary is not empty and onesurvived was false set onesurvived to true
        else
            onesurvived= true;
    }
    inline bool breaktrigger() const {


        //if (star[0]->isempty and star[1]->isempty) return true; //the system is empty
        if ( (broken or onesurvived) and break_at_broken) return true; //system is broken or just one object in the system
        else if (star[0]->isremnant and star[1]->isremnant and break_at_remnant) return true;
        else if ( (!star[0]->isremnant or !star[1]->isremnant)  and break_at_remnant ){
            return false;
        }
        else if(!break_at_remnant){
            return getp(BWorldtime::ID)>=get_tf() or  utilities::areEqual(getp(BWorldtime::ID),get_tf());
        }
        else{
            SevnLogging svlog_tmp;
            svlog_tmp.critical("Star breaktrigger check failed (Maybe Worldtimes are not synchronised?)",__FILE__,__LINE__);
        }


        return false;

    }

    /**
     * Check if a given process is actually ongoing
     * @param id Id of the process to check
     * @return true if the process is ongoing, false otherwise
     */
    inline bool is_process_ongoing(const size_t &id) const {return process[id]->is_process_ongoing();}

    ///Handle processes alarms
    /**
     * Check if the specific process is signaling an special_evolution_alarm
     * @param id Id of the process to check
     * @return true if the special evolution alarm is on, false otherwise
     */
    inline bool processes_alarm(const size_t &id) const {return process[id]->is_special_evolution_alarm_on();}

    /**
     * Check all the processes alarm
     * @return True if at least one of the processes has a special_evolution_alarm switched on, False otherwise
     */
    bool process_alarm_any() const;

    /**
     * Switch off all the processes alarms
     * @return
     */
    inline void reset_all_processes_alarms(){for (auto& proc : process) proc->special_evolution_alarm_switch_off();};

    ///Handle events
    /**
     * Reset all the events, i.e. set to -1 the event member of all the processes
     */
    inline void reset_events(){
        for (auto& proc : process)
            proc->reset_event();
    }
    /**
     * Check if one of the process is signaling an event.
     * Notice that the function use a priority schema to decide which one of the events will be returned.
     * The priority is given by the position of the event in Lookup::EventList, higher is the value higher is the priority
     * The order of the checked process follows the process initialisation in static_main.h
     * @return a process code event
     */
    inline double get_event_from_processes(){

        bool RLOB_triggered=false;
        bool Merger_triggered=false;
        bool Collision_triggered=false;
        bool CE_triggered=false;

        double event=-1.0;
        for (const auto& proc : process){
            double this_event=proc->get_event();
            if (this_event==EventsList::RLOBegin) RLOB_triggered=true;
            else  if (this_event==EventsList::Merger) Merger_triggered=true;
            else if (this_event==EventsList::Collision) Collision_triggered=true;
            else if (this_event==EventsList::CE) CE_triggered=true;
            else if (this_event==EventsList::CE_Merger){CE_triggered=true;Merger_triggered=true;}
            if (this_event>event){
                event=this_event;
            }
        }

        //Check composite events
        if (Collision_triggered and CE_triggered and Merger_triggered) event=EventsList::Collision_CE_Merger;
        else if (Collision_triggered and CE_triggered) event=EventsList::Collision_CE;
        else if (Collision_triggered and Merger_triggered) event=EventsList::Collision_Merger;
        else if (RLOB_triggered and CE_triggered and Merger_triggered) event=EventsList::RLOB_CE_Merger;
        else if (RLOB_triggered and CE_triggered) event=EventsList::RLOB_CE;
        else if (RLOB_triggered and Merger_triggered) event=EventsList::RLOB_Merger;
        else if (CE_triggered and Merger_triggered) event=EventsList::CE_Merger;

        return event;
    }

    ///Outputs

    //TODO, outoput stuff should be removed from the star and binary classes and handeld by a separate class
    /**
     * Check if time to store the outputs
     * @return
     */
    bool isoutputtime();

    inline bool printall() {return print_all_steps;}
    inline bool notprint() {return print_only_end;}
    /**
    * Print the current binary properties.
    * It is supposed that the evolution of this binary is not stopped by an error, therefore
    * the summary is printed in the evolved files
    * It uses the method print_evolved_summary and print_formatted_output of the class IO.h
    */
    void print(){
        io->print_log();
        io->print_evolved_summary(name, rseed, ID);
        io->print_formatted_output(allstates, name, rseed, ID,true);
    }
    /**
    * Print the current binary properties.
    * It is supposed that the evolution of this binary is  stopped by an error, therefore
    * the summary is printed in the failed files
    * It uses the method print_failed_summary and print_formatted_output of the class IO.h
    */

    /**
     * Print the current binary properties.
     * It is supposed that the evolution of this binary is  stopped by an error, therefore
     * the summary is printed in the failed files
     * It uses the method print_failed_summary and print_formatted_output of the class IO.h
     * @param include_in_output  If true flush allstates to the output even if the evolution failed
     */
    void print_failed(bool include_in_output=false){
        io->print_log();
        io->print_failed_summary(name, rseed, ID);
        if (include_in_output)
            io->print_formatted_output(allstates, name, rseed, ID,true);
    }

    /**
     * Print a message to the log file.
     * It is a wrapper of the method log_put of the class IO
     * @param message
     */
    void  inline print_to_log(std::string& message){
        io->log_put(message);
    }

    /**
    * Initialise the two stars in the binary and the other binary properties.
    * This functiom assume thta the Star constructor is called in the standard way from a init param vector with:
    * Mass, Z, Spin, sn_type, t initial, t final, dt print output.
    * @param params A vector of string with the following orderd parameters:
    * Mass1, Z1, Spin1, t_1 initial, Mass2, Z2, Spin2, t_2 initial, sn_type, t final, dt_out, bin_separation, bin_eccentricity
    */
    void init(const std::vector<std::string> &params);

    inline std::vector<double>  get_zams() {

        std::vector<double> zams{star[0]->get_zams(), star[1]->get_zams()};
        return zams;
    };

    inline std::vector<double>  get_Z() {

        std::vector<double> Z{star[0]->get_Z(), star[1]->get_Z()};
        return Z;
    };

    /**
     *  Set the final time of the simulation
     * @param a value to set
     * @param file
     * @param line
     */
    inline void set_tf(const double a, const char* file, const int line) {

        if(std::isnan(a) || std::isinf(a))
            svlog.critical("Tf set to INF or NAN", file, line);
        else if(a < 0.0)
            svlog.critical("Tf cannot be negative", file, line);
        else
            tf = a;
    }

    /**
     * Set the dtout time.
     * @param a value to set
     * @param file
     * @param line
     */
    inline void set_dtout(const double a, const char* file, const int line) {

        if(std::isnan(a) || std::isinf(a))
            svlog.critical("dtout set to INF or NAN", file, line);
        else if(a <0.0)
            svlog.critical("dtout cannot be negative", file, line);
        else
            dtout = a;
    }


    inline void evolution_step_done(){
        last_Timestep = getp_0(BTimestep::ID); //V0 in BTimestep contains the Timestep of the step currently done, V is the time predicted for the new Timestep
        evolution_step_completed=true;
    }
    inline void evolution_guard(const char *file_input = nullptr, int line_input = -1){
        if (!evolution_step_completed)
            svlog.critical("The function recordstate_w_timeupdate can be called only after an evolution step",file_input,line_input,sevnstd::bse_error());
    }
    inline double  get_tf() const { return tf;};

    //TODO We have to take into account that the two stars can have different phases times. For now we just disable the print_per_phase option
    inline double get_dtout()  {
        //if(print_per_phase)
        //    return dtout_phase[(size_t)properties[Phase::ID]->get()];   //TODO it can be or fixed or print all time steps or print all timesteps diveded by 2 (print every 2), or divided by 3 (print every 3) and so on...
        //else
        //    return dtout;
        return dtout;
    }

    inline std::string get_id_name() const {
        std::string mss = "ID: "+utilities::n2s(ID,__FILE__,__LINE__) + " Name: "+utilities::n2s(name,__FILE__,__LINE__);
        return mss;
    }

    /**
     * Check the outcome of the accretion on a compact object
     * @param donorID Id  of the donor star
     * @param accretorID  ID of the accretor star
     * @param DMaccreted  DM accreted on the  acrretor
     * @return utilisties::SNIA_EXPLODE if a SN explosion is triggered, SN_NOT_EXPLODE otherwise
     */
    utilities::sn_explosion check_accretion_on_compact(size_t donorID, size_t accretorID, double DMaccreted);

    /**
     * Call check_accretion_on_compact(size_t donorID, size_t accretorID, double DMaccreted)
     * checking both stars and taking DMaccreted from the processes
     */
    void check_accretion_on_compact(){

        Star *donor    = star[0];
        Star *accretor = star[1];
        double DM_accretor;

        for (size_t i=0; i<2; i++){

            ///Estimate the DM from the processes
            //Consider all the process, NB DM can be negative
            //Cycle over all the processes to get the total DM variation.
            DM_accretor=0;
            for (auto &proc : process){
                DM_accretor += proc->get_var(accretor->get_ID(), Mass::ID);
            }

            //Check only if it is really accreting
            if (DM_accretor>0)
                check_accretion_on_compact(donor->get_ID(), accretor->get_ID(), DM_accretor);


            utilities::swap_stars(donor, accretor);
        }
    }

    /**
     * Check if the binary system is broken or onesurvived is true or empty. In that case apply set_broken to the properties
     * @return true if the system is broken or if onesurvived=true or the system is empty, false otherwise.
     */
    bool check_and_set_broken(){
        //If broken or onesurvived or empty set properties to broken  and return true
        if (broken or onesurvived or empty) {
            for (auto &prop: property)
                prop->set_broken(this);
            return true;
        }

        return false;
    }


    void check_and_set_bavera_xspin();
    void check_if_applied_bavera_xspin();
    double SetXspinBavera(const double, const double) const; // Computes the Xspin following Bavera et al. (2021). Used only for the second BH spin in the case of a binary undergoing WR+BH --> BH+BH.
    ////////////////



protected:

    ///MEMBERS///
    double tf; /*!< End time of the simulation*/
    double dtout; /*!< Time interval for the output*/
    std::vector<double> state;  /*!<Vector containing the  properties of the binary (the one stored in IO::printcolumns_binary). It is filled with value in Binstar::record */
    std::vector<double> combinedstate; /*!<Vector containing the  properties of the stars (see Star::state) and the properties of the binary (see above) */
    std::vector<std::vector<double>> allstates; /*!<  2D Vector storing all the combinedstates. Each time Binstar::record is called combinedstate is pushed to allstates*/
    ////////////////

    ///METHODS///
    ///SYNCS
    /**
     * Set the time step of the next evolution of the two stars to the minimum values between the current timestep (Timestep::V).
     */
    void synchronise_dt_star();

    /**
    * Check if the last single stellar evolution step has been made with the same timestep.
    * If not evolve again the star with the smallest timestep used in the last evolutions (T0)..
    */
    void check_and_sync_sse();


    /**
    * Check if the last single stellar evolution step has been made with the same timestep.
    * If not evolve again the star with the smallest timestep used in the last evolutions (T0)..
    */
    void check_and_sync_bse();

    ///CHECKS AND CORRECT
    /**
     * Check if after the evolution the Mass transfer (mainly Roche Lobe) causes the total Mass to
     * be lower than the core Mass, considering both the He and CO Mass.
     * If the Roche Lobe radius does not include the Core (He or CO) we force the star to lose only the envelope
     * creating a Naked He star, the radius, the accreted and donated mass from the binary evolution are corrected accordingly.
     * If the Roche Lobe radius is inside the He and/or CO core, the stars is allowed to transfer mass from the core
     * and the MCO and/or MHE is updated to be equal to Mtot. If Mtot<MCO and RCO<RL<RHE a nakedCO star is created
     * and the star is forced to lose at most the envelope and He layer from the core.
     *
     */
    void check_nakedHe_or_nakedCO_after_binary_evolution();

    /**
     * Check if after the binary evolution the angular momentum, of the two stars is larger
     * than the critical angular momentum Lcrit=IOmega_crit.
     * In case set L=Lcrit.
     * If L has been updated, call also update_derived_properties_star  to update derived properties that can depend on L
     */
    void check_AngMomSpin_after_binary_evolution();

    /**
    * Function that limit the mass transfer to a given mass value. For example we can limit the mass transfer to the core mass,
    * in that case the maximum amount of mass we can trasnfer is the Envelope mass. We can also limit the mass transfer to the CO core mass.
    * If the amount of mass transfered is larger, the function estimate the correction to exactly loss the right amount of mass and update the properties.
    * The accreted mass is corrected accordingly.
    * Notice to work this function assume that Total Mass < mass_limiting_property_id
    * @param donor  Pointer to the Star where we want to limit the mass transfer
    * @param accretor  Pointer to the Star that can accrete the mass lost from the donor
    * @param mass_limiting_property_id ID of the mass limiting property (usually MHE or MCO).
    */
    void limit_and_correct_mass_transfer_for_donor_from_binary(Star *donor, Star *accretor,
                                                               std::size_t mass_limiting_property_id);

    /**
     * Check if the condition for Quasi Homogeneous Evolution is triggered after RLO mass transfer.
     * This function follows the recipe in Sec. 2.2 of Eldridge&Stanway11 (https://arxiv.org/abs/1109.0288).
     * If the condition are met, the star is flagged as QHE (ami_following_QHE is set to true)
     * Notice this check function has to be used at the end of the binary evolution so that we are sure
     * the evolution is completed and no repetitions can be called.
     * @param accretor Pointer to a star that is assumed the accretor in the RLO mass transfer
     *          (the fucntion checks if the star is really the accretor or if the RLO is happening
     *          using the property dMcumul_RLO).
     * @return true if the condition is triggered, false otherwise
     */
    bool check_and_set_QHE(Star *accretor);

    /**
     * Call the evolve of the derived properties of a given star.
     * The derived properties are properties that depends on other propertis, therefore when the binary
     * evolution chanfe some of the fundamental properties we can call this function to properly set all the other derived properties.
     * @param s  Pointer to the star
     */
    void update_derived_properties_star(Star *s);

    ///RESET
    /**
     * Reset (put to false) the flag about special evolution: mix, comenv, is_swallowed.
     */
    inline void reset_evolution_flags(){mix=comenv=is_swallowed[0]=is_swallowed[1]=false;}

    ////////////////

protected:
    SevnLogging svlog;  /*!< log object to handle output messages*/


private:

    ///MEMBERS///
    vector<Process*> process;  /*!< List with all the pointers to the processes initialised in static_main.h*/
    vector<BinaryProperty*> property; /*!< List with all the pointers to the binary properties initialised in static_main.h*/
    Star *star[2]={nullptr, nullptr}; /*!< Pointers to the stars in the binary system, initiliased to nullptr*/
    IO *io; /*!< Pointer to the IO object*/
    size_t ID; /*!< ID of the binary (it depends on the position of this system in the list of binary*/
    std::vector<std::string> init_params; /*!< vector of string loaded from the list of binary*/
    std::string name;  /*!< Random name of the binary*/
    unsigned long rseed; /*!< Random seed*/
    ////////////////

    ///METHODS///
    utilities::bse_evolution evolve_binary();

    inline void resynch(const double &dt){

        for (auto &prop : property) prop->restore();
        property[BTimestep::ID]->resynch(dt);

    }

    /**
     * Method to clean the heap from all the heap allocated members of this class
     */
    void default_destructor();

    /**
     * Call the constructor of the two stars with parameters params_star1 and params_star2.
     * The parameters are Mass, Z, Spin, sn_type, t initial, t final, dt print output.
     * @param params_star1  Mass, Z, Spin, sn_type, t initial, t final, dt print output of the first star.
     * @param params_star2 Mass, Z, Spin, sn_type, t initial, t final, dt print output of the second star
     * @note The random seeds of the stars are set to the same random seed of the binary.
     */
    void call_stars_constructor(std::vector<std::string> &params_star1, std::vector<std::string> &params_star2);

    /**
     * Transform the parameter in input it the stardard form needed from the Star constructor and then call the stars constructor
     * NB: This is the version that I like more, but it is different with respect the
     * old SEVN V1 formalism (implemented in init_star_legacy)
     * The constructor is called in the standard way from a init param vector with:
     * Mass, Z, Spin, sn_type, t initial, t final, dt print output.
     * @param params A vector of string with the following orderd parameters:
     * Mass1, Z1, Spin1, sn_type1,  t_1 initial, Mass2, Z2, Spin2, sn_type2,  t_2 initial,  bin_separation, bin_eccentricity, t_final, dt_out
     */
    void init_stars(const std::vector<std::string> &params);

    /**
    * Transform the parameter in input in the standard format needed from the Star constructor and then call the stars constructor.
    * The input format is in the old SEVN1 style.
    * The constructor is called in the standard way from a init param vector with:
    * Mass, Z, Spin, sn_type, t initial, t final, dt print output.
    * @param params A vector of string with the following orderd parameters:
    * Mass1, Mass2, Z1, Z2, Spin1, Spin2, bin_separation, bin_eccentricity, t_final, t_ini, tstep,  sn_type1, sn_type2, dt_out
    * NB: tstep not used at the moment
    */
    void init_stars_legacy(const std::vector<std::string> &params);

    /**
     * Initialise the Semimajor checking the value
     * @param a_ini Initial value of the Semimajor axis, it should be larger than 0 or the functions throw a critical error
     */
    inline void set_initial_semimajor(double a_ini){
        if (a_ini<=0){
            svlog.critical("Error on binary initialisation (" + get_id_name() +"): Semimajor is " + utilities::n2s(a_ini, __FILE__, __LINE__) + ": out of range (Semimajor>0)", __FILE__, __LINE__,sevnstd::sevnio_error());
        }
        property[Semimajor::ID]->init(a_ini);
    }

    /**
     * Initialise the Semimajor checking the value
     * @param a_ini string with the initial semimajor axis, if the string cannot be transformed to a number the function throw an io error.
     * a_ini has to be larger than 0 or the functions throw a critical error
     */
    inline void set_initial_semimajor(std::string a_ini){

        double a_inid;
        try{
            a_inid = utilities::s2n<double>(a_ini, __FILE__, __LINE__);
        }
        catch(sevnstd::sevnio_error& e){
            svlog.critical("Impossible to transform the given initial Semimajor axis to a number ("+
            get_id_name()+"): its value is "+a_ini,__FILE__,__LINE__,sevnstd::sevnio_error());
        }
        set_initial_semimajor(a_inid);
    }


    /**
     * Initialise the Eccentricity checking the value
     * @param ecc_ini  Initial value of the eccentricity, it should be lower than 1 and no lower than 0
     */
    inline void set_initial_eccentricity(double ecc_ini){
        if (ecc_ini<0 or ecc_ini>=1){
            svlog.critical("Error on binary initialisation (" + get_id_name() +"): Eccentricity is " + utilities::n2s(ecc_ini, __FILE__, __LINE__) + ": out of range (0<=Eccentricity<1)", __FILE__, __LINE__,sevnstd::sevnio_error());
        }
        property[Eccentricity::ID]->init(ecc_ini);
    }

    /**
     * Initialise the Eccentricity checking the value
     * @param ecc_ini string with the initial eccentricity, if the string cannot be transformed to a number the function throw an io error.
     * The initial value of the eccentricity has to  be lower than 1 and no lower than 0
     */
    inline void set_initial_eccentricity(std::string ecc_ini){

        double ecc_inid;
        try{
            ecc_inid = utilities::s2n<double>(ecc_ini, __FILE__, __LINE__);
        }
        catch(sevnstd::sevnio_error& e){
            svlog.critical("Impossible to transform the given initial Eccentricity  to a number ("+
                           get_id_name()+"): its value is "+ecc_ini,__FILE__,__LINE__,sevnstd::sevnio_error());
        }
        set_initial_eccentricity(ecc_inid);
    }


    /**
     * Initialise the binary parameters.
     * @param params  A vector of string with the following ordered parameters:
     *          if binput_mode is legacy: Mass1, Mass2, Z1, Z2, Spin1, Spin2, bin_separation, bin_eccentricity, t_final, t_ini, tstep  sn_type1, sn_type2, dt_out
     *          if binput_mode is new: Mass1, Z1, Spin1, sn_type1,  t_1 initial, Mass2, Z2, Spin2, sn_type2,  t_2 initial,  bin_separation, bin_eccentricity, t_final, dt_out
     */
    inline void init_binary_properties(const std::vector<std::string> &params){

        int inchoice=inputmapbin.at(io->binput_mode).first;

        std::string a,ecc;

        if (inchoice==InputBinaryOption::_new){
            a=params[10];
            ecc=params[11];
        }
        else if (inchoice==InputBinaryOption::_legacy){
            a=params[6];
            ecc=params[7];
        }

        //Init semimajor
        set_initial_semimajor(a);

        //Init eccentricity
        set_initial_eccentricity(ecc);

    }

    /**
    * In this function we set the parameters that are already set in each star when we call init_stars or init_stars_legacy.
    * These parameters are: break_at_broken, break_at_remnant, print_all_steps and print_per_phase, tf, dtout.
    * The function just check that  the two stars store the same value for the parameters to set. It this is not true it throws a critical error.
    * @note this function should be called after the initialisation of the stars (init_stars or init_stars_legacy).
    */
    void  init_other_params();

    /**
     * Set time when the tf in input is equal to broken
     * @param tf
     */
    void set_break_at_broken(std::string &tf){
        if (tf=="broken"){
            break_at_broken=true;
            tf="14000"; //Set time to Hubble time
        }
        else{
            break_at_broken=false;
        }
    }

    inline void set_rseed(const unsigned long a, const char* file, const int line){

        if(std::isnan(a) || std::isinf(a))
            svlog.critical("Rseed set to INF or NAN", file, line);
        else if (a<=0)
            svlog.critical("Rseed cannot be negative or 0", file, line);
        else
            rseed = a;
    }

    /**
     * Check if params in input has the expected number of parameters as defined in INPUTMAPBIN and INPUTMAPBIN_PARAM (lookup_and_phases.cpp)
     * @param params params to pass to the init function.
     */
    void inline check_init_param(const std::vector<std::string> &params){


        int par_expected_size = inputmapbin.at(io->binput_mode).second; //Number of expected parameters given the binput_mode stored in io.
        //If rseed is provided we have an extra parameter that is the seed.
        if (io->rseed_provided())
            par_expected_size+=1;
        if ((int)params.size()!=par_expected_size){
            std::string par_exp=utilities::n2s(par_expected_size, __FILE__, __LINE__);
            std::string inparam_size=utilities::n2s(params.size(), __FILE__, __LINE__);
            std::string this_ID=utilities::n2s(ID, __FILE__, __LINE__);
            svlog.critical("Number of Binary params needs to be " +par_exp+", it is instead " + inparam_size + " (ID_bin: " + this_ID +"). Maybe you are using inputs containing the ran"
                                                                                                                                       "dom seed, but rseed parameter is false." , __FILE__, __LINE__);
        }
    }
    ////////////////

public:

    ///DENSTRUCTORS
    ~Binstar();

};



#endif //SEVN_BINSTAR_H
