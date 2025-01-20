//
// Created by Giuliano Iorio on 2020-12-26.
//

#ifndef SEVN_EVOLVE_H
#define SEVN_EVOLVE_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <star.h>
#include <binstar.h>
#include <IO.h>
#include <vector>
#include <random>
#include <errhand.h>
#include <sevnlog.h>
using sevnstd::SevnLogging;

//TODO With the new functor implementation, the evolve function can be passed as parameters as used for other options. TO be implemented
//TODO the evolve function should not be put in the SEVN libraries, so we should move it in the main folder

/**
 *
 * IDEA BEHIND THE EVOLVE FUNCTOR
 *
 * The Evolve Functors are used to handle the evolution, the record of properties and the output
 * of both stars and binaries. The evolve function of the classes Star and Binaries just evolves the system
 * of a single adaptive timestep forward, the Evolve fuctor described here have the role to drive the complete evolution
 * stopping it when needed, print in in output the results and handle errors.
 * The Base EvolveFunctor constructor has 4 parameters:
 *      - Reference to Sevnlogging istance: it is used to initialise the protected member _svlog
 *         svlog is used to return logs (debug, info, errors) in the class.
 *         Use always it to return infos to the standard output or raise errors.
 *      - bool record_state: it is used    to initialise the protected member _record_state.
 *          if this parameter is true, all the states record and print has to be disabled
 *      - bool include_failed> it is used to initialse the protected member _include_failed.
 *          if _include_failed is true, we have to print to the output files also for star with failed evolution
 *      - Options *option: pointer to the struct Options, it is used to initialise the protected member _options
 *          this struct containt a vector of doubles and a vector of string. This is an optional member (it set to nullptr by default)
 *          used to add other parameters for a given Functor.
 *
 *
 *  The key method of a Functor is the override of the operators(Binstar& binstar) and operators(Star& star).
 *  This is the function that perform the actual evolution. As a general guide the evolve should contain an infinite
 *  for that call star or binstar evolution.
 *
 *  ** CHECK STALLING **
 *  At beginning of each evolution step the check for stalling should be added through the static method check_stalling(system).
 *  System could be both a Star or a binary, for example
 *      check_stalling(binary) or check_stalling(star).
 *  Notice that the method is static so it could (and should be called) also for custom evolution functions:
 *      EvolveFunctor::check_stalling(binary) or EvolveFunctor::check_stalling(star).
 *  The method returns a sevnerror if the binary or the stars has been evolving for more a time set at runtime through
 *  the propertis check_stalling_time.
 *  The check can be disabled at runtime by using -check_stalling false
 *
 *  After the evolution it should be checked if the breaktrigger condition have been
 *  reached. In that case the for cycle is stopped through a break. Inside the for it should be checked also
 *  if the condition for record of the states has been reached. If this is the case we have to call the method
 *  recordstate_w_timeupdate(). The method records the states and update the property NextOutput. Notice that
 *  if the memmber _record_state is false, the states should be never recorded indepentely from the input
 *  and EvolveFunctor options.
 *  If we want to handle the errors without stopping the whole computation we have to use a try .. catch .. structure.
 *  The evolution is inserted inside the try part while inside the catch we have to set the member  evolution_results=EXIT_FAILURE;
 *  and a brak for the for cycle should be inserted.
 *  At the end of the function before we have to add the mandatory following lines:
 *      final_print(system); //Handle the output printing
 *      return _evolution_result; //Exit_SUccess (default) or EXIT_FAILURE if it has been set in the evolution (if the evolution failed).
 *
 *
 *
 *
 * HOW TO ADD A NEW EVOLVE FUNCTION
 *
 * The evolve function is a Functor that should have two mandatory override operator:
 * - inline int operator() (Binstar& binary)
 * - inline int operator() (Star& star)
 * The constructor need to call the base class constructor (that has parameters svlog and record_state).
 *
 * Implementation steps
 *    Assume we are implementing a new evolve functor EvolveNizzi
 *    1- Create a class derived from the the base class EvolveFunctor
 *          class EvolveNizzi : public EvolveFunctor
 *    2- Implement the class constructor that calls the base class constructor
 *
 *          EvolveNizzi(SevnLogging& svlog, bool record_state=true) : EvolveFunctor(svlog,record_state) {...}
 *          @param svlog, instance of class Svlogging
 *          @param record_state, If true enable the state recording during the evolution (following the parameters in input)
 *
 *    3- Override the pure virtual operator () for both stars and binstars
 *          - inline int operator() (Binstar& binary) override {....}
 *          - inline int operator() (Star& star) override {....}
 *          These are the function that will be actually called to perform the evolution of a star or a binary
 *          NOTICE: if the function is the same for both Star and Binstar it is better to implement
 *          a separate (protected of private) template function that will be called inside the operator.
 *          REMEMBER to add the check for stalling inside the evolution function at the beginning of each evolution step (see above):
 *              check_stalling(binary) or check_stalling(star).
 *          e.g.
 *          protected:
 *              template <typename System>
 *              inline int evolve_common(System& system){...}
 *
 *          and then
 *
 *          public:
 *                  inline int operator() (Binstar& binary) override  { return evolve_common(binary);}
 *                  inline int operator() (Star& star) override  { return evolve_common(star);}
 *
 *
 *          NOTICE: The functions have to return EXIT_SUCCESS if the evolution ended without problems or EXIT_FAILURE or
 *          raise an Error in the other cases
 *
 * WRITE DOCUMENTATION OF THE CLASS
 * The documentation of new EvolveFunctor should follow the schema:
 *      - A brief description of the Evolve purpose
 *      - @record_state_policy: Describe when the states are recorded and how the dtout in input is used  and interpreted
 *      - @break_evolve_policy: Describe when the evolution is halted and how the tf in input is used and interpreted
 *      - @Extra_options:  It the struct Options is used for extra input options, describe these extra parameters.
 *      - @binstar_evolve: If the evolve can be applied to binaries write Yes, otherwise give information (e.g. what kind of error is raised)
 *      - @star_evolve: If the evolve can be applied to stars write Yes, otherwise give information (e.g. what kind of error is raised)
 *      - @catch: Write what kind of errors (if any) is catched during the evolve sequence
 *
 * HOW TO USE IT IN THE MAIN FUNCTION
 *
 * - Create a pointer to  the given Functor
 *      e.g. evolve_utility::EvolveFunctor* evolve_function = new evolve_utility::EvolveNizzi(svlog,true);
 * - Use directly the function ( (*evolve_function)(binstar) ) or use it as parameter in the chunk_dispatcher:
 *      e.g. evolve_utility::chunk_dispatcher(evolve_function,Nchunk,sevnio,stars,true);
 * - If the function is stoared in the heap, remember to delete the pointer after the evolution
 *      e.g. delte evolve_function.
 *
 */



namespace evolve_utility{

    /**
     * Struct to pass Options to the EvolveFunctors
     * It contains two members:
     *      - A vector of doubles (num_parameters)
     *      - A vector of strings (str_paramters)
     */
    struct Options{
        std::vector<double> num_parameters;
        std::vector<std::string> str_parameters;
    };



    /**
     * Abstract Base class for Evolve functions.
     * The derived class needs to call the base constructor with parameters:
     * @param svlog, instance of class Svlogging
     * @param record_state, If true enable the state recording during the evolution (following the parameters in input)
     * @param include_failed, If true include the results of  failed evolutions in the output files
     * @param Options [Optional], Pointer to an istance of the struct evolve_utility::Options containing a vector of doubles and a vector of strings.
     * The struct is used to pass extra parameters to the given Evolve FUnctor (if needed)
     */
    class EvolveFunctor {
    public:
        explicit EvolveFunctor(SevnLogging& svlog, bool record_state=true, bool include_failed=false, Options* options=nullptr)
                : _svlog(svlog), _record_state(record_state), _include_failed(include_failed), _options{options} {
            static bool instance_created=false;

            //Check if this is the first isstance otherwise return an error
            if (instance_created) _svlog.critical("Only one instance of the EvolveFunctor is allowed",__FILE__,__LINE__,sevnstd::sanity_error());
            else instance_created=true;

        };
        virtual ~EvolveFunctor() = default;
        virtual inline int operator() (_UNUSED Binstar& binary){ return EXIT_SUCCESS;}
        virtual inline int operator() (_UNUSED Star& star){ return EXIT_SUCCESS;}
        virtual inline std::string name() { return "Pure Virtual Base EvolveFunctor";}

        /**
         * Check if the system has been evolving for more than a given threshold (standard 20 s)
         * If this is the case throw a sevnerror
         * The check is performed only if the option check_stalling is set to true
         * @TODO it could be better to move the check_stalling parameter from bool to double including the maximum time
         * allowed for the evolution
         * @tparam System An object of the class Star or Binary
         * @param system the system for which we want to check if the evolution is stalling
         * @return EXIT success if the system is not stalling otherwise throw an error
         * @error throw a sevnerror if the system has been evolving for more than the second set as threshold
         */
        template <typename System>
        static inline bool check_stalling(System& system) {
            static long dt_threshold = long(system.get_svpar_num("check_stalling_time")); //Threshold to check for stalling system in seconds
            if (system.get_svpar_bool("check_stalling")){

                thread_local static std::chrono::steady_clock::time_point clock_begin;
                thread_local static size_t _id{0};
                thread_local static std::string _name;



                //Set value and start new clock for a new system
                if (system.get_name()!=_name or  system.get_ID()!=_id){
                    clock_begin = std::chrono::steady_clock::now();
                    _id         = system.get_ID();
                    _name       = system.get_name();
                }

                if (elapsed_time(clock_begin)>dt_threshold){
                    SevnLogging svlog;
                    svlog.critical("System with name " + system.get_name() +
                                   " and ID " + utilities::n2s(system.get_ID(),__FILE__,__LINE__) +
                                   " evolved for more than " + utilities::n2s(dt_threshold,__FILE__,__LINE__) +
                                   "s this is not expected."
                                   "The evolution is skipped and the system flagged ad failed.",
                                   __FILE__,
                                   __LINE__,
                                   sevnstd::sevnerr(""));
                }

            }

            return EXIT_SUCCESS;
        }

    protected:
        SevnLogging& _svlog;
        const bool _record_state;
        const bool _include_failed;
        const Options* _options;

        template <class System>
        inline void final_print(System& system, bool evolution_result){

            if (evolution_result==EXIT_SUCCESS)
                system.print();
            else
                system.print_failed(_include_failed);
        }

        /**
         * Special error message printed to std::cerr.
         * This message is always print do not throw errors and can be used to print a message of a catched error during
         * during a system evolution
         * @param e Instance of Sevnerror
         * @param system_name name of the system (sytem.get_name())
         * @param system_ID ID of the system (system.get_ID())
         */
        template <class System>
        static void catched_error_message(System& system, const sevnstd::sevnerr &e){
            std::string errmedd= "Evolution of system "+system.get_name() + "(ID "+ std::to_string(system.get_ID())+
                                 ") broken by an error: "+e.what();
            std::cerr<<errmedd<<std::endl;
        }

        /**
         * Measure the elapsed time in seconds from a given starting clock
         * @param start_clock, an instance of std::chrono::steady_clock::time_point,
         * the elapsed time is estimated from this point
         * @return elapsed time in seconds (approximated to the nearest integer)
         */
        static inline long elapsed_time(std::chrono::steady_clock::time_point start_clock){
            return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-start_clock).count();
        }
    };

    /**
     * Evolve functor for debug purpose. It does not evolve the stars/binarys and print just the evolved files
     * @record_state_policy No records
     * @break_evolve_policy No break
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve Yes
     * @catch sevnerr
     */
    class EvolveDebug : public EvolveFunctor {
    public:
        explicit EvolveDebug(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveFunctor(svlog,record_state,include_failed) {};
        inline int operator() (Binstar& binary) override  {
            return evolve_debug(binary);
        }
        inline int operator() (Star& star) override {
            return evolve_debug(star);
        }
        inline std::string name() override { return "EvolveDefault";}

    private:

        template <typename System>
        inline int evolve_debug(System& system){

            bool evolution_result = EXIT_SUCCESS; //Reset to true
            final_print(system,evolution_result); //print out all the recorded states for the star
            /********************************/

            return evolution_result;
        }
    };


    /**
     * Default evolve.
     * Same evolve for both stars and binstars (see method evolve_default).
     * @record_state_policy Following the dtout option in input
     * dtout=end record only the first and last state of the binary
     * @break_evolve_policy The evolution is stopped when a break trigger is set to true (depening on the tf option in input).
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve Yes
     * @catch sevnerr
     */
    class EvolveDefault : public EvolveFunctor {
    public:
        explicit EvolveDefault(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveFunctor(svlog,record_state,include_failed) {};
        inline int operator() (Binstar& binary) override  {
            return evolve_default(binary);
        }
        inline int operator() (Star& star) override {
            return evolve_default(star);
        }
        inline std::string name() override { return "EvolveDefault";}

    private:
        template <typename System>
        inline int evolve_default(System& system){

            bool evolution_result = EXIT_SUCCESS; //Reset to true

            if (_record_state) system.recordstate(); //always record the initial state

            for(;;) {
                //TODO Myabe is better to use try..catch.. outside the for?
                try {
                    check_stalling(system); //Check for stalling systems
                    system.evolve();
                    if (system.breaktrigger()) //stopping condition for evolving the star
                        break; //go to evolve the next star

                    if ( (system.isoutputtime() or system.printall()) and _record_state) {
                        system.recordstate_w_timeupdate();
                    }
                }
                catch(sevnstd::sevnerr& e){

                    catched_error_message(system,e);

                    evolution_result=EXIT_FAILURE;
                    break;
                }
            }
            if (_record_state) {
                system.recordstate(); //always record the final state
            }


            //print out all the recorded states for the star
            final_print(system,evolution_result); //print out all the recorded states for the star
            /********************************/

            return evolution_result;
        }
    };

    /**
     * Evolve functor template that can be used to create Functor that stop (and print) at given condition
     * @record_state_policy The state is recorded only the first time the conditions  in the method stop_condition are fulfilled
     * @break_evolve_policy The evolution is stopped when the conditions  in the method stop_condition are fulfilled
     * or following the tf input option.
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve Yes
     * @catch sevnerr
     */
    class EvolveStopCondition : public EvolveFunctor{
    public:
        explicit EvolveStopCondition(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
        : EvolveFunctor(svlog,record_state,include_failed) {};
        virtual inline int operator() (Binstar& binary) override {


            bool evolution_result = EXIT_SUCCESS; //Reset to true

            for(;;) {

                try {
                    check_stalling(binary); //Check for stalling systems
                    ///Evolve
                    binary.evolve();

                    ///Check record state
                    //If special break condition is flagged and _record_state is true:
                    //  record_state, print it to the file and exit.
                    if (stop_condition(binary) and _record_state){
                        binary.recordstate(); //always record the final state
                        break; //Exit
                    }
                    //If special break condition or natural breaktrigger  (controlled by input options) just break
                    else if (stop_condition(binary) or binary.breaktrigger())
                        break;
                }
                catch(sevnstd::sevnerr& e){

                    catched_error_message(binary,e);

                    evolution_result = EXIT_FAILURE;
                    break;

                }
            }
            /********************************/
            final_print(binary,evolution_result); //print out all the recorded states for the star
            return evolution_result;
        }
        virtual inline int operator() (Star& star) override {


            bool evolution_result = EXIT_SUCCESS; //Reset to true

            for(;;) {

                try {
                    check_stalling(star); //Check for stalling systems
                    ///Evolve
                    star.evolve();

                    ///Check record state
                    //If special break condition is flagged and _record_state is true:
                    //  record_state, print it to the file and exit.
                    if (stop_condition(star) and _record_state){
                        star.recordstate(); //always record the final state
                        break; //Exit
                    }
                        //If special break condition or natural breaktrigger  (controlled by input options) just break
                    else if (stop_condition(star) or star.breaktrigger())
                        break;
                }
                catch(sevnstd::sevnerr& e){

                    catched_error_message(star,e);

                    evolution_result = EXIT_FAILURE;
                    break;

                }
            }
            /********************************/
            final_print(star,evolution_result); //print out all the recorded states for the star
            return evolution_result;
        }
        inline std::string name() override { return "EvolveStopCondition";}
    protected:
        virtual inline bool stop_condition(_UNUSED Binstar& binstar){
            return false;
        }

        virtual inline bool stop_condition(_UNUSED Star& star){
            return false;
        }
    };

    /**
     * Evolve functor template that can be used to create Functor that record (and print) only if  given conditions
     * are fulfilled.
     * @record_state_policy The state is recorded only when  the conditions  in the method record_condition are fulfilled
     * and if record_state is set to true. The information to record depends on the dtout option,
     * dtout=end record only the the first and the last time the condition are fulfilled
     * @break_evolve_policy The evolution is stopped when a break trigger is set to true (depening on the tf option in input).
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve Yes
     * @catch sevnerr
     */
    class EvolveRecordCondition : public  EvolveFunctor{
    public:
        explicit EvolveRecordCondition(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveFunctor(svlog,record_state,include_failed) {};
        virtual inline int operator() (Binstar& binary) override {


            bool evolution_result = EXIT_SUCCESS; //Reset to true
            unsigned int Record_counter=0;

            for(;;) {

                try {
                    check_stalling(binary); //Check for stalling systems
                    ///Evolve
                    binary.evolve();

                    //If record condition is true, update the BLC counter and record state (if allowed by other flags)
                    if (record_condition(binary)){
                        Record_counter++;
                        //Record only if it is time (or print only end and this is the start of the Record condition Record_counter=1) and if _record state is true
                        if ( ( (binary.isoutputtime() or binary.printall()) or Record_counter==1 ) and _record_state ) {
                            binary.recordstate_w_timeupdate();
                        }
                    }
                        //If record_condition is false, but Record_counter>0 we just exit from the record condition, so record state (if allowed by other flags) ans reset the Record_counter
                    else if(Record_counter>0){
                        //Always print the end
                        if (_record_state){
                            binary.recordstate_w_timeupdate();
                        }
                        //reset the Record_counter counter
                        Record_counter=0;
                    }

                    if (binary.breaktrigger()) //stopping condition for evolving the star
                        break; //go to evolve the next star

                }
                catch(sevnstd::sevnerr& e){


                    catched_error_message(binary,e);

                    evolution_result = EXIT_FAILURE;
                    break;

                }
            }
            /********************************/
            final_print(binary,evolution_result); //print out all the recorded states for the star
            return evolution_result;
        }
        virtual inline int operator() (Star& star) override {


            bool evolution_result = EXIT_SUCCESS; //Reset to true
            unsigned int Record_counter=0;

            for(;;) {

                try {
                    check_stalling(star); //Check for stalling systems
                    ///Evolve
                    star.evolve();

                    //If record condition is true, update the BLC counter and record state (if allowed by other flags)
                    if (record_condition(star)){
                        Record_counter++;
                        //Record only if it is time (or print only end and this is the start of the Record condition Record_counter=1) and if _record state is true
                        if ( ( (star.isoutputtime() or star.printall()) or Record_counter==1 ) and _record_state ) {
                            star.recordstate_w_timeupdate();
                        }
                    }
                        //If record_condition is false, but Record_counter>0 we just exit from the record condition, so record state (if allowed by other flags) ans reset the Record_counter
                    else if(Record_counter>0){
                        //Always print the end
                        if (_record_state){
                            star.recordstate_w_timeupdate();
                        }
                        //reset the Record_counter counter
                        Record_counter=0;
                    }

                    if (star.breaktrigger()) //stopping condition for evolving the star
                        break; //go to evolve the next star

                }
                catch(sevnstd::sevnerr& e){

                    catched_error_message(star,e);

                    evolution_result = EXIT_FAILURE;
                    break;

                }
            }
            /********************************/
            final_print(star,evolution_result); //print out all the recorded states for the star
            return evolution_result;
        }
        inline std::string name() override { return "EvolveRecordCondition";}
    protected:


        virtual inline bool record_condition(_UNUSED Binstar& binstar){
            return false;
        }

        virtual inline bool record_condition(_UNUSED Star& star){
            return false;
        }
    };


    /**
     * Evolve functor specialised for the analysis of compact object binaries.
     * @record_state_policy The state is recorded only the first time the system becomes a compat object  binary
     * dtout is not used
     * @break_evolve_policy The evolution is stopped when both objects become compact and the stars are still bound
     * or following the tf input option.
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve raise Not implemented error
     * @catch sevnerr
     */
    class EvolveBinaryCompact : public EvolveStopCondition{
    public:
        explicit EvolveBinaryCompact(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
        : EvolveStopCondition(svlog,record_state,include_failed) {};
        inline int operator() (_UNUSED Star& star) override {
            _svlog.critical("Evolve for star is not implemented for functor EvolveBinaryCompact",__FILE__,__LINE__,
                            sevnstd::notimplemented_error());
            return -1;
        }

        inline std::string name() override { return "EvolveBinaryCompact";}

    protected:
        inline bool stop_condition(Binstar& binstar) override {
            return binstar.getstar(0)->amiCompact() and binstar.getstar(1)->amiCompact() and !binstar.broken;
        }

    };

    /**
     * Evolve functor specialised for the analysis of compact object binaries.
     * @record_state_policy The state is recorded only the first time the system becomes a compat object  binary
     * dtout is not used
     * @break_evolve_policy The evolution is stopped when both objects become compact and the stars are still bound
     * or following the tf input option.
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve raise Not implemented error
     * @catch sevnerr
     */
     //TODO Deprecated to be removed, new implementation is EvolveBinaryCompact
    class EvolveBinaryCompactOld : public EvolveFunctor {
    public:
        explicit EvolveBinaryCompactOld(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveFunctor(svlog,record_state,include_failed) {};
        inline int operator() (Binstar& binary) override {

            bool special_break_condition = false;
            bool evolution_result = EXIT_SUCCESS; //Reset to true

            for(;;) {

                try {
                    check_stalling(binary); //Check for stalling systems

                    ///Evolve
                    binary.evolve();

                    ///The special break condition  is activated the first time both stars are compact and the system is not broken
                    special_break_condition  = binary.getstar(0)->amiCompact() and binary.getstar(1)->amiCompact() and !binary.broken;

                    ///Check record state
                    //If special break condition is flagged and _record_state is true:
                    //  record_state, print it to the file and exit.
                    if (special_break_condition and _record_state){
                        binary.recordstate(); //always record the final state
                        break; //Exit
                    }
                    //If special break condition or natural breaktrigger  (controlled by input options) just break
                    else if (special_break_condition or binary.breaktrigger())
                        break;
                }
                catch(sevnstd::sevnerr& e){

                    catched_error_message(binary,e);

                    evolution_result = EXIT_FAILURE;
                    break;

                }
            }
            /********************************/
            final_print(binary,evolution_result); //print out all the recorded states for the star
            return evolution_result;
        }
        inline int operator() (_UNUSED Star& star) override {
            _svlog.critical("Evolve for star is not implemented for functor EvolveBinaryCompact",__FILE__,__LINE__,
                            sevnstd::notimplemented_error());
            return -1;
        }

        inline std::string name() override { return "EvolveBinaryCompact";}

    };

    /**
     * Evolve functor specialised for the analysis of  BH  binaries.
     * During the evolution the only output is the state when the system becomes a BH  binary
     * @record_state_policy The state is recorded only the first time the system becomes a BH  binary
     * dtout is not used
     * @break_evolve_policy The evolution is stopped when both objects become BHs and the stars are still bound
     * or following the tf input option.
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve raise Not implemented error
     * @catch sevnerr
     */
    class EvolveBBH :  public EvolveStopCondition{
    public:
        explicit EvolveBBH(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
        : EvolveStopCondition(svlog,record_state,include_failed) {};
        inline int operator() (_UNUSED Star& star) override {
            _svlog.critical("Evolve for star is not implemented for functor EvolveBinaryCompact",__FILE__,__LINE__,
                            sevnstd::notimplemented_error());
            return -1;
        }

        inline std::string name() override { return "EvolveBBH";}

    protected:
        inline bool stop_condition(Binstar& binstar) override {
            return binstar.getstar(0)->amiBH() and binstar.getstar(1)->amiBH() and !binstar.broken;
        }

    };


    /**
     * Evolve functor specialised for the analysis of stages of binary lifetime with a compact object plus a luminous variable (BLC).
     * @record_state_policy The states start to be recorded when the binary enters the first time in a BLC phase (bound compact + Luminous star)
     * The last record is got the first time the binary exits the BLC phase (e.g. the Luminous star becomes a compact object or a merger occurs).
     * Inside this window the state are recorded following the dtout option. When dtout=end only the first (start of the BLC) and last (end of the BLC)
     * states  are recorded.
     * @break_evolve_policy TThe evolution is broken when a broken trigger is true (depending on the input tf) or after we conclude a BLC phase
     * since there is no way that a new BLC phase can occur.
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve raise Not implemented error
     * @catch sevnerr
     */
    class EvolveBLC : public EvolveFunctor {
    public:
        explicit EvolveBLC(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveFunctor(svlog,record_state,include_failed) {};
        inline int operator() (Binstar& binary) override {

            bool BLC_flag = false, BLC_case1=false, BLC_case2=false;
            //The BLC counter counts the evolution steps the systems has been in the BLC configuration
            //It is increased by one if the BLC_flag is zero, while it is reset to zero when BLC_flag is false
            unsigned int BLC_counter=0;
            bool evolution_result = EXIT_SUCCESS; //Reset to true

            for(;;) {
                check_stalling(binary); //Check for stalling systems

                try {

                    ///Evolve
                    binary.evolve();

                    BLC_case1 = binary.getstar(0)->amiCompact() and !binary.getstar(1)->amiCompact();
                    BLC_case2 = binary.getstar(1)->amiCompact() and !binary.getstar(0)->amiCompact();

                    ///The special break condition  is activated the first time one star is compact and the other not and the system is not broken
                    BLC_flag  = (BLC_case1 or BLC_case2) and !binary.broken;

                    //If BLC_flag is true, update the BLC counter and record state (if allowed by other flags)
                    if (BLC_flag){
                        BLC_counter++;
                        //Record only if it is time (or print only end and this is the start of the BLC conditions BLC_counter=) and if _record state is true
                        if ( ( (binary.isoutputtime() or binary.printall()) or BLC_counter==1 ) and _record_state ) {
                            binary.recordstate_w_timeupdate();
                        }
                    }
                    //If BLC_flag is false, but BLC_counter>0 we are at the end of a BLC phase, so record state (if allowed by other flags) ans reset the BLC counter
                    //FInally break, in fact if this system was in a BLC phase and not it is not, there is no way it will enter again in a BLC phase
                    else if(BLC_counter>0){

                        //Always print the end
                        if (_record_state){
                            binary.recordstate_w_timeupdate();
                        }
                        //reset the BLC counter
                        BLC_counter=0;
                        break;
                    }

                    if (binary.breaktrigger()) //stopping condition for evolving the star
                        break; //go to evolve the next star

                }
                catch(sevnstd::sevnerr& e){

                    catched_error_message(binary,e);

                    evolution_result=EXIT_FAILURE;
                    break;
                }
            }
            /********************************/
            final_print(binary,evolution_result); //print out all the recorded states for the star
            return evolution_result;
        }
        inline int operator() (_UNUSED Star& star) override {
            _svlog.critical("Evolve for star is not implemented for functor EvolveBinaryCompact",__FILE__,__LINE__,
                            sevnstd::notimplemented_error());
            return -1;
        }
        inline std::string name() override { return "EvolveBLC";}


    };




    /**
     * Evolve functor specialised for analysis of the Cygnus X-1 X-ray binaries based on the constrain by https://arxiv.org/pdf/2102.09091.pdf
     * @record_state_policy The state is recorded only the first time the system becomes BH - Main sequence system
     * with some other requirements (see the stop_condition method). dtout is not used
     * @break_evolve_policy The evolution is stopped when the conditions  in the method stop_condition are fulfilled
     * or following the tf input option.
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve raise Not implemented error
     * @catch sevnerr
     */
    class EvolveBCX1 : public EvolveStopCondition {
    public:
        explicit EvolveBCX1(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveStopCondition(svlog,record_state,include_failed) {};
        inline int operator() (_UNUSED Star& star) override {
            _svlog.critical("Evolve for star is not implemented for functor EvolveBinaryCompact",__FILE__,__LINE__,
                            sevnstd::notimplemented_error());
            return -1;
        }

        inline std::string name() override { return "EvolveBCX1";}

    protected:
        inline bool stop_condition(Binstar& binstar) override {
            bool BH_cond = binstar.getstar(0)->amiBH() and binstar.getstar(0)->getp(Mass::ID)>=14 and binstar.getstar(0)->getp(Mass::ID)<=24;
            bool  MS_cond = binstar.getstar(1)->getp(Phase::ID)==1 and binstar.getstar(1)->getp(Mass::ID)>=30 and binstar.getstar(1)->getp(Mass::ID)<=50;
            bool Semimajor_cond = binstar.getp(Semimajor::ID)<=1000000000000.0;

            return BH_cond and MS_cond and Semimajor_cond;
        }

    };

    /**
     * Evolve functor specialised for analysis of the BH-MS objects by https://ui.adsabs.harvard.edu/abs/2022NatAs...6.1085S/abstract
     * @record_state_policy The state is recorded only the first time the system becomes BH - Main sequence system
     * with some other requirements (see the stop_condition method). dtout is not used
     * @break_evolve_policy The evolution is stopped when the conditions  in the method stop_condition are fulfilled
     * or following the tf input option.
     * @Extra_options Not used
     * @binstar_evolve Yes
     * @star_evolve raise Not implemented error
     * @catch sevnerr
     */
    class EvolveBHMSTarantula : public EvolveStopCondition {
    public:
        explicit EvolveBHMSTarantula(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveStopCondition(svlog,record_state,include_failed) {};
        inline int operator() (_UNUSED Star& star) override {
            _svlog.critical("Evolve for star is not implemented for functor EvolveBinaryCompact",__FILE__,__LINE__,
                            sevnstd::notimplemented_error());
            return -1;
        }

        inline std::string name() override { return "EvolveBHMSTarantula";}

    protected:
        inline bool stop_condition(Binstar& binstar) override {

            constexpr double  yr_to_day = 265.24;

            bool BH_cond = binstar.getstar(0)->amiBH() and binstar.getstar(0)->getp(Mass::ID)>=6 and binstar.getstar(0)->getp(Mass::ID)<=16;
            bool MS_cond = binstar.getstar(1)->getp(Phase::ID)==1 and binstar.getstar(1)->getp(Mass::ID)>=15 and binstar.getstar(1)->getp(Mass::ID)<=32;
            bool P_cond  = binstar.getp(Period::ID)*yr_to_day < 100;
            bool e_cond  = binstar.getp(Eccentricity::ID) < 0.3;

            return BH_cond and MS_cond and P_cond and e_cond;
        }

    };

    /**
    * Evolve functor specialised for analysis of BHMS
    * @record_state_policy The state is recorded only when  the following conditions are fulfilled:
    *     - System composed by a  BH  and a star (no limitations on the phase)
    *     - The system is not broken
    * In addition record_state has to be set to true.
    * The information to record depends on the dtout option,
    * dtout=end record only the the first and the last time the condition are fulfilled
    * @break_evolve_policy The evolution is stopped when a break trigger is set to true (depending on the tf option in input).
    * @Extra_options Not used
    * @binstar_evolve Yes
    * @star_evolve No
    * @catch sevnerr
    */
    class EvolveBHMS : public EvolveRecordCondition {
    public:
        explicit EvolveBHMS(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveRecordCondition(svlog,record_state,include_failed) {};
        inline int operator() (_UNUSED Star& star) override {
            _svlog.critical("Evolve for star is not implemented for functor EvolveBinaryCompact",__FILE__,__LINE__,
                            sevnstd::notimplemented_error());
            return -1;
        }

        inline std::string name() override { return "EvolveBHMSTarantula";}

    protected:
        inline bool record_condition(Binstar& binstar) override {
            bool cond_1      = binstar.getstar(0)->amiBH() and (binstar.getstar(1)->getp(Phase::ID)<7);
            bool cond_2      = binstar.getstar(1)->amiBH() and (binstar.getstar(0)->getp(Phase::ID)<7);
            bool binary_cond = !binstar.broken;

            return (cond_1 or cond_2) and binary_cond;
        }

    };


    /**
    * Evolve functor specialised for analysis of X-ray pulsator found in the cluster Westerlund1 (private communication)
    * @record_state_policy The state is recorded only when  the following conditions are fulfilled:
    *     - System composed by a compact object (BH/NS) and a star (no limitations on the phase)
    *     - The system is bound and the period is within 1 - 5  hours (1.14e-4 - 5.7e-4 years)
    * In addition record_state has to be set to true.
    * The information to record depends on the dtout option,
    * dtout=end record only the the first and the last time the condition are fulfilled
    * @break_evolve_policy The evolution is stopped when a break trigger is set to true (depending on the tf option in input).
    * @Extra_options Not used
    * @binstar_evolve Yes
    * @star_evolve No
    * @catch sevnerr
    */
    class EvolveW1 : public EvolveRecordCondition {

    public:
        explicit EvolveW1(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveRecordCondition(svlog,record_state,include_failed) {};
        inline int operator() (_UNUSED Star& star) override {
            _svlog.critical("Evolve for star is not implemented for functor EvolveBinaryCompact",__FILE__,__LINE__,
                            sevnstd::notimplemented_error());
            return -1;
        }

        inline std::string name() override { return "EvolveW1";}

    protected:
        inline bool record_condition(Binstar& binstar) override {
            bool BHNS_cond = (binstar.getstar(0)->amiCompact() and !binstar.getstar(1)->amiCompact()) or (binstar.getstar(1)->amiCompact() and !binstar.getstar(0)->amiCompact());
            bool Period_cond = !binstar.broken and binstar.getp(Period::ID)>=Pmin and binstar.getp(Period::ID)<=Pmax;

            return BHNS_cond and Period_cond;
        }
    private:
        static constexpr double hours_to_yr = 1.0/8766.0; //1 hours to year, considering 1 year=365.25 days
        static constexpr double Pmin = 1*hours_to_yr; //Minimum Period 1h
        static constexpr double Pmax = 5*hours_to_yr; //Maximum Period 5h

    };


    /**
    * Evolve functor specialised for retrive the output of only the stars crossing the instability strip
    * with the right condtions to become RRL pulsators.
    * @record_state_policy The state is recorded only when the method amiRRL of a  star (or at least one star in binary) return true
    * @break_evolve_policy The evolution is stopped when a break trigger is set to true (depending on the tf option in input).
    * @Extra_options Not used
    * @binstar_evolve Yes
    * @star_evolve Yes
    * @catch sevnerr
    */
    class EvolveRRL : public EvolveRecordCondition {

    public:
        explicit EvolveRRL(SevnLogging& svlog, bool record_state=true, bool include_failed=false)
                : EvolveRecordCondition(svlog,record_state,include_failed) {};
        inline std::string name() override { return "EvolveRRL";}

    protected:
        inline bool record_condition(Binstar& binstar) override {
            return binstar.getstar(0)->amiRRL() or binstar.getstar(1)->amiRRL();
        }
        inline bool record_condition(Star& star) override {
            return star.amiRRL() or star.amiRRL();
        }
    };


    /**
     * Evolve a binary system DEPREACATED
     * @param binary Binary system to evolve
     * @param svlog Instance of the Sevenlogging class
     * @param record_state If true record and print the binary states (accordingly with the option parameters)
     * @return EXIT_SUCCESS if the evolution is ended without errors. EXIT_FAILURE is sevnerr is raised.
     * In that case the error is catched and an error message is printed in the stardard error output (but the run is not halted)
     */
    inline int evolve_single(Binstar& binary, SevnLogging& svlog, bool record_state=true){

        if (record_state) binary.recordstate(); //always record the initial state
        for(;;) {

            try {
                EvolveFunctor::check_stalling(binary); //Check stalling
                binary.evolve();


                if (binary.breaktrigger()) //stopping condition for evolving the star
                    break; //go to evolve the next star

                if ( (binary.isoutputtime() or binary.printall()) and record_state ) {
                    binary.recordstate_w_timeupdate();
                }
            }
            catch(sevnstd::sevnerr& e){


                svlog.error("Evolution of binary "+binary.get_name() + "(ID "+ utilities::n2s(binary.get_ID(),__FILE__,__LINE__)+
                            ") broken by an error: "+e.what(),__FILE__,__LINE__,false);

                //std::cerr<< binaries[i].get_name()  <<" " <<e.what()<<std::endl << std::flush;
                return EXIT_FAILURE;
            }

            //utilities::wait();


        }
        if (record_state) {
            binary.recordstate(); //always record the final state
            binary.print(); //print out all the recorded states for the star
        }
        /********************************/

        return EXIT_SUCCESS;
    }

    /**
     * Evolve a star
     * @param star Star system to evolve
     * @param svlog  Instance of the Sevenlogging class
     * @param record_state If true record and print the binary states (accordingly with the option parameters)
     * @return EXIT_SUCCESS if the evolution is ended without errors. EXIT_FAILURE is sevnerr is raised.
     * In that case the error is catched and an error message is printed in the stardard error output (but the run is not halted)
     */
    inline int evolve_single(Star& star, SevnLogging& svlog, bool record_state=true){

        if (record_state) star.recordstate(); //always record the initial state
        for(;;) {

            try {

                EvolveFunctor::check_stalling(star); //Check stalling
                star.evolve();


                if (star.breaktrigger()) //stopping condition for evolving the star
                    break; //go to evolve the next star

                if ( (star.isoutputtime() or star.printall()) and record_state ) {

                    star.recordstate_w_timeupdate();
                }
            }
            catch(sevnstd::sevnerr& e){

                svlog.error("Evolution of binary "+star.get_name() + "(ID "+ utilities::n2s(star.get_ID(),__FILE__,__LINE__)+
                            ") broken by an error: "+e.what(),__FILE__,__LINE__,false);

                return EXIT_FAILURE;
            }

            //utilities::wait();


        }
        if (record_state) {
            star.recordstate(); //always record the final state
            star.print(); //print out all the recorded states for the star
        }
        /********************************/

        return EXIT_SUCCESS;
    }


    /**
     * Evolve a list of stars or binstars
     * @tparam System  It could be Binstar or Star
     * @param evolve_function Pointer to functor derived from  base class EvolveFunctor
     * @param systems Vector containing the systems to evolve.
     * @param sevnio Instance of the IO class (the one linked to the binaries)
     * @param Nevolve Number of systems to evolve (first Nevolve). If -1 evolve all the sytems
     * @return Number of failed evolutions
     */
    template <typename System>
    inline int evolve_list(EvolveFunctor* evolve_function, std::vector<System>& systems, _UNUSED IO& sevnio, int Nevolve=-1){
        SevnLogging svlog;

        if (Nevolve==-1) Nevolve=systems.size();
        unsigned Nfailed =0;

#ifdef _OPENMP
#pragma omp parallel num_threads(sevnio.nthreads)
#endif
        {

#ifdef _OPENMP
#pragma omp for schedule(static)  reduction(+: Nfailed)
#endif
            for (size_t i = 0; i < (size_t)Nevolve; i++) {

                //NOTICE: THE FOLLOWING IS EXTREMELY IMPORTANT
                //mtrand is a static threadlocal variaible, i.e. it is a global variable for each thread
                //it is "seeded" each time we initiliase a system (binary or star).
                // In order to gurantee reproducibility (indepently of the number of threads and the position
                // of the system in the list), before to evolve the star we "re-seed"  mtrand withe the sytem random seed
                //Notice also that it is important that mtrand IS NEVER RE-SEEDED during the stellar/binary evolution
                //(for this reason when we create an auxiliary star we save the current mtrand status, then we initiliase the
                //auxiliary star where mtrand is reseed and finally we restore the old mtrand).
                //Notice that this in principle could be a problem when using SEVN as a library (see below),
                //but it is not atm a problem for the SEVN python wrapper since each time we call evolve,
                //the star is re-initialised.
                //TODO We should surely improve this behaviour (see issue 71 on gitlab: https://gitlab.com/sevn/SEVN/-/issues/71 )
                utilities::mtrand.seed(systems[i].get_rseed());//Set random state for reproducibility, DO NOT MOVE OR REMOVE THIS LINE


                /****** EVOLVE  ******/
                if((*evolve_function)(systems[i])==EXIT_FAILURE)
                    Nfailed++;


            }


        }

        return Nfailed;
    }


    /**
     * Evolve using chunk
     * @tparam T  It could be Binstar or Star
     * @param Nchunk  Number of systems to evolve in each chunk
     * @param sevnio Instance of the IO class (the one linked to the binaries)
     * @param systems Vector containing the systems to evolve.
     * @param record_state  If true record and print the systems states (accordingly with the option parameters)
     * @param progress If true print progress information to the standard output
     * @return  Number of failed evolutions
     * @Note  The vector of systems is passed by reference and it is cleared if not empty yet.
     */
    //[[deprecated("Replaced by functor evolve implementation")]] This attribute can be used only from C++14
    template<typename T>
    inline int chunk_dispatcher(unsigned int Nchunk, IO& sevnio, std::vector<T>& systems, bool record_state=true, bool progress=true){

        ///Preliminary reset
        //If systems is not an empty vector clear it
        if (!systems.empty())
            systems.clear();

        ///Preliminary assignment
        //Max size
        unsigned int Ntot = sevnio.STARS_MATRIX.size();
        //Reset Nchunk if needed
        Nchunk = std::min(Nchunk,Ntot);
        //Reserve space for systems
        systems.reserve(Nchunk);


        ///Cycle
        if (progress)
            std::cout<<"Evolving systems:"<<std::endl;

        unsigned int Ndone=0;
        size_t current_idx=0;
        int Ntodo=0;
        int Nfailed=0;
        while (Ndone<Ntot){


            //Assign Ntodo
            Ntodo = std::min(Ntot-Ndone, Nchunk);

            //Fill vector
            for (size_t i = 0; i < Ntodo; i++){
                systems.emplace_back(&sevnio, sevnio.STARS_MATRIX[current_idx], current_idx);
                current_idx++;
            }

            //MEMORY DEBUG
            //utilities::hardwait("After filling systems",__FILE__,__LINE__);

            //Evolve and update Nfailed
            Nfailed+=evolve_list(systems,sevnio,Ntodo,record_state);

            //MEMORY DEBUG
            //utilities::hardwait("After evolving systems",__FILE__,__LINE__);

            //Clear
            systems.clear();

            //MEMORY DEBUG
            //utilities::hardwait("After clearing systems",__FILE__,__LINE__);

            //Update Ndone
            Ndone+=Ntodo;

            //TODO 1: Here we can estimate time needed to perform a chunk run and estimate the time to the end
            //TODO 2: Maybe Nfailed is not needed because svlog has already static counters to warning, error and critical messages
            if (progress){
                std::cout<<"\r"<<Ndone<<"/"<<Ntot<<" (Nfailed:"<<Nfailed<<")";
                std::cout<<std::flush;
            }



        }

        if (progress)
            std::cout<<std::endl;

        return EXIT_SUCCESS;
    }

    /**
     * Evolve using chunk version with functor
     * @tparam T It could be Binstar or Star
     * @param evolve_function  Pointer to functor derived from  base class EvolveFunctor
     * @param Nchunk Number of systems to evolve in each chunk
     * @param sevnio  Instance of the IO class (the one linked to the binaries)
     * @param systems Vector containing the systems to evolve.
     * @param progress If true print progress information to the standard output
     * @return Number of failed evolutions
     * @Note  The vector of systems is passed by reference and it is cleared if not empty yet.
     */
    template<typename T>
    inline int chunk_dispatcher(EvolveFunctor* evolve_function,unsigned int Nchunk, IO& sevnio, std::vector<T>& systems, bool progress=true){

        ///Preliminary reset
        //If systems is not an empty vector clear it
        if (!systems.empty())
            systems.clear();

        SevnLogging sevnlog;

        ///Preliminary assignment
        //Max size
        unsigned int Ntot = sevnio.STARS_MATRIX.size();
        //Reset Nchunk if needed
        Nchunk = std::min(Nchunk,Ntot);
        //Reserve space for systems
        systems.reserve(Nchunk);

        ///Cycle
        if (progress)
            std::cout<<"Evolving systems:"<<std::endl;

        unsigned int Ndone=0;
        size_t current_idx=0;
        int Ntodo=0;
        int Nfailed=0;
        while (Ndone<Ntot){


            //Assign Ntodo
            Ntodo = std::min(Ntot-Ndone, Nchunk);

            //Fill vector
            for (size_t i = 0; i < (size_t)Ntodo; i++){

                try{
                    systems.emplace_back(&sevnio, sevnio.STARS_MATRIX[current_idx], current_idx);
                }
                catch(sevnstd::sevnio_error& e){ //sevnio error contains initialisation errors
                    sevnio.print_failed_initilisation_summary(current_idx);
                    sevnlog.error("Failed initilisation for System with ID="+utilities::n2s(current_idx,__FILE__,__LINE__)+
                    " with message:\n"+e.what(),__FILE__,__LINE__,sevnio.svpar.get_bool("initerror_stop"));
                }
                current_idx++;


            }

            //MEMORY DEBUG
            //utilities::hardwait("After filling systems",__FILE__,__LINE__);
            //Evolve and update Nfailed
            //Nfailed+=evolve_list(evolve_function,systems,sevnio,Ntodo);
            //We use system.size instead of Ntodo because wome of the systema could be not initliased due to initilisation errors
            Nfailed+=evolve_list(evolve_function,systems,sevnio,systems.size());
            //MEMORY DEBUG
            //utilities::hardwait("After evolving systems",__FILE__,__LINE__);
            //Clear
            systems.clear();

            //MEMORY DEBUG
            //utilities::hardwait("After clearing systems",__FILE__,__LINE__);

            //Update Ndone
            Ndone+=Ntodo;

            //TODO 1: Here we can estimate time needed to perform a chunk run and estimate the time to the end
            //TODO 2: Maybe Nfailed is not needed because svlog has already static counters to warning, error and critical messages
            if (progress){
                std::cout<<"\r"<<Ndone<<"/"<<Ntot<<" (Nfailed:"<<Nfailed<<")";
                std::cout<<std::flush;
            }



        }

        if (progress)
            std::cout<<std::endl;

        return EXIT_SUCCESS;
    }



}


#endif //SEVN_EVOLVE_H
