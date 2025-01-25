//
// Created by mario on 09/02/19.
//

#ifndef SEVN_STELLAR_H
#define SEVN_STELLAR_H

#include <vector>
#include <IO.h>
#include <lookup_and_phases.h>
#include <property.h>
#include <errhand.h>
#include <iostream>
#include <utilities.h>
#include <memory>
#include <params.h>
#include <algorithm>
#include <sevnlog.h>
#include <chrono>
#include <limits>


using sevnstd::SevnLogging;

//using namespace std;

#define FOR4 for(size_t _i = 0; _i < 4; _i++)

class Binstar;
class supernova;
class Staremnant;

//TODO, Leave the class star free from binary method moving all of them to the binary class?
/**
 *  Class  representing a single star in SEVN.
 *  It contains a list of properties, the info about the interpolating tracks.
 *  It has also method to call the single stellar evolution and the print function to returns the properties.
 */
class Star{

private:

    friend class Binstar;


    supernova *SN = nullptr;  /*!< Pointer to the SN model. Notice this will point to dynamically allocated memory, it has to be freed in the destructor  */
    Staremnant *staremnant= nullptr; /*!< Pointer to the remnant class */

    std::unique_ptr<utilities::ListGenerator> dtout_generator=nullptr; /*!< Unique pointer to a list generator to handle dtout */
    inline bool is_dtout_generator_callable(){
        if (dtout_generator!= nullptr){
            if (!dtout_generator->empty()){
                return true;
            }
        }

        return false;
    }

    ///FLAGS
    Lookup::Remnants initialise_as_remnant;
    bool faulty_initialisation=false; /*!< If true this object has not been properly initiliased  */
    bool auxiliary_star=false;  /*!< If true the star is considered auxiliary  */
    bool needsinit, isremnant, isempty, iskicked;
    bool tobe_SNIA=false;  /*!< Flag that can be used to rember that this star is exploding as SN Ia  */
    bool isconaked;  /*!< If true the star is COnaked  */
    bool once_conaked; /*!< If true the star was a nakedCO before to become a remnant*/
    bool break_at_remnant; /*!< Break the evolution when the stars become remnant */
    bool print_all_steps;  /*!< Print al the intermediate steps */
    bool print_per_phase;  /*!< Print each phase, currently disabled */
    bool print_only_end;   /*!< If true do not print intermediate steps */
    bool print_events;   /*!< If true print only the timesteps where events happen */
    bool initialised_on_zams; /*!< Set to true when the tini is zams or cheb for pureHE stars */
    bool force_jump;       /*!< If true force the jump */
    bool ami_on_pureHE_track;  /*!< True if the star is evolving on pureHE tracks */
    bool ami_following_QHE; /*!< True if the star is following a Quasi-homogeneous evolution*/
    bool evolution_step_completed;



    std::vector<std::vector<double>> allstates;
    std::string sntype;


    SevnLogging svlog;

    IO *io;           /*!< Pointer to io class*/
    std::vector<std::string> init_params;
    size_t ID;
    std::string name;

    std::vector<double> state;
    std::vector<std::vector<std::vector<double>*>> tables; //A vector of vector of pointers to a vector<double>
    vector<Property*> properties;
    //TODO transform raw pointers to unique_ptr? (see issue #47)
    //vector<std::unique_ptr<Property>> properties;
    std::vector<double> dtout_phase;
    unsigned long rseed=0;
    /**
     * Set the Rzams
     * //WARNING! Since this use an auxiliary star, this method must no be called inside this special constructor otherwise it will generate
     * a infinite nested call to this function.
     */
    void set_rzams();

    /**
     * Reset the staremnant pointer and cler the associated memory space.
     * The pointer is then set to nullptr
     */
     void reset_staremnant();

     /**
      * Perform all the required actions to safely destroy heap allocated members
      */
     void default_destructor();

protected:

    double tphase[Nphases] = {0.}; //Initialise all the elements to 0

    double t_co_start, t_he_start;
    std::vector<std::vector<double>> tphase_neigh;
    double tf,tini;
    double dtout;
    unsigned int dtout_type = Lookup::DTout_type::_DTUNKNOWN;

    double Z0;  /*!< Initial metallicity  (set only in the Star constructor)  */
    double Z;  /*!< Current metallicity  */
    double mzams0; /*!< Initial zams  (set only in the Star constructor) (actually it is the mass at the PreMS). Notice if we initiliase the star as remnant, this is the initial remant mass not the mass of the progenitor */
    double mzams; /*!< Current zams (actually it is the mass at the PreMS) */
    double rzams0; /*!< initial radius at t=tzams  (set only in the standard Star constructor). Notice if we initiliase the star as remnant, this indicate the radius of a non remnant star with zams=mass */
    double rzams=-1; /*!< raiud at t=tzams for the current interpolating track */
    size_t jtrack_counter=0; /*!< Counter for track jumps  */
    double MCO_max=std::nan(""); /*!< maximum MCO reachable by this star. It is set the first time the star develop a CO core, i.e. pass from the phase 4 to 5  */
    inline void set_MCO_max(){MCO_max= get_Vlast(MCO::ID);}
    double MHE_min=std::nan(""); /*!< minimum MHE reachable by this star. It is set the first time the star develop a MHE core, i.e. pass from the phase 1 to 2  */
    inline void set_MHE_min(){MHE_min= get_Vlast(MHE::ID);}

    double Mtrack[4] = {0.}; //zams masses of the neighbors (per metallicity), initialise all the elements to 0
    double Ztrack[2] = {0.}; //Z of the neighbor tracks, initialise all the elements to 0



    void remnant(){

        for (auto &prop : properties) {
            prop->set_remnant(this);
        }

    }

    void remnant_in_bse(Binstar *b){

        for (auto &prop : properties) {
            prop->set_remnant_in_bse(this,b);
        }

    }

    void turn_into_remnant();

    void crem_transition_to_NS();

    void crem_transition_to_BH();


    /**
     * This function takes care of the initiliasation of the auxiliary stars, i.e. stars that are generated
     * starting from another star. Essentially, we got the initial parameters from another stars initialising a new stars.
     * Notice, before to make the initiliasation we save the state of utilities::mtrand (it is a static threadlocal local varialbe, e.g. global variable for each thread),
     * indeed during the initiliastion
     * the mtrand engin will be reset, restarting the random number generation from the beginning. After the initliasation
     * mtrand is restored to the old state so that the random number generation is not affcted by the creation of a auxiliary star.
     * Without this save/restore, after a change of track mtrand is reset and if this happens before a SN kick you could end with a systems
     * with the same natal kick for both stars
     * chaning some of their parameters, i.e. Mzams,Z,tini, pureHE or not
     * @param s  pointer to an instance of class Star to get the initial parameters
     * @param mzams New Mzams mass in Msun
     * @param Z New Z metallicity
     * @param tini New initial time (string)
     * @param pureHE if true initialise it as a pureHE otherwsie as an H star
     */
    inline void initialise_auxiliary_star(const Star *s, double mzams, double Z, std::string tini, bool pureHE){
        //Default initialiser reset the extern variable mtrand resetting the random seed by default.
        auto saved_mtrand = utilities::mtrand;
        //io and ID are already set now, we have just to define params and call the default_initialiser.
        //params is a vector of string containing:
        //0-mzams; 1-Z; 2-spin; 3-sn_type; 4-tini; 5-tf; 6-dtout;
        init_params = s->get_initparams(true); //This is  safe since the vector assignement op  makes a deep copy.
        change_params(mzams, Z, tini); //Modify init params
        default_initialiser(false,pureHE);
        //Now restore the extern variable mtrand so that its has been reseed inside default initaliser we can instead start
        //again from the old generator (this avoid to restart the number generation after a change of track)
        utilities::mtrand = saved_mtrand;
    }

    /**
     * Wrapper of the function  initialise_auxiliary_star(const Star *s, double mzams, double Z, std::string tini, bool pureHE),
     * where tini is a double instead of a string. Internally it will transform tini to a string and call the origianl function
     * @param s  pointer to an instance of class Star to get the initial parameters
     * @param mzams New Mzams mass in Msun
     * @param Z New Z metallicity
     * @param tini New initial time (double)
     * @param pureHE if true initialise it as a pureHE otherwsie as an H star
     */
    inline void initialise_auxiliary_star(const Star *s, double mzams, double Z, double tini, bool pureHE){
        std::string tini_ss = tini==utilities::NULL_DOUBLE ? utilities::NULL_STR : utilities::n2s(tini,__FILE__,__LINE__);
        initialise_auxiliary_star(s,mzams,Z,tini_ss,pureHE);
    }



public:


    double ttimes[4]; /*!< Contains the  current times of the interpolating tracks (estimated using the plife). It is updated in tracktimes() (called inside
 *
 * @return star() )*/
    size_t pos[4];   /*!< It contains the index of the current position in the tables of the interpolating tracks. The current time is between  i nad i+1 */
    double *times_in[4];   /*!< Pointer to array containing the largest time present in the Time table of the interpolating tracks  that is smaller than ttimes. It is updated in lookuppositions() (called inside evolvestar() ) */
    bool repeatstep; /*!< If true, the evolution is restarted from the last point where repeastep=false. It is set in Timestep::evolve(Star *s) in property.h or in BTimestep::evolve(Binstar *binstar) in BinaryProperty.h */
    bool changedphase; /*!< If true, move star and interpolating track to a new phase. It is set in Timestep::evolve(Star *s) in property.h  */
    double vkick[4];  /*!< Vkick*/
    inline double get_last_fromtable(Lookup::Tables tab_id, unsigned int interpolating_track) const  {
        return tables[tab_id][interpolating_track]->back();
    }



    /**
     *
     */


    /**
     * Standard Constructor.
     * It does nothing.
     */
    Star(){}

    /**
    *  Constructor through input properties (tables and parameters).
    *  @param _io instance of the class IO. Note the instance should have already loaded the tracks tables.
    *  @param  params a vector of string with the stellar paramters (mass, Z, spin, sn_type, tfinal_evolution, t_output_evolution)
    *  @param  _ID integer ID that identifies the star.
    *  @param pureHE If true initialise the star on pureHE tracks. Notice even if false the tracks can end in a pure HE track if explicitly stated in the params
    *  @param _rseed Random seed, if 0 it will geneterated automatically
    */
    Star(IO *_io, std::vector<std::string> &params, size_t &_ID, bool pureHE=false, unsigned long _rseed=0) : io{_io}, init_params{params}, ID{_ID}, rseed{_rseed}{

        //GI 03/02/2023: Based on a issue regarding data leak after failed initialisation, I implement
        //this solution to robustly clean the heap allocated variables after a failed object initilisation.
        //Indeed, in such cases, the object is never build, so the desctructor is not available.
        //So, I have implemented a private method default_destructor that is called from the destructor to clean the heap
        //but can be called also here. So we use a try-catch-rethrow schema to check if the initilisation
        //was successful, otherwise we clean the heap and then rethrow the error.
        //We use this schema also for the other two constructors.

        //Try the initialisation
        try {
            default_initialiser(false, pureHE);
            //TODO, Notice set_rzams and the Inertia set could be in principle put inside the default initiliaser
            //However in this case, each time a Star is initialised, even if just an auxiliary star, two stars are initialised
            //instead (the second one is the star initialised at time zams to get rzmas). This could represent a significant overhead when
            //a lot of auxiliary stars are created. Therefore, we decided to set rzams only when a true star is initialised (with the current constructor)
            //set rzams and rzams0
            set_rzams(); //set rzams
            rzams0 = rzams;
            //Set propertly the inertia that could need to know rzams
            if (!get_svpar_bool("tabuse_inertia") and !amiremnant()) properties[Inertia::ID]->evolve(this);
        }
        //Catch any kind of error, dangerous in general, but we use it just to clean the heap and then
        //the error is rethrown and handled somewhere else.
        catch (...){ //Catch all errors
            default_destructor();  //Clean the heap associated to class members
            throw; //Rethrown the exact same  error
        }
    }

    /******** SPECIAL CONSTRUCTOR FOR AUXILIARY STARS *******************/
    /**
     * Special constructor. It creates a star with a  mass \p mass and metallicity \p Z but with starting time \p tini,
     * and all the other properties that are copied from the star \p s.
     * In practice, it is like we call the classical constructor Star(io, params, ID)
     * mixing the new mass, metallicity and starting age with the other information from \p s.
     * @param s  pointer to a star object, we take all the infos (except for mzams and Z and tini) from this star.
     * @param _ID integer ID that identifies the star.
     * @param mzams  Mass [in Msun] on the ZAMS of the star (actually it is the mass at the beginning of the pre-MS) [default NULL]
     * If default the mzams is the same of the star \p s.
     * @param Z   Metallicity of the star [default NULL]. If default the Z is the same of the star \p s.
     * @param tini Starting age of the stars [default NULL]. Note, this is a string value so that we can use also the string initialisation
     * (e.g. phase init (zams, tms, hsb, cheb, tcheb, thesb) or %phase (e.g. %20:1 for 20% of pfile in the zams).
     *  If default tini is the same of the star \p s.
     * @param pureHE If true initialise the star on pureHE tracks. Notice even if false the tracks can end in a pure HE track if explicitly stated in the params
     * @param _rseed Random seed, if 0 it will geneterated automatically
     */
    Star(const Star* s, size_t _ID, double mzams=utilities::NULL_DOUBLE, double Z=utilities::NULL_DOUBLE, std::string tini=utilities::NULL_STR, bool pureHE=false, unsigned long _rseed=0) : auxiliary_star{true},io{s->io}, ID{_ID}, rseed{_rseed}{
        //NOTICE, if you make changes here remember to change also the constructor below (that just change because tini is a double)
        //TODO Myabe it is better to have a template instead of two constructors

        //Try the initialisation
        try{
            initialise_auxiliary_star(s,mzams,Z,tini,pureHE);
        }
        catch (...){ //Catch all errors
            default_destructor();  //Clean the heap associated to class members
            throw; //Rethrown the exact same  error
        }
    }
    /**
      * Special constructor. It creates a star with a  mass \p mass and metallicity \p Z but with starting time \p tini,
      * and all the other properties that are copied from the star \p s.
      * In practice, it is like we call the classical constructor Star(io, params, ID)
      * mixing the new mass, metallicity and starting age with the other information from \p s.
      * @param s  pointer to a star object, we take all the infos (except for mzams and Z and tini) from this star.
      * @param _ID integer ID that identifies the star.
      * @param mzams  Mass [in Msun] on the ZAMS of the star (actually it is the mass at the beginning of the pre-MS) [default NULL]
      * If default the mzams is the same of the star \p s.
      * @param Z   Metallicity of the star [default NULL]. If default the Z is the same of the star \p s.
      * @param tini Starting age of the stars [default NULL]. It has to be a number.
      * @param pureHE If true initialise the star on pureHE tracks. Notice even if false the tracks can end in a pure HE track if explicitly stated in the params
      * @param _rseed Random seed, if 0 it will geneterated automatically
      */
    Star(const Star* s, size_t _ID, double mzams=utilities::NULL_DOUBLE, double Z=utilities::NULL_DOUBLE, double tini=utilities::NULL_DOUBLE, bool pureHE=false, unsigned long _rseed=0) : auxiliary_star{true},io{s->io}, ID{_ID},rseed{_rseed}{
        //NOTICE, if you make changes here remember to change also the constructor above (that just change because tini is as a string)
        //TODO Myabe it is better to have a template instead of two constructors

        //Try the initialisation
        try {
            initialise_auxiliary_star(s, mzams, Z, tini, pureHE);
        }
        catch (...){ //Catch all errors
            default_destructor();  //Clean the heap associated to class members
            throw; //Rethrown the exact same  error
        }
    }


    inline bool has_been_properly_initiliased(){return !faulty_initialisation;}
    /********************************************************************/


    /****************** Return M, Radius for core and envelope **********************/
    /*
     * Return the mass of the Core, if the star is nakedco Mcore=0
     */
    inline double Mcore(bool old=false){

        const auto& getP = [this,&old](int _ID){
            if (old)
                return this->getp_0(_ID);
            else
                return this->getp(_ID);
        };

        if (utilities::areEqual(getP(MCO::ID),getP(Mass::ID)))
            return 0.0;
        else if (utilities::areEqual(getP(MHE::ID),getP(Mass::ID)))
            return getP(MCO::ID);
        else
            return getP(MHE::ID);
    }
    /*
    * Return the radius of the Core, if the star is nakedco Rcore=0
    */
    inline double Rcore(bool old=false){

        const auto& getP = [this,&old](int _ID){
            if (old)
                return this->getp_0(_ID);
            else
                return this->getp(_ID);
        };

        if (utilities::areEqual(getP(MCO::ID),getP(Mass::ID)))
            return 0.0;
        else if (utilities::areEqual(getP(MHE::ID),getP(Mass::ID)))
            return getP(RCO::ID);
        else
            return getP(RHE::ID);
    }
    /*
     * Return the mass of the envelope, if the star is nakedco Menvelope=Mass
     */
    inline double Menvelope(bool old=false){

        const auto& getP = [this,&old](int _ID){
            if (old)
                return this->getp_0(_ID);
            else
                return this->getp(_ID);
        };

        if (utilities::areEqual(getP(MCO::ID),getP(Mass::ID)))
            return getP(Mass::ID);
        else if (utilities::areEqual(getP(MHE::ID),getP(Mass::ID)))
            return getP(Mass::ID)-getP(MCO::ID);
        else
            return getP(Mass::ID)-getP(MHE::ID);
    }
    /*
     * Return the mass of the envelope, if the star is nakedco Renvelope=Radius
     */
    inline double Renvelope(bool old=false){

        const auto& getP = [this,&old](int _ID){
            if (old)
                return this->getp_0(_ID);
            else
                return this->getp(_ID);
        };

        if (utilities::areEqual(getP(MCO::ID),getP(Mass::ID)))
            return getP(Radius::ID);
        else if (utilities::areEqual(getP(MHE::ID),getP(Mass::ID)))
            return getP(Radius::ID)-getP(RCO::ID);
        else
            return getP(Radius::ID)-getP(RHE::ID);
    }
    /*******************************************************************************/

    utilities::jump_convergence find_track_after_CE_Ebinding(double Ebind, double Min_Mass, double Max_mass, bool pureHE=false);

    inline bool isoutputtime() {
        if (printevents())
            return getp(Event::ID)!=Lookup::EventsList::NoEvent;

        //return(properties[Worldtime::ID]->get() == properties[NextOutput::ID]->get());
        //GI 12/07/2022 more robust alternative
        return utilities::areEqual(properties[Worldtime::ID]->get(),properties[NextOutput::ID]->get());
    }

    inline bool printall() {return print_all_steps;}
    inline bool notprint() {return print_only_end;}
    inline bool printevents() {return print_events;}

    template <class T> void castp(const size_t &id, T classtype){properties[id] = classtype;}

    /**
     * Return the initial parameters, the order of the parameters are the one in the Enum Star::Init_params.
     * Notice if the parameter effective_values is false and there are some parameters that have been set in the
     * command line option (e.f. Z) there will be a mismatch between the parameter returned by this function and the one
     * currently used in the run. Be careful!
     * @param effective_values If false return the init_parameter as it is stores in init_params conserving the
     * possible placeholder and the values from the input file, if true replace the placeholders and input values
     * with the current effective values (Z from get_Z, sntype from get_sntype,
     * tini from get_tini(), tfin from get_tfin(), dtout from get_dtout_original(),spin from get_svpar_str("spin") if != list or
     * from the values stored in init_params)
     * @return A vector of strings
     */
    inline std::vector<std::string> get_initparams(bool effective_values=true) const{

        if (!effective_values)
            return init_params;

        std::vector<std::string> init_params_tmp=init_params;

        for (unsigned int i=0; i<init_params_tmp.size(); i++){

            switch(i) {
                case Star::Init_params::_Mzams:
                    init_params_tmp[i] = utilities::n2s(get_zams(),__FILE__,__LINE__);
                    break;
                case Star::Init_params::_Z:
                    init_params_tmp[i] = utilities::n2s(get_Z(),__FILE__,__LINE__);
                    break;
                case Star::Init_params::_Tini:
                    init_params_tmp[i] = utilities::n2s(get_tini(),__FILE__,__LINE__);
                    break;
                case Star::Init_params::_Tfin:
                    init_params_tmp[i] = utilities::n2s(get_tf(),__FILE__,__LINE__);
                    break;
                case Star::Init_params::_SN_type:
                    init_params_tmp[i] = get_sn_type();
                    break;
                case Star::Init_params::_Dtout:
                    init_params_tmp[i] = utilities::n2s(get_dtout_original(),__FILE__,__LINE__);
                    break;
                case Star::Init_params::_Spin:
                    init_params_tmp[i] = get_svpar_str("spin")=="list" ? init_params[i] : get_svpar_str("spin");
                    break;
                default:
                    break;
            }

        }

        return init_params_tmp;
    }



    std::vector<std::vector<double>> load_auxiliary_table(std::string tab_name) const {
        return (io->load_auxiliary_table(tab_name));
    }

    /**
     * Record the current state and put it inside the 2D vector (allstates) containing all the states.
     */
    void recordstate(){
        for(size_t i = 0; i < io->printcolumns_star.size(); i++){
            //Concerning the Timestep, put the used one not the proposed one
            if ((size_t)io->printIDs_star[i]==Timestep::ID)
                state[i] = getp_0((size_t)io->printIDs_star[i]);
            else
                state[i] = getp((size_t)io->printIDs_star[i]);
        }



        //svlog.debug("State star "+utilities::n2s(ID,__FILE__,__LINE__)+" "+utilities::n2s(state[0],__FILE__,__LINE__)+" "+utilities::n2s(state[1],__FILE__,__LINE__));

        allstates.push_back(state);
    }



    /**
     * NOTICE, this function can be called only after a complete evolution step.
     * Record the current state and put it inside the 2D vector (allstates) containing all the states.
     * Then update the dtout_generator
     * Then Evolve the NextOutput property.
     */
    void recordstate_w_timeupdate() {

        ///Check that the evolution has been done
        evolution_guard(__FILE__,__LINE__);
        //Record state
        recordstate();
        //update_nextoutput();

        while (is_dtout_generator_callable() and dtout_generator->get()<=getp(Worldtime::ID)){
            dtout_generator->next();
        }

        if (is_dtout_generator_callable()){
            properties[NextOutput::ID]->special_evolve(this);
        }

    }



    /**
    * Print the current star properties.
    * It is supposed that the evolution of this stars is not stopped by an error, therefore
    * the summary is printed in the evolved files
    * It uses the method print_evolved_summary and print_formatted_output of the class IO.h
    */
    void inline print(){

        io->print_log();
        io->print_evolved_summary(name, rseed, ID);
        io->print_formatted_output(allstates, name, rseed, ID,false);

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
    * Print the current star properties.
    * It is supposed that the evolution of this stars is  stopped by an error, therefore
    * the summary is printed in the failed files
    * It uses the method print_failed_summary and print_formatted_output of the class IO.h
    * @param include_in_output  If true flush allstates to the output even if the evolution failed
    */
    void print_failed(bool include_in_output=false){
        io->print_log();
        io->print_failed_summary(name, rseed, ID);
        if (include_in_output){
            io->print_formatted_output(allstates, name, ID,false);
        }
    }

    /**
     * Create a log message about SNIa explosion
     * @return A string containing the log message
     */
    std::string log_message_SNIa();


    /**
     * This function return true if it exploded as SNII (remnant is NS or BH)
     * @return true ifcurrent remnant type is NS or BH and getp_0(Phase::id) is not remnant
     */
    inline bool just_exploded(){

        if ( amiremnant() and getp_0(Phase::ID)!=Lookup::Phases::Remnant)
            return true;
        return false;
    }


    ///Flag check

    /**
     * Check if it is currently  a naked Helium star
     * @return True if it is a naked He star, False otherwise
     */
    inline bool aminakedhelium() const {
        return  ami_on_pureHE_track and getp(MHE::ID)!=0;
        //TODO, this function return true also if the star is nakedC0, maybe we can separete and just check  ami_on_pureHE_track when needed and add the condition aminakedco here
        //return ami_on_pureHE_track and getp(MHE::ID)!=0 and !aminackedco();
    }

    /**
     * Check if it is currently following pureHe tables
     * @return True if it is a naked He star, False otherwise
     */
    inline bool amifollowingpureHetrack() const {
        return  ami_on_pureHE_track;
    }

    /**
     * Check if the star is a Wolf-Rayet based on the Eq.2 of Spera+19
     * @return True if the star is a Wolf-Rayet, false otherwise
     */
    inline bool amiWR() const {
        return getp(MHE::ID)>=(1- get_svpar_num("star_tshold_WR_envelope"))*getp(Mass::ID);
    }

    /**
     * Check if the star was a Wolf-Rayet based on the Eq.2 of Spera+19
     * @return True if the star is a Wolf-Rayet, false otherwise
     */
    inline bool amiWR_0() const {
        return getp_0(MHE::ID)>=(1- get_svpar_num("star_tshold_WR_envelope"))*getp_0(Mass::ID);
    }


    inline bool aminakedco() const {
        return isconaked and getp(MCO::ID)!=0;
    }


    inline bool amiempty() const { return  isempty;}

    inline bool amiremnant() const {
        return (int)getp(Phase::ID)==Remnant;
        //TODO, isremnant should be deprecated
        //return isremnant;
    }

    inline bool amiremnant_0() const {
        return (int)getp_0(Phase::ID)==Remnant;
        //TODO, isremnant should be deprecated
        //return isremnant;
    }

    /**
     * Check if the star is a WD
     * @return true if the star is a WD
     */
    inline bool amiWD() const {

        return getp(RemnantType::ID)==Lookup::Remnants::HeWD or
               getp(RemnantType::ID)==Lookup::Remnants::COWD or
               getp(RemnantType::ID)==Lookup::Remnants::ONeWD;;
    }

    /**
     * Check if the star was a WD in the previous step
     * @return true if the star is a WD
     */
    inline bool amiWD_0() const {

        return getp_0(RemnantType::ID)==Lookup::Remnants::HeWD or
               getp_0(RemnantType::ID)==Lookup::Remnants::COWD or
               getp_0(RemnantType::ID)==Lookup::Remnants::ONeWD;;
    }

    /**
     * Check if the star was a NS
     * @return true if the star was a NS
     */
    inline bool amiNS() const {

        return getp(RemnantType::ID)==Lookup::Remnants::NS_ECSN or
                getp(RemnantType::ID)==Lookup::Remnants::NS_CCSN;
    }

    /**
     * Check if the star was a NS
     * @return true if the star was a NS
     */
    inline bool amiNS_0() const {
        return getp_0(RemnantType::ID)==Lookup::Remnants::NS_ECSN or
               getp_0(RemnantType::ID)==Lookup::Remnants::NS_CCSN;
    }

    /**
     * Check if the star is a BH
     * @return true if the star is a BH
     */
    inline bool amiBH() const {
        return getp(RemnantType::ID)==Lookup::Remnants::BH;
    }

    /**
     * Check if the star is a BH
     * @return true if the star is a BH
     */
    inline bool amiBH_0() const {
        return getp_0(RemnantType::ID)==Lookup::Remnants::BH;
    }

    /**
     * Check if the star is a compact remnant (NS or BH)
     * @return true if the star is a BH or a NS
     */
     inline bool amiCompact() const {
        return amiNS() or amiBH();
     }

    /**
    * Check if the star was a compact remnant (NS or BH) in the previous step
    * @return true if the star was a BH or a NS
    */
    inline bool amiCompact_0() const {
        return amiNS_0() or amiBH_0();
    }

    /**
     * Check if the star is a degenerate remnant (WD NS or BH)
     * @return true if the star is a BH or a NS
     */
    inline bool amIdegenerate() const {
        return amiWD() or amiNS() or amiBH();
    }

    /**
     * Check if the star was a degenerate remnant (WD NS or BH) in the previous step
     * @return true if the star was a WD, BH or a NS
     */
    inline bool amIdegenerate_0() const {
        return amiWD_0() or amiNS_0() or amiBH_0();
    }



    /**
     * Check if the stars has a clear core/envelope separation.
     * @return True if the stare has a core+envelope structure, otherwise false
     */
    inline bool haveienvelope() const {
         if (getp(MHE::ID)==0 or (aminakedhelium() and getp(MCO::ID)==0)
             or  aminakedco() or amiremnant())
             return false;

         return true;
     }


     /**
      * Check if the star is following the Quasi homegenous evolution path.
      * E.g. the radius does not evolve during MS and after MS the star becomes a pureHE.
      * @return true or false
      */
    inline bool amifollowingQHE() const {
        return ami_following_QHE;
    }

    /**
     * Check if the star is an auxiliary star
     * @return true or false
     */
    inline bool amiauxiliary() const {
        return auxiliary_star;
    }

    /**
     * Check if the star has the right condition to be a RRL pulsators.
     * The check  is based on Temperature, Luminosity and phase of the star.
     * The temperature range is taken from Karczmarek+17, the check on luminositu and phase
     * depends on the function paramter.
     * @param Lmin minimum luminosity threshold in Lsun (default 1.3 Lsun)
     * @param Lmax  maximum luminosity cap in Lsun (defulat 1.9 Lsun)
     * @param just_He_core  if false consider only the star with phase 4 (CoreHeBurning), otherwise all the stars
     * with an He core (Phase>MainSequence) (default false)
     * @return true if all the conditions are satisfied, false otherwise.
     */
    inline bool amiRRL(double Lmin=1.3, double Lmax=1.9, bool just_He_core=false) const {

        if (aminakedhelium() or amiremnant())
            return false;

        double Tred  = -0.05*std::log10(getp(Luminosity::ID)) + 3.94; //Karczmarek+17, Eq. 4
        double Tblue = -0.05*std::log10(getp(Luminosity::ID)) + 4.00; //Karczmarek+17, Eq. 4
        bool temp_range = getp(Temperature::ID)>=Tred and getp(Temperature::ID)<=Tblue;
        bool lum_range = getp(Luminosity::ID)>=Lmin and  getp(Luminosity::ID)<=Lmax;
        bool Hecore_condition;
        if (just_He_core) Hecore_condition = getp(MHE::ID)>0;
        else Hecore_condition=getp(Phase::ID)==Lookup::Phases::CoreHeBurning;

        return temp_range and lum_range and Hecore_condition;
    }


    /**
     *
     * What material is transferred from this star to other objects through winds
     * or rlo.
     * @return An item of the enum list Lookup::Material
     */
    Material whatamidonating() const;

    /**
     * What material is transferred from this star to other objects through winds
     * or rlo considering the get_0 properties.
     * @return An item of the enum list Lookup::Material
     */
    Material whatamidonating_0() const;

    inline bool breaktrigger() const {

        if(isremnant && break_at_remnant) return true;
        else if (!isremnant && break_at_remnant) return false;
        else if(!break_at_remnant){
            return properties[Worldtime::ID]->get() >= tf or utilities::areEqual(properties[Worldtime::ID]->get(),tf);
        }
        else {
            SevnLogging svlog_tmp;
            svlog_tmp.critical("Star breaktrigger check failed",__FILE__,__LINE__);
        }

        return false;
    }


    int get_bse_phase() const;
    int get_bse_phase_0() const;



    inline unsigned long get_rseed(){
        if (rseed==0)
            svlog.critical("rseed has not been initialised",__FILE__,__LINE__,sevnstd::sse_error());
        return rseed;}
    inline double get_max_zams() const {
        return ami_on_pureHE_track ?  get_svpar_num("max_zams_he") : get_svpar_num("max_zams");
    }
    inline double get_min_zams() const {
        return ami_on_pureHE_track ?  get_svpar_num("min_zams_he") : get_svpar_num("min_zams");
    }
    inline double get_max_Z() const {
        return ami_on_pureHE_track ?  get_svpar_num("max_z_he") : get_svpar_num("max_z");
    }
    inline double get_min_Z() const {
        return ami_on_pureHE_track ?  get_svpar_num("min_z_he") : get_svpar_num("min_z");
    }

    inline double  get_zams() const {return mzams;}
    inline double  get_Z() const  {return Z;}
    /**
     * Get rzams, if it has never been set set it.
     * @return
     */
    inline double get_rzams()  {
        if (rzams<0) set_rzams();
        return rzams;
    }
    inline double  get_zams0() const {
        if (rzams0<0)
            svlog.critical("Trying to get rzams0, but it has never been set.",__FILE__,__LINE__,sevnstd::sanity_error());
        return mzams0;
    }
    inline double  get_Z0() const {return Z0;}
    inline double get_rzams0() const {return rzams0;}
    inline double get_MCO_max() const {return MCO_max;}
    inline double get_MHE_min() const {return MHE_min;}
    inline double  get_tf() const { return tf;}
    inline double get_tini() const { return  tini;}
    inline double get_dtout_original() const {return  dtout;}
    inline size_t  get_jcounter() const { return jtrack_counter;}
    inline const supernova* get_supernova() const {return SN;} //The pointer returned cannot modify what is inside
    double get_staremnant(size_t ID);
    inline  Staremnant* get_staremnant() const {return staremnant;}


    ///Get the value of the parameters store in the svpar class in io
    inline double get_svpar_num(std::string name) const { return io->svpar.get_num(name);}
    inline std::string get_svpar_str(std::string name) const { return io->svpar.get_str(name);}
    inline bool get_svpar_bool(std::string name) const { return io->svpar.get_bool(name);}

//    inline double get_dtout() { return dtout;}; //NB GI: there is another definition of get_dtout that check also if print_per_phase is True. Myabe we have to decide how to manipulate the output directly from dtout

    inline const std::vector<double> & getstate() {return state;} //Here I return a reference, because the state iterators are used  in binstar::record to
    // add the state elements to a bigger vectors. However, without reference this function returns a copy and using insert(..., getstate().begin, getstate().end())
    // raises an errore beacsue the first (that calls begin) and the second elements (that calls end) are different objects.

    inline double getp(const size_t &id) const {return properties[id]->get(this);}
    inline double getp_fk(const size_t &id) const {return properties[id]->get_fk(this);}
    inline double getp_0(const size_t &id) const {return properties[id]->get_0(this);}
    inline double getp_fk0(const size_t &id) const {return properties[id]->get_0_fk(this);}
    inline double get_Vlast(const size_t &id) const {return properties[id]->get_Vlast(this);}
    /**
     * Method that return the current Timestep, i.e.
     *   the value from getp(Timestep::ID) if evolution_step_completed=false, i.e. we are currently evolving the properties
     *   the value from getp_0(Timestep::ID) if evolution_step_completed=false, i.e. if the the evolution step successfully ended
     *   and now getp(Timestep::ID) returns the Timestep predicted for the next step
     * @return current Timestep in Myr
     */
    inline double get_last_Timestep() const {
        return evolution_step_completed ?
                getp_0(Timestep::ID) : getp(Timestep::ID);
    }
    inline unsigned int get_dtout_type(){return dtout_type;}

    //Get and return a string
    inline std::string getps(const size_t &id) const {return utilities::n2s(getp(id),__FILE__,__LINE__);}
    inline std::string getps_fk(const size_t &id) const {return utilities::n2s(getp_fk(id),__FILE__,__LINE__);}
    inline std::string getps_0(const size_t &id) const {return utilities::n2s(getp_0(id),__FILE__,__LINE__);}
    inline std::string getps_fk0(const size_t &id) const {return utilities::n2s(getp_fk0(id),__FILE__,__LINE__);}
    inline std::string get_id_name() const {
        std::string mss = "ID: "+utilities::n2s(ID,__FILE__,__LINE__) + " Name: "+utilities::n2s(name,__FILE__,__LINE__);
        return mss;
    }

    //Get raw V,v values, Notice these functions just return the V,V0,value,value0 withouth processing them, therefore
    //getp_raw(Radius::ID) will return log10(R) since V and V0 in the property radius are stored as log10
    inline double getp_raw(const size_t &id) const {return properties[id]->get_raw(this);}
    inline double getp_fk_raw(const size_t &id) const {return properties[id]->get_fk_raw(this);}
    inline double getp_0_raw(const size_t &id) const {return properties[id]->get_0_raw(this);}
    inline double getp_fk0_raw(const size_t &id) const {return properties[id]->get_0_fk_raw(this);}


    inline bool table_loaded(size_t propID) const {return properties[propID]->are_table_loaded();}



    /**
     * Get the pointer to a given position inside the tables for Properties in table_id
     * for given interpolating_track at position pos.
     * For example get_table_atpos(_Mass,0,0) returns the pointer to the Mass tables for
     * the interpolating track 0 (with Zams_0 and Z_0) at position in the array 0.
     * @param tabled_id Id of the table (see Lookup::Tables)
     * @param interpolating_track Id of the interpolating track
     * @param pos Position inside the table
     * @return Pointer to the value at a given position of the table.
     * If the given table is not loaded return a nullptr.
     */
    inline double* get_table_pointer_atpos(size_t table_id, int interpolating_track, int pos){

        if (tables[table_id][interpolating_track]==nullptr){
            return nullptr;
        } else{
            return &tables[table_id][interpolating_track]->at(pos);
        }

    }

    inline double get_current_tphase() const {return tphase[(int)properties[Phase::ID]->get()];}
    inline double get_next_tphase() const {return tphase[(int)properties[Phase::ID]->get()+1];}
    inline std::string get_phase_str () const { return Lookup::literal((Lookup::Phases)(int)(getp(Phase::ID)));}

    inline size_t  get_ID() const  {return ID;}
    inline std::string & get_name() {return name;}
    inline std::string  get_sn_type() const {return sntype;}

    inline double* get_mtrack() {return &Mtrack[0];}
    inline double* get_ztrack() {return &Ztrack[0];}

    inline const double* get_tphase() const {return &tphase[0];}
    inline double & costart() {return t_co_start;}
    inline double & hestart() {return t_he_start;}
    inline bool & kicked() {return iskicked;}

    /**
     * Pubblic wrapper for the property function copy_V_from
     * @param prop_id Id of the property for which we have to change V copying from another property
     * @param copy_from_id Id of the property from which we have to copy the value of v
     */
    inline void copy_property_from(const size_t &prop_id, const size_t &copy_from_id){
        properties[prop_id]->copy_V_from(properties[copy_from_id]);
    }

    ///Set flags
    inline void set_kicked() {iskicked = true;}
    inline void set_forcejump() {force_jump = true;}
    inline void set_empty(){
        isempty=true;
        set_remnant(); //Set remnant before setting property to empty otherwise Timestep set_empty will raise an error if the star is also nakedco
        for (auto &prop : properties) {
            prop->set_empty(this); //V0=V=nan
        }
        //Now safely delete staremnant
        reset_staremnant();

    }
    inline void set_empty_in_bse(Binstar *b){
        isempty=true;
        set_remnant(); //Set remnant before setting property to empty otherwise Timestep set_empty will raise an error if the star is also nakedco

        for (auto &prop : properties) {
            prop->set_empty_in_bse(this,b); //V0=V=nan
        }
        //Now safely delete staremnant
        reset_staremnant();
    }

    void set_staremnant(Staremnant *remnant);


    inline void set_remnant(){
        //Reset conaked and henaked
        if (aminakedco()){
            isconaked = false; //Now this star is a remnant, reset nakedco
            once_conaked = true; //If the star was a nakedCO before to become a remnant, take memory of this (useful in case of repetition)
        }
        isremnant=true;
    }


    inline void set_COnaked(){


        //Check
        //TODO THis is a sanity check, myabe we can remove it at certain point
        //Notice we check the difference, beacuse since Radius is stored internally in log, there is always some round-off error when compared with RCO and RHE that are not log
        if (std::abs(getp(MCO::ID)-getp(Mass::ID))>get_svpar_num("ev_naked_tshold")){
            svlog.critical("Called set_COnaked, but  MCO ("+getps(MCO::ID)+")" +
                           " have not the same value of Mass (" +getps(Mass::ID)+
                           ").", __FILE__,__LINE__,sevnstd::bse_error());

        }

        if (std::abs(getp(RCO::ID)-getp(Radius::ID))>1e-10){
            properties[RHE::ID]->copy_V_from(properties[RCO::ID]);
            properties[Radius::ID]->copy_V_from(properties[RCO::ID]);
        }


        if (getp(MCO::ID)<=get_svpar_num("sn_co_lower_ecsn")){
            //Notice, the order here is important: We have to 1) Adavance the Localtime 2)Let the Phase change 3)turno into remnant
            //Notice, the Phase is updated only in Phase special evolve the remnant type inside turno_into_remnant
            properties[Localtime::ID]->update_from_binary(this, 1.001*(get_tphase()[Remnant] - getp(Localtime::ID)));
            properties[Phase::ID]->special_evolve(this);
            turn_into_remnant();
            return;
        }
        isconaked=true;
        once_conaked = false; //Reset

    }
    inline void set_QHE(){
        if (getp(Phase::ID)>2 and getp_0(Phase::ID)>2)
            svlog.critical("Trying to set the flag ami_following_QHE when the star is more evolved than the MS phase",
                    __FILE__,__LINE__,sevnstd::bse_error());
        ami_following_QHE = get_svpar_bool("rlo_QHE");
        std::string w = utilities::log_print("QHE", this,getp(Mass::ID),getp(Radius::ID),getp(dMcumul_RLO::ID));
        print_to_log(w);
    }
    inline void reset_QHE(){
        ami_following_QHE = false;
    }

    /**
     * Update MCO_max after a merge with another star.
     * If the two stars already developed a CO core the maximum CO is the sum of the two, otherwise it is just the MCO_max of the star that already developed a CO core
     * If neither of them has a CO core, the MCO_max remains unset (it is nan)/
     * @param merged_star  Pointer to the merged star
     */
    inline void set_MCO_max_aftermerge(Star *merged_star){

        if (!std::isnan(MCO_max) and !std::isnan(merged_star->MCO_max)){
            MCO_max = MCO_max + merged_star->MCO_max;
        }
        else if (std::isnan(MCO_max) and !std::isnan(merged_star->MCO_max)){
            MCO_max = merged_star->MCO_max;
        }
    }

    /**
     * Update MHE_max after a merge with another star.
     * If the two stars already developed a CO core the minimum MHE is the sum of the two, otherwise it is just the MHE_min of the star that already developed a CO core
     * If neither of them has a CO core, the MHE_min remains unset (it is nan)/
     * @param merged_star  Pointer to the merged star
     */
    inline void set_MHE_min_aftermerge(Star *merged_star){

        if (!std::isnan(MHE_min) and !std::isnan(merged_star->MHE_min)){
            MHE_min = MHE_min + merged_star->MHE_min;
        }
        else if (std::isnan(MHE_min) and !std::isnan(merged_star->MHE_min)){
            MHE_min = merged_star->MHE_min;
        }
    }

    /**
     * Set some flags/methods at the end of a successful evolution step, i.e. when all the repetitions are done
     */
    inline void evolution_step_done(){
        evolution_step_completed=true;
    }
    inline void evolution_guard(const char *file_input = nullptr, int line_input = -1){
        if (!evolution_step_completed)
            svlog.critical("The function recordstate_w_timeupdate can be called only after an evolution step",file_input,line_input,sevnstd::sse_error());
    }

    inline bool get_COnaked(){return isconaked;}


    inline double get_dtout()  {


        if(is_dtout_generator_callable())
            return dtout_generator->get()-getp(Worldtime::ID);

        return 1E30;

        //if(print_per_phase)
        //    return dtout_phase[(size_t)properties[Phase::ID]->get()];   //TODO it can be or fixed or print all time steps or print all timesteps diveded by 2 (print every 2), or divided by 3 (print every 3) and so on...
        //else
        //    return dtout;

    }

    inline double get_dtmax(){ //avoid big jumps in the look-up tables

        double max_dt_tables = 1E30;
        double max_dt_phase = 1E30;
        double max_dt_dout=1E30;

        ///1) Max dt is the time to reach the end of the simulation
        //If we are already at the end of the simulation just set a fake max_tend
        double max_tend =  get_tf() - getp(Worldtime::ID);
        if(max_tend==0) return 1E-15;


        ///2) Max_dt_dtout, this is the maximum time to reach the next output time
        if (is_dtout_generator_callable()){
            //Notice, in case the time in dout generate is equal to the Worldtime, we use forecast instead of get
            //because we are going to the next step. We use forecast instead of next, because next is only used in recordstate_w_timeupdate.
            //This is necessary to avoid that we skip a dout due to a repeat that happens exactly at the same time in which dtout_generator->get()==getp(Worldtime::ID).
            //In that cases the Worldtime comes back but not the dtout. Instead the recordstate_w_timeupdate can be called only at the end of the evolution
            //step where we are sure all the necessary repetions have been already done.

            if (dtout_generator->get()==getp(Worldtime::ID) and !std::isnan(dtout_generator->forecast()))
                max_dt_dout = dtout_generator->forecast()-getp(Worldtime::ID);
            else if (dtout_generator->get()>getp(Worldtime::ID))
                max_dt_dout =  dtout_generator->get()-getp(Worldtime::ID);

            if (max_dt_dout<0){
                std::cout<<dtout_generator->get()<<" "<<getp(Worldtime::ID)<<std::endl;
                svlog.critical("Max timestep from dtout check is negative.",__FILE__,__LINE__,sevnstd::sse_error());
            }
        }



        ///3) Now we have two conditions for non remnant and not nakedco stars
        if (!amiremnant() and !aminakedco()){
            ///3a)  Do not jump more than two steps in the interpolating tables
            FOR4 {
                if(tables[_Time][_i]->size() > pos[_i]+2)
                    max_dt_tables = std::min(max_dt_tables, tables[_Time][_i]->at(pos[_i]+2) - tables[_Time][_i]->at(pos[_i])); //+2 means "maximum jump = 2 points in the look-up tables"
                else
                    max_dt_tables = std::min(max_dt_tables, tables[_Time][_i]->at(pos[_i]+1) - tables[_Time][_i]->at(pos[_i]));
            }
            ///3b) Use a dt so that we should have a minimum number of points per phase as defined in input (parameter ts_min_points_per_phase)
            //TODO We have to write in the documentation that the parameter ts_min_point can be overwritten by the above criteria
            max_dt_phase = (get_next_tphase() - get_current_tphase())/get_svpar_num("ts_min_points_per_phase");

        }

        //TODO Write an utility function to the take the minimum of a variadic list
        return std::min(max_tend,std::min(max_dt_tables,std::min(max_dt_phase,max_dt_dout)));
    }



    void evolve_to_tini(){


        /*** NOTE: this is just a "fake" evolution that aim just to move the position on the interpolating tables
         * to the right value at given tini. This evolution is made with a single step and there are not checks that
         * ca repeat the step***/
        svlog.debug("Bringing star to the initial Localtime",__FILE__,__LINE__);

        ///The proposed timestep is equal to the initial time of the star (tini)
        properties[Timestep::ID]->resynch(tini,false);



        //Evolve Localtime and Phase, but not Worldtime since this is just a fake pre-evolution phase to bring
        //the stars to the correct initial position given by tini.
        properties[Localtime::ID]->special_evolve(this);
        properties[Phase::ID]->special_evolve(this);


        ///Set the rigth positions on the interpolating tracks
        tracktimes();
        lookuppositions();

        //Evolve the properties (Essentially we are interested only on the fake evolution)
        //Essentially in this step we make fake evolution putting v to the right value considering the interpolating tracks and v0 is 0
        //Since v0 is 0, Dvalue is automatically set to 0 and as a consequence the evolution of V is just V=v.
        //Notice this is not true for OptionaTable, the properties override the synch method to take into account v is not set, see below.
        for (auto &prop : properties){
            prop->evolve(this);
        }


        //At this points V=v for all the properties and v=v0=0 except for the Timestep that does nothing in evolve but changes V and V0 in special evolve
        //Note if we print now properties like Radius, Luminosity and Inertia their v and v0 are equal to 1, this because
        //internally they are set to 0 as the other properties but they are log so get automatically calculate 10**v and 10**v0.


        //Now propose a conservative next tstep  to be sure to not cross the next phase.
        double proposed_dt = get_dtmax();
        if (get_svpar_num("ts_min_dt")>0)
            proposed_dt = std::max(proposed_dt,get_svpar_num("ts_min_dt"));
        if (get_svpar_num("ts_max_dt")>0)
            proposed_dt = std::min(proposed_dt,get_svpar_num("ts_max_dt"));
        //get_dtmax do not take into account that with this timestep we can directly change phase or even became a remnant
        //With the following we avoid to jump directly to other phases
        proposed_dt = std::min(proposed_dt,0.99*(get_next_tphase()-getp(Localtime::ID)));

        //Prepare Timestep
        //1-Reset all to 0, so that V0=0
        properties[Timestep::ID]->resynch(0.0, true);
        //2-Set the new timestep
        properties[Timestep::ID]->resynch(proposed_dt, false);


        //Now sync (V=V0=v0=v) all the properties to be ready to start the real evolution
        //Notice for the optional property synch is overrided and if the table are not loaded we have v0=v=V0=V.
        //Notice also that some optional properties can be disabled even if loaded, in this case, we cannot set v0=v=V0=V at this stage.
        //but the evolution without table handles this problem setting v0=V0 and v=V after each evolution.
        for (auto &prop : properties) prop->synch(); //init also the real star!!


        //GI 21/07/2022, not needed anymore, since MHE and RHE are properly set on the initialisation
        //If Naked Helium resynch MHE and RHE
        /*
        if (aminakedhelium()){
            copy_property_from(MHE::ID,Mass::ID);
            copy_property_from(RHE::ID,Radius::ID);
        }
        */



        //initialize also the first value for the next output
        properties[NextOutput::ID]->init(get_dtout());
        //Initialize OmegaSpin from Spin (given in input)
        properties[OmegaSpin::ID]->init(properties[Spin::ID]->get()*utilities::omega_crit(getp(Mass::ID),getp(Radius::ID)));
        //Initialize AngMomSpin from OmegaSpin
        properties[AngMomSpin::ID]->init(properties[OmegaSpin::ID]->get()*properties[Inertia::ID]->get());


        svlog.debug("Star " + std::to_string(ID) + " has been evolvod to the inital LocalTime "+std::to_string(properties[Localtime::ID]->get())
                    +". \n"+" Next proposed Timestep is "+std::to_string(properties[Timestep::ID]->get()),__FILE__,__LINE__);
    }

    void evolve_to_tini_as_remnant();


    /**
    * Evolve the star from the current state to the end of simulation (or one of the breaking condition is verified).
    * It each step, it evolves the times and the Phase, then it checks for breaking condition (isempty, isremnant, isconade).
    * If they are not verified it set  the right position in the interpolating tracks and then it evolves all the properties in the list.
    * Finally it checks if the repeatstep  member (it is managed by the TimeStep property), if  is true it reset the properties to the last stage and
    * repeat the evolution with a smaller timestep.
    */
    utilities::evolution evolve();

    /**
     * Handle the losing of the envelope for the star
     * @return EXIT_SUCCESS if the envelope has been removed, EXIT_FAILURE otherwise (e.g. when the stars has not a core-envelope segregation)
     */
    int lose_the_envelope(){

        //Case 1: Star has an H envelope and a He core, it loses the H envelope
        if (!aminakedhelium() && !aminakedco() && getp(MHE::ID)>0){
            //The radius of the star is the radius of the RHE
            properties[Radius::ID]->copy_V_from(properties[RHE::ID]);
            //The new mass is MTOT=MHE
            properties[Mass::ID]->copy_V_from(properties[MHE::ID]);
            //Jump to pureHE tracks
            //jump_to_pureHE_tracks();
            //TODO in sevn1 in common_envelope.cpp::lose_the_envelope the radius is not changed
        }
        //Case 1: Star is a naked Helium, loosing the envelope we remain with a naked CO
        //TODO is it transformed to a WD?
        else if (aminakedhelium() && getp(MCO::ID)>0){
            //The radius of the star is the radius of the RCO, R=RHE=RCO
            properties[Radius::ID]->copy_V_from(properties[RCO::ID]);
            properties[RHE::ID]->copy_V_from(properties[RCO::ID]);

            //The new mass is MCO=MHE=MTOT
            properties[Mass::ID]->copy_V_from(properties[MCO::ID]);
            properties[MHE::ID]->copy_V_from(properties[MCO::ID]);
            //set_COnaked();
            //utilities::hardwait("CO NAKED",get_ID(),__FILE__,__LINE__);

            //TODO in sevn1 in common_envelope.cpp::lose_the_envelope the luminosity is forced to 1e-4
        }
        //Case 3: In all the other cases (nakedHe without CO, naked CO, MS star with only H) do nothing
        else{
            //svlog.pwarning("Star", ID, "does not have any envelope, it cannot lose the envelope.",
            //        "\n Mass: ",getp(Mass::ID),
            //        "\n Mass HE: ",getp(MHE::ID),
            //        "\n Mass CO: ",getp(MCO::ID),
            //        "\n",__FILE__,__LINE__);
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;

    }

    /**
     * Update the stellar properties after binary evolution (if any)
     * @param Property_id  id of the Property to update in the list
     * @param DV  Variation of the property to update
     */
    inline void  update_from_binary(int Property_id, double DV, Binstar* b = nullptr){

        properties[Property_id]->update_from_binary(this, DV,b);

    }


    /**
     * Restore all the properties and handle some flags
     */
    inline void restore(){

        ///Restore all the properties
        for (auto &prop : properties) prop->restore();
        ///Restore other keywords
        //TODO Now isremnant is depreacted, since amiremnant check directly the phase, but the final transition to new formalism has to be completed
        //Check if the star came back from a SN explosion, in case reset isremnant
        if (getp(Phase::ID)!=Lookup::Remnant){
            isremnant=false;
            //Now, if the star is not nakedco but once_conaked is true
            //set again nakedCO. THis can happen when in the last step the star
            //becomes a remnants, but then a repeat is called, so we need to restore conaked.
            //once_conaked is set to true in turn_into_remnant if the star was a nakedCO and
            //it is reset to false in set_COnaked
            if (!aminakedco() and once_conaked)
                set_COnaked();
        }
        //Reset changed phase
        //In case the previous step starts the reduction of the timestep for the next change of phase, reset it since we are repeating
       changedphase = false;

    }

    inline void resynch(const double &dt){
        restore();
        properties[Timestep::ID]->resynch(dt);
    }



    /**
     * Sync the Timestep with the star @p s.
     * This function set the value of timestep of this object to the value stored in the s Timestep.
     * @param s
     */
    inline void sync_with(Star *s){

        properties[Timestep::ID]->resynch(s);

    }

    //set only the Timestep V (not V0).
    inline void sync_with(double dt){

        properties[Timestep::ID]->resynch(dt, false);

    }

    /**
    * Estimate the fraction of life spent in the current phase.
    * @return  (LocalTime-t_start_phase)/(t_start_next_phase - t_start_phase)
    */
    inline double plife() const {
        int currentphase = (int)getp(Phase::ID);
        return ((getp(Localtime::ID) - tphase[currentphase])/(tphase[currentphase+1] - tphase[currentphase]));
    }

    inline void tracktimes(){

        double perc = plife();
        int currentphase = (int)getp(Phase::ID);

        for(int i = 0; i < 4; i++) {

            double time_step = (tphase_neigh[i][currentphase + 1] - tphase_neigh[i][currentphase]);
            ttimes[i] = (tphase_neigh[i][currentphase]  + perc * time_step);

            //TODO THis is an experimental feature, we are artificially changing the plife for this track but of a very small amount
            if (ttimes[i]>=tphase_neigh[i][currentphase + 1]){

                double smallest_difference=std::numeric_limits<double>::epsilon()*std::fabs(tphase_neigh[i][currentphase + 1]);
                ttimes[i]=tphase_neigh[i][currentphase + 1]-smallest_difference;
                std::string w = utilities::log_print("CPLIFE",this,i);
                print_to_log(w);

                if (ttimes[i]>=tphase_neigh[i][currentphase + 1]){
                    svlog.critical("In tracktimes at track " + utilities::n2s(i,__FILE__,__LINE__) +
                                   " the difference in timestep is smaller than the machine precision, therefore this track"
                                   " will use values of the next phase.",__FILE__,__LINE__,sevnstd::sanity_error());
                }
                else{
                    svlog.warning("In tracktimes at track " + utilities::n2s(i,__FILE__,__LINE__) +
                                   " the difference in timestep is smaller than the machine precision, therefore this track"
                                   " will use values of the next phase. We used instead the smallest timestep possible",__FILE__,__LINE__);
                }
            }

            /* HOLD IMPLEMENTATION
            if (ttimes[i]>=tphase_neigh[i][currentphase + 1]){


                    svlog.critical("In tracktimes at track " + utilities::n2s(i,__FILE__,__LINE__) +
                                   " the difference in timestep is smaller than the machine precision, therefore this track"
                                   " will use values of the next phase.",__FILE__,__LINE__,sevnstd::sanity_error());
            }
             */


        }



    }


    //calculate the actual positions on the lookup tables for all the interpolating tracks, given a target time
    inline void lookuppositions() {

        //TODO We have a series of a call to vector.at instead of vector[], this robust against bugs, but maybe at certain point we can test what is the performance gain using directly the operator[]

        //  tables[starparameter::_time][0] is a pointer to a vector...
        // &tables[starparameter::_time][0]->at(0) is a pointer to the first element of the pointed vector



        //unroll
        if(needsinit) { //if the star need to be (re-)initialized use binary search algorithms for tables
            FOR4 pos[_i] = utilities::binary_search(&tables[_Time][_i]->at(0), 0, tables[_Time][_i]->size() - 1, ttimes[_i]);
            needsinit=false; //reset needsinit
        }
        else{

            for(size_t i = 0; i < 4; i++) {



                //Sanity check to track possible errors on pointers
                if (pos[i]>tables[_Time][i]->size())
                    svlog.critical("Position of lookup table pointer for interpolating track " + utilities::n2s(i,__FILE__,__LINE__) +
                    " is " + utilities::n2s(pos[i],__FILE__,__LINE__)+
                            " and it is larger than the table size " + utilities::n2s(tables[_Time][i]->size(),__FILE__,__LINE__),
                            __FILE__,__LINE__,sevnstd::sse_error());

                if (ttimes[i]==tables[_Time][i]->back())
                    svlog.critical("The ttimes of interpolating track " + utilities::n2s(i,__FILE__,__LINE__) +
                   " is equal to the last point in the time tables. This is likely due to time differences smaller than the machine precision.",__FILE__,__LINE__,sevnstd::sanity_error());


                //Time one step ahead of the current pointer position
                double time_right = tables[_Time][i]->at(pos[i] + 1);



                // If ttimes is larger than just one step ahead try going another step ahead and so on
                //TODO maybe we need to check that pos is inside the Time table size
                while (ttimes[i] >= time_right) {
                    pos[i]++;
                    time_right = tables[_Time][i]->at(pos[i] + 1);
                }



                //If for some reason ttimes is lower that current time pointe in the table use binary search
                //TODO the ttimes[i] > tables[_Time][i]->at(pos[i] + 1) seemes pointless here given the while above, check
                if (ttimes[i] < tables[_Time][i]->at(pos[i]) || ttimes[i] > tables[_Time][i]->at(pos[i] + 1)) {
                    FOR4 pos[_i] = utilities::binary_search(&tables[_Time][_i]->at(0), 0, tables[_Time][_i]->size() - 1, ttimes[_i]);
                    break;
                }
            }
        }


        FOR4 times_in[_i] = &tables[_Time][_i]->at(pos[_i]);


    }

    /**
     * Explode as SNI leaving an empty remnant
     */
    void explode_as_SNI();

    /**
     * Set the flat tobe_SNIa to true.
     * This is flag that can be used to delay the SNIa explosion that will be managed elsewhere in the code
     */
     void tobe_exploded_as_SNI();

     /**
      * Get the flag tobe_SNIa
      */
     bool is_exploding_as_SNI() const;

    /**
     * Explode as SNI leaving an empty remnant triggered by BSE
     */
    void explode_as_SNI(_UNUSED Binstar *b);

    /**
     * Check if the remnant has exceeded the Chandrasekhar or VTO limits and call the transition to the new  remnant
     * @return utilities::SNIA_EXPLODE if a SN Ia has been triggered otherwise utilities::SN_NOT_EXPLODE;
     */
     utilities::sn_explosion check_remnant_mass_limits(){

        //If it is not a remnant exit
        if (!amiremnant())
            return utilities::SN_NOT_EXPLODE;

        /* If the Chandrasekhar limit is exceeded for a white dwarf then destroy
        * the white dwarf in a supernova. If the WD is ONe then a neutron star
        * will survive the supernova.
        */
        double Mchandra = get_svpar_num("sn_Mchandra");
        double Mtov     = get_svpar_num("sn_max_ns_mass");  //Tolman–Oppenheimer–Volkoff limit
        double const &Mass = getp(Mass::ID);
        double remnant = getp(RemnantType::ID);



        if (Mass>Mchandra && amiWD()) {
            if (remnant == Remnants::HeWD or remnant == Remnants::COWD) {
                utilities::wait("SN1a triggered by accretion Mchandra limit", __FILE__, __LINE__);
                explode_as_SNI();
                return utilities::SNIA_EXPLODE;
            } else if (remnant == Remnants::ONeWD) {
                //TODO test if the transition is working
                //Check if the mass is larger than the maximum allowed NS mass
                if (Mass>Mtov)
                    svlog.critical("We are transforming a ONeWD to a NS but the final mass (" + utilities::n2s(Mass,__FILE__,__LINE__) +
                                   " Msun) is larger than the maximum mass allowed for a NS (" + utilities::n2s(get_svpar_num("sn_max_ns_mass"),__FILE__,__LINE__) +
                                   " Msun).",__FILE__,__LINE__,sevnstd::bse_error());
                crem_transition_to_NS();
                utilities::wait("We are transforming a WD to a NS (",Mass,
                                " Msun) at Worldtime ",getp(Worldtime::ID), "Myr",__FILE__,__LINE__);
                return utilities::SN_NOT_EXPLODE;
            }
        }
        else if(Mass>Mtov and amiNS()){
            crem_transition_to_BH();
            utilities::wait("We are transforming a NS to a BH (",Mass,
                            " Msun) at Worldtime ",getp(Worldtime::ID), "Myr",__FILE__,__LINE__);
            return utilities::SN_NOT_EXPLODE;
        }



        return utilities::SN_NOT_EXPLODE;
     }

    /**
    * Initialise the star.
    * Mass, Z, Spin, sn_type, t initial, t final, dt print output.
    * @param params A vector of string with the following orderd parameters:
    * Mass, Z, Spin, sn_type, t initial, t final, dt print output.
    */
    void init(std::vector<std::string> &params, bool isareset){

        //Check if dimension is correct. The check is made here since init is publicly exposed (init_stars and init_binary are privates).
        //TODO put this as a constant static parameter?
        int par_expected_size = 7;
        //If rseed is provided we have an extra parameter that is the seed.
        if (io->rseed_provided())
            par_expected_size+=1;

        if ((int)params.size()!=par_expected_size){
            std::string par_exp=utilities::n2s(par_expected_size, __FILE__, __LINE__);
            std::string inparam_size=utilities::n2s(params.size(), __FILE__, __LINE__);
            std::string this_ID=utilities::n2s(ID, __FILE__, __LINE__);
            svlog.critical("Number of Star params needs to be " +par_exp+", it is instead " + inparam_size + " (ID_bin: " + this_ID +")" , __FILE__, __LINE__);
        }



        //If this is called by a reset we don't have to reset properties, random seed and name
        if (!isareset) {

            /*** Load the properties ***/
            if (properties.size() == 0) {
                for (auto Prop : Property::all)
                    //for (size_t i = 0; i < Property::all.size(); i++) {
                    //svlog.debug("Property " + Property::all[i]->name() + " added");
                    properties.push_back(Prop->Instance()); //create instances for all properties
                //TODO transform raw pointers to unique_ptr? (see issue #47)
                //std::unique_ptr<Property> prop_temp(Property::all[i]->Instance());
                //properties.push_back(std::move(prop_temp)); //create instances for all properties
            }


            /*** Resize tables ***/
            state.resize(io->printcolumns_star.size());
            tables.resize(Lookup::_Ntables);
            for (auto &i : tables) i.resize(4);
            tphase_neigh.resize(4);

            //TODO Here happens that if parameter -rseed is in io, the rseed is got
            //from the initial parameters, otherwise we look for the rseed parameter of the function
            //This can be dangerous since one can believe that changing rseed in inpunt is possibile
            //to use a custom seed, but this is not always true, we should change this and create a
            //more clear and robust initialisation of rseed
            if (io->rseed_provided()) {
                set_rseed(inspect_param_rseed(init_params[7]), __FILE__, __LINE__); //7th column is the random seed, if provided
            }
            else if(rseed == 0) set_rseed(utilities::gen_rseed(), __FILE__, __LINE__);


            //initialize the random number generator with the seed (to generate IDs for binary systems)
            utilities::mtrand.seed(rseed);

            //Generate name
            //TODO we should create a special constuctor taking a pointer to binary as parameter, so that stars inside a binary can be created in this way, and the information like the name and the ID of the binary can be know also at the Star initilisation
            name = get_svpar_str("name_prefix")+utilities::random_keygen(&utilities::mtrand);
        }



        //GIQ: I have used the same order of the original implementation. Is the this ordered needed?
        init_1(params);
        //If we are initialising the star as a remnant we don't need to init lookup tables
        //TODO We are to check if this lack of initialisation can create indefinite behaviour somewhere else
        //I think not because V and V0 are handled with set_remnant. In any case I already check that it seems
        //to work even if we remove this if and we let to initialise also for stars that are already remnants,
        //in case the mass of the remnant is lower or higher than the lower or upper zams limits, the interpolating
        // tracks will be the first two or the last two zams.
        if (initialise_as_remnant==NotARemnant){
            init_on_lookup();
        }



        init_2(params);


        svlog.debug("Init Done",__FILE__,__LINE__);


    }

    /**
     * Set the dtout parameter from  the dtout in input
     * @param p String
     * @return Return the value dtout, if 1e30 it means that the dtour has been set with special variables (e.g. print_all_steps=true)
     */
    double inspect_param_dt(std::string p);
    /**
     * Set in Myr from tini. If the input is a number just return the number.
     * Alternatively it can be ne name of a phase (zams,tms,hsb,cheb,tcheb,hesb,end) or a given fraction of phase lifetime,
     * e.g. %80:2.
     * @param p Number or string
     * @param ignore_global if true dot not consider the global input parameter tini
     * @return Time in Myr
     */
    inline double inspect_param_tini(std::string p, bool ignore_global=false){

        ///Check if the global variable is set
        if (!ignore_global){
            if (get_svpar_str("tini")!="list"){
                p=get_svpar_str("tini");
            }
            else if(get_svpar_str("tini")=="list" and p==utilities::PLACEHOLDER){
                svlog.critical("Using placeholder option for tf but the parameter tini is"
                               " set to list",__FILE__,__LINE__,sevnstd::params_error());
            }
        }

        double _tini;

        if(utilities::string_is_number<double>(p)){
            _tini = utilities::s2n<double>(p, __FILE__, __LINE__);
        }
        else{
            ///Transform all to lowercase
            std::transform(p.begin(),p.end(),p.begin(),[](unsigned char c){ return std::tolower(c);});//using algorithm + lambda function

            ///TSTART
            if (p=="zams") {
                if (ami_on_pureHE_track)
                    svlog.critical("Using tini=zams when initialising a pureHE star. Maybe you are using a global tini overwriting list values?",__FILE__,__LINE__,sevnstd::sevnio_error(""));
                initialised_on_zams = true;
                _tini = tphase[MainSequence];
            }
            else if (p=="tms" or p=="tams"){
                if (ami_on_pureHE_track)
                    svlog.critical("Using tini=tams when initialising a pureHE star. Maybe you are using a global tini overwriting list values?",__FILE__,__LINE__,sevnstd::sevnio_error(""));
                _tini = tphase[TerminalMainSequence];
            }
            else if (p=="hsb" or p=="shb"){
                if (ami_on_pureHE_track)
                    svlog.critical("Using tini=shb when initialising a pureHE star. Maybe you are using a global tini overwriting list values?",__FILE__,__LINE__,sevnstd::sevnio_error(""));
                _tini = tphase[ShellHBurning];
            }
            else if (p=="cheb"){
                if (ami_on_pureHE_track)
                    initialised_on_zams = true;
                _tini = tphase[CoreHeBurning];
            }
            else if (p=="tcheb"){
                _tini = tphase[TerminalCoreHeBurning];
            }
            else if (p=="hesb" or p=="sheb"){
                _tini = tphase[ShellHeBurning];
            }
                //TODO ATM setting the stars to remant put bad values in the semimajor axis
                //else if (p.find("rem") != std::string::npos || p.find("REM") != std::string::npos){
                //   return tphase[Remnant];
                //}
            else if (    p.find('%')  == 0) {
                std::string token = p.substr(1);

                std::stringstream ss(token);
                int phase;
                double perc;
                double tstart;

                if (!(ss >> perc)) svlog.critical("Error while reading percentage of initial time", __FILE__, __LINE__,sevnstd::sevnio_error());
                if (ss.peek() == ':') ss.ignore();
                if (!(ss >> phase)) svlog.critical("Error while reading phase for the initial time", __FILE__, __LINE__,sevnstd::sevnio_error());


                if (phase > Remnant)
                    svlog.critical("You want to evolve a star starting at " + utilities::n2s(perc, __FILE__, __LINE__) +
                                   " percent of its life as a remnant. This does not make any sense. Please provide meaningful values",
                                   __FILE__, __LINE__);

                if (phase == Remnant) {
                    tstart = (perc/100) * tphase[phase];
                } else {
                    double time_begin = tphase[phase];
                    double time_end = tphase[phase + 1];
                    tstart = (perc/100) * (time_end - time_begin) + time_begin;
                }

                if (tstart>tphase[Remnant])
                    svlog.critical("You want to evolve a start starting at time "+utilities::n2s(tstart, __FILE__, __LINE__)+
                                   " Gyr that is larger than the star life time ["+utilities::n2s(tphase[Remnant], __FILE__, __LINE__)+ " Gyr].",__FILE__,__LINE__, sevnstd::sevnio_error());

                _tini = tstart;
            }
            else
                svlog.critical("Option "+p+" not implemented for tini",__FILE__,__LINE__,sevnstd::notimplemented_error());

        }


        if (_tini==0.)
            _tini = utilities::TINY;
        else if (_tini>=tphase[Remnant] and initialise_as_remnant==NotARemnant)
            svlog.critical("Tini ("+ utilities::n2s(_tini,__FILE__,__LINE__)+" Myr) is larger than the star lifetime ("+utilities::n2s(tphase[Remnant],__FILE__,__LINE__)+" Myr)");

        return _tini;

    }
    /**
     * Set time in Myr for tfin
     * @param p  number or string end
     * @return  Time in Myr
     */
    inline double inspect_param_tf(std::string p){

        ///Check if the global variable is set
        if (get_svpar_str("tf")!="list"){
            p=get_svpar_str("tf");
        }
        else if(get_svpar_str("tf")=="list" and p==utilities::PLACEHOLDER){
            svlog.critical("Using placeholder option for tf but the parameter tf is"
                           " set to list",__FILE__,__LINE__,sevnstd::params_error());
        }

        if(utilities::string_is_number<double>(p)){
            return utilities::s2n<double>(p, __FILE__, __LINE__);
        }
        else if(p=="end") {
            break_at_remnant = true;
            //return tphase[Remnant];
            //GI, 07/09/22, If tf=end just put a large value since
            //the evolution will be stopped when the object is a remnant
            //This is already what SEVN does for the binary evolution,
            //If tf end it sets tf=utilities::LARGE for both stars.
            //In this way, the evolution of timestep is also more robust.
            return utilities::LARGE;
        }
        else
            svlog.critical("Option "+p+" not implemented for tfin",__FILE__,__LINE__,sevnstd::notimplemented_error());

        return EXIT_FAILURE;
    }
    /**
     * Set the zams mass of the star in Msun
     * @param p zams in Msun, if HE is added at the end the star will be initiliased as a pureHE
     * @return zams in Msun
     */
    inline double inspect_param_mass(std::string p){

        std::size_t pos;
        double mass;


        if ( (pos=p.find("BH"))!=std::string::npos  ){

            std::string token = p.substr(0,pos);
            initialise_as_remnant=Lookup::Remnants::BH;
            mass = utilities::s2n<double>(token, __FILE__, __LINE__);
        }
        else if ( (pos=p.find("NS"))!=std::string::npos  ){
            std::string token = p.substr(0,pos);
            initialise_as_remnant=Lookup::Remnants::NS_CCSN;
            mass = utilities::s2n<double>(token, __FILE__, __LINE__);
        }
        else if ( (pos=p.find("NSEC"))!=std::string::npos  ){
            std::string token = p.substr(0,pos);
            initialise_as_remnant=Lookup::Remnants::NS_ECSN;
            mass = utilities::s2n<double>(token, __FILE__, __LINE__);
        }
        else if ( (pos=p.find("HEWD"))!=std::string::npos  ){
            std::string token = p.substr(0,pos);
            initialise_as_remnant=Lookup::Remnants::HeWD;
            mass = utilities::s2n<double>(token, __FILE__, __LINE__);
        }
        else if ( (pos=p.find("COWD"))!=std::string::npos  ){
            std::string token = p.substr(0,pos);
            initialise_as_remnant=Lookup::Remnants::COWD;
            mass = utilities::s2n<double>(token, __FILE__, __LINE__);
        }
        else if ( (pos=p.find("ONEWD"))!=std::string::npos  ){
            std::string token = p.substr(0,pos);
            initialise_as_remnant=Lookup::Remnants::ONeWD;
            mass = utilities::s2n<double>(token, __FILE__, __LINE__);
        }
        else if( (pos=p.find("HE"))!=std::string::npos ){
            std::string token = p.substr(0,pos);
            ami_on_pureHE_track=true;
            mass = utilities::s2n<double>(token, __FILE__, __LINE__);
        }
        else if(utilities::string_is_number<double>(p)){
            mass = utilities::s2n<double>(p, __FILE__, __LINE__);
        }
        else
            return EXIT_FAILURE;

        if (mass==0){
            initialise_as_remnant=Lookup::Remnants::Empty;
        }

        return mass;

    }
    /**
     * Set the Z of the star
     * @param p zams can be list or a number
     * @param ignore_global if true dot not consider the global input parameter tini
     * @return If list return the Z listed in the input file otherwise use the one set in the parameter Z
     */
    inline double inspect_param_Z(std::string p, bool ignore_global=false){
        ///Check if the global variable is set


        if (!ignore_global){
            if (get_svpar_str("Z")!="list"){
                p=get_svpar_str("Z");
            }
            else if(get_svpar_str("Z")=="list" and p==utilities::PLACEHOLDER){
                svlog.critical("Using placeholder option for Z but the parameter Z is"
                               " set to list",__FILE__,__LINE__,sevnstd::params_error());
            }
        }


        return utilities::s2n<double>(p, __FILE__, __LINE__);
    }
    /**
     * Set the rseed for the star
     * @param p rseed
     * @return rseed as unsigned long
     */
    inline unsigned long inspect_param_rseed(std::string p){

        if(utilities::string_is_number<unsigned long>(p)){

            unsigned long _rseed = utilities::s2n<unsigned long>(p, __FILE__, __LINE__);
            if (_rseed <=0)
                svlog.critical("Rseed in input () "+p+" cannot be negative or zero",__FILE__,__LINE__,sevnstd::sevnio_error());
            else
                return _rseed;
        }
        else
            svlog.critical("Rseed in input () "+p+" cannot be intepreted as integer",__FILE__,__LINE__,sevnstd::sevnio_error());

        return 1;
    }

    inline std::string inspect_param_sntype(std::string p){
        ///Check if the global variable is set
        if (get_svpar_str("snmode")!="list"){
            p=get_svpar_str("snmode");
        }
        else if(get_svpar_str("snmode")=="list" and p==utilities::PLACEHOLDER){
            svlog.critical("Using placeholder option for sntype but the parameter snmode is"
                           " set to list",__FILE__,__LINE__,sevnstd::params_error());
        }

        return p;
    }

    inline double inspect_param_spin(std::string p){

            ///Check if the global variable is set
            if (get_svpar_str("spin")!="list"){
                p=get_svpar_str("spin");
            }
            else if(get_svpar_str("spin")=="list" and p==utilities::PLACEHOLDER){
                svlog.critical("Using placeholder option for sntype but the parameter spin is"
                               " set to list",__FILE__,__LINE__,sevnstd::params_error());
            }


        return utilities::s2n<double>(p, __FILE__, __LINE__);
    }


    ///CHANGE OF TRACKS METHODS

    /**
     * Find new interpolating track given the new conditions on mass of the star.
     * This function checks if a track change is needed, in case it calls the appropriate change track functions
     * @param is_a_merger if true, this function is called to match the total mass after a merger
     * @return EXIT_SUCCESS or EXIT_FAILURE
     */
    utilities::jump_convergence find_new_track(bool is_merger=false);

    /**
     * Function to call after a merger, it set the proper flags and call find_new_track
     * @return
     */
    inline utilities::jump_convergence find_new_track_after_merger(){
        set_forcejump();
        return  find_new_track(true); //true is to set is_merging=true, this enable to search the mass within MIN_ZAMS - MAX_ZAMS
    }


    utilities::jump_convergence jump_to_pureHE_tracks();

    utilities::jump_convergence jump_to_normal_tracks();
    ////////////////////////


protected:

    /**
     * Enum linked to the init_param list
     */
    enum Init_params {
        _Mzams = 0,
        _Z,
        _Spin,
        _SN_type,
        _Tini,
        _Tfin,
        _Dtout
    };




    inline void change_params(double mzams=utilities::NULL_DOUBLE, double Z=utilities::NULL_DOUBLE, std::string tini=utilities::NULL_STR){

        //TODO This implementation  relies on the fact that params order remain constant. Use a more robust implementation with set?

        //io and ID are already set now, we have just to define params and call the default_initialiser.
        //params is a vector of string containing:
        //0-mzams; 1-Z; 2-spin; 3-sn_type; 4-tini; 5-tf; 6-dtout;
        if (mzams!=utilities::NULL_DOUBLE)  init_params[0] = utilities::n2s(mzams,__FILE__,__LINE__,20);
        if (Z!=utilities::NULL_DOUBLE)      init_params[1] = utilities::n2s(Z,__FILE__,__LINE__,20);
        if (tini!=utilities::NULL_STR)      init_params[4] = tini;
    }

    inline void change_params(double mzams=utilities::NULL_DOUBLE, double Z=utilities::NULL_DOUBLE, double tini=utilities::NULL_DOUBLE){

        //TODO This implementation  relies on the fact that params order remain constant. Use a more robust implementation with set?

        //io and ID are already set now, we have just to define params and call the default_initialiser.
        //params is a vector of string containing:
        //0-mzams; 1-Z; 2-spin; 3-sn_type; 4-tini; 5-tf; 6-dtout;
        if (mzams!=utilities::NULL_DOUBLE)     init_params[0] = utilities::n2s(mzams,__FILE__,__LINE__,20);
        if (Z!=utilities::NULL_DOUBLE)         init_params[1] = utilities::n2s(Z,__FILE__,__LINE__,20);
        if (tini!=utilities::NULL_DOUBLE)      init_params[4] = utilities::n2s(tini,__FILE__,__LINE__,20);
    }

    /**
     * Initialise the class to the default
     * @param isareset if true avoid to reinitialise unnecessary staff
     * @param pureHE if true initialise the stars on the pureHE track
     * @param _rseed Random seed, if 0 it will generated automatically
     */
    void default_initialiser(bool isareset = false, bool pureHE = false);


    /**
     * Reset the current star to a new state.
     * This is similar to the call of a constructor but without creating a new object.
     * Notice the reset cannot change tracks from normal to pureHE and viceversa
     * @param mzams New zams Mass of the star in Msun, default[Old values]
     * @param Z  New Z  of the star, default[Old values]
     * @param tini New initial age of the star in Myr, default[Old values]
     */
    inline void reset(double mzams=utilities::NULL_DOUBLE, double Z=utilities::NULL_DOUBLE, std::string tini=utilities::NULL_STR){

        change_params(mzams, Z, tini); //Modify init params

        //reset all the V,V0,value and value0 to 0. This is necessary for a number of reason.
        //E.g. the mass cannot increase in SSE, so it can happens the the value of the mass is restored to the old one
        //in correct_interpolation_errors and correct_interpolation_errors_real in property.h
        for (auto& prop:properties)
            prop->reset();

        default_initialiser(true,ami_on_pureHE_track); //Initialise with the new values

        jtrack_counter=0; //reset jump tracks counter
    }

    /**
     * Reset the current star to a new state.
     * This is similar to the call of a constructor but without creating a new object
     * @param mzams New zams Mass of the star in Msun, default[Old values]
     * @param Z  New Z  of the star, default[Old values]
     * @param tini New initial age of the star in Myr, default[Old values]
     */
    inline void reset(double mzams=utilities::NULL_DOUBLE, double Z=utilities::NULL_DOUBLE, double tini=utilities::NULL_DOUBLE){
        std::string tini_str = utilities::n2s(tini,__FILE__,__LINE__);
        reset(mzams, Z, tini_str);
    }

    ///CHANGE OF TRACKS METHODS
    //TODO, make find_mass_linear and find_mass_bisection as an utility function in utility.h?
    /**
     * Iterative process to find the zams that match the mass defined in property_id at given plife and Phase, using a linear interpolation.
     * We start with a low and high guess for mzams. If one of them has a relative error smaller than best_rel_err return the given zams, otherwise
     * start an iterative approach, where a the line that connect low and high are used to predict a new zams. If this one is not enough to reach the convergence,
     * the new zams takes tha place of the low or high limit and the iteration continues.
     * @param mlow Initial low guess for mzams
     * @param mhigh  Initial high guess for mzams
     * @param MAX_MASS  Maximum allowed value of mzams
     * @param MIN_MASS  Minimum allowe value of mazams
     * @param plife  Percentage of life at given phase (see below)
     * @param Phase  Phase
     * @param pureHE If true use the tracks from the pureHE tables
     * @param best_rel_err  Best relative error  considering the property_id (see below)
     * @param convergence  0-utilities::JUMP_CONVERGE 1-utilities::JUMP 2-utilities::NO_JUMP
     * @param property_id  Property to match, it can be Mass, MHE, MCO
     * @param z_look Metallicity to use to find the zams, if equal to utilities::NULL_double it is the same Z of the calling star
     * @return The value of the Mzams for which we obtain the best match. To interpret the match check the convergence (see above)
     */
    double find_mass_linear(double mlow, double mhigh, const double MAX_MASS, const double MIN_MASS, const double plife, const size_t Phase, const bool pureHE, double & best_rel_err, utilities::jump_convergence& convergence, const size_t property_id=Mass::ID, const double z_look=utilities::NULL_DOUBLE);

    /**
     * Iterative process to find the zams that match the mass defined in property_id at given plife and Phase, using a bisection method.
     * @param mlow Initial low guess for mzams
     * @param mhigh  Initial high guess for mzams
     * @param MAX_MASS  Maximum allowed value of mzams
     * @param MIN_MASS  Minimum allowe value of mazams
     * @param plife  Percentage of life at given phase (see below)
     * @param Phase  Phase
     * @param pureHE If true use the tracks from the pureHE tables
     * @param best_rel_err  Best relative error  considering the property_id (see below)
     * @param convergence  0-utilities::JUMP_CONVERGE 1-utilities::JUMP 2-utilities::NO_JUMP
     * @param property_id  Property to match, it can be Mass, MHE, MCO
     * @param z_look Metallicity to use to find the zams, if equal to utilities::NULL_double it is the same Z of the calling star
     * @return The value of the Mzams for which we obtain the best match. To interpret the match check the convergence (see above)
     */
    double find_mass_bisection(double mlow, double mhigh, const double MAX_MASS, const double MIN_MASS, const double plife, const size_t Phase, const bool pureHE, double & best_rel_err, utilities::jump_convergence& convergence, const size_t property_id=Mass::ID, const double z_look=utilities::NULL_DOUBLE);

    /**
     * Jump to a new track, matching the total mass (used for MS and pureMS stars)
     * @param best_rel_err  Best relative error of the match
     * @param is_a_merger if true, this function is called to match the total mass after a merger
     * @return 0-utilities::JUMP_CONVERGE 1-utilities::JUMP 2-utilities::NO_JUMP
     */
    utilities::jump_convergence match_M(double & best_rel_err, bool is_merger=false);
    /**
     * Jump to a new track matching the innermost core Mass at the same Phase and plife.
     * It uses a bisection method to find the matching track.
     * @param best_rel_err, reference to a double, it stores the best relative error reached in the matching
     * @param is_a_merger if true, this function is called to match the total mass after a merger
     * @return 0-utilities::JUMP_CONVERGE 1-utilities::JUMP 2-utilities::NO_JUMP
     */
    utilities::jump_convergence match_core(double & best_rel_err, _UNUSED bool is_merger=false);

    /**
     * Specialised function that finds new interpolating track for stars without He or CO core (Phase 1)
     * @param best_rel_err, reference to a double, it stores the best relative error reached in the matching.
     * @return 0-utilities::JUMP_CONVERGE 1-utilities::JUMP 2-utilities::NO_JUMP
     */
    utilities::jump_convergence match_HE_and_binding(double Ebind, double & best_rel_err, double Min_Mass=utilities::NULL_DOUBLE, double Max_Mass=utilities::NULL_DOUBLE, bool pureHe=false);
    /**
     * Set the interpolating tracks info from the one stored in the input star  \p s.
     * The function modify the pointer to the interpolating tracks tables and update the properties weights
     * @param s star
     * @return 0-utilities::JUMP_CONVERGE 1-utilities::JUMP 2-utilities::NO_JUMP
     */
    utilities::jump_convergence jump_tracks(Star *s);

    /**
     * Add 1 to the counter of the tracks changes.
     */
    void update_jtrack_counter(){jtrack_counter++;}

    /**
     * Take into account the correction due to errors on evolution (E.g. Mhe>Mtot ecc)
     * So far we use only the correct_interpolation_errors_real from properties
     */
    inline void correct_evolution_error(){

        for (auto& p : properties)
            p->correct_interpolation_errors_real(this);

    }
    ////////////////////////

    /**
     * Core function to translate SEVN phase to BSE phases.
     * Notice, we assume that nakedHe BSE is equivalent to pureHe in SEVN (i.e., star on pureHe tracks)
     * @param old if true use getp_0 properties
     * @return an integer representing the BSE phase (see Hurley+02)
     */
    int _main_get_bse_phase(bool old=false) const;

    inline std::string log_mess_COnaked(){

        std::string w = utilities::log_print("CONAKED",this,getp_0(Mass::ID),getp_0(MHE::ID),getp_0(MCO::ID),getp_0(Radius::ID),getp_0(RCO::ID),getp(MCO::ID),getp(RCO::ID),getp(Phase::ID));

        return w;
    }

    inline std::string log_mess_HEnaked(utilities::jump_convergence outcome){

        std::string w = utilities::log_print("HENAKED",this,getp_0(Mass::ID),getp_0(MHE::ID),getp_0(MCO::ID),getp_0(Radius::ID),getp(Radius::ID),getp(Phase::ID),get_zams0(),get_zams(), outcome);

        return w;
    }

    inline std::string log_mess_jumped(_UNUSED utilities::jump_convergence outcome){

        //TODO juste temporarly, is nothing wrong just to reduce the amount of logprinting
        //std::string w = utilities::log_print("NEWTRACK",this,getp_0(Mass::ID),getp_0(MHE::ID),getp_0(MCO::ID),getp_0(Radius::ID),getp(Mass::ID),getp(MHE::ID),getp(MCO::ID),getp(Radius::ID),getp(Phase::ID),get_zams0(),get_zams(), outcome);
        std::string w="";

        return w;
    }

private:


    inline void set_mzams(const double a, const char* file, const int line) {

        const double MAX_ZAMS=get_max_zams();
        const double MIN_ZAMS=get_min_zams();

        if(std::isnan(a) or std::isinf(a))
            svlog.critical("Error on star initialisation (" + get_id_name() +"): MZAMS set to INF or NAN", file, line, sevnstd::sevnio_error());
        else if(a <= 0.0 and initialise_as_remnant!=Empty)
            svlog.critical("Error on star initialisation (" + get_id_name() +"): MZAMS cannot be zero and not initiliased as an Empty remnant, current value " +
                                   utilities::n2s(a, __FILE__, __LINE__), file, line,sevnstd::sevnio_error());
        else if(a < 0.0)
            svlog.critical("Error on star initialisation (" + get_id_name() +"): MZAMS cannot be  negative, current value " +
                           utilities::n2s(a, __FILE__, __LINE__), file, line,sevnstd::sevnio_error());
        else if ( (a > MAX_ZAMS or a < MIN_ZAMS) and initialise_as_remnant==NotARemnant){
            svlog.critical("Error on star initialisation (" + get_id_name() +"): MZAMS is " + utilities::n2s(a, __FILE__, __LINE__) + ": out of range", file, line,sevnstd::sevnio_error());
        }
        else
            mzams = a;
    }

    inline void set_Z(const double a, const char* file, const int line) {

        const double MAX_Z= get_max_Z();
        const double MIN_Z= get_min_Z();


        if(std::isnan(a) || std::isinf(a))
            svlog.critical("Error on star initialisation (" + get_id_name() +"): Z set to INF or NAN", file, line,sevnstd::sevnio_error());
        else if(a <= 0.0)
            svlog.critical("Error on star initialisation (" + get_id_name() +"): Z cannot be negative", file, line,sevnstd::sevnio_error());
        else if (a > MAX_Z || a < MIN_Z)
            svlog.critical("Error on star initialisation (" + get_id_name() +"): Z is " + utilities::n2s(a, __FILE__, __LINE__) + ": out of range", file, line,sevnstd::sevnio_error());
        else
            Z = a;
    }

    inline void set_sntype(const std::string &a, _UNUSED const char* file, _UNUSED const int line) {
            sntype = a;
    }

    inline void set_tf(const double a, const char* file, const int line) {

        if(std::isnan(a) || std::isinf(a))
            svlog.critical("Error on star initialisation (" + get_id_name() +"): Tf set to INF or NAN", file, line,sevnstd::sevnio_error());
        else if(a < 0.0)
            svlog.critical("Error on star initialisation (" + get_id_name() +"): Tf cannot be negative", file, line,sevnstd::sevnio_error());
        else
            tf = a;
    }

    inline void set_dtout(const double a, const char* file, const int line) {


        if(std::isnan(a) || std::isinf(a))
            svlog.critical("Error on star initialisation (" + get_id_name() +"): dtout set to INF or NAN", file, line,sevnstd::sevnio_error());
        else if(a <0.0)
            svlog.critical("Error on star initialisation (" + get_id_name() +"): dtout cannot be negative", file, line,sevnstd::sevnio_error());
        else
            dtout = a;
    }

    inline void set_rseed(const unsigned long a, const char* file, const int line){

        if(std::isnan(a) || std::isinf(a))
            svlog.critical("Error on star initialisation (" + get_id_name() +"): Rseed set to INF or NAN", file, line,sevnstd::sevnio_error());
        else if (a<=0)
            svlog.critical("Error on star initialisation (" + get_id_name() +"): Rseed cannot be negative or 0", file, line,sevnstd::sevnio_error());
        else
            rseed = a;
    }

    /**
    * Initialise the Mzams, Z, Spin, sn_type
    * @param params A vector of string with the following orderd parameters:
    * Mass, Z, Spin, sn_type, t initial, t final, dt print output.
    */
    void init_1(std::vector<std::string> &params){


        bool ignore_global = amiauxiliary();

        set_mzams(inspect_param_mass(params[0]), __FILE__, __LINE__);
        set_Z(inspect_param_Z(params[1],ignore_global), __FILE__, __LINE__);
        set_sntype(inspect_param_sntype(params[3]), __FILE__, __LINE__);
    }


    /**
    * Initialise the LocalTime, tf and dtout
    * @param params A vector of string with the following orderd parameters:
    * Mass, Z, Spin, sn_type, t initial, t final, dt print output.
    */
    void init_2(std::vector<std::string> &params){

        bool ignore_global = amiauxiliary();

        properties[Spin::ID]->init(inspect_param_spin(params[Init_params::_Spin]));

        //properties[Localtime::ID]->init(inspect_param(params[4]));
        tini=inspect_param_tini(params[4],ignore_global);
        //if (tini==0.) tini=utilities::TINY;

        properties[Localtime::ID]->init(0);
        properties[Timestep::ID]->resynch(tini,false);

        set_tf(inspect_param_tf(params[5]), __FILE__, __LINE__);
        dtout=inspect_param_dt(params[6]);
        set_dtout(dtout, __FILE__, __LINE__);

    }

    /**
     * Set the 4 interpolating tracks and set the various variables to look into them.
     * NB: The info here are used by the properties to set their interpolating weights.
     * The weights can be different between properties. E.g. for the Mass the weights
     * are linear, while for Luminosity and Radius the weights are estimated as the Log difference (between the interpolating Mass tracks).
     */
    void init_on_lookup() {

        //TODO this function is still a bit messy, there are a lot of prints and portion that are not executed. It should be cleaned a bit.

        auto & _Z = ami_on_pureHE_track ? io->Z_HE : io->Z;
        auto & _allzams = ami_on_pureHE_track ? io->allzams_HE : io->allzams;
        auto & _tables = ami_on_pureHE_track ? io->tables_HE : io->tables;

        auto *Z_id = new size_t [2];

        /*************************************************************/
        /********** Z and zams of the 4 interpolating tracks *********/
        /*************************************************************/
        //By default we assume that the index of the interpolating tracks
        // i and i+1 (for z) and j, j+1 (for Mass) are always such that Z[i+1]>Z[i] and M[i+1]>M[i]
        // M[i][j+1]>M[i][j] and [i+1][j+1]>M[i+1][j].
        // The only problem is at the zams and Z upper limit because i+1 and j+1 will be out of index.
        // We solve this problem setting i and j has the minimum between the value obtained from the binary search and the dimension of the array -2.
        // Therefore when Z and/or Zams in input are coincident with the upper limit, i (j) corresponds to the second-last element of the tables and i+1 (j+1) to the last

        ///Find Zlow and Zhigh
        Z_id[0] = std::min(utilities::binary_search(&_Z[0], 0, _Z.size()-1, Z),_Z.size()-2); //search the metallicity of the star (Z)
        Z_id[1] = Z_id[0] + 1;

        //interpolating metallicity values
        Ztrack[0] = _Z[Z_id[0]];
        Ztrack[1] = _Z[Z_id[1]];

        ///Find neighbour ZAMSs at zlow and zhigh
        size_t *zams_id = new size_t [4];

        for(size_t i = 0, k = 0; i < 2; i++, k+=2){
            size_t dim = _allzams[Z_id[i]].size();
            zams_id[k] =  std::min(utilities::binary_search(&_allzams[Z_id[i]][0], 0, dim-1, mzams), dim-2); //search for the ZAMS of the star (mzams)
            zams_id[k+1] = zams_id[k] + 1;
        }


        //interpolating zams values
        Mtrack[0] = _allzams[Z_id[0]][zams_id[0]];
        Mtrack[1] = _allzams[Z_id[0]][zams_id[1]];
        Mtrack[2] = _allzams[Z_id[1]][zams_id[2]];
        Mtrack[3] = _allzams[Z_id[1]][zams_id[3]];
        /*************************************************************/


        //get the right table pointers
        for(size_t i = 0; i < 2; i++) {
            size_t zid = Z_id[i];

            for(size_t k = 0; k < 2; k++) {
                size_t mid = zams_id[k + 2*i];

                for (int j = 0; j < _Ntables; j++) {
                    if (_tables[j][zid].empty()){
                        tables[j][k + 2*i] = nullptr; //Tables have not be loaded, use a nullptr
                    } else{
                        tables[j][k + 2*i] = &_tables[j][zid][mid]; //points to the j-th table, at the zid-th metallicity, at the mid-th zams
                    }
                    //utilities::hardwait("A",j,tables[j][k + 2*i],__FILE__,__LINE__);
                    //star is always init at t=0, then we call the evolve function to evolve it at the desired time (read from the param file)
                }
            }
        }

        //std::cout<<" Star pointers set "<<tables[_Mass][0]->at(0)<<"   "<<tables[_Mass][1]->at(0)<<std::endl;
        //std::cout<<" Star pointers set "<<tables[_Mass][2]->at(0)<<"   "<<tables[_Mass][3]->at(0)<<std::endl;
        //exit(1);


        //We use linear weights that are functions of mzams to calculate the evolution of all the stellar
        //evolution parameters.

        //1st point: maybe, we should use weights that use luminosities to calculate luminosities, radii to calculate
        //radii and so on... Example: if mzams=60 and the two interpolating tracks are at mzams_1 = 61 and mzams_2 = 59
        //the weights are both 0.5, therefore 0.5*61 + 0.5*59 = 60 (correct). But this is not true for luminosities,
        //for example, on the MS, L = C*M^3.5, which is not linear in masses!!

        //2nd point: linear weights works well on the ZAMS because Mass vs ZAMS at t=tzams is, obviously, a linear relation.
        //Still, the relation Mass vs ZAMS at t!=tzams is NOT a linear relation therefore linear weights are wrong.
        //If the mass-grid of the look-up tables is dense enough we can assume that the relation Mass vs ZAMS is locally linear.
        //But, how dense that should be to obtain a good approximation? The same, of course, holds for the other physical
        //stellar parameters.




        for(size_t i = 0; i < 2; i++) {
            size_t zid = Z_id[i];

            for (size_t k = 0; k < 2; k++) {
                size_t mid = zams_id[k + 2 * i];

                //tphase_neigh[k + 2 * i].resize(_tables[_Phase][zid][mid].size()/2 + 1); //+1 = supernova time.. starparameter::Remnant
                tphase_neigh[k + 2 * i].resize(Nphases + 1); //+1 = supernova time.. starparameter::Remnant

                //format file: time phase time phase time phase


                //check if there are too many phases in the look-up tables
                if(_tables[_Phase][zid][mid].size()/2. + 1 > Nphases)
                    svlog.critical("Too many phases specified in the look-up tables. Please check your tables", __FILE__, __LINE__);

                //initialize all tphase_neigh to -1
                for (double &j : tphase_neigh[k + 2 * i])
                    j = -1;


                for (size_t t = 0; t < _tables[_Phase][zid][mid].size(); t+=2) {
                    double value_t = _tables[_Phase][zid][mid][t];
                    int    value_p = (int) _tables[_Phase][zid][mid][t+1];



                    if(tphase_neigh[k + 2 * i][value_p] != -1)
                        svlog.critical("Phase time has already been set. Do you have duplicated in the look-up tables?", __FILE__, __LINE__);
                    else
                        tphase_neigh[k + 2 * i][value_p] = value_t;

                }

                //set also remnant time, which is not included in the look-up tables
                size_t last_point = _tables[_Time][zid][mid].size();
                tphase_neigh[k + 2 * i][Remnant] = _tables[_Time][zid][mid][last_point-1];


#if 0
                bool heset = false;
                for (size_t j = 0; j < io->tables[_MCO][zid][mid].size(); j++) {

                    if (io->tables[_MCO][zid][mid][j] != 0.0) {
                        t_neigh_co[k + 2 * i] = io->tables[_Time][zid][mid][j - 1];
                        break;
                    }

                    if (io->tables[_MHE][zid][mid][j] != 0.0 && !heset) {
                        t_neigh_he[k + 2 * i] = io->tables[_Time][zid][mid][j - 1];
                        heset = true;
                    }
                }
#endif
            }
        }

#if 0
        //check if some phase is missing and alert the user
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < starparameter::Nphases; j++) {
            if (tphase_neigh[i][j] == -1) {
                //TODO print a message that the following phase is missing...
            }
        }
    }
#endif
        /*
        std::cout<<" No weird phases... printout phases"<<std::endl;
        for(int i = 0; i < 4; i++) {
            for (int j = 0; j < Nphases; j++) {
                //std::cout<<tphase_neigh[i][j]<<"   ";
            }
            //std::cout<<std::endl;
        }
         */

        /*
        for(int i = 0; i < 4; i++) {
            std::cout<<" t(he_start) = "<<tphase_neigh[i][TerminalMainSequence]<<"   t(co_start) = "<<tphase_neigh[i][TerminalCoreHeBurning]<<std::endl;
        }
        */

        //cout<<" nphases = "<<Nphases<<endl;



        //now I am ready to set all the weights
        for (auto &prop : properties) prop->set_w(this);

        //Get weights from Phases
        double *massw = properties[Phase::ID]->get_wm();
        double *metalw = properties[Phase::ID]->get_wz();
        //initialize tphase vector of the star using the tphase_neigh vectors
        for (int j = 0; j < Nphases; j++){
            if(tphase_neigh[0][j] != -1 && tphase_neigh[1][j] != -1 && tphase_neigh[2][j] != -1 && tphase_neigh[3][j] != -1) {
                double tphase_zlow = tphase_neigh[0][j] * massw[0] + tphase_neigh[1][j] * massw[1];
                double tphase_zhigh = tphase_neigh[2][j] * massw[2] + tphase_neigh[3][j] * massw[3];
                tphase[j] = metalw[0] * tphase_zlow + metalw[1] * tphase_zhigh;
                //cout<<" pesi = "<<metalw[0]<<"   "<<metalw[1]<<endl;
                //cout<<massw[0]<<"   "<<massw[1]<<"   "<<massw[2]<<"   "<<massw[3]<<endl;
            }
            else
                tphase[j] = -1;
        }


        /*
        std::cout<<" time stellar phases "<<std::endl;
        for (int j = 0; j < Nphases; j++){
            if(tphase[j] != -1)
                //std::cout<<tphase[j]<<"   ";
        }
        std::cout<<std::endl;
        */




        //calculate also t_he_start and t_co_start for the subphases
        double t_he_zlow = tphase_neigh[0][TerminalMainSequence] * massw[0] + tphase_neigh[1][TerminalMainSequence] * massw[1];
        double t_he_zhigh = tphase_neigh[2][TerminalMainSequence] * massw[2] + tphase_neigh[3][TerminalMainSequence] * massw[3];
        t_he_start = metalw[0] * t_he_zlow + metalw[1] * t_he_zhigh;

        double t_co_zlow = tphase_neigh[0][TerminalCoreHeBurning] * massw[0] + tphase_neigh[1][TerminalCoreHeBurning] * massw[1];
        double t_co_zhigh = tphase_neigh[2][TerminalCoreHeBurning] * massw[2] + tphase_neigh[3][TerminalCoreHeBurning] * massw[3];
        t_co_start = metalw[0] * t_co_zlow + metalw[1] * t_co_zhigh;


        /*** Set the pointer of the interpolatinc tracks of the properties to the 0 position of the right table ***/
        for (auto &prop : properties) {
            svlog.pdebug(" Set pointer for property/TabID: ",prop->name(),"/",prop->TabID(),__FILE__,__LINE__);
            prop->set_refpointers(this);
        }



        svlog.pdebug(" done pointers ",__FILE__,__LINE__);

        delete [] Z_id;
        delete [] zams_id;


    }

public:

    /*** Define the copy/move constructor and assignement + destructor ***/

    //TODO the following should be enabled if we start using unique_ptr

    /**
     * Destructor
     * use default
     */
    ~Star();


    /**
     * Copy constructor
     * Disable copy constructor
     */
    Star (const Star&) = delete;


    /**
     * Copy assignement
     * Disable copy
     * @return
     */
    Star& operator= (const Star&) = delete;


    /**
     * Move constuctor
     * Use default
     * @param other
     */
    Star (Star && other) = default;


    /**
     * Move assignement
     * Use default
     * @param other
     */
    Star& operator=(Star&& other) = default;



};

//TODO This class should be used to generate auxiliary stars that are used on as temporary object (e.g. during changing of tracks)
class Star_auxiliary : public Star{

};

class Empty : public Star{

};


//OVERLOAD OSTREAM FOR CLASS
/**
 * Overlead of the ostream opeator << to make a custome print of the class star.
 * @param os ostream operator
 * @param s  object of class star
 * @return an ostream with the star infos (Name, Mzams, Z, sn_type, tfinale, dt out)
 */
inline std::ostream& operator<< (std::ostream &os,  Star &s){
    os <<  "Name: " << s.get_name() << ", M: " << s.get_zams() << ", Z: " << s.get_Z() << ", sn_type: " << s.get_sn_type();
    os << ", tf: " << s.get_tf() << ", dtout: " << s.get_dtout_original();
    return os;
}





#endif //SEVN_STAR_H
