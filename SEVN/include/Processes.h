//
// Created by spera on 13/02/19.
//

/**
 * When Adding a process, remember that the process should usually estimate:
 * BSE
 * - DA: Variation of Semimajor axis
 * - DE: Variation of Eccentricity
 * SSE
 * - DM: Variation of Mass
 * - DAngMomSpin: Variation of stellar angular momentum spin
 *
 *
 * Notice that the various processes should take into account that during the RLO,
 * the radius should be replaced by an effective radius equal to max(RL,Rc).
 *
 */


#ifndef SEVN_PROCESSES_H
#define SEVN_PROCESSES_H


#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <random>
#include <property.h>
#include <BinaryProperty.h>
#include <Orbit.h>
#include <IO.h>
#include <lookup_and_phases.h>

#include <sevnlog.h>
using sevnstd::SevnLogging;

class Star;
class Binstar;
class Orbital_change;

class Process {

public:
    Process(){
        size++;
        VB.resize(BinaryProperty::all.size());
        VS.resize(2);

        VS[0].resize(Property::all.size());
        VS[1].resize(Property::all.size());

        set_V_to_0();
    }

    virtual ~ Process();

    static std::vector<Process*> all;
    typedef std::map<std::string,size_t> _PrintMap;
    static _PrintMap PrintMap;

    virtual Process * Instance(_UNUSED IO *_io){ return nullptr; }
    virtual inline std::string name(){return "Property (generic)";}
    virtual int evolve(_UNUSED Binstar *binstar){return EXIT_SUCCESS;}
    virtual int special_evolve(_UNUSED Binstar *binstar){return 0;} //0 Means system is not broken


    double  get(const size_t &id){ return VB[id]; }
    double  get_var(const size_t &starID, const size_t &propID){ return VS[starID][propID]; }

    std::string & get_msg() {return message;}

    inline void restore(){set_V_to_0();}

    //Special functions to  modify the VB matrix
    /**
     * Correct the VB cell containing the Semimajor variation
     * @param factor correction factor the new DA will be DA*factor.
     * @return EXIT_SUCCESS
     */
    int modify_SemimajorDV_by_a_factor(double factor){
        VB[Semimajor::ID]*=factor;
        return EXIT_SUCCESS;
    }
    /**
     * Correct the VB cell containing the Eccentricity variation
     * @param factor correction factor the new DE will be DE*factor.
     * @return EXIT_SUCCESS
     */
    int modify_EccentricityDV_by_a_factor(double factor){
        VB[Eccentricity::ID]*=factor;
        return EXIT_SUCCESS;
    }

    /**
     * Check if this process is changing the mass
     * @return false if the DV regarding the Mass of both stars is currently set to 0, true otherwise
     */
    bool is_mass_transfer_happening(){
        if (VS[0][Mass::ID]==0 && VS[1][Mass::ID]==0)
            return false;
        else
            return true;
    }

    /**
     * Check if the current process is ongoing
     * @return true if the process is ongoing, false otherwise
     */
    virtual bool is_process_ongoing() const {return orb_change->is_process_ongoing();};

    ///Events handling
    inline void set_event(double code){event_code=code;}
    inline double get_event(){return event_code;}
    inline void reset_event(){event_code=-1;}

    ///Three methods to handle the private member special_evolution_alarm
    inline void special_evolution_alarm_switch_on(){special_evolution_alarm=true;}
    inline void special_evolution_alarm_switch_off(){special_evolution_alarm=false;}
    inline bool is_special_evolution_alarm_on() const {return special_evolution_alarm;}

    ///Public interface to allow the user to get information without actually evolving the systems.
    ///NOTICE: to be self-consistent the protected and private methods (such DA and DE should use these methods)
    ///when changing the properties of stars and binaries.

    /**
     * This function should be an interface for the user to return the amount of mass accreted on the accretor
     * given the binary, the donor and the available amount of mass DM.
     * NOTICE: internally the accrete mass should use this method to be self-consistent
     * @param DM amount of mass available to be accreted in Msun.
     * @param donor Pointer to the Star that is donating the mass
     * @param accretor  Pointer to the star that is receiving the mass
     * @param binstar Pointer to the binary system hosting the donor and the accretore
     * @return the amount of accreted mass
     */
    virtual double estimate_accreted_mass(_UNUSED  double DM,  _UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binstar) const {return 0.0;}

private:
    static size_t size;
    std::vector<double> VB; /**< values of the properties of the binary system (e.g. processes change eccentricity, semimajor...) */
    std::vector<std::vector<double>> VS; /**< variations of the single star parameters due to binary stellar evolution processes (e.g. mass, radius...) */
    std::string message;
    double event_code=-1; //
    bool special_evolution_alarm=false; /**< Member used to check if a special evolution  is required */

    //std::string log_message_core

protected:

    //Notice the pointer are set to nullptr because now the Collision process don't initialise them
    //so when the Process destructor is called it could be possible that it will try to delete a fake memory
    //address. With orb_change=nullptr this will not happen and we can avoid SEGFAULT ERRORS
    Star *donor= nullptr;
    Star *accretor= nullptr;
    SevnLogging svlog;
    Orbital_change *orb_change= nullptr;

    void Register(Process *_p, size_t *id, const std::string &_name){
        Process::all.push_back(_p);
        *id = Process::size - 1;
        Process::PrintMap.insert(std::pair<std::string,size_t>(_name, *id));
        svlog.debug("Binary process "+ name() + " registered" + " (Nproperties: " + utilities::n2s(all.size(),__FILE__,__LINE__) + ")");

    }
    std::uniform_real_distribution<double> _uniform_real;


    void set(const size_t &id, const double &value){ VB[id] = value; }
    void set_var(const size_t &starID, const size_t &propID, const double &value){ VS[starID][propID] = value; }
    void set_msg(const std::string &str) {message = str;}

    inline void set_V_to_0(){

        for(size_t i = 0; i < VB.size(); i++)
            VB[i] = 0.0;

        for(size_t i = 0; i < 2; i++){
            for(size_t k = 0; k < VS[i].size(); k++){
                VS[i][k] = 0.0;
            }
        }
    }

};


class MaccretionProcess : public  Process{
public:

    /**
     * This function estimate the increment or decrement of NS Bmag due do the accretion of material
     * @param s Pointer to the star that is accreting
     * @param dM Mass accreted in Msun
     * @return THe difference in Bmag after accretion
     */
    virtual double NS_DBmag_accretion  (_UNUSED Binstar *b, _UNUSED Star *s, _UNUSED double DM) const;
    /**
     * This function estimate the increment or decrement of NS angular velocity OmegaRem due do the accretion of material
     * @param s Pointer to the star that is accreting
     * @param dM Mass accreted in Msun
     * @return THe difference in OmegaRem after accretion
     */
    virtual double NS_DOmegaRem_accretion  (_UNUSED Binstar *b, _UNUSED Star *s, _UNUSED double DM) const;

protected:

    /**
     * Function in which the mass accretion is handled.
     * This function should take care of changing all the single stellar evolution properties (not only the mass)
     * The change of binary parameters such as Semimajor axis and Eccentricity are instead changed inside evolve
     * using the methods DA and DE.
     * @param binstar Pointer to the binary star
     * @return EXIT SUCCESS
     */
    virtual int accrete_mass(_UNUSED Binstar* binstar){return EXIT_SUCCESS;}

    /**
     *It is a wrapper of all the necessray function calls to set the DV of NS properties after mass accretion
     * @param s Pointer to the ns star
     * @param dM Mass accreted in Msun
     */
    virtual void handle_NS_massaccretion(_UNUSED Binstar *b, _UNUSED Star *s, _UNUSED double DM);


};

class RocheLobe : public MaccretionProcess{

public:

    //The RocheLobe in a eccentric orbit
    //https://iopscience.iop.org/article/10.1086/513736/pdf
    RocheLobe(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }

        if (_io!= nullptr){
            std::string orb_change_name=Lookup::rlmap_name.at(Lookup::rlmap.at(_io->RL_mode));
            orb_change = Orbital_change_RL::Instance(orb_change_name);
            if(orb_change== nullptr) svlog.critical("RL mode " + _io->RL_mode + " not available", __FILE__, __LINE__);
            //TODO io should be passed directly to Instance
            orb_change->set_options(_io);
        }

    }

    static size_t ID;
    static RocheLobe _rochelobe;

    RocheLobe *Instance(IO *_io) override {
        return (new RocheLobe(_io, false));
    }

    inline std::string name() override { return "RocheLobe"; }

    int evolve(Binstar *binstar) override;

    /**
     * Roche Lobe speciale evolve: dynamic_swallowing.
     * The star is entirely swallowed through the RLO, but we don't have a mix, rather just a certain amount of mass
     * can be accreted on the other star.
     * @param binstar Pointer to the binary
     * @return EXIT_SUCCES
     */
    int special_evolve(Binstar *binstar) override;

    /**
     * Build the log message at the beginning of the RLO
     * @param binstar pointer to the binary class
     * @param q  donor to accretor mass-ration
     * @param qcrit  critical mass ratio
     * @param swap_donor_accretor
     *          if false:  use the default order for printing the star properties: first the donor, then the accretor
     *          if true:   use the reverse order: first the accretor, then the donor
     * @return  a string containing the log message
     */
    std::string log_message_start(Binstar *binstar, double q, double qcrit, bool swap_donor_accretor=false);

    /**
     * Build the log message at the end of the RLO
     * @param binstar pointer to the binary class
     * @param swap_donor_accretor
     *          if false:  use the default order for printing the star properties: first the donor, then the accretor
     *          if true:   use the reverse order: first the accretor, then the donor
     * @return  a string containing the log message
     */
    std::string log_message_end(Binstar *binstar, bool swap_donor_accretor=false);

    static std::string log_message_swallowed(Binstar *binstar, Star *swallowed, Star *other,double DM_accreted);
    std::string log_message_swallowed(Binstar *binstar, double DM_accreted);

protected:

    bool RLO_last_step=false;

    /**
     * Base method to print the log message
     * @param binstar pointer to the binary class
     * @param swap_donor_accretor
     *          if false:  use the default order for printing the star properties: first the donor, then the accretor
     *          if true:   use the reverse order: first the accretor, then the donor
     * @return the base log message
     */
    std::string _log_message(Binstar *binstar, bool swap_donor_accretor=false);

    inline int reset_DM_global(){DM_global_0=DM_global_1=0.0; return EXIT_SUCCESS;}

    /**
     *  Check if during the RLO the donor and accretor have swapped.
     *  This can happens when, for example, the donor becomes a nakedHe so that R<RL in the same phase
     *  in which the other starts trigger a R>RL.
     *  Since the start and the end of a RLO is checked only during the binary evolution, the processes recognises that the RLO is still ongoing and
     *  does not trigger the procedure to setup the end of the RLO.
     *  With this method we check if the donor/accretor swapped (using the member id_donor and RLOswapped).
     *  Then, we set properly the flags the id_donor and RLOswapped and RLO_last_step (to false, to let the evolve to properly set everything)
     * @param binstar pointer to the binary class
     * @return true if the donor/accretor swapped, false otherwise
     */
    bool check_and_handle_swapping_condition(Binstar *binstar);

    int reset_dMcumul_RLO_in_stars(Binstar *binstar);

private:

    double DM_global_0=0.0, DM_global_1=0.0;
    double id_donor=-1.;
    bool RLOswapped=false; /*!< If true the donor/accretor swapped during the last timestep without exting from the RLO*/

};

class CommonEnvelope : public MaccretionProcess{

public:
    CommonEnvelope(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }

        if (_io!= nullptr){
            std::string orb_change_name=Lookup::cemap_name.at(Lookup::cemap.at(_io->CE_mode));
            orb_change = Orbital_change_CE::Instance(orb_change_name);
            if(orb_change== nullptr) svlog.critical("CE mode " + _io->CE_mode + " not available", __FILE__, __LINE__);
        }


    }

    static size_t ID;
    static CommonEnvelope _commonenvelope;

    CommonEnvelope *Instance(IO *_io) override {
        return (new CommonEnvelope(_io, false));
    }

    inline std::string name() override { return "CommonEnvelope"; }

    int special_evolve(Binstar *binstar) override;

    static std::string log_mess(Binstar *binstar, Star *primary, Star *secondary);

    /**
     * Swallowed after CE, assumint always that the swallowed star is the primary
     * @param binstar pointer to the binary
     * @return A string log message
     */
    std::string log_message_swallowed(Binstar *binstar);


protected:

    int init(Binstar *binstar);

    /**
     * Initialise the primary and secondary attribute with the stars in the binary system.
     * The primary is the star in the system that has a core and is overfilling the Roche Lobe.
     * If both stars satisfied these conditions, the primary is the star that is more overfilling the Roche Lobe (R1/RL1>R2/RL2).
     * @param binstar Pointer to the binary system
     * @return EXIT_SUCCESS or throw an error.
     *
     * @throws sevenstd::bse_error Thrown if  neither of the two stars has a core and is overfilling the Roche Lobe (The system should not start a CE).
     *
     */
    int whoisprimary(Binstar *binstar);

    /**
     * Estimate the final mass after CE coalescence using a simple precription:
     * Mf = (Mc1 + Mc2) + K_NCE * (M2_NCE)  + k_CE * (M1_CE + M2_CE),
     * where:
     *      - Mc1, Mc2: Masses of the cores (He + CO)
     *      - K_NCE: fraction [0,1] of the (non core) mass of the secondary M2_NCE not participating to the CE (e.g. MS) reteained after the CE
     *      - K_CE: fraction [0,1] of the (non core) mass of the primary M1_CE and secondary M2_CE  participating to the CE reteained after the CE
     * @param primary Pointer to the primary star, NB this is the star starting the CE, it has to have a core.
     * @param secondary Pointer to the secondary star
     * @return Final total mass after CE
     */
    double final_mass_after_coalescence();

    /**
     * Let the two stars coalesce after the common envelope.
     * @param binstar Pointer to the binary system
     * @return A string with the log message from MIX::log_message
     */
    std::string main_coalesce(Binstar *binstar);
    /**
     * Handling the coalesce after the common envelope using the SEVN2 formalism.
     * @param binstar Pointer to the binary system
     * @return A string with the log message from MIX::log_message
     * */
    std::string coalesce(Binstar *binstar);
    /**
     * Handling the coalesce after coalesce after the common envelope using the SEVN1 binding energy approach
     * @param binstar Pointer to the binary system
     * @return A string with the log message from MIX::log_message
     */
    std::string coalesce_with_bindingEnergy(Binstar *binstar);
    /**
     * Remove the envelope after the CE phase.
     * @param binstar Pointer to the binary system.
     * @return EXIT_SUCCESS
     */
    int lose_the_envelope(Binstar *binstar);

    inline double hurley_rm_x(double Z){

        double zeta  = std::log(Z/0.02);

        double ret = 0.30406 + 0.0805*zeta + 0.0897*zeta*zeta + 0.0878*zeta*zeta*zeta + 0.0222*zeta*zeta*zeta*zeta;


        return ret;

    }

    double hurley_final_mass(const double Ebindf,  const double toll=1e-3, const double maxit=100){

        //Check
        if (!isinisialised)
            svlog.critical("CE process has not been initialised",__FILE__,__LINE__,sevnstd::ce_error());

        double XX = 1 + hurley_rm_x(primary->get_Z());


        double MC22=secondary->getp(MHE::ID);
        double MC1=M_core_primary;
        double M1=primary->getp(Mass::ID);
        double M2=secondary->getp(Mass::ID);
        double DELY, DERI;
        int Niter = 0;

        double Eratio = Ebindf/Ebind_ini;

        double CONST = (pow(M1+M2,XX))*(M1-MC1+M2-MC22)*Eratio;

        double MF=std::max(MC1+MC22, (M1+M2)*pow((Eratio),(1/XX)));
        DELY = 10*toll*MF;


        while(std::abs(DELY/MF)>toll && Niter<maxit){

            DELY = pow(MF,XX)*(MF-MC1-MC22) - CONST;

            DERI = pow(MF,(XX-1))*((1+XX)*MF - XX*(MC1+MC22));

            MF -= DELY/DERI;

            Niter++;



        }



        return MF;

    }



private:

    //Common variables
    Star * primary;
    Star * secondary;

    //CE lambda, alpha formalism
    double lambda; //Concentraction factor to define the envelope binding energy
    double alpha; //Fraction of orbital energy used to expel the envelope

    //Orbital
    double a_fin, ecc_fin;
    //Mass and radii
    double M_env_primary, M_core_primary;
    double R_core_primary;
    double M_env_secondary, M_core_secondary;
    double R_core_secondary;
    double RL_primary_final, RL_secondary_final; //Roche Lobe of primary and secondary after CE
    //Energies
    double Ebind_ini; // Initial envelope binding energy
    double Eorb_ini, Eorb_fin; //Initial and final orbital energy of the cores

    //Initilisation check
    bool isinisialised=false;
};

class SNKicks : public Process{

public:
    SNKicks(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }

        if (_io!= nullptr){
            std::string orb_change_name=Lookup::snkmap_name.at(Lookup::snkmap.at(_io->SNK_mode));
            orb_change = Orbital_change_SNKicks::Instance(orb_change_name);
            if(orb_change== nullptr) svlog.critical("RL mode " + _io->RL_mode + " not available", __FILE__, __LINE__);
        }

    }

    static size_t ID;
    static SNKicks _snkicks;

    SNKicks *Instance(IO *_io) override {
        return (new SNKicks(_io, false));
    }

    inline std::string name() override { return "SNKick"; }

    /**
     * This special evolve it is called before the other processes. If it returns 1 it means that the bianry has been broken an the
     * other processes are not considered.
     * @param binstar  Pointer to the binary
     * @return 1 if the binary has been broken by the kick, 0 if not.
     */
    int special_evolve(Binstar *binstar) override;

    static std::string log_message(Binstar *binstar,double a_fin,double e_fin, double cos_nu, double vcom);


};

class GWrad : public Process{

public:
    GWrad(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }

        if (_io!= nullptr){
            std::string orb_change_name=Lookup::gwmap_name.at(Lookup::gwmap.at(_io->GW_mode));
            orb_change = Orbital_change_GW::Instance(orb_change_name);
            if(orb_change== nullptr) svlog.critical("GW mode " + _io->GW_mode + " not available", __FILE__, __LINE__);
        }


    }

    static size_t ID;
    static GWrad _gwrad;

    GWrad *Instance(IO *_io) override {
        return (new GWrad(_io, false));
    }

    inline std::string name() override { return "GWrad"; }

    int evolve(Binstar *binstar) override;

private:

    bool first_time_GW = true; //This true until the GW orbital energy loss starts to be estimated (see inside evolve)

};

class Tides : public Process{

public:
	Tides(IO * _io= nullptr, bool reg = true) {
		if (reg) {
			Register(this, &ID, name());
		}

		if (_io!= nullptr){
            std::string orb_change_name=Lookup::tidesmap_name.at(Lookup::tidesmap.at(_io->tides_mode));
            orb_change = Orbital_change_Tides::Instance(orb_change_name);
			if(orb_change== nullptr) svlog.critical("Tides mode " + _io->tides_mode + " not available", __FILE__, __LINE__);
		}
	}

	static size_t ID;
	static Tides _tides;

	Tides *Instance(IO * _io) override {
		return (new Tides(_io,false));
	}

	inline std::string name() override { return "Tides"; }

	int evolve(Binstar *binstar) override;
};

class Mix : public  Process{

public:

    Mix(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }

        if (_io!= nullptr){
            std::string orb_change_name=Lookup::mixmap_name.at(Lookup::mixmap.at(_io->mix_mode));
            orb_change = Orbital_change_Mix::Instance(orb_change_name);
            if(orb_change== nullptr) svlog.critical("Mix mode " + _io->mix_mode + " not available", __FILE__, __LINE__);
        }


    }

    static size_t ID;
    static Mix _mix;

    Mix *Instance(IO *_io) override {
        return (new Mix(_io, false));
    }

    inline std::string name() override { return "Mix"; }

    int special_evolve(Binstar *binstar) override;

    static std::string log_message(Binstar *binstar, Star *accretor, Star *donor);


};








#endif //SEVN_PROCESSES_H
