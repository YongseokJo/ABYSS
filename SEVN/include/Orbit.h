//
// Created by spera on 21/02/19.
//

#ifndef SEVN_ORBIT_H
#define SEVN_ORBIT_H

#include <lookup_and_phases.h>
#include <iostream>
#include <vector>
#include <errhand.h>
#include <memory>
#include <sevnlog.h>
#include <cmath>
#include <map>
#include <utilities.h>
#include <star.h>
#include <MTstability.h>

#define _UNUSED __attribute__ ((unused))

/*
#ifdef HAVE_STD__MAKE_UNIQUE
    using std::make_unique;
#else
    using utilities::make_unique;
#endif
*/

using sevnstd::SevnLogging;

class Binstar;
class Process;
class Star;
struct Qcrit;



class Orbital_change  {


public:

    Orbital_change(){}

    virtual ~Orbital_change(){
        //std::cout<<"Orb "<<  this << std::endl;
    }

    Orbital_change (const Orbital_change&) = delete;

    Orbital_change& operator= (const Orbital_change&) = delete;

    Orbital_change( _UNUSED std::string argv){

    }

    virtual inline std::string name(){return "Orbital change (generic)";}

    /**
     * Function that set all the variable needed in the estimate of DA, DE, DM
     * @param b
     */
    virtual void init(_UNUSED Binstar *b) {}
    virtual double DA(_UNUSED Binstar *b, _UNUSED int procID){ return 0;}
    virtual double DE(_UNUSED Binstar *b, _UNUSED int procID){ return 0;}
    virtual double DM(_UNUSED Binstar *b, _UNUSED int procID, _UNUSED int starID){ return 0;}
    virtual double DAngMomSpin(_UNUSED Binstar *b, _UNUSED int procID, _UNUSED int starID){ return 0;}


     ///Other function that return useful information
     /**
      * Check if the process is really ongoing, i.e. it is really changing properties.
      * @return
      */
     virtual inline bool is_process_ongoing(){return false;}


    /**
     * Function to be called (if needed) in the special_evolve function of Processes.
     * Each Process can have its particular special evolve. See the proper documentation.
     * @param b Pointer to the binary.
     * @return
     */
    virtual int speciale_evolve(_UNUSED Binstar *b){ return EXIT_SUCCESS;}

    virtual void set_options(_UNUSED IO* _io){return;}

protected:
    SevnLogging svlog;
};


/**********Tides*****************/
class Orbital_change_Tides : public Orbital_change {

public:

	Orbital_change_Tides(){

	}

	virtual ~Orbital_change_Tides(){
        //std::cout<<"Tid "<<  this << std::endl;
    }

	virtual Orbital_change_Tides *instance(){ return nullptr;}

	static Orbital_change_Tides  *Instance(const std::string &name);

	inline std::string name() override {return "Orbital change Tides";}


protected:

	void Register(Orbital_change_Tides *ptr, const std::string &_name) {
		GetStaticMap().insert(std::make_pair(_name, ptr));
		GetUsed().resize(GetStaticMap().size());
		GetUsed()[GetUsed().size()-1] = 0;
	}

    inline static double f1(double e2) {
        return 1+e2*(15.5+e2*(31.875+e2*(11.5625+e2*0.390625)));
    }
    inline static double f2(double e2) {
        return 1.0+e2*(7.5+e2*(5.625+e2*0.3125));
    }
    inline static double f3(double e2) {
        return 1.0+e2*(3.75+e2*(1.875+e2*7.8125e-02));
    }
    inline static double f4(double e2) {
        return 1+e2*(1.5+e2*0.125);
    }
    inline static double f5(double e2) {
        return 1.0+e2*(3.+e2*0.375);
    }
    inline static double pseudosync_spin(double e2, double omegaorb) {
        double syncspin = f2(e2) / (f5(e2) *  pow(1.0 - e2, 1.5)) * omegaorb;
        return syncspin;
    }

    /// Override this in derived class
    inline virtual double get_tconv(_UNUSED Star *star_wtides, _UNUSED Binstar *b) {return 1.0;}

    /**
     * Compute kt in convective envelopes
     * @param Mass Total mass of the star [Msun]
     * @param Menv_cnv mass of the convective envelope [Msun]
     * @param t_cnv Convective turnover time [yr]
     * @param tide_freq Spin frequency [yr^-1]
     * @return Return the kt in yr^-1
     */
    double compute_kt_zahn_conv(double Mass, double Menv_cnv, double t_cnv, double tide_freq);

    /** Computes kt in convective envelopes in star_wtides
      *  Output is in yr^-1
      */
    virtual double compute_kt_zahn_conv(Star *star_wtides, Binstar *binstar);




    /** Computes kt in radiative envelopes in star_wtides induced by star_pert
     *  Output is in yr^-1
     */
    double compute_kt_zahn_rad(Star *star_wtides, Star *star_pert, Binstar *binstar);
    double compute_tconv_rasio96(Star *star_wtides, Binstar *b);
    double compute_tconv_pwheel18(Star *star_wtides, Binstar *b);
    //inline double get_tconv_grom_tracks(Star *star_wtides) {return star_wtides->getp(Tconv::ID) }

    /** Returns da/dt from the tides induced from star_pert on star_wtides
     *  Eq 9 Hut 81
     */
    double dadt_hut(Star *star_wtides, Star *star_pert, Binstar *binstar, double kt, double spin);

    /** Returns de/dt from the tides induced from star_wtides on star_pert
      * Eq 10 Hut 81
     */
    double dedt_hut(Star *star_wtides, Star *star_pert, Binstar *binstar, double kt, double spin);

    /** Returns dspin/dt from the tides induced from star_wtides on star_pert
      * Eq 11 Hut 81
     */
    double dspindt_hut(Star *star_wtides, Star *star_pert, Binstar *binstar, double kt, double spin);

    /**
     * Mass of the convective envelope  in Msun
     * @param star_wtides Pointer to the star
     */
    virtual double Mcnv(Star *star_wtides) const;

    /**
     * Depth of the convective envelope  in Rsun, the depth is calculated from the surface torward the inner parts
     * @param star_wtides Pointer to the star
     */
    virtual double Dcnv(Star *star_wtides) const;

private:

	static std::map<std::string, Orbital_change_Tides *> & GetStaticMap(){
		static std::map<std::string, Orbital_change_Tides *> _locmap;
		return _locmap;
	}

	static std::vector<int> & GetUsed(){
		static std::vector<int> _used;
		return _used;
	}
};

class disabled_tides : public Orbital_change_Tides {

public:

    disabled_tides(bool reg = true){
        if (reg) Register(this, name());
    };

    disabled_tides * instance()  override { return  new disabled_tides(false);}

    static disabled_tides _disabled_tides;

    inline std::string name() override { return Lookup::tidesmap_name.at(Lookup::TidesMode::_Tdisabled);}

    inline bool is_process_ongoing() override {return false;}

};

class Tides_simple : public Orbital_change_Tides {

public:

	Tides_simple(bool reg = true){
        if (reg){
            Register(this, name());
        }
	};

	~Tides_simple(){
        //std::cout<<"T simple "<<  this << std::endl;
	}

	Tides_simple * instance()  override {return new Tides_simple(false);}

	static Tides_simple _Tides_simple;

    void init(_UNUSED Binstar *binstar) override;
    
    inline std::string name() override { return Lookup::tidesmap_name.at(Lookup::TidesMode::_Tsimple);}

	double DA(_UNUSED Binstar *binstar, _UNUSED int procID) override;

	double DE(_UNUSED Binstar *binstar, _UNUSED int procID) override;

    double DAngMomSpin(_UNUSED Binstar *b, _UNUSED int procID, _UNUSED int starID) override;


protected:

    /// Here we use tconv_rasio96
    inline double get_tconv(Star *star_wtides, Binstar *b) override { return compute_tconv_rasio96(star_wtides,b); }

    double kt[2];
    double syncspin;

    double _DA=0.,_DE=0.;
    double _DLspin[2]={0.,0.};

};

class Tides_simple_notab : public Tides_simple{

public:
    Tides_simple_notab(bool reg = true){
        if (reg){
            Register(this, name());
        }
    };

    ~Tides_simple_notab(){
        //std::cout<<"T simple "<<  this << std::endl;
    }

    Tides_simple_notab * instance()  override {return new Tides_simple_notab(false);}

    static Tides_simple_notab _Tides_simple_notab;

    inline std::string name() override { return Lookup::tidesmap_name.at(Lookup::TidesMode::_Tsimple_notab);}

protected:

    /**
     * Mass of the convective envelope  in Msun
     * @param star_wtides Pointer to the star
     */
    double Mcnv(Star *star_wtides) const override;

    /**
     * Depth of the convective envelope  in Rsun, the depth is calculated from the surface torward the inner parts
     * @param star_wtides Pointer to the star
     */
    double Dcnv(Star *star_wtides) const override;

};

/***************************/

/**********GW*****************/
class Orbital_change_GW : public Orbital_change{

public:
    Orbital_change_GW(){}

    virtual ~Orbital_change_GW(){
        //std::cout<<"Orb GW "<<  this << std::endl;

    }

    virtual Orbital_change_GW *instance(){ return nullptr;}

    static Orbital_change_GW  *Instance(const std::string &name);

    inline std::string name() override {return "Orbital change GW";}


protected:

    void Register(Orbital_change_GW *ptr, const std::string &_name) {
        GetStaticMap().insert(std::make_pair(_name, ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }

private:

    static std::map<std::string, Orbital_change_GW *> & GetStaticMap(){
        static std::map<std::string, Orbital_change_GW *> _locmap;
        return _locmap;
    }

    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }

};

class disabled_gw : public Orbital_change_GW {

public:

    disabled_gw(bool reg = true){
        if (reg){
            Register(this, name());
        }
    };

    disabled_gw * instance()  override {return new disabled_gw(false);}

    static disabled_gw _disabled_gw;

    inline std::string name() override { return Lookup::gwmap_name.at(Lookup::GWMode::_GWdisabled); }

    inline bool is_process_ongoing() override {return false;}

};

class Peters_gw : public Orbital_change_GW{

    Peters_gw(bool reg = true){
        if (reg){
            Register(this, name());
        }


    };

    ~Peters_gw(){
        //std::cout<<"Orb GW Peters "<<  this << std::endl;

    }

    Peters_gw * instance()  override {return new Peters_gw(false);}

    static Peters_gw _peters_gw;

    inline std::string name() override { return Lookup::gwmap_name.at(Lookup::GWMode::_GWPeters); }

    /**
     * Change in semimajor axis obtained using Eq. 5.6 in Peters64
     * @param b  binary sistem
     * @param procID  Id of the process calling the function
     * @return Variation of semimajor axis in the time interval stored in BTimestep
     */
    double DA(_UNUSED Binstar *b, _UNUSED int procID) override;

    /**
     * Change in eccentricity using Eq. 5.7 in Peters64
     * @param b  binary sistem
     * @param procID  Id of the process calling the function
     * @return Variation of eccentricity in the time interval stored in BTimestep
     */
    double DE(_UNUSED Binstar *b, _UNUSED int procID) override;


protected:

    //G^3/c^5 in sun units
    const double G3_over_c5=utilities::G3_over_c5;



    //dat/dt/a =  -K_a*G3_over_c5*m1*m2*(m1+m2)/(a^4 * (1-e^2)^7/2) * (1+c2_a*e^2+c4_a*e^4) (Eq. 5.6 Peter64)
    const double K_a=64./5., c2_a=73./24., c4_a=37./96.;
    // det/dt/e =  -K_e*G3_over_c5*m1*m2*(m1+m2)/(a^4 * (1-e^2)^5/2) * (1+c2_e*e^2) (Eq. 5.7 Peter64)
    const double K_e=304./15., c2_e=121./304.;

    /**
     * Core function both for estimate of edot/e and adot/a (see below)
     * @param M1  Mass of the first star [Msun]
     * @param M2  Mass of the second star [Msun]
     * @param a semi-major axis  [Rsun]
     * @return G^3/c^5 * (M1*M2)*(M1+M2)/a^4
     */
    inline double function_core(double M1, double M2, double a){

        double a4=a*a*a*a;

        return G3_over_c5*M1*M2*(M1+M2)/a4;

    }


    /**
     * Estimate of the eccentricity time derivative over current eccentricity from Eq. 5.7 in Peters64
     * @param M1  Mass of the first star [Msun]
     * @param M2  Mass of the second star [Msun]
     * @param a   semi-major axis  [Rsun]
     * @param e   eccentricity
     * @return  edot/e = -304/15 (G^3/c^5 \p M1* \p M2* (\p M1 +\p M2))/a^4 (\p 1 +  121/24 * \p e^2 )/(1-e*e)^(5/2)
     */
    inline double Peters_GW_de(double M1, double M2, double a, double e){

        double cost=K_e*function_core(M1, M2, a);
        double e2=e*e;

        double num=(1+c2_e*e2);
        double den=(1-e2);
        den=std::sqrt(den*den*den*den*den); //(1-e^2)^(5/2)

        return -cost*num/den;

    }


    /**
     * Estimate of semimajor axis time derivative over current semimajor axis from Eq. 5.6 in Peters64
     * @param M1  Mass of the first star
     * @param M2  Mass of the second star
     * @param a   semi-major axis  [Rsun]
     * @param e   eccentricity
     * @return adot/a = -64/5 (G^3/c^5 \p M1* \p M2* (\p M1 +\p M2))/a^4 (\p 1 +  73/24 * \p e^2 + 37/96 * \p e^4 )/(1-e*e)^(7/2)
     */
    inline double Peters_GW_da(double M1, double M2, double a, double e){

        double cost=K_a*function_core(M1, M2, a);
        double e2=e*e;
        double e4=e2*e2;
        double den=(1-e2);
        den=std::sqrt(den*den*den*den*den*den*den); //(1-e^2)^(7/2)

        double num=(1+c2_a*e2+c4_a*e4);

        return -cost*num/den;

    }




};
/***************************/

/**********Roche Lobe*****************/
class Orbital_change_RL : public Orbital_change{

public:
    Orbital_change_RL(){}

    virtual ~Orbital_change_RL(){
        //std::cout<<"Orb RL "<<  this << std::endl;
        delete mt_stability;
        mt_stability= nullptr;
    }

    virtual Orbital_change_RL *instance(){ return nullptr;}

    static Orbital_change_RL  *Instance(const std::string &name);

    void set_options(IO* _io) override {

        //Set MTstability pointer
        std::string mtstability_option = _io->svpar.get_str("rlo_stability");
        mt_stability = MTstability::Instance(mtstability_option);
        if (mt_stability == nullptr)
            svlog.critical("Unknown mtstability option: [" + mtstability_option + "]", __FILE__, __LINE__,
                           sevnstd::sevnio_error());
    }

    inline std::string name() override {return "Orbital change RL";}

    /**
     * Init function to be called at the beginning of the Process::evolve
     * @param b Pointer to the binary
     */
    void init(_UNUSED Binstar *b) override {};

    int speciale_evolve(_UNUSED Binstar *b) override;

    inline double  get_DM(int ID){
        if (ID!=0 && ID!=1)
            svlog.critical("Stars ID in binary can be only 0 or 1, current value "+
                           utilities::n2s(ID,__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::rl_error());

        return _DM[ID];
    };

    double DM(_UNUSED Binstar *b, _UNUSED int procID, _UNUSED int starID) override {
        return get_DM(starID);
    }

    inline bool is_swallowed(int starID){
        return _is_swallowed[starID];
    }

    inline bool is_process_ongoing() override {return _is_RLO_happening;}
    inline bool is_colliding()  {return _is_colliding;}


    virtual inline double get_mt_value() const {svlog.critical("get_q not implemented for class Orbital_change_RL",__FILE__,__LINE__,sevnstd::notimplemented_error());
        return  0.0;}

    virtual inline double get_mt_tshold() const {svlog.critical("get_qcrit not implemented for class Orbital_change_RL",__FILE__,__LINE__,sevnstd::notimplemented_error());
        return  0.0;}

    bool novae=false;
    bool super_eddington=false;




protected:
    void Register(Orbital_change_RL *ptr, const std::string &_name){
        GetStaticMap().insert(std::make_pair(_name, ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }

    MTstability *mt_stability;

    /**
     * Actual RocheLobe Overflow process
     * @param donor Pointer to the star that is overflowing the RocheLobe
     * @param accretor Pointer to the star that is accreting mass
     * @param b Pointer to the binary
     */
    virtual void RLO(_UNUSED Star *donor,_UNUSED Star *accretor,_UNUSED Binstar *b){};

    /**
     * General function to  be called to correct the accreted mass for whatever reason.
     * In this general version it just returns the same amount of DM_accreted in input.
     * The idea is that specialised version of this process could need to add special
     * correction without change the rest of the process.
     * @param DM_accreted mass accreted by the accretor (before this correction) in Msun
     * @param DM_donor mass donated by the donor in Msun
     * @param donor Pointer to the star that is overflowing the RocheLobe
     * @param accretor Pointer to the star that is accreting mass
     * @param b Pointer to the binary
     */
    virtual double correct_accreted_mass(double DM_accreted, _UNUSED double DM_donor,
                                       _UNUSED Star *donor,_UNUSED Star *accretor,_UNUSED Binstar *b) {
        return DM_accreted;
    }

    /**
     * Roche Lobe radius normalised over the star radius
     * @param star pointer to the star for which we are estimateting fRL
     * @param b pointer to the binary
     * @return R_rl/R_star
     */
    virtual double fRL_radius(_UNUSED Star* star, _UNUSED Binstar *b){svlog.critical("fRL not implemented for this RLO option ("+name()+")",__FILE__,__LINE__,sevnstd::rl_error()); return -1;}

    inline int  set_DM(double DM,int starID){
        if (starID!=0 && starID!=1) {
            svlog.critical("Stars ID in binary can be only 0 or 1, current value " +
                           utilities::n2s(starID, __FILE__, __LINE__), __FILE__, __LINE__, sevnstd::rl_error());
        }
        _DM[starID] = DM;
        return EXIT_SUCCESS;
    };
    inline int  reset_DM(){
        _DM[0]=_DM[1]=0;
    return EXIT_SUCCESS;
    }

    inline int set_swallowed(int starID){
        if (starID!=0 && starID!=1) {
            svlog.critical("Stars ID in binary can be only 0 or 1, current value " +
                           utilities::n2s(starID, __FILE__, __LINE__), __FILE__, __LINE__, sevnstd::rl_error());
        }
        _is_swallowed[starID]=true;
        return EXIT_SUCCESS;
    }
    inline int reset_swallowed(int starID){
        _is_swallowed[starID]=false;
        return EXIT_SUCCESS;
    }
    inline int reset_swallowed(){
        _is_swallowed[0]=_is_swallowed[1]=false;
        return EXIT_SUCCESS;
    }

    inline void set_mix(){mix= true;}
    inline void unset_mix(){mix= false;}
    inline bool get_mix(){ return  mix;}
    inline void set_comenv(){comenv= true;}
    inline void unset_comenv(){comenv= false;}
    inline bool get_comenv(){ return  comenv;}
    inline void set_is_RLO_happening(){_is_RLO_happening = true;}
    inline void unset_is_RLO_happening(){_is_RLO_happening = false;}
    inline void set_is_colliding(){_is_colliding = true;}
    inline void unset_is_colliding(){_is_colliding = false;}

    ///Mass transfer
    /**
     * Function that handle the situation in which the mass transfer is dynamic unstable but the stars are not mixing or doing a ce.
     * @param b  Pointer to the binary
     * @return EXIT_SUCCESS;
     */
    virtual int dynamic_swallowing(_UNUSED Binstar *b){return EXIT_SUCCESS;}


    ///q stuff
    /**
     * Simply estimate the mass fraction of the donor over the accretor
     * @param donor  Star that is donating mass through the RLO
     * @param accretor  Star that is accreting mass
     * @return mass_fraction
     */
    virtual double estimate_q(Star *donor, Star *accretor);

    ///Utilities
    /**
     * Set binary mix or common evelope to true depending on the stellar bse phases in case of
     * both stars are overfilling the RL
     * @param b Pointer to the binary
     * @return EXIT_SUCCESS
     */
    virtual int outcome_double_RLO(_UNUSED Binstar *b){svlog.critical("outcome_double_RLO not implemented for this RLO option ("+name()+")",__FILE__,__LINE__,sevnstd::rl_error()); return -1;}
    /**
     * Set binary mix or common evelope to true depending on the stellar bse phases in case of
     * collision at the pericenter
     * @param b Pointer to the binary
     * @return EXIT_SUCCESS
     */
    virtual int outcome_collision(_UNUSED Binstar *b){svlog.critical("outcome_collision not implemented for this RLO option ("+name()+")",__FILE__,__LINE__,sevnstd::rl_error()); return -1;}
    /**
     * Set binary mix or common evelope to true depending on the stellar bse phases in case of
     * collision at the pericenter
     * @param b Pointer to the binary
     * @return EXIT_SUCCESS
     */
    virtual int outcome_collision(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *b){svlog.critical("outcome_collision not implemented for this RLO option ("+name()+")",__FILE__,__LINE__,sevnstd::rl_error()); return -1;}

    /**
     * Estimate the H fraction based on the functional form .... (found in SEVN1)
     * @param s Pointer to the star
     * @return
     */
    virtual  double Hfrac(Star *s);

    /**
     * Estimate the rate of mass transfer considering the Eddington limit
     * @param donor  Star that is donating mass through the RLO
     * @param accretor  Star that is accreting mass
     * @return
     */
    virtual double dMdt_eddington(Star *donor, Star *accretor, Binstar *b);

    /**
     * Check if the star can be considered giant or giant-like in the BSE type classification
     * @param star_type_bse BSE star type
     * @return true if the type is between 2 and 9 (inclusive) but not 7, false otherwise
     */
    inline bool is_giant_like(int star_type_bse){

        if (star_type_bse>=2 && star_type_bse<=9 && star_type_bse!=7)
            return true;
        else
            return false;

        return false;
    }


    ///Mass transfer
    /**
     * Rate of dynamic  mass transfer
     * @param s pointer to the Star that is losing mass
     * @param b pointer to the binary
     * @return dm/dt or thrown an error
     *
     * @throws rl_error if the method is not implemented
     *
     */
    virtual double dynamic_dmdt(_UNUSED Star *s,_UNUSED Binstar *b){svlog.critical("dynamic_dmdmt not implemented for this RLO option ("+name()+")",__FILE__,__LINE__,sevnstd::rl_error()); return -1;};
    /**
     * Rate of nuclear  mass transfer
     * @param s pointer to the Star that is losing mass
     * @param b pointer to the binary
     * @return dm/dt or thrown an error
     *
     * @throws rl_error if the method is not implemented
     */
    virtual double nuclear_dmdt(_UNUSED Star *s,_UNUSED Binstar *b){svlog.critical("nuclear_dmdt not implemented for this RLO option ("+name()+")",__FILE__,__LINE__,sevnstd::rl_error()); return -1;}
    /**
     * Rate of nuclear  mass transfer
     * @param s pointer to the Star that is losing mass
     * @param b pointer to the binary
     * @return dm/dt or thrown an error
     *
     * @throws rl_error if the method is not implemented
     */
    virtual double thermal_dmdt(_UNUSED Star *s,_UNUSED Binstar *b){svlog.critical("thermal_dmdt not implemented for this RLO option ("+name()+")",__FILE__,__LINE__,sevnstd::rl_error()); return -1;}

    /**
     *
     * @param donor
     * @param accretor
     * @param b
     * @return
     */
    virtual double thermal_nuclear_accreted_mass(_UNUSED double DM_donor, _UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *b){svlog.critical("accreted mass not implemented for this RLO option ("+name()+")",__FILE__,__LINE__,sevnstd::rl_error()); return -1;}

    ///time scales
    /**
     * Kelvin-Helmholtz time  scale  from the modified classical expression  (Eq. 61, Hurley+2002)
     * Directly taken from SEVN1
     * @param s Pointer to star
     * @return Kelvin-Helmholtz time  scale in yr
     */
    virtual double kelvin_helmotz_tscale(Star *s);

    /**
     * Dynamic time  scale  (Eq. 63, Hurley+2002)
     * Directly taken from SEVN1
     * @param s Pointer to star
     * @return Dynamic time  scale in yr
     */
    virtual double dynamic_tscale(Star *s);

    ///Checks
    /**
     * Check that for compact remnants  (WD,NS,BH) the accretor is as compact or more compact with respect to the donor.
     * @param donor_type_bse  BSE type of the donor
     * @param accretor_type_bse  BSE type of the accretor
     * @return EXIT_SUCCESS or throw an error
     *
     * @throws sevnstd::rl_error
     */
    virtual int check_compact_accretor(int donor_type_bse, int accretor_type_bse){


        if (donor_type_bse>=10 && donor_type_bse<=12 && accretor_type_bse<10)
            svlog.critical("RLO with a donor with bse type "+utilities::n2s(donor_type_bse,__FILE__,__LINE__)+
                           " and accretor with bse type "+utilities::n2s(accretor_type_bse,__FILE__,__LINE__)+" is not expected!",__FILE__,__LINE__,sevnstd::rl_error());
        else if (donor_type_bse>=13 && accretor_type_bse<donor_type_bse)
            svlog.critical("RLO with a donor with bse type "+utilities::n2s(donor_type_bse,__FILE__,__LINE__)+
            " and accretor with bse type "+utilities::n2s(accretor_type_bse,__FILE__,__LINE__)+" is not expected!",__FILE__,__LINE__,sevnstd::rl_error());

        return EXIT_SUCCESS;

    }

private:


    static std::map<std::string, Orbital_change_RL *> & GetStaticMap(){
        static std::map<std::string, Orbital_change_RL *> _locmap;
        return _locmap;
    }
    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }
    double _DM[2]={0,0};
    bool _is_swallowed[2]={false, false}; /*!< Array of bool corresponding to the ID of stars in binary. If true, the star has been swallowed through the RLO*/
    bool mix = false;
    bool comenv = false;
    bool _is_RLO_happening=false;
    bool _is_colliding=false;




};

class disabled_rl :public Orbital_change_RL {

public:
    disabled_rl(bool reg = true){
        if (reg){
            Register(this, name());
        }
    }

    ~disabled_rl(){
        //std::cout<<"Orb RL dis "<<  this << std::endl;

    }

    disabled_rl * instance() override { return new disabled_rl(false);}

    static disabled_rl _disabled_rl;

    inline std::string name() override { return Lookup::rlmap_name.at(Lookup::RLMode ::_RLdisabled); }

    inline bool is_process_ongoing() override {return false;}



};

class Hurley_rl : public Orbital_change_RL {

public:

    Hurley_rl(bool reg=true){
        if (reg){
            Register(this, name());
        }
        //std::cout<<"Make Orb RL H "<<  this << std::endl;


    }

    Hurley_rl* instance() override  { return new Hurley_rl(false);}

    ~Hurley_rl(){
        //std::cout<<"Orb RL H "<<  this << std::endl;

    }

    static Hurley_rl _Hurley_rl;

    inline std::string name() override { return Lookup::rlmap_name.at(Lookup::RLMode ::_RLHurley); }

    void init(_UNUSED Binstar *b) override;

    double DA(_UNUSED Binstar *b, _UNUSED int procID) override;

    double DE(_UNUSED Binstar *b, _UNUSED int procID) override;

    double DAngMomSpin(_UNUSED Binstar *b, _UNUSED int procID, _UNUSED int starID) override;


protected:

    /**
     * Set the common variable used by DM, DE, DA assuming a given donor and accretor star
     * @param donor_star Pointer to the star that is filling the RL
     * @param accretor_star  Pointer to the star that is accreting mass
     * @param b  Pointer to the binary
     */
    void init_common_variable(Star *donor_star, Star *accretor_star, Binstar *b);

    ///AUX
    /**
     * Check if both stars are filling the RL
     * The outcome of this two events can be a Mix or a CE.
     * @param b  Pointer to the binary
     * @return true if one of two conditions is fulfilled
     */
    bool check_doubleRLO(Binstar *b);

    /**
     * Check if both stars are colliding at the periastron.
     * The outcome of this two events can be a Mix or a CE.
     * @param b  Pointer to the binary
     * @return true if one of two conditions is fulfilled
     */
    bool check_collision_at_periastron(Star *donor, Star *accretor, Binstar *b);

    /**
     * Set binary mix or common evelope to true depending on the stellar bse phases in case of
     * both stars are overfilling the RL
     * @param b Pointer to the binary
     * @return EXIT_SUCCESS
     */
    int outcome_double_RLO(_UNUSED Binstar *b) override;
    /**
     * Set binary mix or common evelope to true depending on the stellar bse phases in case of
     * collision at the pericenter. In that version the check is made in a symmetric way (both the donor and accretor have the same role)
     * @param donor_star Pointer to the star that is filling the RL
     * @param accretor_star  Pointer to the star that is accreting mass
     * @param b Pointer to the binary
     * @return EXIT_SUCCESS
     */
    int outcome_collision(_UNUSED Binstar *b) override;

    /**
     * Set binary mix or common evelope to true depending on the stellar bse phases in case of
     * collision at the pericenter. In that version the check is made in a asymmetric way, the donor is the star filling the Roche Lobe.
     * @param b Pointer to the binary
     * @return EXIT_SUCCESS
     */
    int outcome_collision(Star *donor, Star *accretor, _UNUSED Binstar *b) override;


    /**
     * Actual RocheLobe Overflow process as reported in Hurely+2002 and BSE
     * @param donor Pointer to the star that is overflowing the RocheLobe
     * @param accretor Pointer to the star that is accreting mass
     * @param b Pointer to the binary
     *
     * @note: Sometime the code does not strictly follows what reported in Hurley2002 paper, this because BSE code has been updated.
     * The current implementation is directly ported from SEVN1 with a modification made by MM on the fraction of the mass accreted that is now
     * tuned by the parameter rlo_f_mass_accreted.
     */
    void RLO(_UNUSED Binstar *b);

    virtual void set_orbital(_UNUSED Binstar *b);

    ///Mass transfer

    /**
     * Rate of dynamic  mass transfer from Eq. 62 in Hurley+2002
     * @param s pointer to the Star that is losing mass
     * @param b pointer to the binary
     * @return dm/dt
     *
     */
    double dynamic_dmdt(_UNUSED Star *s,_UNUSED Binstar *b) override;
    /**
     * Rate of dynamic  mass transfer from Eq. 60 in Hurley+2002
     * @param s pointer to the Star that is losing mass
     * @param b pointer to the binary
     * @return dm/dt
     *
     */
    double nuclear_dmdt(_UNUSED Star *s,_UNUSED Binstar *b) override;
    /**
     * Rate of dynamic  mass transfer from Eq. 58 and Eq.59 in Hurley+2002
     * @param s pointer to the Star that is losing mass
     * @param b pointer to the binary
     * @return dm/dt
     *
     * @note the implementation is a bit different with respect to what is reported in Hurley+2002, but it comes directly from BSE
     *
     */
    double thermal_dmdt(_UNUSED Star *s,_UNUSED Binstar *b) override;
    /**
     * Function that handle the situation in which the mass transfer is dynamic unstable but the stars are not mixing or doing a ce.
     * Note in that case we use the current stellar properties (after single and binary evolution) not the properties at time 0.
     * @param b  Pointer to the binary
     * @return EXIT_SUCCESS;
     */
    int dynamic_swallowing(_UNUSED Binstar *b) override;

    /**
     * Handle the mass accreted on the companion on thermal or nuclear timescale
     * It takes into account also special cases such as accretionon compact objects.
     * The basic accreted mass is simply DM_accreted=fMT*DM_donor (fmt is the runtime parameter rlo_f_mass_accreted)
     * @param DM_donor  Total amount of mass lost by the donor star [Msun]
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param b Pointer to the binary
     * @return Final amount of mass accreted on the accretor
     */
    double thermal_nuclear_accreted_mass(_UNUSED double DM_donor, _UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *b) override;


    /**
     * Roche Lobe radius normalised over the star radius
     * @param star pointer to the star for which we are estimateting fRL
     * @param b pointer to the binary
     * @return R_rl/R_star
     */
    double fRL_radius(_UNUSED Star* star,_UNUSED Binstar *b) override;


    /**
     *
     * @return
     */

    Star* get_donor(){ return  donor;}
    Star* get_accretor(){ return  donor;}
    double get_rldonor(){ return rl1;}
    void set_DA(double DA){_DA=DA;}
    void set_DE(double DE){_DE=DE;}

    inline double get_mt_value() const override {return  mt_value;}

    inline double get_mt_tshold() const override {return  mt_tshold;}

    /***Define common variables ***/
    //Stars
    Star *donor;
    Star *accretor;
    //BSE type (1 for donor, 2 for accretor)
    int star_type1 = utilities::NULL_INT;
    int star_type2 = utilities::NULL_INT;
    //Radii
    double frl1    = utilities::NULL_DOUBLE; // RL over R1
    double rl1     = utilities::NULL_DOUBLE; //Radius of the Roche Lobe for the donor
    bool is_RL_in_core = false; //Is the RL inside the Core?
    ///mt stuff - NOTICE, the interpretation of these two values depend on the mt_stability option
    ///if the option is from the qcrit family, mt_value=q and mt_tshold=qcrit,
    ///if the option is from the zeta family, mt_value=zeta_adiabatic and mt_tshold=zeta_RL
    double mt_value = utilities::NULL_DOUBLE;  //Mdonor/Maccretore
    double mt_tshold = utilities::NULL_DOUBLE; //qcrit for dynamic mass transfer

    ///Accretion efficency
    double f_MT = utilities::NULL_DOUBLE; //Mass accretion efficenciecy
    double dmdt_edd = utilities::NULL_DOUBLE; //Eddingtont rate of the accretion on the accretor

    ///Times
    double dt    = utilities::NULL_DOUBLE;//Equations are in yr
    double tkh1  = utilities::NULL_DOUBLE;//KH timescale of the star filling the RL
    double tdyn1 = utilities::NULL_DOUBLE; //Dynamical timescale of the star filling the RL
    double taum  = utilities::NULL_DOUBLE; //Time scale for dynamic Mass Loss through RL //TODO Why? Why not just tdyn1 as in Hurley, Time scale for Mass loss through
    double tkh2  = utilities::NULL_DOUBLE; //KH timescale of the accretor

    //Orbital
    double _DA = 0;
    double _DE = 0;

    //Spin
    double _DLspin[2]={0.,0.};


    /********************************/



};

class Hurley_rl_bse : public Hurley_rl{
public:

    Hurley_rl_bse(bool reg=true){
        if (reg){
            Register(this, name());
        }
    }

    Hurley_rl_bse* instance() override  { return new Hurley_rl_bse(false);}

    ~Hurley_rl_bse(){
    }

    static Hurley_rl_bse _Hurley_rl_bse;

    inline std::string name() override {return Lookup::rlmap_name.at(Lookup::RLMode ::_RLHurleyBSE);}

protected:
    /**
     * Handle the mass accreted on the companion on thermal or nuclear timescale
     * It takes into account also special cases such as accretionon compact objects.
     *
     * @param DM_donor  Total amount of mass lost by the donor star [Msun]
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param b Pointer to the binary
     * @return Final amount of mass accreted on the accretor
     */
    double thermal_nuclear_accreted_mass(_UNUSED double DM_donor, _UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *b) override;


    /**
     * Function that handle the situation in which the mass transfer is dynamic unstable but the stars are not mixing or doing a ce.
     * @param b  Pointer to the binary
     * @return EXIT_SUCCESS;
     */
     int dynamic_swallowing(_UNUSED Binstar *b) override;

};

class Hurley_rl_propeller : public  Hurley_rl{

public:

    Hurley_rl_propeller(bool reg=true){
        if (reg){
            Register(this, name());
        }
    }

    Hurley_rl_propeller* instance() override  { return new Hurley_rl_propeller(false);}

    inline std::string name() override { return Lookup::rlmap_name.at(Lookup::RLMode ::_RLHurleyprop); }

    static Hurley_rl_propeller _RLHurleyprop;

protected:

    /**
     *  Disable accreation on NS for propeller mechanism (Campana+18, https://www.aanda.org/articles/aa/pdf/2018/02/aa30769-17.pdf)
     * @param DM_accreted mass accreted by the accretor (before this correction) in Msun
     * @param DM_donor mass donated by the donor in Msun
     * @param donor Pointer to the star that is overflowing the RocheLobe
     * @param accretor Pointer to the star that is accreting mass
     * @param b Pointer to the binary
     */
    double correct_accreted_mass(double DM_accreted, _UNUSED double DM_donor,
                                         _UNUSED Star *donor,_UNUSED Star *accretor,_UNUSED Binstar *b)  override;

};


class Hurley_mod_rl : public Hurley_rl {

public:

    Hurley_mod_rl(bool reg=true){
        if (reg){
            Register(this, name());
        }
    }

    Hurley_mod_rl* instance() override  { return new Hurley_mod_rl(false);}

    inline std::string name() override { return Lookup::rlmap_name.at(Lookup::RLMode ::_RLHurleymod); }

    static Hurley_mod_rl _Hurley_mod_rl;

protected:

    void set_orbital(_UNUSED Binstar *b) override;

};
/*************************************/

/**********Mix*****************/
class Orbital_change_Mix : public Orbital_change{

public:

    Orbital_change_Mix(){}

    virtual ~Orbital_change_Mix(){
        //std::cout<<"Orbital_change_Mix "<<  this << std::endl;

    }

    virtual Orbital_change_Mix *instance(){ return nullptr;}

    static Orbital_change_Mix  *Instance(const std::string &name);

    inline std::string name() override {return "Orbital change Mix";}

protected:

    void Register(Orbital_change_Mix *ptr, const std::string &_name) {
        GetStaticMap().insert(std::make_pair(_name, ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }


private:

    static std::map<std::string, Orbital_change_Mix *> & GetStaticMap(){
        static std::map<std::string, Orbital_change_Mix *> _locmap;
        return _locmap;
    }

    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }

};

class disabled_mix :  public Orbital_change_Mix{

public:
    disabled_mix(bool reg = true){
        if (reg){
            Register(this, name());
        }
    }

    disabled_mix * instance() override { return new disabled_mix(false);}

    static disabled_mix _disabled_mix;

    inline std::string name() override { return Lookup::mixmap_name.at(Lookup::MixMode::_Mixdisabled); }

    inline bool is_process_ongoing() override {return false;}


};

class simple_mix :  public Orbital_change_Mix{

public:
    simple_mix(bool reg = true){
        if (reg){
            Register(this, name());
        }
    }

    ~simple_mix(){
        //std::cout<<"simple mix "<<  this << std::endl;

    }

    simple_mix * instance() override { return new simple_mix(false);}

    static simple_mix _simple_mix;

    inline std::string name() override { return Lookup::mixmap_name.at(Lookup::MixMode::_Mixsimple); }



};

/**********SN kicks*****************/
class Orbital_change_SNKicks : public Orbital_change {

public:

    Orbital_change_SNKicks(){}

    virtual ~Orbital_change_SNKicks(){
        //std::cout<<"orb sn "<<  this << std::endl;

    }

    virtual Orbital_change_SNKicks *instance(){ return nullptr;}

    static Orbital_change_SNKicks  *Instance(const std::string &name);

    inline std::string name() override {return "Orbital change SNKicks";}

    inline double get_cos_nu() const { return cos_nu;}
    inline double get_vcom() const { return vcom;}

protected:

    void Register(Orbital_change_SNKicks *ptr, const std::string &_name) {
        GetStaticMap().insert(std::make_pair(_name, ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }

    inline void set_cos_nu(double _cos_nu)  {  cos_nu=_cos_nu; return;}
    inline void set_vcom(double _vcom)  {  vcom=_vcom; return;}

private:

    static std::map<std::string, Orbital_change_SNKicks *> & GetStaticMap(){
        static std::map<std::string, Orbital_change_SNKicks *> _locmap;
        return _locmap;
    }

    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }

    double cos_nu; //Angle between old and new ang mom
    double vcom; //New centre-of-mas velocity

};

/***
 * Implement disabled SNkicks, i.e. does not apply the kicks and does not change orbital parameters
 */
class disabled_SNKicks :  public Orbital_change_SNKicks{

public:
    disabled_SNKicks(bool reg = true){
        if (reg){
            Register(this, name());
        }
    }

    disabled_SNKicks * instance() override { return new disabled_SNKicks(false);}

    static disabled_SNKicks _disabled_SNKicks;

    inline std::string name() override { return Lookup::snkmap_name.at(Lookup::SNKickMode::_SNKickdisabled); }

    inline bool is_process_ongoing() override {return false;}

};

class Hurley_SNKicks : public Orbital_change_SNKicks{

public:
    Hurley_SNKicks(bool reg = true){
        if (reg){
            Register(this, name());
        }
    }

    ~Hurley_SNKicks(){
        //std::cout<<"orb sn  hu "<<  this << std::endl;

    }

    Hurley_SNKicks * instance() override { return new Hurley_SNKicks(false);}

    static Hurley_SNKicks _Hurley_SNKicks;

    inline std::string name() override { return Lookup::snkmap_name.at(Lookup::SNKickMode::_SNKickHurley); }

    void init(_UNUSED Binstar *b) override;

    double DA(_UNUSED Binstar *b, _UNUSED int procID) override;

    double DE(_UNUSED Binstar *b, _UNUSED int procID) override;

protected:

    /**
     * Estimate the final orbital paratmers after a Sn kick using the Hurley+02 formalism (see Appendix A1)
     * @param sn Pointer to the star that is exploding as SN
     * @param other Pointer to the other star
     * @param binstar Pointer to the binary system
     * @param a_fin Variable where to store the final semimajor axis
     * @param ecc_fin Variable where to store the final eccentricity
     * @return EXIT_SUCCESS
     */
    int kick_star(Star *sn, Star *other, Binstar *binstar, double &a_fin, double &ecc_fin);

    //double Vorb(double Mass1, double Mass2, );

    const double G = utilities::G;

private:

    std::uniform_real_distribution<double> _uniform_real;

    inline void set_DA(double DA){_DA=DA;}
    inline void set_DE(double DE){_DE=DE;}

    //Orbital
    double _DA = 0;
    double _DE = 0;
};
/*************************************/

/**********Common Envelope*****************/
class Orbital_change_CE : public  Orbital_change {

public:

    Orbital_change_CE(){}

    virtual ~Orbital_change_CE(){
        //std::cout<<"orb ce "<<  this << std::endl;

    }

    virtual Orbital_change_CE *instance(){ return nullptr;}

    static Orbital_change_CE  *Instance(const std::string &name);

    inline std::string name() override {return "Orbital change CE";}

protected:

    void Register(Orbital_change_CE *ptr, const std::string &_name) {
        GetStaticMap().insert(std::make_pair(_name, ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }

private:

    static std::map<std::string, Orbital_change_CE *> & GetStaticMap(){
        static std::map<std::string, Orbital_change_CE *> _locmap;
        return _locmap;
    }

    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }
};

class energy_CE : public  Orbital_change_CE{
public:
    energy_CE(bool reg = true){
        if (reg){
            Register(this, name());
        }
    }

    ~energy_CE(){
        //std::cout<<"orb CE  energy "<<  this << std::endl;

    }

    energy_CE * instance() override { return new energy_CE(false);}

    static energy_CE _energy_CE;

    inline std::string name() override { return Lookup::cemap_name.at(Lookup::CEMode::_CEEnergy); }

};

class disabled_CE : public  Orbital_change_CE{
public:
    disabled_CE(bool reg = true){
        if (reg){
            Register(this, name());
        }
    }

    ~disabled_CE(){
        //std::cout<<"orb CE  disabled "<<  this << std::endl;

    }

    disabled_CE * instance() override { return new disabled_CE(false);}

    static disabled_CE _disabled_CE;

    inline std::string name() override { return Lookup::cemap_name.at(Lookup::CEMode::_CEEnergy); }

};
/*************************************/


#endif //SEVN_ORBIT_H
