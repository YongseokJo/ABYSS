//
// Created by spera on 13/02/19.
//

#ifndef SEVN_BINPROP_H
#define SEVN_BINPROP_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <random>
#include <utilities.h>
#include <lookup_and_phases.h>

#include <sevnlog.h>
using sevnstd::SevnLogging;


#define _UNUSED __attribute__ ((unused))

class Binstar;
class Star;

class BinaryProperty{

public:
    BinaryProperty(){
        size++;
    }

    virtual ~BinaryProperty(){}

    static std::vector<BinaryProperty*> all;
    typedef std::map<std::string,size_t> _PrintMap;
    static _PrintMap PrintMap;

    virtual BinaryProperty * Instance(){ return nullptr; }
    virtual inline std::string name(){return "Binary property (generic)";}
    virtual inline std::string units()  {return "UNKNOWN";}

    virtual inline int getID() const {return -1;}

    virtual void resynch(_UNUSED const double &dt, _UNUSED bool set0=true) {}

    virtual inline void check_boundaries(_UNUSED double *val, _UNUSED Binstar *binstar) {}

    virtual void init(_UNUSED const double par){
        set(par);
        set_0(par);
    };

    virtual void init_derived(_UNUSED Binstar *binstar){};

    virtual int evolve(_UNUSED Binstar *binstar);

    virtual int special_evolve(_UNUSED Binstar *binstar){
        svlog.info(" my special evolve ");
        return 0;
    } //if you have some properties that you want to evolve outside the main cycle

    //Just return the s
    double get_variation(Binstar *binstar) const;


    virtual int set_broken(_UNUSED Binstar *binstar){
        V0=V=std::nan("");
        return EXIT_SUCCESS;
    }

    virtual inline double get(_UNUSED const Binstar* b=nullptr) {return V;}
    inline double get_0(_UNUSED const Binstar* b=nullptr) {return V0;}
    virtual inline void restore(){ V = V0;}

private:
    static size_t size;



protected:
    void Register(BinaryProperty *_p, size_t *id, const std::string &_name){
        BinaryProperty::all.push_back(_p);
        *id = BinaryProperty::size - 1;
        BinaryProperty::PrintMap.insert(std::pair<std::string,size_t>(_name, *id));
        svlog.debug("Binary property "+ name() + " registered" + " (Nproperties: " + utilities::n2s(all.size(),__FILE__,__LINE__) + ")");
    }

    std::uniform_real_distribution<double> _uniform_real;

    bool isbad(const double &val) {
        return (std::isnan(val) || std::isinf(val));
    }

    virtual void set(const double &a) { V = a;}
    virtual void set_0(const double &a) { V0 = a;}
    SevnLogging svlog;
    double V, V0;
};

///Fundamental
class Eccentricity : public BinaryProperty{

public:
    Eccentricity(bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Eccentricity _eccentricity;

    Eccentricity *Instance() override {
        return (new Eccentricity(false));
    }

    inline std::string name() override { return "Eccentricity"; }
    inline std::string units()  override {return "";}
    inline int getID() const override  {return (int)ID;}

    void init(_UNUSED const double par) override {

        if (par<0 or par>=1){
            svlog.critical("Error in eccentricity initialisation, it should be between 0 and 1(excluded), you use e="
            +utilities::n2s(par,__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::sevnio_error(""));
        }

        BinaryProperty::init(par);
    }


    int evolve(_UNUSED Binstar *binstar) override;


    inline void check_boundaries(double *val,  _UNUSED Binstar *binstar) override;

};

class Semimajor : public BinaryProperty{

public:
    Semimajor(bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Semimajor _semimajor;

    Semimajor *Instance() override {
        return (new Semimajor(false));
    }

    inline std::string name() override { return "Semimajor"; }
    inline std::string units() override  {return "Rsun";}

    inline int getID() const override {return (int)ID;}

    void init(_UNUSED const double par) override {

        if (par<=0){
            svlog.critical("Error in Semimajor initialisation, it should be larger than 0, you use Semimajor="
                           +utilities::n2s(par,__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::sevnio_error(""));
        }

        BinaryProperty::init(par);
    }


    /**
     * Alternative to classic check_boundary.
     * If a<=0 it set a to the sum of Radii of the single stars
     * @param val Value to check
     * @param binstar Pointer to the binstar
     */
    void check_boundaries(double *val, Binstar *binstar) override;




};



///Derived
/**
 * Base class for properties that are derived from other properties
 */

class Derived_Property_Binary : public BinaryProperty{


    void init_derived(_UNUSED Binstar *binstar) override {
        evolve(binstar); //Set V
        set_0(get()); //Set V_0
    };

    void init(_UNUSED const double par) override {
            svlog.critical("Cannot initialise the derived Binary property"
                           +name(),__FILE__,__LINE__,sevnstd::sevnio_error(""));
    };


    virtual inline bool amiderived() { return  true;}

};

class dadt : public Derived_Property_Binary{

public:

    dadt(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    inline std::string name() override { return "dSemimajordt"; }
    inline std::string units() override  {return "Rsun/Myr";}


    static size_t ID;
    static dadt _dadt;

    dadt * Instance() override {
        return (new dadt(false));
    }

    inline int getID() const override {return (int)ID;}
    int evolve(_UNUSED Binstar *binstar) override;

};

class dedt : public Derived_Property_Binary{

public:

    dedt(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    inline std::string name() override { return "dEccentricitydt"; }
    inline std::string units() override  {return "/Myr";}


    static size_t ID;
    static dedt _dedt;

    dedt * Instance() override {
        return (new dedt(false));
    }

    int evolve(_UNUSED Binstar *binstar) override;
    inline int getID() const override {return (int)ID;}

};

class AngMom : public Derived_Property_Binary{

public:

    AngMom(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    inline std::string name() override { return "AngMom"; }
    inline std::string units() override  {return "MsunRsun^2";}

    static size_t ID;
    static AngMom _angmom;

    AngMom * Instance() override {
        return (new AngMom(false));
    }

    int evolve(_UNUSED Binstar *binstar) override;
    inline int getID() const override {return (int)ID;}

protected:

    /**
     * Estimate the orbital angular momentum from orbital properties
     * @param a semi_major axis [Rsun]
     * @param e eccentricity
     * @param M1  Mass of the first star [Msun]
     * @param M2  Mass of the second star  [Msun]
     * @return  Orbital angular momentum in Rsun^2 Msun yr^-1
     */
    inline double AngMom_from_orbit(double a, double e, double M1, double M2){
        double G = utilities::G;
        return  (M1*M2)*std::sqrt( (a*(1-e*e)*G)/(M1+M2));
    }


};

class Period : public Derived_Property_Binary{

public:

    Period(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    inline std::string name() override { return "Period"; }
    inline std::string units() override  {return "yr";}

    static size_t ID;
    static Period _period;

    Period * Instance() override {
        return (new Period(false));
    }

    int evolve(_UNUSED Binstar *binstar) override;
    inline int getID() const override {return (int)ID;}

protected:

    /**
     * Estimate the binary period (Kepler Third Law)
     * @param a Semimajor axis in Rsun
     * @param Mtot  Total mass of the binary (M1+M2) in Msun
     * @return Period in years
     */
    inline double calc_period(const double & a, const double & Mtot){
        return 2.0*M_PI*std::sqrt( (a*a*a) / (utilities::G * (Mtot)) ); //returns the period in years
    }


};

class GWtime : public Derived_Property_Binary{

public:

    GWtime(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    GWtime *Instance() override {
        return (new GWtime(false));
    }


    inline std::string name() override { return "GWtime"; }
    inline std::string units() override  {return "Myr";}


    static size_t ID;
    static GWtime _gwtime;

    inline int getID() const override {return (int)ID;}

    int evolve(_UNUSED Binstar *binstar) override;

protected:

    const double time_scaling = utilities::yr_to_Myr*0.01953125/utilities::G3_over_c5; //tGW=time_scaling*f(a,e,m1,m2) [Myr]
    //time scaling=5/256 * (c^5/G^3), see evolve for f.

};

class _RL : public Derived_Property_Binary{

public:

    _RL(){}

    _RL * Instance() override { return nullptr; }

    inline std::string name() override { return "_RL";}
    inline std::string units() override  {return "Rsun";}

    /**
     * Estimate the RL follwing the Eggleton formalism (Eq. 53 Hurley02)
     * @param primary  Pointer to the primary star (the star for which we are estimating the RL)
     * @param secondary  Pointer to the secondary star.
     * @param b  Pointer to the binary object
     * @return  Roche Lobe radius in Rsun
     */
    virtual double RL_Eg(Star* primary, Star* secondary, Binstar* b);


    int evolve(_UNUSED Binstar *binstar) override { return 0; };

};

class RL0 : public _RL{

public:
    RL0(bool reg=true){
        if(reg){
            Register(this, &ID, name());
        }
    }

    RL0 *Instance() override{
        return (new RL0(false));
    }

    inline std::string name() override { return "RL0";}
    static size_t ID;
    static RL0 _rl0;

    inline int getID() const override { return (int)ID;}

    int evolve(_UNUSED Binstar *binstar) override;


};

class RL1 : public _RL{

public:
    RL1(bool reg=true){
        if(reg){
            Register(this, &ID, name());
        }
    }

    RL1 *Instance() override{
        return (new RL1(false));
    }

    inline std::string name() override { return "RL1";}
    static size_t ID;
    static RL1 _rl1;

    inline int getID() const override { return (int)ID;}

    int evolve(_UNUSED Binstar *binstar) override;


};


//Time handling
class BWorldtime : public BinaryProperty {

public:
    BWorldtime(bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }
        set_0(0.0);
        set(0.0);
    }

    BWorldtime *Instance() override {
        return (new BWorldtime(false));
    }

    inline std::string name() override { return "BWorldtime"; }
    inline std::string units() override  {return "Myr";}

    static size_t ID;
    static BWorldtime _bworldtime;

    inline int getID() const override {return (int)ID;}

    int evolve(_UNUSED Binstar *binstar) override {return 0;}

    int special_evolve(_UNUSED Binstar *binstar) override;

    int set_broken(_UNUSED Binstar *binstar) override { return EXIT_SUCCESS;}

};

class BTimestep : public BinaryProperty{

public:
    BTimestep(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
        set(1.0e30); //Initial timestep set to a large value
        set_0(0.0); //Initial last timestep set to 0

    }

    BTimestep * Instance() override {
        return (new BTimestep(false));
    }

    inline std::string name() override { return "BTimestep"; }
    inline std::string units() override  {return "Myr";}

    static size_t ID;
    static BTimestep _btimestep;


    double max_timestep(_UNUSED Binstar *s);

    inline bool check_repeat(_UNUSED Binstar *binstar, const double &V0, const double &V, const double &derivative);

    inline bool check_repeat_eccentricity(_UNUSED Binstar *binstar, const double &m0, const double &m, const double &derivative);



    /**
     * Check if the Mass transfer during the RLO is too much with respect to the maximum variation
     * @param binstar  Pointer to the binary system
     * @return true if the system has to repeat the evolution with a new timestep, false otherwise
     */
    bool check_repeat_DM_RLO(_UNUSED Binstar *binstar);

    /**
     * Check if the SpinUp or Spindown of the NS Spin has varied too much with respect to the maximum variation
     * @param binstar Pointer to the binary system
     * @return true if the system has to repeat the evolution with a new timestep, false otherwise
     */
    bool check_repeat_OmegaRem_NS(_UNUSED Binstar *binstar);

    /**
     * Check the variation of stellar rotations due to binary processes
     * @param binstar Pointer to the binary system
     * @return true if the system has to repeat the evolution with a new timestep, false otherwise
     */
    bool check_repeat_stellar_rotation(_UNUSED Binstar *binstar);

    int evolve(_UNUSED Binstar *binstar) override;

    inline int getID() const override {return (int)ID;}

    void resynch(const double &dt, bool set0) override {
        set(dt);
        if (set0) set_0(get());
    }

    int set_broken(_UNUSED Binstar *binstar) override {
        set(max_timestep(binstar));
        return EXIT_SUCCESS;
    }

    int special_evolve(_UNUSED Binstar *binstar) override {
        set_0(get());
        set(max_timestep(binstar));
        return 0;
    } //if you have some properties that you want to evolve outside the main cycle


protected:

    /**
     * Check if the proposed dt is larger than max_dt or lower than min_dt.
     * In that cases force dt to be one of the limits.
     * @param dt  reference to the timestep to analyse
     * @param s  pointer to the star
     */
    void check_dt_limits(double &dt, Binstar *binstar);

    double tiny_dt(Binstar *binstar);
};

class BEvent : public BinaryProperty{

public:

    BEvent(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
        init(Lookup::EventsList::NoEvent); //Initiliase to no event
    }

    inline std::string name() override { return "BEvent"; }
    inline std::string units() override  {return "";}

    static size_t ID;
    static BEvent _bevent;

    BEvent * Instance() override {
        return (new BEvent(false));
    }

    int evolve(_UNUSED Binstar *binstar) override{set((double)Lookup::EventsList::NoEvent); return EXIT_SUCCESS;}

    int special_evolve(_UNUSED Binstar *binstar) override;


};



/**
 * Derived class for properties that are derived from other properties
 * JIT means Just in Time, these are properties that do not need to be
 * evolved at each time step, but just to be computed when needed.
 */
class BJIT_Property : public Derived_Property_Binary{

public:

    inline std::string name()  override { return "Generic JIT Property"; }

    //Special evolve
    int evolve(_UNUSED Binstar *binstar) override {
        evolve_number++;
        return EXIT_SUCCESS;
    }


    inline double get(_UNUSED const Binstar* b=nullptr) override {
        //If evolve_number>last_evolve_number, we have to obtain a new estimate for lambda
        if(new_estimate_needed()){
            specialised_set(b);
        }
        return V;
    }

    inline void restore() override {evolve_number++;}


protected:


    virtual void specialised_set(_UNUSED const Binstar* b=nullptr){
        return;
    }

    bool new_estimate_needed(){
        return evolve_number>last_evolve_number;
    }

    void set(const double &a) override  {
        V=a;
        last_evolve_number=evolve_number;
    }



private:


    //The following  numbers check if an evolution of a binary evolution has been called
    //from the last time this property has been estimated.
    //evolve_number increases by 1 each time evolve or update_from_binary has been called.
    //last_evolve_number stores the value of evolve_number when getp is called.
    //When calling getp it checks if the two numbers are equal, if they are the property is just last stored one.
    unsigned int last_evolve_number=0, evolve_number=0;

};





#endif //SEVN_BINPROP_H
