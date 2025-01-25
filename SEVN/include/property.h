//
// Created by mario on 09/02/19.
//
/*
Adding a new Property

        As an example we are adding a very simple property the Schwarzwald radius.

* Go to Property.h and add the new class for the property

class Rs : public Property{

};

* Add the public static members ID and obejct instance.

static size_t ID;
static Rs _rs;

* Add the constructor  (it is standard)

Rs(bool reg = true){
    if(reg) {
        Register(this, &ID, name());
    }
}

* Override the basic public  method Instance

Rs * Instance() override {
return (new Rs(false));
}

* Override the basic public method name

inline std::string name() const override { return "Rs"; }

* Override the public method evolve, this is the core of the property, it should be defined in property.h, but written in propery.cpp because it depends on the star s that has been only pre-called here (to avoid cross reference).

void evolve(_UNUSED Star *s) override;

* In property.cpp  write the evolve function

void Rs::evolve( Star *s){

    set_0(get());

    double rs=2*utilities::G_over_c2*s->getp(Mass::ID); //rs=2GM/c^2

    set(rs);
}
NOTE: remember to write set_0(get) at the beginning.

In that case is simple we are just setting rs=2*GM/c^2.


                                                     * If needed Override the public method set_remant. How this property change when the star become a remnant? The default set_remant in the class Property set all the values V,V0,value,value0,Dvalue to 0. In our case we don’t want this, rather since Rs, depends just on the mass, we use exactly the same set_remannt of mass

void set_remnant(_UNUSED Star *s) override {
V0=V;
value0=value;
}

* If needed Override the public method changed_track.  How this property change when we jump to a new track? The default method in Property does nothing. Since this a property derived from the real value V of the Mass this does not change during the jump we don’t need to override.


* Other function to override
* Special_evolve, if they ave to evolve outside the cycle see ecc. Etc.
* set_w, only for interpolating property
* get_ecc only if internally the function is for example in log (see Luminosity)

Add in static_main.h

        * Go in static_main.h and add the static members
* Add the the ID  size_t Rs::ID;
* Add the static instance
        Rs Rs::_rs;

NOTE: The position in which we define the static instance is importante, because this is the order they will have during the evolution.
Therefore NEVER write something after Timestep, and for derived property NEVER write  before the property the are derived from , follow the comments in static_main.h
*/

#ifndef SEVN_PROPERTY_H
#define SEVN_PROPERTY_H

#include <sevnlog.h>
#include <iostream>
#include <vector>
#include <errhand.h>
#include <lookup_and_phases.h>
#include <cmath>
#include <sevnlog.h>
#include <utilities.h>
#include <lambda_base.h>
#include <lambda_klencki21.h>

using namespace Lookup;

using std::cout;
using std::endl;
using std::vector;
using sevnstd::SevnLogging;


#define _UNUSED __attribute__ ((unused))
#define FOR4 for(size_t _i = 0; _i < 4; _i++)

class Star;
class Staremnant;
class Binstar;

class Property{

public:

    typedef std::map<std::string,size_t> _PrintMap;
    static _PrintMap PrintMap;

    Property(){
        _size++;
        value = value0 = Dvalue = V = V0 = 0.0; //NB DO NOT CHANGE THIS
    }


    virtual ~Property(){

    }

    virtual Property * Instance() = 0;

    virtual inline std::string name() const {return "Property (generic)";}

    virtual inline std::string units(){return "UNKOWN";}


    /**
     * Check if this is a derived. A property is defined dervied if it is  estimated by the
     * combination of other properties and it is not the derivative of a property.
     * @return true if it is derived, false if it is not (overrided in the derived properties)
     */
    virtual inline bool amiderived() { return  false;}

    virtual void resynch(_UNUSED const double &dt, _UNUSED bool set0=true) {}

    virtual void resynch(_UNUSED Star *s) {}

    virtual void set_empty(_UNUSED Star *s) {
        V0=V=std::nan("");
        value=value0=0;
        Dvalue=0;
    }
    /**
     * Set empty method to be used during BSE evolution
     * @param s pointer to Star
     * @param b pointer to Binary
     */
    virtual void set_empty_in_bse(_UNUSED Star *s, _UNUSED Binstar *b) { set_empty(s);}
    virtual void evolve_empty(_UNUSED Star *s) {
        set_empty(s);
    }

    virtual void set_remnant(_UNUSED Star *s){

        V0 = V;
        value0 = value;

        V=value=0;
        Dvalue=0;
    }
    /**
     * Set remnant method to be used during BSE evolution
     * @param s pointer to Star
     * @param b pointer to Binary
     */
    virtual void set_remnant_in_bse(_UNUSED Star *s, _UNUSED Binstar *b) { set_remnant(s);}
    virtual void evolve_remnant(_UNUSED Star *s){V0=V;}
    virtual void evolve_nakedco(_UNUSED Star *s){V0=V;}

    //Main evolve function... the sub-(virtual)functions can be adapted for each Property
    virtual void evolve(_UNUSED Star *s){

        evolve_fake(s);
        correct_interpolation_errors(s);
        update_variation();



        evolve_real();
        correct_interpolation_errors_real(s);

        safety_check();

    }



    virtual size_t TabID() { return 0;}

    virtual void special_evolve(_UNUSED Star *s){
        //cout<<" my special evolve "<<endl;
    } //if you have some properties that you want to evolve outside the main cycle


    virtual void evolve_fake(Star *s);


    virtual void init(_UNUSED const double &a) {svlog.critical("Init not allowed for this property",__FILE__,__LINE__,sevnstd::sevnerr(""));}

    virtual void evolve_real(){

        if(value0 == 0.0 && value != 0.0)
            V = value;
        else {
            V0 = V;
            V = V0 + Dvalue * V0;
        }



        //svlog.pdebug("V",V,"V0",V0,"value",value,"value0",value0,"Dvalue",Dvalue,__FILE__,__LINE__);
        //svlog.pdebug("V",get(),"V0",get_0(),"value",get_fk(),"value0",get_0_fk(),"Dvalue",Dvalue,__FILE__,__LINE__);

    }

    virtual void update_from_binary(_UNUSED Star* s, const double &DV, _UNUSED  Binstar *b= nullptr) {

        if (std::isnan(DV) or std::isinf(DV))
            svlog.critical("Update from Binary of property " +
                           name()+" is nan or infinite",__FILE__,__LINE__);

        V = V + DV;

    };

    /** This function can be called to update the value of V in derived proprties (without setting V0)
     * It is empty for all the others. This is needed since after a binary evolution, the derived properties have to be updated
     * without setting V0 (same as update_from_binary).
     */
    virtual void update_derived(_UNUSED Star *s){}

    virtual void copy_V_from(Property *p){

        if (typeid(*p).name() == typeid(*this).name())
            set(p->get());
        else
            svlog.critical("Copying from property " + p->name() + " to property " + this->name() + " is not allowed",__FILE__,__LINE__,sevnstd::sevnerr());
    }

    virtual inline void update_variation() {

        if (value0 != 0.0)
            Dvalue = (value - value0) / value0;
        else
            Dvalue = 0.0;
    }

    virtual void correct_interpolation_errors(_UNUSED Star *s){}
    virtual void correct_interpolation_errors_real(_UNUSED Star *s){}

    virtual void restore(){
        V = V0;
        value = value0;
    }

    virtual void synch(){
        V = V0 = value0 = value;
    }

    /**
     * Reset to 0
     */
    virtual inline void reset(){V=V0=value=value0=0;}

    /**
     * Handle the modification of the property after a track change.
     * It does nothing by default.
     * @param s  The star for which we are changing the tracks.
     * @param s_jtrack  The star that we are using to get info of the new tracks.
     */
    virtual inline void changed_track(_UNUSED Star* s, _UNUSED Star* s_jtrack){};

    ///Here Star is optional, because it is needed only by property that have to estimate V on the fly.
    ///Notice that const Star*s means that inside the function we can only call method that we are sure will not
    ///change the Star s.
    virtual inline double get_fk(_UNUSED const Star* s=nullptr) {return value;}
    virtual inline double get_0_fk(_UNUSED const Star* s=nullptr) {return value0;}
    virtual inline double get(_UNUSED const Star* s=nullptr) {return V;}
    virtual inline double get_0(_UNUSED const Star* s=nullptr) {return V0;}

    //Raw get value, with these functions we are sure that the return is always the raw value of V,V0, value and value0
    //They are not virtual and cannot be changed.
    inline double get_fk_raw(_UNUSED const Star* s=NULL) const {return value;}
    inline double get_0_fk_raw(_UNUSED const Star* s=NULL) const {return value0;}
    inline double get_raw(_UNUSED const Star* s=NULL) const {return V;}
    inline double get_0_raw(_UNUSED const Star* s=NULL) const {return V0;}

    /**
     * Get the value of the property expected at the end of the evolution (before becoming a remnant).
     * It is based on the last values of the tables.
     * @param s Pointer to the star
     * @return The value of the property at the end of the evolution estimated as the weighted mean of the interpolating track values
     */
    virtual inline double get_Vlast(_UNUSED const Star* s) const {svlog.critical("get_vend not implemented for property"+
    name()+".",__FILE__,__LINE__,sevnstd::notimplemented_error()); return 0.;}

    double * get_wm() {return &wM[0];}
    double * get_wz() {return &wZ[0];}

    virtual void set_w(_UNUSED Star *s){}

    //global properties
    static inline size_t size() {return _size;}
    static vector<Property*> all;

    virtual bool are_table_loaded()  const  {return false;}


    //virtual inline void set_refpointers(_UNUSED std::vector<std::vector<std::vector<double>*>> &tables){return;}
    virtual inline void set_refpointers(_UNUSED Star *s){return;}


    double *val_ref[4];
    double *val_in[4];


protected:
    static size_t _size;

    double value, value0, Dvalue;
    double V, V0;
    vector<double> VBIN;

    double interpolating_values[4];
    double wM[4], wZ[2];

    //TODO GI changed the set functions to virtual, since as is the case of get some properties works internally in log, but here we it seesm that we assumed always a non log value, to be checked
    void set_fk(const double &a) { value = a;}
    void set_0_fk(const double &a) { value0 = a;}
    virtual void set(const double &a) { V = a;}
    virtual void set_0(const double &a) { V0 = a;}
    void set_VBIN(const size_t &id, const double &a) {VBIN[id] = a;}

    virtual inline void safety_check(){};


protected:
    //TODO, name is not needed because it can be linked directly to the pointer p->name()
    inline void Register(Property *_p, size_t *id, const std::string &_name){
        Property::all.push_back(_p);
        *id = Property::size() - 1;
        Property::PrintMap.insert(std::pair<std::string,size_t>(_name, *id));
        svlog.debug("Stellar property " + name() + " registered (Nproperties: " + utilities::n2s(all.size(),__FILE__,__LINE__) + ")");
    }
    SevnLogging svlog;
};

/** Timeing **/
/**
 * Base class for Properties about the times
 */
class Time_object : virtual public Property{

    inline std::string units() override {return "Myr";}

    /**
     * Disable copy
     * @param p
     */
    inline void copy_V_from(_UNUSED Property *p) override {

        svlog.error("Copy to the Time property " + this->name() + " is not allowed",__FILE__,__LINE__,true,sevnstd::sevnerr());
    }

    void evolve_remnant(_UNUSED Star *s) override {}
    void evolve_nakedco(_UNUSED Star *s) override {}
    virtual void set_empty(_UNUSED Star *s) override {}
    virtual void evolve_empty(_UNUSED Star *s) override {}


};

class Localtime :  public Time_object {

    public:
        Localtime(bool reg = true) {
            if (reg) {
                Register(this, &ID, name());
            }
        }

        static size_t ID;
        static Localtime _localtime;

        Localtime *Instance() override {
            return (new Localtime(false));
        }

        void synch() override {};

        inline void init(const double &a) override {V = V0 = value = value0 = a;}

        inline std::string name() const override { return "Localtime"; }

        //this property evolves outside the main loop
        inline void evolve(_UNUSED Star *s) override {};

        void set_remnant(_UNUSED Star *s) override {};

        void evolve_remnant(_UNUSED Star *s) override {}

        void special_evolve(_UNUSED Star *s) override;

        inline void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;

private:

    };

class Worldtime :  public Time_object {

public:
    Worldtime(bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }
        set(0.0); //worltime always initialized to ZERO... every time a new instance of this class is created!!
    }

    static size_t ID;
    static Worldtime _worldtime;

    Worldtime *Instance() override {
        return (new Worldtime(false));
    }

    inline std::string name() const override { return "Worldtime"; }
    void synch() override {};

    inline void set_remnant  (_UNUSED Star *s) override {};


    //this property evolves outside the main loop
    void evolve(_UNUSED Star *s) override {}
    void special_evolve(_UNUSED Star *s) override;
    void evolve_remnant(_UNUSED Star *s) override {}

private:

};

class Timestep :  public Time_object{
public:
    Timestep(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
        set(0.0); //Timestep initialized to 0.0 by default
    }

    Timestep * Instance() override {
        return (new Timestep(false));
    }

    inline std::string name() const override { return "Timestep"; }

    //TODO Dependent on the Compact remnant subtype
    double timestep_remnant(Star *s);
    void set_remnant(_UNUSED Star *s) override;
    /**
     * The only difference with set_remnant is that we do not change T0 because we are inside BSE
     * @param s Pointer to Star
     * @param b Pointer to Binary
     */
    void set_remnant_in_bse(_UNUSED Star *s, _UNUSED Binstar *b) override;
    void evolve_remnant(_UNUSED Star *s) override {set_remnant(s);};
    void evolve_nakedco(_UNUSED Star *s) override {set_remnant(s);};
    void set_empty(_UNUSED Star *s) override {set_remnant(s);};
    /**
     * The only difference with set_empty is that we do not change T0 because we are inside BSE
     * @param s Pointer to Star
     * @param b Pointer to Binary
     */
    void set_empty_in_bse(_UNUSED Star *s, _UNUSED Binstar *b) override;
    void evolve_empty(_UNUSED Star *s) override {set_remnant(s);};



    static size_t ID;
    static Timestep _timestep;


    /**
     * This function check if the current evolution time step has caused a too fast evolution of a property with values @p m0 (before the evolution), @p m (after the evolution)
     * and time derivative (@p derivative). If the property is related to some subphase (e.g. MHE, MCO) we should also set the initial time  of this phase with the parameter @p tstart.
     * If the property has not subphases (e.g. the Radius) @p tstart should be set to 0.
     * At the beginning of the function the new_dt is set to the max dt allowed (Delta_t phase/starparameter::min_points_per_phase).
     * Then if tstart is 0 the function check if the relative variation (@p m -@p m0)/ @p m is lower that 2 times  starparameter::maximum_variation (default 0.05).
     * If it is lower return false, else set a new dt to half of the current dt and set repeat to true.
     * if tstart is not  0 there are different options:
     *      if m0=0 and m>starparameter::maximum_variation:
     *             set a new_dt as 0.5*(Localtime - tstart) + (tstart - Localtime_0)  (time from Localtime_0 to start + half of Localtime-tstart)
     *             if the realtive variaiton  ((Localtime + new_dt)-tstart)/tstart <1e-2 do not repeat
     *             else repeat with the new dt (return immediatly).
     *      if m0>2 and relative difference (m-m0)/m0 < 2*starparameter::maximum_variation not needed to repeat
     *      else set repeat to true  and set new_dt=starparameter::maximum_variation*(@p m/@p derivative)
     *      if m02<2 and absolute difference (m-m0) <4*starparameter::maximum_variation not needed to repeat
     *      else set repeat to true and new_dt=2.0*starparameter::maximum_variation*(@p m/@p derivative)
     *
     * Finally there is a final check, if the new_dt is larger than 0.8*dt return False.
     *
     *
     * Example of usage
     * @code
     * check_repeat(s, s->getp_0(MHE::ID), s->getp(MHE::ID), s->getp(dMHEdt::ID), s->hestart());
     * check_repeat(s, s->getp_0(Radius::ID), s->getp(Radius::ID), s->getp(dRdt::ID), 0.0);
     * @endcode
     * @param s  Object of class star
     * @param m0 Value to check (before the evolution)
     * @param m  Value to check (after the evolution)
     * @param derivative  derivative of the value to check
     * @param tstart Starting time of a given subphase if we are checking the core mass of that subphase. E.g. s->hestart(). If we are checking a varaible without subphase tstart
     * should be set to 0.
     * @return Return true if the evolve step should be reset and repeated with a new smaller time step. The new time step is selected inside the
     * function with set(new_dt).
     * @note This function is called inside Timestep::evolve. Notice that the new_dt is only used for the repeated step. If this function return false
     * the proper new time step is estimated in  Timestep::evolve.
     */
    inline bool check_repeat(_UNUSED Star*s, const double &m0, const double &m, const double &derivative, const double tstart);
    //TODO To discuss: Since m0 and m should be related to the same properties can we just use the property ID as input?
    // Moreover the derivative could be a child class of the property. In this case all the infor are got by the class star.
    // However we are not really using the star here.
    //TODO This function does too much stuff. It both check and assigne a new timestep. I think these two operations should be separated or
    // at least we should change the name in something more informative.

    /**
     * This function check if the star is close to become a pureHE or nakedCO (the difference between Mass and MHE or MHE and CO or
     * Mass and CO are smaller than ev_naked_tshold)
     * @param s
     * @return
     */
    inline bool check_repeat_almost_naked(_UNUSED Star*s);

    //TODO mass should has already been evolved!!!
    /**
     *  Evolve the Timestep, the new value V will be the dt to use in the next iteration.
     *  The function passes through the following steps:
     *      - Check if the MHE, MCO or the Radius have evolved too much (comparing with  starparameter::maximum_variation). If this happened the function
     *         set repeat to True and set a new dt (see Timestep::checkrepeat)
     *      - If not needed to repeat determine next dt by imposing a maximum relative variation of stellar mass and radius.
     *      - Check if the dt should be changed because we are close to phase change (e.g. formation of the He core, CO core, remnant, ...)
     *              - Set reference_time to 0 if Localtime + new_dt<time_next_phase else reference_time=time_next_phase
     *              - if reference_time is 0, new_dt=min(new_dt, 0.999*(reference_time - Localtime) ) else new_dt=min(new_dt, 1.001 *(reference_time -Localtime).
     *      - Finally set Timestep::V to new_dt
     * @param s  Object of class star
     */
    void evolve(_UNUSED Star *s) override;

    void synch() override {};

    void resynch(const double &dt, bool set0) override {
        set(dt);
        if (set0) set_0(get());
    }

    void resynch(Star *s) override;

    inline void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;

protected:

    /**
     * Check if the proposed dt is larger than max_dt or lower than min_dt.
     * In that cases force dt to be one of the limits.
     * @param dt  reference to the timestep to analyse
     * @param s  pointer to the star
     */
    void check_dt_limits(double &dt, Star *s);

    /**
     * Set or reset the appropriate flag on Star after a check and repeat
     * @param s Pointer to the star s
     */
    void handle_star_check_repeat(Star *s);

private:
    bool checked_almost_naked=false; //Used  in check_repeat_almost_naked
};

class NextOutput : public Time_object {

public:
    NextOutput(bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }
        set(1.0e30); //nextoutput is always initialized to a very large number (to avoid multiple prints of the first snapshots)
    }

    static size_t ID;
    static NextOutput _nextoutput;

    NextOutput *Instance() override {
        return (new NextOutput(false));
    }

    inline std::string name() const override { return "NextOutput"; }

    void init(const double &a) {V = V0 = value = value0 = a;}

    //this property evolves outside the main loop
    inline void evolve(_UNUSED Star *s) override {}
    void special_evolve(Star *s) override;

    inline void set_remnant  (_UNUSED Star *s) override {};

    inline void evolve_remnant  (_UNUSED Star *s) override {};
    void evolve_nakedco(_UNUSED Star *s) override {};

    inline void restore() override {}




private:

};

/******************************************
 *  Table Properties
 *****************************************/
class TableProperty  : virtual public Property{
public:

    virtual size_t TabID() const = 0; //Each subclass of TableProperties has to implement this

    void set_refpointers(_UNUSED Star *s) override;

    bool are_table_loaded()  const override {return table_loaded;}

    void set_w(_UNUSED Star *s) override;


protected:

    virtual void set_wZ(_UNUSED Star *s);
    virtual void set_wM(_UNUSED Star *s);
    //Option for set_wM
    /**
     * Simple linear interpolation in Mass.
     *
     * @param s
     */
    void set_wM_linear(_UNUSED Star *s);
    /**
     * Linear interpolation with an additional weighting correction
     * @param s
     */
    void set_wM_rational(_UNUSED Star *s);
    /**
     * Log weights
     * @param s
     */
    void set_wM_log(_UNUSED Star *s);


    bool table_loaded=false;
};


/**** Classes about radii ***/

/**
 * Base class for properties concerning the Radii
 */
class R_object : virtual public TableProperty{


    inline std::string units() override {return "Rsun";}


    /**
     * Allow copy from/to R_obejct properties
     * @param p
     */
    inline void copy_V_from(Property *p) override {

        if (dynamic_cast<const R_object*>(p) != nullptr){
            set(p->get());
        }
        else
            svlog.error("Copying from property " + p->name() + " to property " + this->name() + " is not allowed",__FILE__,__LINE__,true,sevnstd::sevnerr());
    }


};

class Radius : public R_object{

public:
    Radius(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    size_t TabID() const override {return _Radius;}

    inline std::string name() const override { return "Radius"; }

    void set_remnant(_UNUSED Star *s) override;
    void evolve_remnant(_UNUSED Star *s) override;

    static size_t ID;
    static Radius _radius;

    Radius * Instance() override {
        return (new Radius(false));
    }

    void evolve(_UNUSED Star *s) override;
    void evolve_fake(Star *s) override;
    void update_from_binary(_UNUSED Star* s, const double &DV, _UNUSED  Binstar *b= nullptr) override;



    inline double get_fk(_UNUSED const Star* s=NULL) override {return pow(10.0,value);}
    inline double get_0_fk(_UNUSED const Star* s=NULL) override {return pow(10.0,value0);}
    inline double get(_UNUSED const Star* s=NULL) override {return pow(10.0, V);}
    inline double get_0(_UNUSED const Star* s=NULL) override {return pow(10.0, V0);}

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;

protected:

    void set(const double &a) override { V = std::log10(a);}
    void set_0(const double &a) override { V0 = std::log10(a);}
    void set_wM(_UNUSED Star *s) override;

    void evolve_fake_linear(Star *s);


    virtual inline void safety_check(){
        //TODO Use V<log10(TINY), could be faster.
        if (get()<utilities::TINY)
            svlog.critical("The property "+name()+" becomes extremely small. Something is seriously broken.",__FILE__,__LINE__);
    }

};


/**** Classes about Mass ***/
/**
 * Base class for Properties concerning masses
 */
class Mass_obejct : virtual public TableProperty{

    inline std::string units() override {return "Msun";}

    /**
     * Allow copy from/to Mass_obejct properties
     * @param p
     */
    inline void copy_V_from(Property *p) override {

        if (dynamic_cast<const Mass_obejct*>(p) != nullptr){
            set(p->get());
        }
        else
            svlog.error("Copying from property " + p->name() + " to property " + this->name() + " is not allowed",__FILE__,__LINE__,true,sevnstd::sevnerr());
    }

    double get_Vlast(const Star* s) const override;

public:
    void update_from_binary(Star* s, const double &DV, _UNUSED  Binstar *b) override;


};

class Mass : public Mass_obejct{

public:
    Mass(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Mass _mass;

    size_t TabID() const override {return _Mass;}

    Mass * Instance() override {
        return (new Mass(false));
    }

    inline std::string name() const override { return "Mass"; }

    void set_remnant(Star *s) override;

    void correct_interpolation_errors(_UNUSED Star *s) override;

    void correct_interpolation_errors_real(_UNUSED Star *s) override;

    inline void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;

protected:

    virtual inline void safety_check(){
        if (get()<utilities::TINY)
            svlog.critical("The property "+name()+" becomes extremely small. Something is seriously broken.",__FILE__,__LINE__);
    }

};

class MHE : public Mass_obejct {

public:
    MHE(bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }
    }

    size_t TabID() const override { return _MHE; }

    static size_t ID;
    static MHE _masshe;

    MHE *Instance() override {
        return (new MHE(false));
    }

    inline std::string name() const override { return "MHE"; }


    inline void evolve(_UNUSED Star *s);

    void correct_interpolation_errors(_UNUSED Star *s) override;
    void correct_interpolation_errors_real(_UNUSED Star *s) override;

    void changed_track(_UNUSED Star *s, _UNUSED Star *s_jtrack) override;

};

class MCO : public Mass_obejct{

public:
    MCO(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    size_t TabID() const override {return _MCO;}

    static size_t ID;
    static MCO _massco;

    MCO * Instance() override {
        return (new MCO(false));
    }

    inline std::string name() const override { return "MCO"; }




    inline void evolve(_UNUSED Star *s);


    void correct_interpolation_errors(_UNUSED Star *s) override;

    void correct_interpolation_errors_real(_UNUSED Star *s) override;

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


};

/*** Other Table properties ***/
class Phase : public TableProperty{

public:
    Phase(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Phase _phase;

    size_t TabID() const override {return _Phase;}

    Phase * Instance() override {
        return (new Phase(false));
    }

    inline std::string name() const override { return "Phase"; }
    void synch() override {};
    void set_remnant(_UNUSED Star *s) override {};
    void set_empty(_UNUSED Star *s) override {
        V0 = V = Lookup::Phases::Remnant;
    }
    void evolve_empty(_UNUSED Star *s) override {
        set_empty(s);
    }

    //this property evolves outside the main loop
    void special_evolve(_UNUSED Star *s) override;

    inline void evolve(_UNUSED Star *s) override{};
    void evolve_remnant(_UNUSED Star *s) override {}
    inline void changed_track(_UNUSED Star* s, _UNUSED Star* s_jtrack) override;

protected:

    void set_wM(Star *s) override;

};
class Luminosity : public TableProperty{

public:
    Luminosity(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Luminosity _luminosity;

    size_t TabID() const override {return _Lumi;}

    inline std::string name() const override { return "Luminosity"; }

    void set_remnant(_UNUSED Star *s) override;
    void evolve_remnant(_UNUSED Star *s) override;

    Luminosity * Instance() override {
        return (new Luminosity(false));
    }

    inline double get_fk(_UNUSED const  Star* s=NULL) override {return pow(10.0,value);}
    inline double get_0_fk(_UNUSED const  Star* s=NULL) override {return pow(10.0,value0);}
    inline double get(_UNUSED const  Star* s=NULL) override {return pow(10.0, V);}
    inline double get_0(_UNUSED const  Star* s=NULL) override {return pow(10.0, V0);}
    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;

protected:

    void set_wM(Star *s) override;


};


/*** Optional tables ***/
class OptionalTableProperty  : virtual public TableProperty{

public:

    void set_refpointers(_UNUSED Star *s) override;

    void evolve(Star *s) override {

        if (table_loaded){
            Property::evolve(s);
        } else{
            evolve_without_table(s);
        }

    }
    //NOTICE: VEry important! It tables are not loaded this property is de-facto as a derived property
    inline bool amiderived() override { return  !table_loaded;}
    /** If not following a table, these are derived properties,
     * so they have to follow the update_derived function same as derived properties
     */
    void update_derived(Star *s) override;
    virtual void synch() override{

        if (table_loaded)
            TableProperty::synch();
        else
            value = value0 = V0 = V;
    }

protected:

    virtual inline void evolve_without_table(_UNUSED Star *s){return;}

    inline void synch_v_value_without_table(){
        value0=V0;
        value=V;
    }
};

///Inertia
class Inertia : public OptionalTableProperty{

public:
    Inertia(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Inertia  _inertia;

    size_t TabID() const override {return _Inertia;}

    Inertia * Instance() override {
        return (new Inertia(false));
    }


    void set_refpointers(_UNUSED Star *s) override;

    void set_remnant(_UNUSED Star *s) override;
    void evolve_remnant(_UNUSED Star *s) override;
    void evolve_nakedco(_UNUSED Star *s) override;

    /** If we are following the tracks, it could happen that after a change of tracks the Inertia evolve to 0
     * in this case we just set a  minimum value of 0.1
     * @param s Pointer to the star
     */
     void correct_interpolation_errors_real(_UNUSED Star *s) override;


    inline std::string name() const override { return "Inertia"; }

    inline double get_fk(_UNUSED const  Star* s=NULL) override {return pow(10.0,value);}
    inline double get_0_fk(_UNUSED const  Star* s=NULL) override {return pow(10.0,value0);}
    inline double get(_UNUSED const  Star* s=NULL) override {return pow(10.0, V);}
    inline double get_0(_UNUSED const  Star* s=NULL) override {return pow(10.0, V0);}
    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;

    /**
     * Estimate the Inertia considering an homogeneous sphere considering the presence of a cavity
     * @param Mass total mass of the sphere
     * @param Outer_radius  total radius of the sphere
     * @param Inner_radius sphere of the inner cavity
     * @return the inertia in Msun Rsun^2
     */
    double estimate_Inertia_homogeneous_sphere(double Mass,double Outer_radius, double Inner_radius=0.);

protected:

    inline void evolve_without_table(_UNUSED Star *s) override;

    /**
     * Estimate the log10 of the Inertia assuming a uniform sphere
     * @param s Pointer to star
     * @return  log10 (2/5 MR^2)
     */
    double estimate_logInertia(_UNUSED Star *s);
    /**
     * Estimate the Inertia of a star assuming that the star is an homogeneous sphere
     * @param s Pointer to star
     * @return the inertia in Msun Rsun^2
     */
    double estimate_Inertia_homogeneous_sphere(_UNUSED Star *s);
    /**
     * Estimate the Inertia of a star assuming that the star is composed by two elements, an inner homogeneous sphere (the core)
     * and an outer homogeneous sphere with a cavity (the envelope).
     * @param s  Pointer to star
     * @return the inertia in Msun Rsun^2
     */
    double estimate_Inertia_homogeneous_sphere_wcore(_UNUSED Star *s);
    /**
     * Estimate the inertia following the formalism in Hurley+00 (Eq. 109)
     * @param s Pointer to star
     * @return the inertia in Msun Rsun^2
     */
    double estimate_Inertia_Hurley(_UNUSED Star *s);
    /**
     * Estimate the inertia following the formalism in DeMink+13 (Appendix A, A1,A2,A3)
     * @param s Pointer to the star
     * @return the inertia in Msun Rsun^2
     */
    double estimate_Inertia_DeMink(_UNUSED Star *s);

    void set_wM(Star *s) override;

    void check_and_set_rzams(_UNUSED Star *s);

private:

    double (Inertia::*inertia_func)(_UNUSED Star* s) = nullptr; //Pointer to the inertia func, set in set_refpointers
    double spin_core;

};

///Core radii Rhe, Rco
class CoreRadius : public OptionalTableProperty,  public R_object{

protected:

    /**
     * Estimate  Rcore using the equation 78 in Hurley+02, the equation has been slightly reshaped
     * so that the normalisation is the value of the radius at Mcore= 1 Msun
     * @param Mcore Mass of the core in Msun
     * @param scale_factor  value of the radius at Mcore= 1 Msun
     * @return estimate of the Radius in Rsun at given Mcore Mass.
     */
    virtual  inline double estimate_Rcore(double Mcore, double scale_factor){

        double norm = scale_factor *  1.1685; //1.1685 is the inverse value of the function at M=1 not considering the scale factor
        double M3 = Mcore*Mcore*Mcore;

        return norm*std::pow(Mcore,4.6) / (  Mcore*M3 + 0.162*M3 + 0.0065   );

    }

};

class RHE :  public CoreRadius {

public:
    RHE(bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static RHE _rhe;

    size_t TabID() const override {return _RHE;}

    RHE *Instance() override {
        return (new RHE(false));
    }

    void set_refpointers(_UNUSED Star *s) override;

    inline std::string name() const override { return "RHE"; }


    inline void evolve(_UNUSED Star *s);


    void correct_interpolation_errors(_UNUSED Star *s) override;
    void correct_interpolation_errors_real(_UNUSED Star *s) override;

    inline void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


protected:

    void evolve_without_table(_UNUSED Star *s) override;

    double estimate_Rcore(_UNUSED Star *s);

};

class RCO : public CoreRadius {

public:
    RCO(bool reg = true) {
        if (reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static RCO _rco;

    size_t TabID() const override {return _RCO;}

    RCO *Instance() override {
        return (new RCO(false));
    }


    void set_refpointers(_UNUSED Star *s) override;


    inline std::string name() const override { return "RCO"; }


    inline void evolve(_UNUSED Star *s);


    void correct_interpolation_errors(_UNUSED Star *s) override;
    void correct_interpolation_errors_real(_UNUSED Star *s) override;
    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


protected:

    void evolve_without_table(_UNUSED Star *s) override;
    double estimate_Rcore(_UNUSED Star *s);

};

///Chemical Surface composition
class SurfaceAbundanceTable : public OptionalTableProperty{
public:

    void set_refpointers(_UNUSED Star *s) override;

};

/** This property stores the H fraction on the star surface
 */
class Hsup : public SurfaceAbundanceTable{

public:
    Hsup(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Hsup _hsup;

    size_t TabID() const override {return _Hsup;}

    Hsup * Instance() override {
        return (new Hsup(false));
    }

    inline std::string name() const override { return "Hsup"; }

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


};
/** This property stores the He fraction on the star surface
 */
class HEsup : public SurfaceAbundanceTable{

public:
    HEsup(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static HEsup _hesup;

    size_t TabID() const override {return _HEsup;}

    HEsup * Instance() override {
        return (new HEsup(false));
    }

    inline std::string name() const override { return "HEsup"; }

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


};
/** This property stores the C fraction on the star surface
 */
class Csup : public SurfaceAbundanceTable{

public:
    Csup(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Csup _csup;

    size_t TabID() const override {return _Csup;}

    Csup * Instance() override {
        return (new Csup(false));
    }

    inline std::string name() const override { return "Csup"; }

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


};
/** This property stores the N fraction on the star surface
 */
class Nsup : public SurfaceAbundanceTable{

public:
    Nsup(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Nsup _nsup;

    size_t TabID() const override {return _Nsup;}

    Nsup * Instance() override {
        return (new Nsup(false));
    }

    inline std::string name() const override { return "Nsup"; }

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


};
/** This property stores the O fraction on the star surface
 */
class Osup : public SurfaceAbundanceTable{

public:
    Osup(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Osup _osup;

    size_t TabID() const override {return _Osup;}

    Osup * Instance() override {
        return (new Osup(false));
    }

    inline std::string name() const override { return "Osup"; }

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


};

/// Convective envelope properties ****/
class ConvectiveTable : public OptionalTableProperty{

public:
    void set_refpointers(_UNUSED Star *s) override;

};

/** This property stores the Mass fraction of the convective envelope.
 * In order to estimate the final mass it has to be multiplied by the total Mass.
 */
class Qconv : public ConvectiveTable{

public:
    Qconv(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Qconv _qconv;

    size_t TabID() const override {return _Qconv;}

    Qconv * Instance()  override {
        return (new Qconv(false));
    }

    inline std::string name() const override { return "Qconv"; }

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


protected:

    /**
     * If Qconv table is not loaded, consider all the envelope as convective
     * except for nakedHelium and main sequence star
     * @param s
     */
    void evolve_without_table(_UNUSED Star *s) override;

    double estimate_Qconv(_UNUSED Star *s);

};
/** This property stores the T(?) scale  of the convective envelope
 */
class Tconv : public ConvectiveTable{

public:
    Tconv(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Tconv _tconv;

    size_t TabID() const override {return _Tconv;}

    Tconv * Instance() override {
        return (new Tconv(false));
    }

    inline std::string name() const override { return "Tconv"; }

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


};
/** This property stores the Depth of the convective envelope normalised to the
 * stellar radius.
 * @Note it has to be multiplied with the stellar radius.
 */
class Depthconv : public ConvectiveTable{

public:
    Depthconv(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }

    }

    static size_t ID;
    static Depthconv _depthconv;

    size_t TabID() const override {return _Depthconv;}

    Depthconv * Instance() override {
        return (new Depthconv(false));
    }

    inline std::string name() const override { return "Depthconv"; } //to multiply with the stellar radius!

    void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;
protected:

    /**
     * If DepthConv table is not loaded, consider all the envelope as convective
     * except for nakedHelium and main sequence star
     * @param s
     */
     void evolve_without_table(_UNUSED Star *s) override;

     double estimate_Dconv(_UNUSED Star *s);

};
/******************************************/

/*** Other non tab properties ***/

class Bmag : public Property{

public:
    Bmag(bool reg=true){
        if (reg){
            Register(this,&ID,name());
        }
        V=V0=value=value0=0;
    }
    static size_t ID;
    static Bmag _bmag;
    inline std::string name() const override {return "Bmag";}
    virtual inline std::string units(){return "Gauss";}

    Bmag *Instance() override{
        return new Bmag(false);
    }

    void evolve(_UNUSED Star *s) override;

    void set_remnant(_UNUSED Star *s) override;
    void evolve_remnant(_UNUSED Star *s) override;

};

class OmegaRem : public Property{
public:
    OmegaRem(bool reg=true){
        if (reg){
            Register(this,&ID,name());
        }
        V=V0=value=value0=std::nan("");
    }
    static size_t ID;
    static OmegaRem _omegarem;
    inline std::string name() const override {return "OmegaRem";}
    virtual inline std::string units(){return "s^-1";}

    OmegaRem *Instance() override{
        return new OmegaRem(false);
    }

    //Do nothing when evolve
    void evolve(_UNUSED Star *s) override {};

    //Here do stuff
    void set_remnant(_UNUSED Star *s) override;
    void evolve_remnant(_UNUSED Star *s) override;

};

class RemnantType : public Property{

public:
    RemnantType(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
        set(Lookup::Remnants::NotARemnant); //initialize to not-a-remnant
    }

    static size_t ID;
    static RemnantType _remnanttype;


    RemnantType * Instance() override {
        return (new RemnantType(false));
    }

    inline std::string name() const override { return "RemnantType"; }
    void synch() override {};

    //this property evolves outside the main loop
    void set_remnant(_UNUSED Star *s) override;

    void set_empty(_UNUSED Star *s) override;
    void evolve_empty(_UNUSED Star *s) override {set_empty(s);}

    inline void evolve(_UNUSED Star *s) override;

};

/*!  This property stores the DM (positive or negative) accumulated in the stars from binary processeses.
 *  Note it does not store the mass lost from Wind because it is a single stellar process and it is already included in the track interpolation.
 *  However it stores the mass accreted from the Wind or RLO of the companion star and/or the mass lost in RLO.
 *  It is reset to 0   everytime the star jump to  a new track
 *  In the SSE it remains always set to 0.
 * */
class dMcumul_binary : public Property{

public:


    dMcumul_binary(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
            value = value0 = Dvalue = V = V0 = 0.0; //The Property constructor already set this so to 0. But just to be sure for future changes.
        }

    }

    dMcumul_binary * Instance() override {
        return (new dMcumul_binary(false));
    }

    inline std::string name() const override { return "dMcumul_binary"; }


    //Static definitions
    static size_t ID;
    static dMcumul_binary _dMcumul_binary;


    inline void evolve(_UNUSED Star *s) override {V0=V;}
    inline void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override{V=V0=value=value0=0;}



};

/*!  This property stores the DM (positive or negative) accumulated during an episode of RLO.
 *  It is reset to 0 everytime the RLO episode stops
 *  In the SSE it remains always set to 0.
 * */
class dMcumul_RLO : public Property{
public:

    dMcumul_RLO(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
            value = value0 = Dvalue = V = V0 = 0.0; //The Property constructor already set this so to 0. But just to be sure for future changes.
        }

    }

    dMcumul_RLO * Instance() override {
        return (new dMcumul_RLO(false));
    }

    inline std::string name() const override { return "dMcumul_RLO"; }



    //Static definitions
    static size_t ID;
    static dMcumul_RLO _dMcumul_RLO;

    inline void evolve(_UNUSED Star *s) override {V0=V;}
    //changed_track do nothing same as the base changed_track

    void update_from_binary(_UNUSED Star* s, const double &DV, _UNUSED  Binstar *b) override;


};


/*!  This property stores the DM  accumulated due to the mass wind accretion.
 *  In the SSE it remains always set to 0.
 * */
class dMcumulacc_wind : public Property{
public:

    dMcumulacc_wind(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
            value = value0 = Dvalue = V = V0 = 0.0; //The Property constructor already set this so to 0. But just to be sure for future changes.
        }

    }

    dMcumulacc_wind * Instance() override {
        return (new dMcumulacc_wind(false));
    }

    inline std::string name() const override { return "dMcumulacc_wind"; }



    //Static definitions
    static size_t ID;
    static dMcumulacc_wind _dMcumulacc_wind;

    inline void evolve(_UNUSED Star *s) override {V0=V;}
    //changed_track do nothing same as the base changed_track
};

/**
 * Stellar angular momentum
 * This is the fundamental property for stellar rotation.
 * This is the one that is modified by the sse and bse.
 * All the others are just derived: OmegaSpin and Spin
 */
class AngMomSpin : public Property{
public:
    AngMomSpin(bool reg=true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    inline std::string name() const override { return "AngMomSpin"; }
    inline std::string units() override {return "Msun Rsun^2 yr^-1";}

    //Angmom is set at the star initialisation (through spin)
    inline void init(const double &a) override {V = V0 = value = value0 = a;}
    static size_t ID;
    static AngMomSpin _angmomspin;

    AngMomSpin * Instance() override {
        return (new AngMomSpin(false));
    }

    void evolve(_UNUSED Star *s) override;
    /**
     *  We consider that angular momentum is conserved, however we have to check
     * if the new angmom is below the critical L, otherwise set is to Lcrit
     * @param s  Pointer to the star that is changing track
     * @param s_jtrack Pointer to the auxiliary stars used to change track (not used)
     */
    void changed_track(_UNUSED Star* s, _UNUSED Star* s_jtrack) override;


protected:
    double evolve_angmom(_UNUSED Star *s);
};





/******************************************
 *  Derived Properties
 *****************************************/

/**
 * Base class for properties that are estimated as the derivative of basic properties
 */
class Derivative_Property : virtual public Property{

    /**
     * Disable copy
     * @param p
     */
    inline void copy_V_from(_UNUSED Property *p) override {

            svlog.error("Copy to the Derived property " + this->name() + " is not allowed",__FILE__,__LINE__,true,sevnstd::sevnerr());
    }
};

class dMdt : public Derivative_Property{

public:

    dMdt(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    inline std::string name() const override { return "dMdt"; }
    inline std::string units() override {return "Msun/Myr";}


    static size_t ID;
    static dMdt _dmdt;

    dMdt * Instance() override {
        return (new dMdt(false));
    }

    //void synch() override {};

    void evolve(_UNUSED Star *s) override;

};

class dMHEdt : public Derivative_Property{

public:

    dMHEdt(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static dMHEdt _dmhedt;

    inline std::string name() const override { return "dMHEdt"; }
    inline std::string units() override {return "Msun/Myr";}


    dMHEdt * Instance() override {
        return (new dMHEdt(false));
    }

    //void synch() override {};

    void evolve(_UNUSED Star *s) override;

};

class dMCOdt : public Derivative_Property{

public:

    dMCOdt(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static dMCOdt _dmcodt;

    inline std::string name() const override { return "dMCOdt"; }
    inline std::string units() override {return "Msun/Myr";}


    dMCOdt * Instance() override {
        return (new dMCOdt(false));
    }

    //void synch() override {};

    void evolve(_UNUSED Star *s) override;

};

class dRdt : public Derivative_Property{

public:

    dRdt(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    dRdt * Instance() override {
        return (new dRdt(false));
    }

    inline std::string name() const override { return "dRdt"; }
    inline std::string units() override {return "Rsun/Myr";}


    static size_t ID;
    static dRdt _drdt;

    //void synch() override {};

    void evolve(_UNUSED Star *s) override;

};

/*** Derived properties ***/

/**
 * Base class for properties that are derived from other properties
 */
class Derived_Property : virtual public Property{

    /**
     * Disable copy
     * @param p
     */
    inline void copy_V_from(_UNUSED Property *p) override {

        svlog.error("Copy to the Derived property " + this->name() + " is not allowed",__FILE__,__LINE__,true,sevnstd::sevnerr());
    }

    virtual inline bool amiderived() { return  true;}

    //A bit different synch for derived, since value and value0 are always 0.
    void synch() override {
        value=value0=V0=V;
    };



    //For the derived property when set_remnant or evolve_remnant is called just run evolve to set the proper value
    void set_remnant(_UNUSED Star *s){ evolve(s);}
    void evolve_remnant(_UNUSED Star *s){ evolve(s);}
    void evolve_nakedco(_UNUSED Star *s){ evolve(s);}
    /** This function can be called to update the value of V in derived proprties (without setting V0)
     * It is empty for all the others. This is needed since after a binary evolution, the derived properties have to be updated
     * without setting V0 (same as update_from_binary).
     */
     void update_derived(Star *s) override;
};
class Temperature : public Derived_Property{

public:

    Temperature(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    Temperature * Instance() override {
        return (new Temperature(false));
    }

    inline std::string name() const override { return "Temperature"; }
    inline std::string units() override {return "K";}

    static size_t ID;
    static Temperature _temperature;

    //TODO mass should has already been evolved!!!
    void evolve(_UNUSED Star *s) override;

    inline void changed_track (_UNUSED Star* s, _UNUSED Star* s_jtrack) override;

};

class Rs : public Derived_Property{

public:

    static size_t ID;
    static Rs _rs;

    Rs(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    Rs * Instance() override {
        return (new Rs(false));
    }
    inline std::string name() const override { return "Rs"; }
    inline std::string units() override {return "Rsun";}

    //TODO mass should has already been evolved!!!
    void evolve(_UNUSED Star *s) override;

};

/**
 * OmegaSpin is derived from AngmomSpin and the Inertia,
 * it represents the angular velocity of the star
 */
class OmegaSpin : public Derived_Property{

public:

    OmegaSpin(bool reg=true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    inline std::string name() const override { return "OmegaSpin"; }
    inline std::string units() override {return "yr^-1";}

    //Spin is set at the star initialisation
    inline void init(const double &a) override {V = V0 = value = value0 = a;}
    static size_t ID;
    static OmegaSpin _omegaspin;

    OmegaSpin * Instance() override {
        return (new OmegaSpin(false));
    }

    void evolve(_UNUSED Star *s) override;
    void changed_track(Star* s, _UNUSED Star* s_jtrack) override { evolve(s);}

};

/**
 * Derived properties representing the value OmegaSpin/OmegaCrit
 */
class Spin : public Derived_Property{

public:
    static size_t ID;
    static Spin _spin;

    Spin(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }

    Spin * Instance() override {
        return (new Spin(false));
    }
    inline std::string name() const override { return "Spin"; }
    inline std::string units() override {return "";}

    //Classica synch here, since we initialise it with init
    void synch() override {
        V=V0=value=value0;
    };

    //Spin is set at the star initialisation
    inline void init(const double &a) override {V = V0 = value = value0 = a;}

    void evolve(_UNUSED Star *s) override;
    void changed_track(Star* s, _UNUSED Star* s_jtrack) override { evolve(s);}


    /**
     * Estimate Spin (i.e. the ratio between OmegaSpin and Omegacrit) from the OmegaSpin (i.e. the angular velocity)
     * and the stellar properties
     * @param OmegaSpin Angular velocity (yr^-1)
     * @param Mass Mass of the star in Msun
     * @param Radius Radius of the star in Rsun
     * @return Return of the ration between OmegaSpin and Omegacrit
     */
    static double Spin_from_OmegaSpin(double OmegaSpin, double Mass, double Radius);

};

/**
 * Black hole spin,
 * it is classified as a derived property since it could be updated during the BSE
 */
class Xspin : public Derived_Property{

public:

    Xspin(bool reg = true){
        if(reg){
            Register(this, &ID, name());
        }
        V=V0=value=value0=std::nan("");
    }

    static size_t ID;
    static Xspin _xspin;

    Xspin * Instance() override {
        return (new Xspin(false));
    }

    inline std::string name() const override {return "Xspin";}
    inline void evolve (_UNUSED Star *s) override {};
    void set_remnant(_UNUSED Star *s) override;
    void evolve_remnant(_UNUSED Star *s) override;

};


/******************************************
 *  JIT Properties
 *****************************************/

/**
 * Derived class for properties that are derived from other properties
 * JIT means Just in Time, these are properties that do not need to be
 * evolved at each time step, but just to be computed when needed.
 */
class JIT_Property : public Derived_Property{

public:

    inline std::string name() const override { return "Generic JIT Property"; }

    //Special evolve
    void evolve(_UNUSED Star *s) override {
        evolve_number++;
    }
    void update_from_binary(_UNUSED Star* s, _UNUSED const double &DV, _UNUSED  Binstar *b) override{
        evolve_number++;
    }

    //Do nothing here, empty is already handled in get
    void set_empty(_UNUSED Star *s) override {}
    void evolve_empty(_UNUSED Star *s) override {}

    void set_remnant(_UNUSED Star *s) override {evolve(s);}
    void evolve_remnant(_UNUSED Star *s) override {evolve(s);}
    void evolve_nakedco(_UNUSED Star *s) override {evolve(s);}
    /**Not needed to update JIT properties after binary evolution, they are estimated on request**/
    void update_derived(_UNUSED Star *s) override{};

    virtual void restore() override  {}
    virtual void synch() override {}


    //Disable calling some get methods (only get is important for JIT properties)
    inline double get_fk(_UNUSED const Star* s=NULL) override {
        svlog.critical("get_fk for Property "+name()+" should not be called",__FILE__,__LINE__,sevnstd::notimplemented_error());
        return -1.0;
    }
    inline double get_0_fk(_UNUSED const Star* s=NULL) override {
        svlog.critical("get_0_fk for Property "+name()+" should not be called",__FILE__,__LINE__,sevnstd::notimplemented_error());
        return -1.0;
    }
    inline double get_0(_UNUSED const Star* s=NULL) override {
        svlog.critical("get_0 for Property "+name()+" should not be called",__FILE__,__LINE__,sevnstd::notimplemented_error());
        return -1.0;
    }


protected:

    bool new_estimate_needed(){

        return evolve_number>last_evolve_number;
    }

    void set(const double &a) override  {
        V = a;
        last_evolve_number=evolve_number;
    }
    //virtual void set_0(const double &a) { V0 = a;}



private:


    //The following  numbers check if an evolution of a binary evolution has been called
    //from the last time this property has been estimated.
    //evolve_number increases by 1 each time evolve or update_from_binary has been called.
    //last_evolve_number stores the value of evolve_number when getp is called.
    //When calling getp it checks if the two numbers are equal, if they are the property is just last stored one.
    unsigned int last_evolve_number=0, evolve_number=0;

};

class Lambda : public JIT_Property{

public:

    static size_t ID;
    static Lambda _lambda;

    Lambda(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    Lambda * Instance() override {
        return (new Lambda(false));
    }
    inline std::string name() const override { return "Lambda"; }
    inline std::string units() override {return "Number";}

    void set_remnant(_UNUSED Star *s) override {V=std::nan("");}
    void evolve_remnant(_UNUSED Star *s) override {V=std::nan("");}
    void evolve_nakedco(_UNUSED Star *s) override {V=std::nan("");}

    /**
     * Estimate the lambda for a given star.
     * The binding energy of an envelope is Ebind = Integrate^Mtot_Mcore dM - G M(r)/r + alpha_th*U,
     * where U is the internal energy and alpha_th a fraction of internal energy used to unbound the envelope.
     * For simplicity we can write Ebind_s = -G Mtot Menv / lambda R, where lambda is the factor for which Ebind_s=Ebind.
     * Usually, the lambda is taken between two extrema: lambda_g obtained considering only the gravitational energy (alpha_th=0)
     * and lambda_b obtained consider a fully efficiency of the internal energy in unbounding the envelope (alpha_th=1).
     * The implementation  depends on the Common envelope models and common envelope parameters
     * @param star Star for which we want to estimate the binding energy
     * @return
     */
    double get(_UNUSED const Star* s=NULL) override;


protected:
    /**
     * Estimate the lambda for a given star.
     * The binding energy of an envelope is Ebind = Integrate^Mtot_Mcore dM - G M(r)/r + alpha_th*U,
     * where U is the internal energy and alpha_th a fraction of internal energy used to unbound the envelope.
     * For simplicity we can write Ebind_s = -G Mtot Menv / lambda R, where lambda is the factor for which Ebind_s=Ebind.
     * Usually, the lambda is taken between two extrema: lambda_g obtained considering only the gravitational energy (alpha_th=0)
     * and lambda_b obtained consider a fully efficiency of the internal energy in unbounding the envelope (alpha_th=1).
     * The implementation  depends on the Common envelope models and common envelope parameters
     * @param star Star for which we want to estimate the binding energy
     * @return
     */
    double estimate_lambda (const Star *star);

private:
    /**
     * Estimate lambda following exactly the implementation of lambda in Mobse
     * The mass fraction of the convective envelope over the whole envelope will be used
     * @param star Star for which we want to estimate the binding energy
     * @param whole_cenv If true conside that all the envelope is convective
     * @return
     */
    double estimate_lambda_BSE  (const Star *star, bool whole_cenv=false, bool m0_as_hurley=false);

    /**
     * Estimate lambda from the Parsec tracks.
     * This is an experimental function based on the first results of Mario and her student.
     * At the moment we use a very super preliminary results: it seems that expect for the very last
     * stellar phase, the lambda estimated from Parsec tracks is a scaled (0.1) version of Claeys et al. 2014 (Appendix A)
     * @param star Star for which we want to estimate the binding energy
     * @return
     */
    double estimate_lambda_Parsec  (const Star *star);

    /**
     * Estimate lambda following Claeys+14 (Appendix A, https://www.aanda.org/articles/aa/pdf/2014/03/aa22714-13.pdf)
     * The mass of the convective envelope in Msun will be used
     * @param star Star for which we want to estimate the binding energy
     * @param whole_cenv If true conside that all the envelope is convective
     * @return
     */
    double estimate_lambda_Claeys14(const Star *star, bool whole_cenv=false);

    /**
     * Estimate lambda following Izzard+04 (Appendix E, http://personal.ph.surrey.ac.uk/~ri0005/doc/thesis/thesis.pdf),
     * it is similar to Claeys+14, but  lambda1 is  multiplied by 2 before to add the correction for the ionization terms
     * The mass of the convective envelope in Msun will be used
     * @param star Star for which we want to estimate the binding energy
     * @param whole_cenv If true conside that all the envelope is convective
     * @return
     */
    double estimate_lambda_Izzard04(const Star *star, bool whole_cenv=false);

    /**
     * Estimate lambda using tables from Klencki+21
     * The Lambda fit comes from Appendix A in Klencki+21 (https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..54K/abstract).
     * Lambda is thus estimated as:
     *      log10(lambda) = a*x^3 + b*x^2 + c*x + d
     *      where x=log10(R)
     * The value of a,b,c,d depends on the Mass (Mzams), Z (metallicity) and Radius.
     * @param star Pointer to the star
     * @param interpolate if true, estimate lambda interpolating between lambdas for tracks in Klencki+21 tables
     * if false, use a quantised formalism, i.e. for a Mzams and a Z not present in the table, the
     * closed track with Mtrack>Mzams and Ztrack>Z. If Mzams<10, Mtracks=10 is used, for Z/Zsun<0.01, Z/Zsun=0.01 is used.
     * @return Lambda
     */
    double estimate_lambda_Klencki21(const Star *star, bool interpolate=false);

    /**
     * Estimate Lambda using fitting equations from the work of Xi&Liu10 (https://iopscience.iop.org/article/10.1088/0004-637X/716/1/114) +
     * its errata corrige (https://iopscience.iop.org/article/10.1088/0004-637X/722/2/1985).
     * The implementation is taken directly from the COMPAS code https://github.com/TeamCOMPAS/COMPAS that is based
     * on the STARTRACK implementation (see Dominik+12, https://ui.adsabs.harvard.edu/abs/2012ApJ...759...52D/abstract).
     * Lambda is thus estimated as:
     *
     *      y = a + b1*x + b2*x^2 + b3+*x^3 + b4*x^4 + b5*x^5
     *
     *      where y is in general lambda (gravitational or considering also the thermal energy), sometime it is 1/lambda
     *      or sometime it is log lambda.
     *      x in general is the Radius, but sometime it is the mass of the envelope.
     * @param star  Pointer to the star
     * @param interpolate if true, estimate lambda interpolating between lambdas for tracks in Xu&Lu10
     * if false, use a quantised formalism, i.e. for a Mzams not in  the table use, the
     * closed track with Mtrack>Mzams and Ztrack>Z.
     * @return Lambda
     */
    double estimate_lambda_Nanjing(const Star *star, bool interpolate=false);

    /** Cached values for BSE-like MOBSE-like lambdas **/
    double Mzams_cachedM=0, Zmet_cachedM=0, M0_cached=0; //Cached values for the function get_M0_BSE
    double Mzams_cachedR=0, Zmet_cachedR=0, Rzams_cached=0; //Cached values for the function get_Rzams

    /** Lambda Klencki **/
    bool first_call=true; /*!Used in some lambda estimate (e.g. Klencki) to check if it is needed to initialise stuff*/
    std::unique_ptr<Lambda_Base> lambda_base_ptr; /*!Pointer to Lambda Klencki class*/

    /**
     * Estimate the bse-like parameter M0 (it is the mass at the end of the MS)
     * We estimate it creating an auxilixiary star with the same Mzams and Z and estimate the mass at the tms.
     * Since this value depends only on Mzams and Z, the results is cached in the member variable M0_cached.
     * @param s Pointer to the star
     * @return M0 in Msun
     */
    double get_M0_BSE(const Star *s);
    /**
     * Estimate the Radius at ZAMS
     * We estimate it creating an auxilixiary star with the same Mzams and Z.
     * Since this value depends only on Mzams and Z, the result is cached in the member variable Rzams_cached.
     * @param s Pointer to the star
     * @return Radius at MZAMS in Rsun
     */
    double get_Rzams(const Star *s);



};

///Binding energy of the envelope
class Ebind : public JIT_Property{

public:
    static size_t ID;
    static Ebind _ebind;

    Ebind(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    Ebind * Instance() override {
        return (new Ebind(false));
    }
    inline std::string name() const override { return "Ebind"; }
    inline std::string units() override {return "Msun^2/R/G";}

    void set_remnant(_UNUSED Star *s) override {V=0.;}
    void evolve_remnant(_UNUSED Star *s) override {V=0.;}
    void evolve_nakedco(_UNUSED Star *s) override {V=0.;}

    /**
     * Ebind = -(1/lambda)*M*Menv/R
     * @param s
     * @return
     */
    double get(_UNUSED const Star* s=NULL) override;

};

///BSE Type
class PhaseBSE : public JIT_Property{

public:
    static size_t ID;
    static PhaseBSE _phasebse;

    PhaseBSE(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    PhaseBSE * Instance() override {
        return (new PhaseBSE(false));
    }
    inline std::string name() const override { return "PhaseBSE"; }
    inline std::string units() override {return "";}

    double get(_UNUSED const Star* s=NULL) override;
};

//zams
class Zams : public JIT_Property{

public:
    static size_t ID;
    static Zams _zams;

    Zams(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    Zams * Instance() override {
        return (new Zams(false));
    }
    inline std::string name() const override { return "Zams"; }
    inline std::string units() override {return "Msun";}

    void set_remnant(_UNUSED Star *s) override {}
    void evolve_remnant(_UNUSED Star *s) override {}
    void evolve_nakedco(_UNUSED Star *s) override {}

    double get(_UNUSED const Star* s=NULL) override;
};

//Z
class Zmet : public JIT_Property{

public:
    static size_t ID;
    static Zmet _zmet;

    Zmet(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    Zmet * Instance() override {
        return (new Zmet(false));
    }
    inline std::string name() const override { return "Zmet"; }
    inline std::string units() override {return "";}

    void set_remnant(_UNUSED Star *s) override {}
    void evolve_remnant(_UNUSED Star *s) override {}
    void evolve_nakedco(_UNUSED Star *s) override {}

    double get(_UNUSED const Star* s=NULL) override;


};

//Event
class Event : public JIT_Property{
public:
    static size_t ID;
    static Event _event;

    Event(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
        set(Lookup::EventsList::NoEvent);
    }

    Event * Instance() override {
        return (new Event(false));
    }
    inline std::string name() const override { return "Event"; }
    double get(_UNUSED const Star* s=NULL) override;

private:
    bool is_qhe_set = false;
};

/**
 * JIT property to estimate the current accretion/loss rate on a object.
 * If it is positive it is the accretion rate, if it is negative it is the mass loss rate.
 * The units are Msun/Myr
 */
class dMRLOdt : public JIT_Property{

public:
    static size_t ID;
    static dMRLOdt _dmorlodt;
    dMRLOdt(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    dMRLOdt * Instance() override {
        return (new dMRLOdt(false));
    }

    inline std::string name() const override { return "dMRLOdt"; }
    inline std::string units() override {return "Msun/Myr";}

    double get(_UNUSED const Star* s=NULL) override;

};

/**Property to store the rate of mass accreted through winds
 *
 */
class dMaccwinddt : public JIT_Property{

public:
    static size_t ID;
    static dMaccwinddt _dmaccwinddt;
    dMaccwinddt(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    dMaccwinddt * Instance() override {
        return (new dMaccwinddt(false));
    }

    inline std::string name() const override { return "dMaccwinddt"; }
    inline std::string units() override {return "Msun/Myr";}

    double get(_UNUSED const Star* s=NULL) override;
};




//alpha angle
//TODO I put this property just to retrieve quickly this property for the Cecilia's work, but this information is already present
//in the logfile and the property should be removed.
/**
 * This property is used  to retrieve the angle between the magnetic axis and the the rotaton axis in NS. It is nan in all the other case.
 * It does not evolve but just get its value when set_remnant is called and the remnant is a NS.
 */
class NSsalpha : public JIT_Property{

    static size_t ID;
    static NSsalpha _nssalpha;
    NSsalpha(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    NSsalpha * Instance() override {
        return (new NSsalpha(false));
    }

    inline std::string name() const override { return "NSsalpha"; }
    inline std::string units() override {return "";}

    void evolve(_UNUSED Star *s) override {}
    void set_remnant(_UNUSED Star *s) override;

    double get(_UNUSED const Star* s=NULL) override{return sinalpha;}

protected:

    double sinalpha=std::nan("");
};


/**
 * This property is used  to retrieve the Plife on the star in the current phase
 * If the star is a remnant it is se to 0
 */
class Plife : public JIT_Property{

public:
    static size_t ID;
    static Plife _plife;
    explicit Plife(bool reg = true){
        if(reg) {
            Register(this, &ID, name());
        }
    }
    Plife * Instance() override {
        return (new Plife(false));
    }

    inline std::string name() const override { return "Plife"; }
    inline std::string units() override {return "";}

    void set_remnant(_UNUSED Star *s) override {V=0.;}
    double get(_UNUSED const Star* s=NULL) override;



};



#endif //SEVN_PROPERTY_H

/******************************************/