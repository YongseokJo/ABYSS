//
// Created by iorio on 7/22/22.
//

#ifndef SEVN_CIRCULARISATION_H
#define SEVN_CIRCULARISATION_H

#include <Processes.h>

//Forward declaration
class Star;
class Binstar;

/**
 *  Circularisation process
 *
 *  This process handles the circularisation at the onset of RLO
 */
class Circularisation : public  Process {
public:
    Circularisation(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            //First register to include this Property in static array
            Register(this, &ID, name());
        }
    }
    virtual ~Circularisation(){};

    static size_t ID;
    static Circularisation _circularisation;
    inline std::string name() override { return "Circularisation"; }

    Circularisation *Instance(_UNUSED IO *_io) override;
    virtual Circularisation *instance() {return nullptr;}

protected:
    //this a prototype for the future refactoring of the whole Process classes
    static std::map<std::string, Circularisation *> & GetStaticMap(){
        static std::map<std::string, Circularisation *> _locmap;
        return _locmap;
    }
    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }
    static void Register_specific(Circularisation *ptr) {
        GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }


};




/**
 * disabled option
 * no circularisation at the onset of RLO
 */
class CircularisationDisabled : public Circularisation{

public:

    CircularisationDisabled(_UNUSED IO *_io= nullptr, bool reg = true) : Circularisation(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "disabled"; }
    static CircularisationDisabled _circularisationdisabled;
    CircularisationDisabled* instance() override {
        return (new CircularisationDisabled(nullptr,false));
    }
    int evolve(_UNUSED Binstar *binstar) override;
};

/**
 * Template class for standard circularisation, i.e. change impulsively the value of a and e
 */
class StandardCircularisation : public Circularisation{

public:
    StandardCircularisation(_UNUSED IO *_io= nullptr, bool _UNUSED = true) : Circularisation(nullptr, false) {}

    /**
     * The standard evolve just check if the circularisation conditions have been reached
     * and set a couple of flag.
     * @param binstar Pointer to the binary
     * @return EXIT SUCCESS
     */
     int evolve(_UNUSED Binstar *binstar) override;


    /**
     * Special evolve, just call circularise and disable a and e check in the adaptive timestep
     * @param binstar Pointer to the binary
     * @return 0 (it measn that the binary is not broken)
     */
    int special_evolve(_UNUSED Binstar *binstar) override;

protected:

    bool check_condition=false; /*!< auxiliary variable used to check when the circularisation has to be activated*/


    /**
     * Actual method to estimate the change on a due to the circularisation and set the variation.
     * It is a pure/virtual method, so it has to be implemented in all the other derived classes
     * @param binstar Pointer to the binary
     * @return EXIT_SUCCESS
     */
    virtual int circularise(_UNUSED Binstar *binstar) = 0;

    /**
     * Check the triggering condition of the circularisation
     * @param binstar Pointer to the binary
     * @return True if the RLO is ongoing
     */
    virtual bool check_activation(Binstar *binstar) const;


    //Log message
    /**
     * Standard log-message for circularisation. It just reportes the old and new orbital parameters
     * @param binstar Pointer to the binary
     * @param a_new semimajor axis after circularisation
     * @param e_new eccentricity a
     * @return
     */
    virtual std::string log_message_circ(Binstar *binstar, double a_new, double e_new) const;

};

/**
 * periastron option
 * Circularisation at periastron:
 * anew=a_old(1-e_old)
 *
 * After the circularisation the binary angular momentum decreases of a factor
 * (Jnew-Jold)/Jold = 1/sqrt(1+e) -1, it ranges from 5% for small eccentriciteis (e=0.01) to ~29% for e->1
 *
 */
class CircularisationPeriastron : public StandardCircularisation{

public:

    CircularisationPeriastron(_UNUSED IO *_io= nullptr, bool reg = true) : StandardCircularisation(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "periastron"; }
    static CircularisationPeriastron _circularisationperiastron;
    CircularisationPeriastron* instance() override {
        return (new CircularisationPeriastron(nullptr,false));
    }

protected:

    int circularise(_UNUSED Binstar *binstar) override;

};

/**
 * periastron_full option
 * Circularisation at periastron with a weaker condition
 *
 * The circularisation is same as in CircularisationPeriastron what change is the activation trigger
 * in CircularisationPeriastron, a RLO has to be currently ongoing while here the circularisation is activated
 * if a RLO could be triggered at the periastron assuming a circular orbit with a=r_periastron.
 * In practice we check that the stellar radius is larger than the RL_periastron estimated as r_periastron/a * RL_eggleton.
 */
class CircularisationPeriastronFull : public CircularisationPeriastron{

public:

    CircularisationPeriastronFull(_UNUSED IO *_io= nullptr, bool reg = true) : CircularisationPeriastron(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "periastron_full"; }
    static CircularisationPeriastronFull _circularisationperiastronfull;
    CircularisationPeriastronFull* instance() override {
        return (new CircularisationPeriastronFull(nullptr,false));
    }

protected:

    /**
     * Check the triggering condition of the circularisation
     * Check if the RLO triggering conditions (R>RLO) are satisfied at the periastron.
     * The RLO at periastron is estimated assuming a circular orbit with a=R_periastron
     * @param binstar Pointer to the binary
     * @return True if the RLO conditions (R>RLO) are satisfied at the periastron, false otherwise
     */
     bool check_activation(Binstar *binstar) const override;

};

/**
 * angmom option
 * The binary is circularised conserving the angular momentum:
 * anew=a_old(1-e_old*e_old)
 */
class CircularisationAngMom : public StandardCircularisation{

public:

    CircularisationAngMom(_UNUSED IO *_io= nullptr, bool reg = true) : StandardCircularisation(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "angmom"; }
    static CircularisationAngMom _circularisationangmom;
    CircularisationAngMom* instance() override {
        return (new CircularisationAngMom(nullptr,false));
    }

protected:

    int circularise(_UNUSED Binstar *binstar) override;

};

/**
 * semimajor option
 * The binary is circularised conserving the semimajor axis, so just putting e=0:
 * anew=a_old
 *
 * After the circularisation the binary angular momentum increases of a factor
 * (Jnew-Jold)/Jold = 1/sqrt(1-e^2) -1,  from ~0.5% for small eccentricities (e=0.01) to >60% for large eccentricities (e>0.8)
 */
class CircularisationSemimajor : public StandardCircularisation{

public:

    CircularisationSemimajor(_UNUSED IO *_io= nullptr, bool reg = true) : StandardCircularisation(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "semimajor"; }
    static CircularisationSemimajor _circularisationsemimajor;
    CircularisationSemimajor* instance() override {
        return (new CircularisationSemimajor(nullptr,false));
    }

protected:

    int circularise(_UNUSED Binstar *binstar) override;

};




#endif //SEVN_CIRCULARISATION_H
