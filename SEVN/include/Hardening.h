//
// Created by iorio on 3/30/22.
//

#ifndef SEVN_HARDENING_H
#define SEVN_HARDENING_H

#include <IO.h>
#include <Processes.h>


//Forward declaration
class Star;
class Binstar;

/**
 * Hardenig via three-body encounters that reduce the semi-major axis
 * and increase the binary eccentricity. This process should be enabled
 * if we want to simulate binaries inside dense stellar environments.
 * The prescription here adopted follows Mapelli+2021 and Heggie1975.
 */
class Hardening : public  Process{

public:

    Hardening(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            //First register to include this Property in static array
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Hardening _hardening;
    inline std::string name() override { return "Hardening"; }

    Hardening *Instance(_UNUSED IO *_io) override;
    virtual Hardening *instance() {return nullptr;}

    bool is_process_ongoing() const override {return false;};
    static std::string log_message(_UNUSED Binstar *binstar);


protected:

    //this a prototype for the future refactoring of the whole Process classes
    static std::map<std::string, Hardening *> & GetStaticMap(){
        static std::map<std::string, Hardening *> _locmap;
        return _locmap;
    }
    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }
    static void Register_specific(Hardening *ptr) {
        GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }


};

/**
 * If Hardenig is disabled, we do not include the dynamical effects of the
 * host environment that may affect the binary evolution.
 */
class HardeningDisabled : public Hardening{

public:

    HardeningDisabled(_UNUSED IO *_io= nullptr, bool reg = true) : Hardening(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "disabled"; }
    static HardeningDisabled _hardeningdisabled;
    HardeningDisabled* instance() override {
        return (new HardeningDisabled(nullptr,false));
    }
    int evolve(_UNUSED Binstar *binstar) override;
};

/**
 * This Hardening prescription follows eq.11 and eq.13 in Mapelli+2021
 */
class HardeningFastCluster : public Hardening{

public:

    HardeningFastCluster(_UNUSED IO *_io= nullptr, bool reg = true) : Hardening(nullptr, false){
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "fastcluster"; }
    static HardeningFastCluster _hardeningfastcluster;
    HardeningFastCluster* instance() override {
        return (new HardeningFastCluster(nullptr,false));
    }
    int evolve(_UNUSED Binstar *binstar) override;

protected:
    /**
     * Estimate if the binary can be considered an hard binary
     * @param binstar Pointer to the binary
     * @return true if |Ebinary|>Ekin_average, where
     *      Ebinary orbital energy of the binary
     *      Ekin_average=0.5*sigma^2*m_average (sigma velocity dispersion of the environment, m_average average mass of the perturbers)
     */
    bool is_hard_binary(_UNUSED Binstar *binstar);

};



#endif //SEVN_HARDENING_H
