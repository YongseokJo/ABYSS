//
// Created by iorio on 5/21/24.
//

#ifndef SEVN_WINDACCRETION_H
#define SEVN_WINDACCRETION_H

#include <IO.h>
#include <Processes.h>
#include <binstar.h>
#include <star.h>

//Forward declaration
class Star;
class Binstar;


class Windaccretion : public MaccretionProcess {

public:
    Windaccretion(_UNUSED IO* _io = nullptr, bool reg = true) {
        if (reg) {
            //First register to include this Property in static array
            Register(this, &ID, name());
        }
        betaw = 0.125;
        alphaw = 1.5;
        muw = 1.0;

        if (_io!= nullptr){
            alphaw = _io->svpar.get_num("w_alpha");
            betaw  = _io->svpar.get_num("w_beta");
        }
    }

    static size_t ID;
    static Windaccretion _windaccretion;
    inline std::string name() override { return "Windaccretion"; }

    Windaccretion* Instance(_UNUSED IO* _io) override;
    virtual Windaccretion* instance() { return nullptr; }


    ///----------

    static constexpr int NO_WINDACCRETION = 0;
    virtual double DA(_UNUSED Binstar* b, _UNUSED int procID) { return 0; }
    virtual double DE(_UNUSED Binstar* b, _UNUSED int procID) { return 0; }
    inline bool is_process_ongoing() const override { return false; }

    ///-------------

protected:

    static std::map<std::string, Windaccretion*>& GetStaticMap() {
        static std::map<std::string, Windaccretion*> _locmap;
        return _locmap;
    }
    static std::vector<int>& GetUsed() {
        static std::vector<int> _used;
        return _used;
    }
    static void Register_specific(Windaccretion* ptr) {
        GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size() - 1] = 0;
    }

    double betaw; /*!< Wind escape velocity parameter (Eq. 9, Hurley+02) */
    double alphaw; /*!< Bondi-Hoyle accretion parameter (Eq. 6, Hurley+02) */
    double muw; /*!< Angular momentum  transfer efficiency  (Eq. 11, Hurley+02) */

};

class WindaccretionHurley : public Windaccretion {

public:

    WindaccretionHurley(_UNUSED IO* _io = nullptr, bool reg = true) : Windaccretion(nullptr, false) {
        if (reg) {
            Register_specific(this);
        }
    }

    inline std::string name() override { return "hurley"; }
    static WindaccretionHurley _windaccretionhurley;

    WindaccretionHurley* instance() override {
        return (new WindaccretionHurley(nullptr, false));
    }

    int evolve(Binstar* binstar) override;

    double estimate_accreted_mass(_UNUSED  double DM,  _UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binstar) const  override;
    double DA(_UNUSED Binstar* b, _UNUSED int procID) override;
    double DE(_UNUSED Binstar* b, _UNUSED int procID) override;
    inline bool is_process_ongoing() const override { return true; }

    /**
     * The accretion efficiency of the Hurley implementation diverges for orbital velocity
     * similar to the wind velocity. In this case in the Hurley paper a constant efficiency is used.
     * Here, we create this method flexible enough to be used in possible derived classes
     * @param donor Pointer to the wind donor star
     * @param accretor  Pointer to the wind accretor
     * @param binstar  Pointer to the binary
     * @return maximum accretion efficiency
     */
    virtual double estimate_maximum_accretion_efficiency ( _UNUSED Star *donor,
                                                   _UNUSED Star *accretor,
                                                   _UNUSED Binstar *binstar) const {
        return 0.8;  //From Hurley maximum possible mass fraction  that can be accreated
    }


protected:

    int accrete_mass(Binstar* binstar) override;


};

/**
 * This formalism is a modified version of the Hurley formalism.
 * The only change is about the maximum accretion efficiency.
 * In the Hurley+02 version it is set to 0.8, here instead we use
 * the maximum accretion efficiency derived by Tejeda&Toala24 (https://arxiv.org/pdf/2411.01755)
 * in the limit of high orbital velocity (in which vorb>=vwind).
 */
class WindaccretionHurleyTT24mod : public WindaccretionHurley {

public:

    WindaccretionHurleyTT24mod(_UNUSED IO* _io = nullptr, bool reg = true) : WindaccretionHurley(nullptr, false) {
        if (reg) {
            Register_specific(this);
        }
    }

    inline std::string name() override { return "hurleyTT24mod"; }
    static WindaccretionHurleyTT24mod _windaccretionhurleytt24mod;

    WindaccretionHurleyTT24mod* instance() override {
        return (new WindaccretionHurleyTT24mod(nullptr, false));
    }

    /**
     * We use the results from the paper Tejeda&Toala24 (https://arxiv.org/pdf/2411.01755).
     * They developed a new analytic prescriptions for the wind accretions that works also
     * when the orbital and wind velocity are similar. In this case the mass accretion efficiency
     * reach a plateau at q^2 where q is m_accretore/Mtot. Their formalism is not directly implemented
     * in this class, we just use their results to include a more physically motivated maximum accretion.
     * Overall, we will still have a larger efficiency in the transition regime (between vorb>vwind and vorb<<vwind),
     * @param donor Pointer to the wind donor star
     * @param accretor  Pointer to the wind accretor
     * @param binstar  Pointer to the binary
     * @return
     */
    double estimate_maximum_accretion_efficiency ( _UNUSED Star *donor,
                                                           _UNUSED Star *accretor,
                                                           _UNUSED Binstar *binstar) const override;

};

class WindaccretionDisabled : public Windaccretion {

public:

    WindaccretionDisabled(_UNUSED IO* _io = nullptr, bool reg = true) : Windaccretion(nullptr, false) {
        if (reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "disabled"; }
    static WindaccretionDisabled _windaccretiondisabled;
    WindaccretionDisabled* instance() override {
        return (new WindaccretionDisabled(nullptr, false));
    }
    int evolve(_UNUSED Binstar* binstar) override;

    inline bool is_process_ongoing() const override { return false; }
};



#endif //SEVN_WINDACCRETION_H
