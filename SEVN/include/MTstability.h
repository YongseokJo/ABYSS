//
// Created by Giuliano Iorio on 17/11/2021.
//
// Classes to handle the MT stability
//
//

/**
 *
 * The class MTstability is used to estimate the stability of the RLO mass transfer.
 * The idea behind the class is that the instability is "in general" checked through the condition:
 *  value > tshold (or < tshold).
 *   For example in a classical qcrit implementation, q > qcrit is the condition for unstable mass transfer.
 *   The key methods are:
 *   - bool mt_unstable(Star *donor, Star *accretor, Binstar * binary): This is the method that
 *   returns true if the condition of instability is satisfied.
 *   - bool get(Star *donor, Star *accretor, Binstar * binary): This method returns the value used to check the
 *   stability, e.g. for the classical qcrit formalism, get returns q=Mdonor/Maccretor
 *   - bool get_tshold(Star *donor, Star *accretor, Binstar * binary): this method returns the threshold value used in
 *   the comparison, e.g. for the classical qcrit formalism, get_tshold returns qcrit.
 *
 *   NOTICE: WE ASSUME that MTstability and the inherited classes are used during the BInary processes RLO, therefore
 *   the stellar properties ARE THE OBTAINED WITH GETP_0 (properties at the previous sse step), while the binary properties
 *   are always OBTAINED WITH GETP
 *
 */

/**
 *
 * HOW TO ADD A NEW MTstability IMPLEMENTATION
 *
 *
 *
 * 1- Define the new class, public inherited from MTstability, or another class (e.g. MT_Qcrit)
 * e.g.  class Qcrit_Hurley : public MT_Qcrit
 *
 * 2- Implement the constructor
 * e.g.
 *     Qcrit_Hurley(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }
 *
 * 3- Override the method name, it returns the name of this MT_stability option
 * NOTICE this will be the same name used to choice this MT formalism through the parameter
 * rlo_stability.
 * e.g.
 *       inline std::string name() override { return "qcrit_hurley"; }
 * this particular MT formalism can be called with -rlo_stability qcrit_hurley.
 *
 * 4- Define a static class member that is an instance  of the same class we are defining
 * e.g. static Qcrit_Hurley _qcrit_hurley;
 *
 * 5- Define the method instance that return a new a pointer to a dinamically allocated instance of the class
 * e.g.
 *
 *     Qcrit_Hurley* instance() {
        return (new Qcrit_Hurley(false));
    }
    NOTICE: remember to use false when calling the constructor here
 *
 * 6- If the class is derived directly from MT_stability, you have to define the following public pure virtual classes
 *  - mt_unstable(Star *donor, Star *accretor, _UNUSED Binstar *binary)  const override;
 *  - double get_tshold(Star *donor, Star *accretor, _UNUSED Binstar *binary) const override;
 *  - double get(Star *donor, Star *accretor, _UNUSED Binstar *binary) const override;
 *  If the class is derived from a subclass of MT_stability check what are the methods to be overrided.
 *  e.g. Qcrit_Hurely is derived from MT_Qcrit that already overidde mt_unstable, get and get_tshold, however
 *  there are two  virtual protected methods:
 *      - double q(Star *donor, Star *accretor) const;
 *      - virtual double qcrit(Star *donor, Star *accretor) const = 0;
 *  In particular qcrit is a pure virtual method that has to be implemented in the derived classes.
 *
 *  7- Add the static class member in general/static_main.h
 *  e.g.     Qcrit_Hurley Qcrit_Hurley::_qcrit_hurley;
 *
 * *   NOTICE: WE ASSUME that MTstability and the inherited classes are used during the BInary processes RLO, therefore
 *   the stellar properties ARE THE OBTAINED WITH GETP_0 (properties at the previous sse step), while the binary properties
 *   are always OBTAINED WITH GETP
 */



#ifndef SEVN_MTSTABILITY_H
#define SEVN_MTSTABILITY_H

#include <string>
#include <map>
#include <utility>
#include <vector>

#include <sevnlog.h>

#define set_prop(star,prop,old) old ? star->getp_0(prop::ID) :  star->getp(prop::ID)


class Star;
class Binstar;
using sevnstd::SevnLogging;

/**
 * Base pure virtual class
 */
class MTstability{

public:
    virtual inline std::string name() = 0;
    static MTstability *Instance(std::string const &name);
    virtual MTstability *instance(){ svlog.critical("This is supposed to be pure virtual, do you have an uninitialized module?", __FILE__, __LINE__); return nullptr; }

    /**
     * Check if the mass transfer is stable or not
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return true is the mass transfer is unstable, false otherwise
     */
    virtual bool mt_unstable(Star *donor, Star *accretor, _UNUSED Binstar *binary) const = 0;

    virtual ~MTstability(){
    }

    /**
     * Return the threshold  value used to make the comparison for mass stability
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return Return the threshold  value used to make the comparison for mass stability
     */
    virtual double get_tshold(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binary) const = 0;

    /**
     * Return the value used to make the comparison for mass stability
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return Return the value used to make the comparison for mass stability
     */
    virtual double get(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binary) const = 0;

protected:

    void Register(MTstability *ptr, const std::string &_name) {
        GetStaticMap().insert(std::make_pair(_name, ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }

    SevnLogging svlog;

private:

    static std::map<std::string, MTstability *> & GetStaticMap(){
        static std::map<std::string, MTstability *> _locmap;
        return _locmap;
    }

    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }

};


class MT_Stable : public  MTstability{

public:
    MT_Stable(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "stable"; }

    /**
     * Check if mass transfer is unstable
     * NOTICE: in this option mass transfer is always stable, therefore this method returns always false
     * @param donor  Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return Just a large value (1E30)
     */
     bool mt_unstable(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binary)  const override {return false;}

    /**
     * Return the threshold  value used to make the comparison for mass stability
     * NOTICE: In this option all the mass transfer are stable, so this method is not used
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return Just a large value (1E30)
     */
    double get_tshold(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binary) const override {return 1E30;}

    /**
     * Return the value used to make the comparison for mass stability
     * NOTICE: in this case the mass transfer is always stable, so wr are not checking any property.
     * This method returns just 0
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return 0
     */
    double get(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binary) const override {return 0.;}

    static MT_Stable _mt_stable;

    MT_Stable* instance() override {
        return (new MT_Stable(false));
    }

};

class MT_UnStable : public  MTstability{

public:
    MT_UnStable(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "unstable"; }

    /**
     * Check if mass transfer is unstable
     * NOTICE: in this option mass transfer is always unstable, therefore this method returns always true
     * @param donor  Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     */
    bool mt_unstable(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binary)  const override {return true;}

    /**
     * Return the threshold  value used to make the comparison for mass stability
     * NOTICE: In this option all the mass transfer are always unstable, so this method is not used
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return Just a negative large value
     */
    double get_tshold(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binary) const override {return -1E30;}

    /**
     * Return the value used to make the comparison for mass stability
     * NOTICE: in this case the mass transfer is always stable, so wr are not checking any property.
     * This method returns just 0
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return 0
     */
    double get(_UNUSED Star *donor, _UNUSED Star *accretor, _UNUSED Binstar *binary) const override {return 0.;}

    static MT_UnStable _mt_unstable;

    MT_UnStable* instance() override {
        return (new MT_UnStable(false));
    }

};


/**
 * MT stability implementation based on mass_ratio (q) and critical mass ratio (qc).
 * The instability is check as q > qc.
 * qc can be estimated with different formalism (they have to be implemented)
 */
class MT_Qcrit : public MTstability{

public:

    bool mt_unstable(Star *donor, Star *accretor, _UNUSED Binstar *binary)  const override;

    /**
     * Return qcrit, the critical mass-ratio, if q (get)>qcrit (get_tshold), the mst is unstable
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return qcrit
     */
    double get_tshold(Star *donor, Star *accretor, _UNUSED Binstar *binary) const override {return qcrit(donor,accretor);}
    /**
     * Return q, the critical mass-ratio, if q (get)>qcrit (get_tshold), the mst is unstable
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return qcrit
     */
     double get(Star *donor, Star *accretor, _UNUSED Binstar *binary) const override {return q(donor,accretor);}

protected:

    /**
     * Mass ratio in the previuous step
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return  Mass ratio between donor and accretor in the previous step
     */
    virtual double q(Star *donor, Star *accretor) const;

    /**
     * Critical mass ratio
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit in the previous step
     */
    virtual double qcrit(Star *donor, Star *accretor) const = 0;

};

/**
 * Qcrit stability based on the Hurley+02 formalism
 */
class Qcrit_Hurley : public MT_Qcrit{

public:

    Qcrit_Hurley(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "qcrit_hurley"; }

    static Qcrit_Hurley _qcrit_hurley;

    Qcrit_Hurley* instance() override {
        return (new Qcrit_Hurley(false));
    }

protected:

    /**
     * Qcrit from Hurley+02, the estimate is based on BSE phase,
     * NOTICE that we estimate the BSE phase in the previous step (using getp_0)
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit
     */
    double qcrit(Star *donor, Star *accretor) const override;

    /**
     * Qcrit for giants (BSE phase 3,5,6) following
     * qc = (1.67-zpars(7)+2.0*pow(mcore1/m1),5.))/2.13 Eq.57 in Hurely+02
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit for bse phase 3,5,6
     */
    virtual double qcrit_giant(Star *donor, _UNUSED Star *accretor) const;

};

/**
 * Same as Qcrit_Hurley but with qcrit for giants estimated using equations from Webbink (1988)
 */
class Qcrit_Hurley_Webbink : public Qcrit_Hurley{

public:

    Qcrit_Hurley_Webbink(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "qcrit_hurley_webbink"; }

    static Qcrit_Hurley_Webbink _qcrit_hurley_webbink;

    Qcrit_Hurley_Webbink* instance()  override {
        return (new Qcrit_Hurley_Webbink(false));
    }

protected:

    /**
     * Qcrit for giants (BSE phase 3,5,6) following
     * qc = 0.362 + 1.0/(3.0*(1.0 - Mcore_fraction)) After Eq. 57 (not numbered) in Hurley+02
     * The equation comes from Webbink (1988) and it has been obtained from models
     * for condensed polytropes.
     * @param donor
     * @param accretor
     * @return
     */
    double qcrit_giant(Star *donor, _UNUSED Star *accretor) const override;

};

/**
 * Same qcrit as Hurley_Webbink but with a special treatment when the accretor is a BH coming from
 * Shao+21 paper (https://arxiv.org/pdf/2107.03565.pdf).
 */
class Qcrit_Hurley_Webbink_Shao : public Qcrit_Hurley_Webbink {

public:

    Qcrit_Hurley_Webbink_Shao(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "qcrit_hurley_webbink_shao"; }

    static Qcrit_Hurley_Webbink_Shao _qcrit_hurley_webbink_shao;

    Qcrit_Hurley_Webbink_Shao* instance() override {
        return (new Qcrit_Hurley_Webbink_Shao(false));
    }

protected:

    /**
     * Same qcrit from Hurley+02 with Webbink 1988 for giants,
     * but with a special treatment accretion on a BH following the Shao+21 paper.
     * In the Shao+21 paper (https://arxiv.org/pdf/2107.03565.pdf) there is a special treatment
     * for the accretion on a BH. Therefore if the accretor is not a BH we use the classical qc from Hurley.
     * In Shao+21 there is not really a constant value for qc, but there are a number of criteria that  we
     * have to take into account to estimate the stability of the mass transfer, then if the mass trasnfer is stable we
     * use a very large value for qcrit otherwise a very small one.
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit
     */
    double qcrit(Star *donor, Star *accretor) const override;

};

/**
 * qcrit prescriptions taken from Neijssel+2020, Sec. 2.3
 *
 */
class Qcrit_COSMIC_Neijssel : public Qcrit_Hurley_Webbink{

public:

    Qcrit_COSMIC_Neijssel(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "qcrit_cosmic_neijssel"; }

    static Qcrit_COSMIC_Neijssel _qcrit_cosmic_neijssel;

    Qcrit_COSMIC_Neijssel* instance() override {
        return (new Qcrit_COSMIC_Neijssel(false));
    }

protected:
    /**
     * qcrit prescriptions taken from Neijssel+2020, section 2.3  as implemented in COSMIC
     * This implementation is taken directly from COSMIC (https://github.com/COSMIC-PopSynth/COSMIC/blob/develop/cosmic/src/evolv2.f
     * In the documentation they report:
     * "We convert from radial response to qcrit for MS and HG,
     * which assumes conservative mass transfer,
     * Stable MT is always assumed for stripped stars,
     * Assume standard qcrit from BSE for kstar>=10"
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit
     */
    double qcrit(Star *donor, Star *accretor) const override;

};

/**
 * Qcrit implementation in Cosmic based on Claeys+04
 */
class Qcrit_COSMIC_Claeys: public Qcrit_Hurley_Webbink{

public:

    Qcrit_COSMIC_Claeys(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "qcrit_cosmic_claeys"; }

    static Qcrit_COSMIC_Claeys _qcrit_cosmic_claeys;

    Qcrit_COSMIC_Claeys* instance() override {
        return (new Qcrit_COSMIC_Claeys(false));
    }

protected:
    /**
     * qcrit prescriptions implemented in COSMIC (evolv2.f, row 2060), following Claeys+14
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit
     */
    double qcrit(Star *donor, Star *accretor) const override;


};

/**
 * Qcrit implementation based on Belczynski+2008
 */
class Qcrit_StarTrack : public MT_Qcrit{

public:

    Qcrit_StarTrack(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "qcrit_startrack"; }

    static Qcrit_StarTrack _qcrit_startrack;

    Qcrit_StarTrack* instance() override {
        return (new Qcrit_StarTrack(false));
    }

protected:

    /**
     * Qcrit from Belczynski+08, the estimate is based on BSE phase,
     * NOTICE that we estimate the BSE phase in the previous step (using getp_0)
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit
     */
    double qcrit(Star *donor, _UNUSED  Star *accretor) const override;

};

/**
 * Test qcrit where MS, HG are always stable
 */
class Qcrit_Radiative_Stable : public Qcrit_Hurley_Webbink{
public:

    Qcrit_Radiative_Stable(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "qcrit_radiative_stable"; }

    static Qcrit_Radiative_Stable _qcrit_radiative_stable;

    Qcrit_Radiative_Stable* instance() override {
        return (new Qcrit_Radiative_Stable(false));
    }

protected:
    /**
     * Just a wrapper of the Qcrit_Hurley_Webbink qcrit with
     * the difference that qc for MS (1), HG(2) and pureHE stars(7-9) is set to 1000
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit
     */
    double qcrit(Star *donor, Star *accretor) const override;

};

/**
 * Test qcrit where MS, HG and pureHe stars are always stable
 */
class Qcrit_HRadiative_Stable : public Qcrit_Hurley_Webbink{
public:

    Qcrit_HRadiative_Stable(bool  reg = true){
        if (reg){
            Register(this, name());
        }
    }

    inline std::string name() override { return "qcrit_Hradiative_stable"; }

    static Qcrit_HRadiative_Stable _qcrit_hradiative_stable;

    Qcrit_HRadiative_Stable* instance() override {
        return (new Qcrit_HRadiative_Stable(false));
    }

protected:
    /**
     * Just a wrapper of the Qcrit_Hurley_Webbink qcrit with
     * the difference that qc for MS (1), HG(2) is set to 1000
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit
     */
    double qcrit(Star *donor, Star *accretor) const override;

};

class MT_Zeta : public MTstability{

public:

    bool mt_unstable(Star *donor, Star *accretor, _UNUSED Binstar *binary)  const override;

    /**
     * Retrun qcrit, the critical mass-ratio, if q (get)>qcrit (get_tshold), the mst is unstable
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return qcrit
     */
    virtual double get_tshold(Star *donor, Star *accretor, _UNUSED Binstar *binary) const {return qcrit(donor,accretor);}
    /**
     * Retrun q, the critical mass-ratio, if q (get)>qcrit (get_tshold), the mst is unstable
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @param binary  Pointer to the binary
     * @return qcrit
     */
    virtual double get(Star *donor, Star *accretor, _UNUSED Binstar *binary) const {return q(donor,accretor);}

protected:

    /**
     * Mass ratio in the previuous step
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return  Mass ratio between donor and accretor in the previous step
     */
    virtual double q(Star *donor, Star *accretor) const;

    /**
     * Critical mass ratio
     * @param donor Pointer to the donor star
     * @param accretor  Pointer to the accretor star
     * @return qcrit in the previous step
     */
    virtual double qcrit(Star *donor, Star *accretor) const = 0;

};



#endif //SEVN_MTSTABILITY_H

