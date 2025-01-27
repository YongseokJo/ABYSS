/*
Created by spera on 30/12/18.
Module for SN models, it includes:
 - Class supernova: Base class for SN models
 - Class PisnON:  auxiliary  class-option used to enable pair instability at compile time
 - Class PisnOFF: auxiliary class-option used to disable pair instability at compile time
 - Class NeutrinoMassLossOFF: auxiliary class-option used to disable neutrino mass loss at compile time
 - Class NSfromGau: auxiliary class-option used to force to draw Neutron Star masses from a Gaussian
*/

#ifndef SEVN_SUPERNOVA_H
#define SEVN_SUPERNOVA_H

/**
 * This header contains the class Supernova
 * Supernova is used by the Star class when a stellar remnant has to be created.
 * Notice that this class (and its inherited models) handles also the WD formation even if it is not a SN event.
 *
 *   WHO USE THE CLASS?
 * The main method that is called from the star class is main(Star *s), called when a remnant is created.
 * The other methods that could be used in the Star or Binary class are:
 *      -  void initialise_remnant(Star *s, double Mass_remnant, Lookup::Remnants Remnant_type), called to initiliase
 *        a remnant directly from the initial conditions
 *      -   void explosion_SNI(Star *s) or void explosion_SNI(Binstar *b), called when a SNIa is triggered
 *
 *   WHAT IS SET AT THE CLASS INITIALISATION
 *   In the Supernova constructor the class members are initialised (see below) and two important objects are created
 *   on the heap:
 *      -Kick model: associated to the pointer kick, this is the kick model that will be called during the remnant formation
 *      -PairInstability model: associated to the pointer pairinstability, this is the pairinstability model that will be called during the remnant formation
 *      -NeutrinoMassLoss model: associated to the pointer neutrinomassloss, this is the neutrino mass loss model that will be called during the remnant formation
 *
 *   WHAT ARE THE IMPORTANT CLASS MEMBERS TO SET:
 *   The key class member to set are:
 *      - Mremnant: mass of the remnant, this is used to set the property of the stellar remnant and usually also to set
 *                  the remnant_type. It can be updated more than once, for example applying pair instability corrections
 *      - remnant_type: type of the remnant, it has to be one of the item in the enum Lookup::Remnants, this is a key
 *                  value used to set the stellar remnant object
 *      - fallback_frac: fallback_frac, fraction of the final mass that return on the collapsing core, it is used
 *                       in some Kick model
 *      - pisn_correction: fraction between the remnant mass before and after the pisn_correction (currently it is not used
 *                          outside the class)
 *      - Mejected: preSN mass - Mass of the remnant
 *      - M_neutrinos: mass lost through neutrinos (it is set by the method  Mass_corrections_after_explosion)
 *      - sn_type: type of the SN, an item of the enum Lookup::SNExplosionType (it is set through the pisn correction)
 *      - Average_ejected, Average_remnant: value containing the average mass of the remnant and average mass of the ejected
 *                                          given a certain populatin. They are used in certain kick models (e.g. Unified)
 *                                          Usually they are set in the class Constructor
 *
 *
 *
 *   HOW THE CLASS WORK?
 *   When a SN (or a general remnant formation) has to be called, the main method is main.
 *   When main is called:
 *          - Information of the current stellar  masses are stored (Mass, MHE, MCO)
 *          - The method remnant_properties(Star *s) is called,
 *             this call set all the remnant properties and the remnant object
 *          - Some additional set (e.g. set isremnant to true) calling s->set_remnant();
 *          - Kicks are applied calling kick->apply(s), kick is a pointer to a Kick model
 *          - If allowed push information in the logs
 *
 *   Most of the setting are inside the method remnant_properties(Star *s):
 *          - First, it applies the preSN correction for pair instability, calling set_preSN_Masses(s),
 *            this method set the member struct Mass_preSN, it depends on the pointer pairinstability and the associated model
 *            NOTICE: from this point onward the mass of the struct Mass_preSN has to be used instead of the mass from the getp
 *            method of the star class. To retrive the preSN mass use the method getpreSN_Mass(x) where x is Mass::ID, MHE::ID or MCO::ID
 *          - Depending the preSN mass of the CO one of the following methods are called:
 *              - ECSN to create an electron capture supernova
 *              - WD to from a white dwarf
 *              - CCexplosion to trigger a CC SN explosion (it can be also a PISN or PPISN depending on the pairinstability correction)
 *              In these functions the key members are set
 *           - Finally, the method set_staremnant(s) is called:
 *              - Depengin of the remnant_type, a remnant object is created and set to the Star pointer staremnant
 *
 *   HOW PISN AND NEUTRINO MASS LOSS ARE TAKEN INTO ACCOUNT?
 *          The correction for PISN and Mass loss are applied through the method  Mass_corrections_after_explosion(const double mass, Star *s)
 *          In the method two other methods are called:
 *              - pisn(const double mass, Star *s)
 *                     pisn is a pure virtual method, its implementation is made in two auxiliary inherited class
 *                          -PisnON: here the pisn just use the pointer pairinstability to call apply_afterSN(s, mass);
 *                                   return the updated  Mremnant plus  and set pisn_correction and sn_type.
 *                          -PisnOFF: here the pisn are disabled, therefore the Mrenant is not corrected
 *                                    pisn_correction is set to 1 and sn_type to Lookup::SNExplosionType::CoreCollapse
 *                           NOTICE: PisnON e PisnOFF also set the propert preSN mass, PisnON apply the pairinstability
 *                                  method pairinstability->apply_beforeSN(s), while PisnOFF just the set the preSN mass
 *                                  to the current stellar masses (Mass, MHE, MCO)
 *                           NOTICE2: Each new Supernova model must inherit one of the two auxiliary class depending if the
 *                                    model has to be corrected for pisn (PisnON) or not (PisnOFF, for example because
 *                                    the model already take into accoun the pisn in the estimate of the final remnant mass).
 *
 *              - neutrino_mass_loss(const double mass, Star *s)
 *                      Similar to pisn, in its default implementation
 *
 *
 *
**/

/**
 * ADD NEW SUPERNOVA MODELS
 * A SN class model need to have the following feature:
 * - Be a class derived from class supernova or one of its child
 * - If is derived directly from Supernova, it is necessary do inherit also from one between:
 *      -PisnON: to enable pisn correction
 *      -PisnOFF: to disable pisn correction
 * - If the neutrino mass loss has to be disabled, inherit from
 *      -NeutrinoMassLossOFF
 * - Override of the method class* instance(Star *s) is mandatory.
 * - Override of the method inline std::string name() is mandatory.
 * - Override of the method void explosion(Star *s) is mandatory (it is a pure virtual function in the base class).
 *
 * Other function can be overridden, but usually is not necessary (see description in the Supernova class):
 * - set_remnant_type_after_explosion(Star *s, double mremnant, Lookup::SNExplosionType sntype)
 * - virtual void CCexplosion(Star *s);
 * - virtual void ECSN(Star *s);
 * - virtual void WDformation(Star *s);
 *
 * All the other methods should not be modified (even if declared virtual)
 *
 * Implementation steps
 * Assume we are implementing a new SN class  Dummy, we suggest to
 * create a new header file called dummy.h in include/star/procs/supernova/
 * and the related source file dummy.cpp in src/star/procs/supernova/
 * Remember to add #include <dummy.h> in the static_main.h and in dummy.cpp
 * Since it is possible that you will use the class star, we suggest to forward declare it in dummy.h
 * class Star;
 * and then #include <star.h> in dummy.cpp, remember also to add #include <supernova.h> in dummy.h
 *
 * 1- Declare the class:
 *      the class has to be inherited from the base class Supernova of another supernova inherited class.
 *      if inherited from supernova, it has to inherit also from  one of the two auxiliary PisnON or PisnOFF  depending if the
 *      model has to be corrected for pisn or not. It can also inherit from NeutrinoMassLossOff to disable neutrino mass loss
 *      and from NSfromGau if one wants to force the NS mass to be drawn from a Gaussian.
 *      For example assume that for Dummy we want to enable pisn correction, disable neutrino mass loss and force NS mass drawn from a a Gaussian,
 *      the class declaration will be:
 *      class Dummy : virtual public supernova, PisnON, NeutrinoMassLossOFF, NSfromGau{  };
 *
 * 2: Define and Implement the class constructor
 *      in dummy.h:
 *          Dummy(Star *s = nullptr);
 *      in dummy.cpp:
 *          dummy::dummy(Star *s) : supernova(s){
 *          if (s==nullptr)
 *              Register(this);
 *          else{
 *              ... *Initialise class member  parameters, for example Average_ejected and Average_remnant
 *              - set Average_ejected and Average_remnant with  set_Average_for_Unified(s,default_Average_Mremnant,default_Average_Mejected)
 *              where default_Average_Mremnant and default_Average_Mejected should be given.
 *              If a negative value is given, we are de-facto disabling the Unified kick model, since
 *              SEVN will return an error if unified kick model is chosen   and default_Average_Mremnatn or  default_Average_Mejected is negative or zero
 *          }
 *        }
 *
 * 3: initialise the class static object (dummy.h):
 *      a)
 *      static Dummy _dummy;
 *          *the standard is to use as name of the member the name of the class (all lowercase)
 *          *with a an underscore as prefix.
 *      b) remember to include it in the static_main.h header
 *          Dummy Dummy::_dummy;
 *          *remember to include #include <dummy.h> in  static_main.h
 *
 * 4: override the method instance (dummy.h):
 *      Dummy* instance(Star *s){
 *          return (new Dummy(s));
 *       }
 *
 * 5: override the method name (dummy.h)
 *      inline std::string name() override {return "dummy";}
 *      ** NOTICE: the chosen name is the one that will be used to set the given SN model
 *                 at runtime through the parameter -snmode
 *
 *      ** NOTICE: the chosen name is the one that will be used to set the SN model
 *                 in the input binary list.
 *
 *
 * 6: override, define (supernova.h) and implement (supernova.cpp) the method void explosion(Star *s)
 *
 *      void explosion(Star*){
 *          *Here is mandatory to assign:
 *              *fallback_frac
 *              *Mremnant
 *           **DO NOT set other member unless you know exactly what you are doing and how this influence the other
 *           method calls
 *            * NOTICE: The fallback_frac member is  used in some Kick models  to correct the natal kick. If it is not set, in the kick
 *              class it is assumed to be =0 (no correction).
 *      }
 *
 *      For example assume that in our dummy model, our explosion just set the Mremnant to the minimum between 30 and preSN MCO
 *      and fallback_frac=0
 *      void explosion(Star *s){
 *          const double& MCOpreSN = get_preSN_Mass(MCO::ID);
 *          Mremnant = std::min(MCOpreSN,30);
 *          fallback_frac=0;
 *      }
 *
 *  7: If really needed override   other methods:
 *      - set_remnant_type_after_explosion(Star *s, double mremnant, Lookup::SNExplosionType sntype)
 *      - virtual void CCexplosion(Star *s);
 *      - virtual void ECSN(Star *s);
 *       - virtual void WDformation(Star *s);
 *       In order to override them, follow what is written in description HOW THE CLASS WORK? and in the methods
 *       documentation. Give a look to other models that override them (e.g. DeathMatrix).
 *
 *  8: Add documentation
 *      Try to document as much as possible
 *      In particular, use this template documentation for the class
 *      **
 *      * Brief model description
 *      * BH/NS: how BH and NS are set
 *      * NS mass: whether the NS mass is estimate in model or it is assumed
 *      * max NS mass: what is the max NS mass (if there is not a special limit just write from SEVN option sn_max_ns_mass)
 *      * ECSN: whether the ECSN method has been modified or not (if not, write from base supernova class)
 *      * WD: whether the WD formation has been modified or not (if not, write from base supernova class)
 *      * fallback: the fallback is used in the Hobbs (and derived) sn models.
 *                  Write here if fallback is estimated or just set to some value or unset.
 *                  If it is unset it is set on-the-fly in the SN model.
 *      * ppisn_correction: write applied if the class inherit from PisnON, not applied if inherit from PisnOFF
 *      * neutrino_correction: applied or not applied if the class inherit from NeutrinoMassLossOFF
 *      *
 *
 *      for example for the class Dummy
 *
 *      **
 *      * My dummy model
 *      * BH/NS: agnostic, based on the max NS mass
 *      * NS mass: drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 *                  but always larger than 1.1 Msun
 *      * max NS mass: from SEVN option sn_max_ns_mass
 *      * ECSN: from base supernova class
 *      * WD: from base supernova class
 *      * fallback: set to 0
 *      * ppisn_correction: applied
 *      * neutrino_correction: not applied
 *      class Dummy : virtual public supernova, PisnON, NeutrinoMassLossOFF, NSfromGau{
 **/




#include <string>
#include <map>
#include <vector>
#include <sevnlog.h>
#include <lookup_and_phases.h>
#include <kicks.h>
#include <random>
#include <neutrinomassloss.h>

class Star;
class PairInstability;
using sevnstd::SevnLogging;

#define _UNUSED __attribute__ ((unused))


class supernova {

/***** For Developers that want to add new supernova models: ******/
///Remember that each new proper supernova model has to have a public  static class istance to define in static_main.h
/// e.g. static MySNmodel _mysnmodel; (in the class) and MySNmodel MySNmodel::_mysnmodel; in static_main.h
/// Methods that must  be customized in new models.
/// !!!!NOTICE!!!!!: when overriding this function, remember to not use the s->getp to retrieve the preSN mass of the star
/// use instead the method get_preSN_Mass(x) where x can be Mass::ID, MHE::ID, MCO::ID
public:
    /**Name of the Supernova model, this will be use at runtime to select the  given model**/
    virtual inline std::string name() const { return "GenericSupernova"; }
    /**
     * Create an instance of the class allocated on the heap
     * @return a pointer to a real instance of the class that is allocated on the  heap
     */
    virtual supernova *instance(_UNUSED Star *s) = 0; //Pure virtual

    /**
     * CoreCollapse SN explosion, this method must to:
     * - set the remnant mass (member Mremnant)
     * - if allowed by the model set the fallback_fraction (fallback_frac)
     * NOTICE: DOT NOT  apply pisn or neutrino mass loss correction, DO NOT set remnant_type,
     * these will be automatically set by other methods in CCexplosion.
     * NOTICE: when overriding this function, remember to not use the s->getp to retrieve the preSN mass of the star,
     * use instead the method get_preSN_Mass(x) where x can be Mass::ID, MHE::ID, MCO::ID
     * @param s Pointer to the exploding star
     */
    virtual void explosion(Star *s) = 0;

/// Methods that may  be customized in new models
/// NOTICE-1: for most of the cases, in order to add a new model it is enough to just override the explosion method,
/// the methods below can be overridden in rare some special cases.

/// NOTICE-2: the pisn and neutrino mass loss corrections are handled through the virtual  method pisn and neutrino_mass_loss,
/// in principle, these methods can be overrided, however DON'T DO IT! The options for pisn and neutrino mass loss hare selected
/// through the proper classes PairInstability and NeutrinoMassLoss, create a new model of the classes if you want to change the current formalism.
/// DO NOT override the class,  If you want to disable the corrections. If you want to disable the pisn at-compile time for a given model
/// just inherit from the class PisnOFF, if you want to disable the neutrino mass loss inherit from NeutrinoMassLossOFF,
/// e.g. class MySNmodel : virtual public supernova, PisnOFF, NeutrinoMassLossOFF{
/// If you wan to disable the corrections at runtime just choose the option disabled, e.g. -sn_pairinstability disabled and/or -sn_neutrinomloss disabled
/// NOTICE: each new supernova model must inherit either from PisnON (to enable pisn correction) or from PisnOFF (to disable the pisn correction)

/// !!!!NOTICE-IMPORTANT!!!!!: when overriding this function, remember to not use the s->getp to retrieve the preSN mass of the star
/// use instead the method get_preSN_Mass(x) where x can be Mass::ID, MHE::ID, MCO::ID
protected:

    /**
     * Set the remnant type, i.e. set the member remnant_type based on the remnant Mass and the sntype.
     * The default version is quite simple, if Mass<Max Mass Neutron star (sn_max_ns_mass), the remnant is a NS_CCSN,
     * otherwise it is a BH. If Mremnant is 0, the remnant type is Empty. If sntype is Electron Capture the remnant
     * is NS_ECSN
     * @param s pointer to the exploding star
     * @param mremnant Final mass of the remnant
     * @param sntype SN type
     * @return the remnant type just set
     */
    virtual Lookup::Remnants  set_remnant_type_after_explosion(Star *s, double mremnant, Lookup::SNExplosionType sntype);

    /**
     * Handle the Core Collapse SN explosion, it calls the explosion method that is overridden in the SN models +
     * take into account pisn corrections and neutrino mass loss.
     * - Estimate and set the remnant mass (member Mremnant)
     * - Estimate the corrections for pisn and neutrino mass loss and set the corrected remnant mass
     *   (set members Mremnant, fallback_frac, M_neutrinos, pisn_correction, sn_type)
     * - Set the remnant type  (set remnant_type and/or sn_type)
     * @param s Pointer to the exploding star
     */
    virtual void CCexplosion(Star *s);

    /**
     * Handle the Electron Capture SN, it has to
     * - Estimate and set the remnant mass (member Mremnant)
     * - Estimate the corrections for pisn and neutrino mass loss and set the corrected remnant m
     *   (set members Mremnant, fallback_frac, M_neutrinos, sn_type)
     *   NOTICE: in the default supernova::ECSN only the neutrino mass loss is taken into account
     * - Set the remnant type  (set remnant_type and/or sn_type)  Lookup::SNExplosionType::ElectronCapture
     * NOTICE: when overriding this function, remember to not use the s->getp to retrieve the preSN mass of the star,
     * use instead the method get_preSN_Mass(x) where x can be Mass::ID, MHE::ID, MCO::ID
     * @param s Pointer to the exploding star
     */
    virtual void ECSN(Star *s);

    /**
     * Handle the WD formation, it has to
     * - Estimate and set the remnant mass (member Mremnant)
     * - Estimate and set the remnant_type (member remnant_type)
     * - Estimate and set the SN type (member sn_type) that should be Lookup::SNExplosionType::WDformation
     * NOTICE: when overriding this function, remember to not use the s->getp to retrieve the preSN mass of the star,
     * use instead the method get_preSN_Mass(x) where x can be Mass::ID, MHE::ID, MCO::ID
     * @param s Pointer to the exploding star
     */
    virtual void WDformation(Star *s);





    /*********************************************************************************************************/
    /***** STOP: DO NOT CHANGE THE REST OF THE CODE IF YOU ARE NOT FULLY AWARE OF WHAT  YOU ARE DOING ******/
    /*********************************************************************************************************/
public:
    //TODO Since Star is passed to the constructor, it can be stored as a pointer and the Star* s can be removed
    //in the other function calls
    supernova(Star *s = nullptr);

    virtual ~supernova();

    static supernova *Instance(std::string const &name, Star *s);


    void main(Star *s);

    /**
     * Method used to initialise directy a remnant without passing through an explosion
     * @param s Pointer to the star, used to set remnant properties
     * @param Mass_remnant  Mass of the remnant
     * @param Remnant_type Type of the Remnant, it should be an elemnt of the enum Lookup::Remnant
     */
    void initialise_remnant(Star *s, double Mass_remnant, Lookup::Remnants Remnant_type);


    inline double get_fallback_frac() const { return fallback_frac; }

    inline double get_pisncorrection() const { return pisn_correction; }

    inline double get_M_neutrinos() const { return M_neutrinos; }

    inline double get_Mremnant() const { return Mremnant; }

    /*** Return the Mass of the remnant before the pisn and neutrino mass loss correction ***/
    inline double get_Mremnant_preCorrection() const {
        return Mremnant_preCorrection<0 ? Mremnant : Mremnant_preCorrection;
    }

    inline double get_AverageRemnant() const { return Average_remnant; }
    inline double get_AverageEjected() const { return Average_ejected; }

    inline double get_Mejected() const { return Mejected; }
    inline double get_remnant_type() const { return remnant_type; }


    /**
     * Trigger a SNIa explosion
     * @param s Pointer to star
     */
    virtual void explosion_SNI(Star *s);

    /**
     * Trigger a SNIa explosion during a BSE evolution
     * @param s Pointer to star
     * @param b Pointer to binary
     */
    virtual void explosion_SNI(Star *s, Binstar *b);

    /**
     * Get the preSN mass from the object
     * @param massid if of the property, can be only Mass::ID, MHE::ID, MCO::ID
     * @return the value of the preSN property
     */
    double  get_preSN_Mass(size_t massid) const;


protected:
    double fallback_frac, M_neutrinos;
    double pisn_correction;
    double Mremnant, Mejected, Mremnant_preCorrection;
    double Average_ejected, Average_remnant;
    Lookup::Remnants remnant_type;
    Lookup::SNExplosionType sn_type=Lookup::SNExplosionType::Unknown; //initialised to Unknown

    void Register(supernova *ptr) {
        //Register only if the name is not already in the map,
        //this is necessary to avoid multiple registration of the same model when the constructor is called from
        //inherited classes
        if (GetStaticMap().find(ptr->name())==GetStaticMap().end()){
            GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
            GetUsed().resize(GetStaticMap().size());
            GetUsed()[GetUsed().size() - 1] = 0;
        }
    }


    /**
     * Estimate PISN and neutrino mass loss correction
     * @param mass Mremnant mass to correct
     * @param s Pointer to the exploding star
     * @return the final remnant mass after correction for pisn and neutrino mass loss
     */
    double Mass_corrections_after_explosion(const double mass, Star *s);

    /**
     * Apply pisn correction
     * @param mass Mass of the remnant before the pair isntability correction
     * @param s  Pointer to the exploding star to which we want to apply the correction
     * @return Mass after the PISN correction, it also set the sn_type
     */
    virtual double pisn(const double mass, Star *s) = 0;

    const PairInstability* get_pairinstability() const {return  pairinstability;}

    /**
     * Estimate the neutrino mass loss correction
     * @param mass
     * @param mass Mass of the remnant before the pair isntability correction
     * @param s  Pointer to the exploding star to which we want to apply the correction
     * @return Mass after the neutrino mass loss correction
     */
    virtual double neutrino_mass_loss(const double mass, Star *s) const;

    /**
     * Set the member Average_ejected and Average_remnant. These are used in the Unified Sn model
     * @param s Pointer to the star
     * @param default_Average_Mremnant  default Average_Mremnant for this SN model
     * @param default_Average_Mejected  default Average_Mejected for this SN model
     */
    void set_Average_for_Unified(Star* s, double default_Average_Mremnant, double default_Average_Mejected);


    SevnLogging svlog;

    /**
     * Set the pre SN mass container
     */
    virtual int set_preSN_Masses(Star *s) = 0;

    int _set_preSN_Masses(Star *s, bool allow_pisn);



private:
    Kicks *kick;
    PairInstability *pairinstability;
    NeutrinoMassLoss *neutrinomassloss;
    //std::unique_ptr<Kicks> *kick;

    //preSN Mass container
    utilities::MassContainer* Mass_preSN=nullptr;

    static std::map<std::string, supernova *> &GetStaticMap() {
        static std::map<std::string, supernova *> _locmap;
        return _locmap;
    }

    static std::vector<int> &GetUsed() {
        static std::vector<int> _used;
        return _used;
    }

    void remnant_properties(Star *s);

    /**
     * Use this function to properly initiliase the pointer to the object
     * It uses the protected member Mremnant
     * @param s Pointer to the star
     */
    void set_staremnant(Star *s);



};


/**
Auxiliary class to set classes for which we have to apply pair instability
*/
class PisnON : virtual public supernova{
public:
    PisnON(){}
protected:
    /**
    * Set the pre SN mass container
    */
    int set_preSN_Masses(Star *s) override;

    double pisn(const double mass, Star *s) override;

};

/**
Auxiliary class to set classes for which we don't have to apply pair instability correction
*/
class PisnOFF : virtual public supernova{
public:
    PisnOFF(){}
protected:
    /**
    * Set the pre SN mass container
    */
    int set_preSN_Masses(Star *s) override;

    double pisn(const double mass, Star *s) override;

};

/**
 * Auxiliary class to switch off neutrino mass loss
 */
class NeutrinoMassLossOFF : virtual public supernova{
public:
    NeutrinoMassLossOFF(){}
    NeutrinoMassLossOFF(_UNUSED Star *s){}
protected:
    double neutrino_mass_loss(const double mass, _UNUSED Star *s)  const override {return  mass;}
};

/**
Auxiliary class to handle the random generation of NS masses
*/
class NSfromGau : virtual public supernova {
public:
    NSfromGau(){}
    void ECSN(Star *s) override;
    double get_NS_mass(Star *s);

protected:
    std::normal_distribution<double> normal_dist{0.,1.};
    inline double generate_random_gau(double mean, double std){
      return normal_dist(utilities::mtrand)*std + mean;
    }

    /**
     * Handle the Core Collapse SN explosion, it calls the explosion method that is overridden in the SN models +
     * take into account pisn corrections and neutrino mass loss.
     * In this overridden method, it cals the Base class explotion and if the remnant is a neutron star
     * the mass is re-set and drawn from a Gaussian
     * @param s Pointer to the exploding star
     */
     void CCexplosion(Star *s) override;

};






#endif //SEVN_SUPERNOVA_H
