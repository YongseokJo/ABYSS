/*
 Created by iorio on 7/25/22.
Module for Pair instability models, it includes:
- Struct PISNreturn: a simple struct containing a double and a Lookup::SNExplosionType (enum) to handle some function returns
- Class PairInstability:  Base class for PairInstability models
- Class PIDisabled: Simple PairInstability model, just disable the pair instability correction
*/

#ifndef SEVN_PAIRINSTABILITY_H
#define SEVN_PAIRINSTABILITY_H

/**
* This header contains the class PairInstability.
* They are used to handle the pair instability supernova.
* This class is auxiliary to the Supernova class, each supernova instance has a PairInstability instance
* This class handles the PairInstability SN through two methods:
*   - apply_beforeSN: correction to apply before the proper explosion to set the preSN masses (e.g. remove the envelope)
*                     it returns a struct (utilities::MassContainer) containing the member Mass, MHE, MCO representing the preSN masses
*   - apply_afterSN: correction to apply after the SN explosion, it is usually a correction to the current remnant mass
*                     it returns a struct (PISNreturn) containing the final remnant mass Mremnant_after_pisn and the SN type
**/

/**
*
*  ******* ADD A NEW PAIR INSTABILITY MODEL  **********
* Assume we are implementing a new PairInstability class  PIIorio22, we suggest to
* create a new header file called piiorio22.h in include/star/procs/pairinstability/
*  and the related source file piiorio22.cpp in src/star/procs/pairinstability/
* Remember to add #include <piiorio22.h> in the static_main.h and in piiorio22.cpp
* Since it is possible that you will use the class star, we suggest to forward declare it in piiorio22.h
* class Star;
* and then #include <star.h> in piiorio22.cpp, remember also to add #include <pairinstability.h> in piiorio22.h
*
* 1- Define the new class, public inherited from  PairInstability
* e.g.  class PIIorio22 : public PairInstability
*
* 2- Build the specific constructor
*     PIIorio22(bool reg = true) : PairInstability(false){
*     if (reg) Register(true);
*     }
*     IMPORTANT: the Base constructor has to be called with reg=false
*
* 3- Override the name method
* inline std::string name() override { return "iorio22"; }
* NOTICE: the returned name will be the one used to register the option and the one that will be used
* to select the given option at runtime, e.g.  -sn_pairinstability iorio22
*
* 4- Define the static instance
* static PIIorio22 _piiorio22;
*
* 4b- Add the static instance in static_main.h
*      PIIorio22 PIIorio22::_piiorio22;
* Notice: In this way as soon the static_main.h is loaded the class is added to the available pair instability models
*
* 5- Override the method instance
*     PIIorio22* instance() override {
*            return (new PIIorio22(false));
*       }
* Notice if very important to use reg=false when calling the constructor
*
* 6-  Override the method apply_beforeSN (if needed)
*     utilities::MassContainer apply_beforeSN(_UNUSED Star* s) const;
*
 *     The method has to return a MassContainer object filled with the preSN values of Mass,MHE,MCO
*     (see the method documentation for further details
*     If not overridden, the default method just returns the MassContainer object with
*          -  MassContainer::Mass = s->getp(Mass::ID)
*          -  MassContainer::MHE = s->getp(MHE::ID)
*          -  MassContainer::MCO = s->getp(MCO::ID)
*
* 7-  Override the method apply_afterSN (if needed)
*     PISNreturn apply_afterSN(_UNUSED Star* s, _UNUSED double mremnant) const override;
*
*     The method has to return a PISNreturn object filled with the Mass of the remnant after the PISN correction Mremnant_after_pisn
*     and the SN type (see the method documentation for further details)
*     If not overridden, the default method just returns the PISNreturn object with
*           - Mremnant_after_pisn=mremnant
*           - SNType=Lookup::SNExplosionType::CoreCollapse
*
*  8: Add documentation
*      Try to document as much as possible.
*      In particular, use this template documentation for the class
*
*      **
*      * Brief model description (report also links and citations to the relevant papers)
*      * Modification with respect to the standard:
*      *   - apply_beforeSN: (if same as standard, just write no modifications)
*      *   - apply_afterSN:  (if same as standard, just write no modifications)
*      *
*      e.g.
*
*      **
*      * Simple prescription presented in Iorio+22 (in prep.) based on the Farmer+19 results
*      * see: Farmer+19: https://ui.adsabs.harvard.edu/abs/2019ApJ...887...53F/abstract
*      * Modification with respect to the standard:
*      *   - apply_beforeSN: modify the preSN mass
*      *         - MCO<38 Msun or MHE>135 Msun -> not modifications
*      *         - MCO>=38 Msun and (MHE<1=35 Msun)-> Mass=MHE=MCO
*      *   - apply_afterSN:
*      *         - MCO<38 Msun or MHE>135 Msun, no mass correction ->Core Collapse
*      *         - 38 Msun<=MCO<=60 Msun, no mass correction (but apply_beforeSN already set Mass=MHE=MCO) -> PPISN
*      *         - MCO>60 Msun and (MHE<135 Msun), set Mremnant to 0 -> PISN
*
**/



#include <map>
#include <utilities.h>
#include <vector>
#include <sevnlog.h>
#include <lookup_and_phases.h>

class Star;
using sevnstd::SevnLogging;

/**
 * Auxiliary struct for return of the main function apply
 */
struct PISNreturn{
    const double Mremnant_after_pisn; //Mass of the remnant after the correction
    const Lookup::SNExplosionType SNType; //Type of SNexplosion considering pairinstability
};


/**
 * Class to the handle the pulsation pair instability SN
 */
class PairInstability {
/***** For Developers that want to add new PairInstability models: ******/
///Remember that each new proper PairInstability model has to have a public  static class instance to define in static_main.h
/// e.g. static MyModel _mymodel; (in the class) and MyModel MyModel::_mymodel; in static_main.h
/// Methods that must  be customized in new models
public:
    /**Name of the PairInstability model, this will be use at runtime to select the  given model**/
    virtual inline std::string name() const {return "PairInstability";}

    /**
     * Create an instance of the class allocated on the heap
     * @return a pointer to a real instance of the class that is allocated on the  heap
     */
    virtual PairInstability *instance() = 0; //Pure virtual

    /// Methods that can  be customized in new models:
    ///NOTICE these are the two methods that actually estimate the pisn corrections.
    /// By default they do nothing (indeed they are not changed in the disabled model)
    /// You may want to change one or both of them
public:
    /**
     * Main function to call to apply the pulsation pair instability before the SN explosion
     * The function return a structure containing the modified preSN stellar mass(es) (total mass, MHE, MCO)
     * due to the pulsation pair instability.
     * This function should be called before estimating the final mass in the explosion method in the Supernova classes
     * The method of the base class just copy the Mass,MHE and MCO
     * NOTICE: the default method does nothing, just initialise the utilities::MassContainer with the current stellar mass properties
     * @param s Pointer to the exploding star
     * @return an Instance of the struct PPISNMass with only three members:
     *          - Mass_after_ppi, Mass of the star before the explosion due to the pulsation pair instability
     *          - MHE_after_ppi:  Mass of the HE core of the star before the explosion due to the pulsation pair instability
     *          - MCO_after_ppi:  Mass of the CO core of the star before the explosion due to the pulsation pair instability
     */
    virtual  utilities::MassContainer apply_beforeSN(_UNUSED Star* s) const;


    /**
     * Main function to call to apply the pair instability after the SN explosion
     * This function has to be called after the remnant mass has been defined.
     * It apply a correction to the remnant mass (if needed) and set the SNtype
     * The method of the base class just set Mremnant_after_pi=mremnant and SNType=CoreCollapse
     * NOTICE: the default method does nothing, just initialise the PISNreturn with:
     *          - Mremnant_after_pi=mremnant
     *          - SNType=Lookup::SNExplosionType::CoreCollapse
     * @param s Pointer to the exploding star to which we want to apply the correction
     * @param mremnant Mass of the remnant before the pair instability correction
     * @return an Instance of the struct PISNeturn with only two members:
     *          - Mremnant_after_pi, Mass of the remnant after the application of pi
     *          - SNType: SN explosion type an object of the enum Lookup::SNExplosionType (ElectronCapture, CoreCollapse, PPISN,PISN)
     */
    virtual  PISNreturn apply_afterSN(_UNUSED Star* s, _UNUSED double mremnant) const;



    /*********************************************************************************************************/
    /***** STOP: DO NOT CHANGE THE REST OF THE CODE IF YOU ARE NOT FULLY AWARE OF WHAT  YOU ARE DOING ******/
    /*********************************************************************************************************/

public:

    //Constructor - Destructor
    explicit PairInstability(_UNUSED bool reg=true){};
    virtual ~PairInstability() = default;

    //Create instance
    /**
     * Generic method to call an instance based on name
     * @param name name of the option
     * @return a pointer to a real instance of the class that is allocated on the  heap
     */
    static PairInstability *Instance(std::string const &name);

    //Remove cp and mv, we don't need to copy or mv, therefore to avoid unexpected behaviour we just delete them
    PairInstability (const PairInstability&) = delete; //copy constructor
    PairInstability& operator= (const PairInstability&) = delete; //copy assignment
    PairInstability (PairInstability&& other) = delete; //mv constructor
    PairInstability& operator= (PairInstability&& other) = delete; //mv assignement


protected:
    /**
     * Register function, to be called inside the option constructor to register the given option
     */
    static void Register(PairInstability *ptr);
    SevnLogging svlog;

private:

    ///Option handling
    //TODO We are using this Register-StaticMap formalism a lot, maybe we have to create a base template class containing
    //all these methods (GetStaticMap,Register,GetUsed) to have only one definition of these functions
    static std::map<std::string, PairInstability *> &GetStaticMap() {
        //First time call the function create the map, since c++11 this is thread safe
        static std::map<std::string, PairInstability *> _locmap;
        return _locmap;
    }
    ///TODO we are never using this method, we should remove it
    static std::vector<int> &GetUsed() {
        static std::vector<int> _used;
        return _used;
    }
};

/**
 * disabled option
 * no correction applied
 */
class PIDisabled : public  PairInstability{
public:
    //Ctor
    explicit PIDisabled(bool reg=true);

    //Static instance
    static PIDisabled _pidisabled;

    //name
    inline std::string name() const override { return "disabled";}

    //return an heap allocated instance
    PIDisabled* instance() override {
        return (new PIDisabled(false));
    }

};




#endif //SEVN_PAIRINSTABILITY_H
