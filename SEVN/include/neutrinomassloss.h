/*
 Created by iorio on 7/25/22.
Module for Neutrino mass loss models, it includes:
- Class NeutrinoMassLoss:  Base class for NeutrinoMassLoss models
- Class NMLDisabled: Simple NeutrinoMassLoss model, just disable the neutrino mass loss correction
*/

#ifndef SEVN_NEUTRINOMASSLOSS_H
#define SEVN_NEUTRINOMASSLOSS_H

/**
 * This header contains the class NeutrinoMassLoss.
 * They are used to handle the neutrino mass loss during a SN explosion
 * This class is auxiliary to the Supernova class, each supernova instance has a NeutrinoMassLoss instance
 * The class has only one method that is called in the Supernova class :
 *   - double apply(const double mremnant, Star *s): it get the (pre-correction) remnant mass mremnant and the pointer to the exploding star
 *   as input and return the mass of the remnant after the neutrino mass loss
**/


/**
*
*  ******* ADD A NEW NeutrinoMassLoss MODEL  **********
* Assume we are implementing a new NeutrinoMassLoss class  NMLLattimer89, we suggest to
* create a new header file called nmllattimer89.h in include/star/procs/neutrinomassloss/
*  and the related source file nmllattimer89.cpp in src/star/procs/neutrinomassloss/
* Remember to add #include <nmllattimer89.h> in the static_main.h and in nmllattimer89.cpp
* Since it is possible that you will use the class star, we suggest to forward declare it in nmllattimer89.h
* class Star;
* and then #include <star.h> in nmllattimer89.cpp, remember also to add #include <neutrinomassloss.h> in nmllattimer89.h
*
* 1- Define the new class, public inherited from  PairInstability
* e.g.  class NMLLattimer89 : public NeutrinoMassLoss
*
* 2- Build the specific constructor
*     NMLLattimer89(bool reg = true) : PairInstability(false){
*     if (reg) Register(true);
*     }
*     IMPORTANT: the Base constructor has to be called with reg=false
*
* 3- Override the name method
* inline std::string name() override { return "lattimer89"; }
* NOTICE: the returned name will be the one used to register the option and the one that will be used
* to select the given option at runtime, e.g.  -sn_neutrinomloss lattimer89
*
* 4- Define the static instance
* static NMLLattimer89 _nmllattimer89;
*
* 4b- Add the static instance in static_main.h
*      NMLLattimer89 NMLLattimer89::_nmllattimer89;
* Notice: In this way as soon the static_main.h is loaded the class is added to the available pair instability models
*
* 5- Override the method instance
*     NMLLattimer89* instance() override {
*            return (new NMLLattimer89(false));
*       }
* Notice if very important to use reg=false when calling the constructor
*
* 6-  Override the method apply
*     double apply(_UNUSED const double mremnant, _UNUSED Star *s) const;
*     it get the (pre-correction) remnant mass mremnant and the pointer to the exploding star
*     as input and return the mass of the remnant after the neutrino mass loss
*
*  8: Add documentation
*      Try to document as much as possible.
**/

#include <sevnlog.h>
#include <map>
#include <vector>
#define _UNUSED __attribute__ ((unused))

class Star;

using sevnstd::SevnLogging;

/**
 * Base class for NeutrinoMassLoss
 */
class NeutrinoMassLoss{
/***** For Developers that want to add new NeutrinoMassLoss models: ******/
///Remember that each new proper NeutrinoMassLoss model has to have a public  static class instance to define in static_main.h
/// e.g. static MyModel _mymodel; (in the class) and MyModel MyModel::_mymodel; in static_main.h
/// Methods that must  be overridded in new models
public:
    /**Name of the NeutrinoMassLoss model, this will be use at runtime to select the  given model**/
    virtual inline std::string name() const {return "NeutrinoMassLoss";}

    /**
     * Create an instance of the class allocated on the heap
     * @return a pointer to a real instance of the class that is allocated on the  heap
     */
    virtual NeutrinoMassLoss *instance() = 0; //Pure virtual

    /**
     * Estimate the neutrino mass loss correction
     * @param mremnant Mass of the remnant before the pair instability correction
     * @param s  Pointer to the exploding star to which we want to apply the correction
     * @return Mass after the neutrino mass loss correction
     */
    virtual double apply(_UNUSED const double mremnant, _UNUSED Star *s) const = 0; //Pure virtual


    /*********************************************************************************************************/
    /***** STOP: DO NOT CHANGE THE REST OF THE CODE IF YOU ARE NOT FULLY AWARE OF WHAT  YOU ARE DOING ******/
    /*********************************************************************************************************/

public:
    //Constructor - Destructor
    explicit NeutrinoMassLoss(_UNUSED bool reg=true){};
    virtual ~NeutrinoMassLoss() = default;

    //Create instance
    /**
     * Generic method to call an instance based on name
     * @param name name of the option
     * @return a pointer to a real instance of the class that is allocated on the  heap
     */
    static NeutrinoMassLoss *Instance(std::string const &name);

    //Remove cp and mv, we don't need to copy or mv, therefore to avoid unexpected behaviour we just delete them
    NeutrinoMassLoss (const NeutrinoMassLoss&) = delete; //copy constructor
    NeutrinoMassLoss& operator= (const NeutrinoMassLoss&) = delete; //copy assignment
    NeutrinoMassLoss (NeutrinoMassLoss&& other) = delete; //mv constructor
    NeutrinoMassLoss& operator= (NeutrinoMassLoss&& other) = delete; //mv assignement

protected:
    /**
     * Register function, to be called inside the option constructor to register the given option
     */
    static void Register(NeutrinoMassLoss *ptr);
    SevnLogging svlog;

private:

    ///Option handling
    //TODO We are using this Register-StaticMap formalism a lot, maybe we have to create a base template class containing
    //all these methods (GetStaticMap,Register,GetUsed) to have only one definition of these functions
    static std::map<std::string, NeutrinoMassLoss *> &GetStaticMap() {
        //First time call the function create the map, since c++11 this is thread safe
        static std::map<std::string, NeutrinoMassLoss *> _locmap;
        return _locmap;
    }
    ///TODO we are never using this method, we should remove it
    static std::vector<int> &GetUsed() {
        static std::vector<int> _used;
        return _used;
    }
};

/**
 * option disabled
 *
 * Do nothing, do not apply any correction for neutrino mass loss
 */
class NMLDisabled : public NeutrinoMassLoss{
public:
    //Ctor
    NMLDisabled(bool reg=true);

    //name
    inline std::string name() const override { return "disabled";}

    //Static instance
    static NMLDisabled _nmldisabled;

    //return an heap allocated instance
    NMLDisabled* instance() override {
        return (new NMLDisabled(false));
    }

    /**
     * Estimate the neutrino mass loss correction, in this disabled version do nothing, just return mremnant
     * @param mremnant Mass of the remnant before the pair instability correction
     * @param s  Pointer to the exploding star to which we want to apply the correction
     * @return mremnant
     */
     double apply(const double mremnant, _UNUSED Star *s) const override {return mremnant;}
};



#endif //SEVN_NEUTRINOMASSLOSS_H
