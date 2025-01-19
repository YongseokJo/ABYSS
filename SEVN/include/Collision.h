//
// Created by iorio on 3/1/22.
//

/**
 * The class Kollision represents a prototype for the future refactoring of the processes.
 * The goal is to remove the orb_changer classes and let the Processes to handle all the BSE changes.
 *
 * For each new implementation we have to implement:
 *  A- A generic interface for the Process
 *  B- Specific inherited classes for each option-formalism of the process
 *
 *
 * ****** ADD THE PROCESS INTERFACE    **********
 *
 * In order to implement the interface, we follow the Kollision Process, for example assume
 * we want to add the process MagneticMassLoss,
 *
 * 1- Define the new class, public inherited from Process,
 * e.g.  class MagneticMassLoss : public Process
 *
 * 2- Implement the constructor
 *
 *     MagneticMassLoss(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            //First register to include this Property in static array
            Register(this, &ID, name());
        }
    }
 *
 * 3- Define the static ID and instance
    static size_t ID;
    static MagneticMassLoss _magneticmassloss;
 *  Notice: as default the name of the static instance is in lower case letter
 *
 *  3b- Add the ID and the instance in the static_main.h
 *  size_t MagneticMassLoss::ID;
 *  MagneticMassLoss MagneticMassLoss::_magneticmassloss;
 *
 *  4- Override the name
 *  inline std::string name() override { return "MagneticMassLoss"; }
 *
 *  5- Initialise the function Instance
 *  MagneticMassLoss *Instance(_UNUSED IO *_io) override;
 *  and define it in the .cpp file
 *
 *  MagneticMassLoss *MagneticMassLoss::Instance(_UNUSED IO *_io) {

    auto it = GetStaticMap().find(_io->MMloss);
    if (it != GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    return it == GetStaticMap().end() ? nullptr
                                      : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}
 * Notice the _io->MMloss is the name of the variable in io (that has to be added)
 * to choose the formalism of the process.
 *
 * 6- Define the virtual method instance returning a nullptr
 * virtual MagneticMassLoss *instance() {return nullptr;}
 * Notice: this method will be overridden in the specific process option
 *
 * 7- Add the protected method GetStaticMap (used to select the formalism at runtime)
 *     static std::map<std::string, MagneticMassLoss *> & GetStaticMap(){
        static std::map<std::string, MagneticMassLoss *> _locmap;
        return _locmap;
    }
 *
 * 8- Add the protected method GetUsed
 *     static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }
 *
 * 9- Add the protected method  Register_specific (used to register specific implementation of the processes)
 *     static void Register_specific(MagneticMassLoss *ptr) {
        GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }
 *
 *
 *  ******* ADD THE SPECIFIC PROCESS FORMALISM    **********
 *
 * To add a specific process formalism, it is possibile to follow the simple KollisionDisabled formalism
 * For example assume we want to add a MagneticMassLoss specific formalism called simple, where we just remove 10% of the mass of a star
 *
 * 1- Define the new class, public inherited from the Process interface,
 * e.g.  class MagneticMassLossSimple : public MagneticMassLoss
 *
 * 2- Build the specific constructor
 *     MagneticMassLossSimple(_UNUSED IO *_io= nullptr, bool reg = true) : MagneticMassLoss(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }
 * Notice: here we are using Register_specific to register the specific formalism,
 * Notice: the parent Process interface constructor is called with reg=false (very important)
 *
 * 3- Override the name method
 * inline std::string name() override { return "simple"; }
 * NOTICE: the returned name will the the one used to register the option and the one will be used
 * to select the given option at runtime, e.g.  -mmlossmode simple
 *
 * 4- Define the static instance
 * static MagneticMassLossSimple _magnetimasslosssimple;
 * Notice: In this case the ID is not needed since only the interface will
 * be used to build the vector of processes
 *
 * 4b- Add the static instance in static_main.h
 * MagneticMassLossSimple MagneticMassLossSimple::_magnetimasslosssimple;
 * Notice: In this way as soon the static_main.h is loaded the process formalism is loaded
 *
 * 5- Override the method instance
 *     MagneticMassLossSimple* instance() override {
        return (new MagneticMassLossSimple(nullptr,false));
    }
 * Notice if very important to use reg=false when calling the constructor
 *
 * 6- Finally, override the method evolve
 * int evolve(_UNUSED Binstar *binstar) override;
 *
 * Notice: this is the method where the effect of the process on the Binary and Stellar evolution are taken into account
 *
 * For example in our case we want to remove 10% of the mass of the star, a possible implementation is
 *
 * int evolve(_UNUSED Binstar *binstar){

      set_var(binstar->getstar(0), Mass::ID, 0.1*binstar->getstar(0)->getp_0(Mass::ID);
      set_var(binstar->getstar(1), Mass::ID, 0.1*binstar->getstar(1)->getp_0(Mass::ID);
  }
 * Notice: In the function set_var(starID, propertyID, variation) is called to set in the
 * proper table the variation of the Mass given by the process, the actual property changes happen
 * elsewhere (inside binary_evolve). DO NOT CHANGE directly the stellar and binary properties there.
 * If  you want to change the binary properties use set instead of set_var
 *
 */

#ifndef SEVN_COLLISION_H
#define SEVN_COLLISION_H

#include <IO.h>
#include <Processes.h>
#include <binstar.h>
#include <star.h>

//Forward declaration
class Star;
class Binstar;


class Kollision : public  Process{

public:

    Kollision(_UNUSED IO *_io= nullptr, bool reg = true) {
        if (reg) {
            //First register to include this Property in static array
            Register(this, &ID, name());
        }
    }

    static size_t ID;
    static Kollision _kollision;
    inline std::string name() override { return "Collision"; }

    Kollision *Instance(_UNUSED IO *_io) override;
    virtual Kollision *instance() {return nullptr;}

    bool is_process_ongoing() const override {return false;};



    static std::string log_message(_UNUSED Binstar *binstar);

    typedef unsigned int collision_outcome;
    static constexpr int NO_COLLISION=0;
    static constexpr int COLLISION_CE=1;
    static constexpr int COLLISION_MIX=2;

protected:

    ///Notice this stuff are here and used only for the process Kollision
    //this a prototype for the future refactoring of the whole Process classes
    static std::map<std::string, Kollision *> & GetStaticMap(){
        static std::map<std::string, Kollision *> _locmap;
        return _locmap;
    }
    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }
    static void Register_specific(Kollision *ptr) {
        GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }


};


class KollisionDisabled : public Kollision{

public:

    KollisionDisabled(_UNUSED IO *_io= nullptr, bool reg = true) : Kollision(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "disabled"; }
    static KollisionDisabled _kollisiondisabled;
    KollisionDisabled* instance() override {
        return (new KollisionDisabled(nullptr,false));
    }
    int evolve(_UNUSED Binstar *binstar) override;
};

class KollisionHurley : public Kollision{

public:

    KollisionHurley(_UNUSED IO *_io= nullptr, bool reg = true) : Kollision(nullptr, false) {
        if(reg) {
            //Second register to handle the different options
            Register_specific(this);
        }
    }

    inline std::string name() override { return "hurley"; }
    static KollisionHurley _kollisionhurley;
    KollisionHurley* instance() override {
        return (new KollisionHurley(nullptr,false));
    }
    int evolve(_UNUSED Binstar *binstar) override;

    static bool check_collision(Binstar *b);
    static collision_outcome outcome_collision(Binstar *b);

};





#endif //SEVN_COLLISION_H
