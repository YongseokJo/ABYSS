//
// Created by iorio on 7/22/22.
//

#include <Circularisation.h>
#include <binstar.h>
#include <star.h>
#include <functional>
//Abstract class
Circularisation *Circularisation::Instance(IO *_io) {

    auto it = GetStaticMap().find(_io->CIRC_mode);

    //TODO New way to automatically check the options, this is hardcoded now, but it should becomes
    //a method in the main property class to just be called here
    //If the name is in the list update the counter
    if (it != GetStaticMap().end()){
        GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    }
    //otherwise return an error
    else{

        std::string available_options;
        for (auto const& option : GetStaticMap()) available_options+=option.first+", ";

        svlog.critical("Option "+_io->CIRC_mode+" not available for process "+name()+". The available options are: \n"+available_options,
                       __FILE__,__LINE__,sevnstd::notimplemented_error(" "));
    }


    return (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

//disabled
int CircularisationDisabled::evolve(_UNUSED Binstar *binstar) {
    return EXIT_SUCCESS;
}


//Abstract class for standard circularisation
int StandardCircularisation::evolve(Binstar *binstar) {



    //Notice: in principle we  could add a condition on the Eccentricity
    //such as if Eccentricity=0, just return. However, it is safer to not add it
    //Indeed some process (e.g. the Hardening) can increase the eccentricity during the RLO
    //so a RLO can start on a circular orbit and becomes eccentric later. Using a check
    //on the eccentricity will produce a circularisation not at the onset of RLO, but slightly later.
    //We move the check of eccentricity=0 in the special evolve methods

    //Now check which condition we have to check
    bool _check_activation=check_activation(binstar);

    //First case, the activation conditions are satisfied ongoing, but check condition is false, so  we have to switch on the alarm
    if (!check_condition and _check_activation){
        check_condition = true; //Switch on check condition, because this is the onset of the activation
        special_evolution_alarm_switch_on(); //Switch on special evolution alarm
    }
    //the activation conditions are still satisfied but the  circularisation has been already applied, so just reset the alarm
    else if(check_condition and _check_activation){
        special_evolution_alarm_switch_off(); //Special evolution not needed
    }
    //Activation conditions are not satisfied, so just reset the check_condition
    else{
        check_condition=false;
    }

    return EXIT_SUCCESS;
}

int StandardCircularisation::special_evolve(Binstar *binstar) {

    if(binstar->getp(Eccentricity::ID)>0){
        circularise(binstar);
        //Disable binary check due to the impulsive variation of a and e
        binstar->disable_a_check=binstar->disable_e_check=true;
    }


    return 0; //0 Means system is not broken
}


std::string StandardCircularisation::log_message_circ(Binstar *binstar, double a_new, double e_new) const {
    std::string w = utilities::log_print("CIRC", binstar,
                                         binstar->getp(Semimajor::ID),binstar->getp(Eccentricity::ID),
                                         a_new, e_new);

    return w;
}


bool StandardCircularisation::check_activation(Binstar *binstar) const {
    //The conditions are: a RLO is ongoing and the Eccentricity>0
    return binstar->is_process_ongoing(RocheLobe::ID) and binstar->getp(Eccentricity::ID)>0;
}

///Circularisation Periastron

/**
 *  Circularisation at periastron: anew=a_old(1-e_old)
 * @param binstar pointer to the binary
 * @return EXIT_SUCCESS
 */
int CircularisationPeriastron::circularise(Binstar *binstar) {
    //Start from here
    double a_old=binstar->getp(Semimajor::ID);
    double e_old=binstar->getp(Eccentricity::ID);
    //Assume circularisation at the pericentre
    double a_new=a_old*(1- e_old);

    set(Semimajor::ID, a_new-a_old);
    set(Eccentricity::ID, -e_old);

    //Print to log
    std::string w = log_message_circ(binstar,a_new,0.);
    binstar->print_to_log(w);

    return EXIT_SUCCESS;
}

///Circularisation PeriastronFull
/**
 * Check the triggering condition of the circularisation
 * Check if the RLO triggering conditions (R>RLO) are satisfied at the periastron.
 * The RLO at periastron is estimated assuming a circular orbit with a=R_periastron
 * @param binstar Pointer to the binary
 * @return True if the RLO conditions (R>RLO) are satisfied at the periastron, false otherwise
 */
bool CircularisationPeriastronFull::check_activation(Binstar *binstar) const {
    //We consider the RL_pericentre= r_pericentre f(q) where f(q) is from the Eggleton equation  RL= a f(q)
    //therefore, since r_pericentre = a(1-e), RL_pericentre = a(1-e) / a * RL_eggleton = (1-e) * RL_eggleton

    //Necessary condition if Eccentricy is 0, just return false
    if (binstar->getp(Eccentricity::ID)==0){
        return false;
    }

    //Periastron radius
    double e1 = 1-binstar->getp(Eccentricity::ID);

    //Estimate the RL radii assuming a circular orbit at the pericentre
    double RL0_peri = e1*binstar->getp(RL0::ID);
    double RL1_peri = e1*binstar->getp(RL1::ID);

    Star* star0=binstar->getstar(0);
    Star* star1=binstar->getstar(1);

    if ( (star0->getp_0(Radius::ID)>RL0_peri and !star0->amiremnant()) or
         (star1->getp_0(Radius::ID)>RL1_peri and !star1->amiremnant())){
        return true;
    }

    return false;
}

///Circularisation Angmom

/**
 * Circularisation conserving the binary angular momentum: anew=a_old(1-e_old*e_old)
 * @param binstar pointer to the binary
 * @return EXIT_SUCCESS
 */
int CircularisationAngMom::circularise(Binstar *binstar) {
    //Start from here
    double a_old=binstar->getp(Semimajor::ID);
    double e_old=binstar->getp(Eccentricity::ID);
    //Assume circularisation conserving angular momentum
    double a_new=a_old*(1- e_old*e_old);

    set(Semimajor::ID, a_new-a_old);
    set(Eccentricity::ID, -e_old);

    //Print to log
    std::string w = log_message_circ(binstar,a_new,0.);
    binstar->print_to_log(w);

    return EXIT_SUCCESS;
}

/**
 * Circularisation conserving the semimajor axis: anew=a_old
 * @param binstar pointer to the binary
 * @return EXIT_SUCCESS
 */
int CircularisationSemimajor::circularise(Binstar *binstar) {
    //Start from here
    set(Semimajor::ID, 0.);
    set(Eccentricity::ID, -binstar->getp(Eccentricity::ID));

    //Print to log
    std::string w = log_message_circ(binstar,binstar->getp(Semimajor::ID),0.);
    binstar->print_to_log(w);

    return EXIT_SUCCESS;
}

