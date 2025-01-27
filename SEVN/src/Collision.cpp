//
// Created by iorio on 3/1/22.
//

#include <Collision.h>
#include <binstar.h>
#include <star.h>
#include <lookup_and_phases.h>

Kollision *Kollision::Instance(_UNUSED IO *_io) {

    auto it = GetStaticMap().find(_io->COLL_mode);

    //If the name is in the list update the counter
    if (it != GetStaticMap().end()){
        GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    }
        //otherwise return an error
    else{

        std::string available_options;
        for (auto const& option : GetStaticMap()) available_options+=option.first+", ";

        svlog.critical("Option "+_io->COLL_mode+" not available for process "+name()+". The available options are: \n"+available_options,
                       __FILE__,__LINE__,sevnstd::notimplemented_error(" "));
    }



    return  (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

std::string Kollision::log_message(Binstar *binstar) {



    const double &rlo_primary = binstar->getp(RL0::ID);
    const double &rlo_secondary = binstar->getp(RL1::ID);
    const Star *primary = binstar->getstar(0);
    const Star *secondary = binstar->getstar(1);



    std::string w = utilities::log_print("COLLISION",binstar,
                                         primary->get_ID(),primary->getp_0(Mass::ID),primary->getp_0(Radius::ID),primary->getp_0(Phase::ID),
                                         secondary->get_ID(),secondary->getp_0(Mass::ID),secondary->getp_0(Radius::ID),secondary->getp_0(Phase::ID),
                                         binstar->getp(Semimajor::ID),binstar->getp(Eccentricity::ID),rlo_primary,rlo_secondary);

    return w;
}

int KollisionDisabled::evolve(_UNUSED Binstar *binstar) {
    return Kollision::NO_COLLISION;
}

Kollision::collision_outcome KollisionHurley::outcome_collision(Binstar *b)  {

    const int &pbse = b->getstar(0)->get_bse_phase_0();
    const int &sbse = b->getstar(1)->get_bse_phase_0();

    //TODO, Here we should put that if the optimistic CE is on, we could trigger a CE also for phase=2
    //In addition is myabe better to move that check between CE and MIX directly in the CE process
    if ( ( (pbse>=3) and (pbse<=9) and (pbse!=7) ) or ( (sbse>=3) and (sbse<=9) and (sbse!=7) ) )
        return Kollision::COLLISION_CE;
    else
        return Kollision::COLLISION_MIX;

}


bool KollisionHurley::check_collision(Binstar *b) {

    const double &R1 = b->getstar(0)->getp_0(Radius::ID);
    const double &R2 = b->getstar(1)->getp_0(Radius::ID);
    double Rperiastron = b->getp(Semimajor::ID)*(1-b->getp(Eccentricity::ID));

    if ( R1+R2>Rperiastron )
        return true;

    return false;
}

int KollisionHurley::evolve(Binstar *b) {

    if (check_collision(b)){

        collision_outcome outcome=outcome_collision(b);

        if(outcome==Kollision::COLLISION_CE){
            b->comenv = true;
        } else if (outcome==Kollision::COLLISION_MIX) {
            b->mix = true;
        } else{
            svlog.critical("In KollisionHurley::evolve, check collision is true but outcome collision returns NO_COLLISION",__FILE__,__LINE__,
                           sevnstd::sanity_error(""));
        }
        //Set events and log message
        set_event((double)Lookup::EventsList::Collision);
        auto log_mess = Kollision::log_message(b);
        b->print_to_log(log_mess);
    }

    return EXIT_SUCCESS;

}

