//
// Created by iorio on 3/30/22.
//

#include <Hardening.h>
#include <binstar.h>
#include <star.h>

int HardeningDisabled::evolve(_UNUSED Binstar *binstar) {
    return EXIT_SUCCESS;
}

Hardening *Hardening::Instance(IO *_io) {
    auto it = GetStaticMap().find(_io->HARD_mode);
    //If the name is in the list update the counter
    if (it != GetStaticMap().end()){
        GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    }
    //otherwise return an error
    else{

        std::string available_options;
        for (auto const& option : GetStaticMap()) available_options+=option.first+", ";

        svlog.critical("Option "+_io->HARD_mode+" not available for process "+name()+". The available options are: \n"+available_options,
                       __FILE__,__LINE__,sevnstd::notimplemented_error(" "));
    }


    return  (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

bool HardeningFastCluster::is_hard_binary(Binstar *binstar) {

    double sigma =  binstar->get_svpar_num("hard_sigma")*utilities::kms_to_RSunyr;
    const double &a = binstar->getp(Semimajor::ID);

    //THis is the mean kinetic energy in the external environment (e.g. stellar cluster) in Msun*Rsun^2/yr^2
    //It is initialised just the first time it is called. Since C+11 this is thread safe.
    //Since this value is constant and it just depends on general parameters it is thread safe anyway
    //The alternative is to make it a private class member
    const static double Ekin_extern_average  = 0.5*sigma*sigma*binstar->get_svpar_num("hard_mass_average");
    //This is the orbital energy (absolute value) of the binary in Msun*Rsun^2/yr^2, NOTICE: the massed are taken at t0 because the binary processes are
    //taken into account after the single evolution of the stars but we want to use the properties at the beginning of the timestep.
    //TODO In princple Ebinary could be a JIT Binary property of a method of the binary class.
    // 1- Do we need Ebinary elsewhere? In this case we can think to defined a method
    // 2- I am againist the inclusion in the Binary properties, it can be simply estimated and we cannot add all
    //  the derived quantities in the properties since each new property means slower updates (for loop over all the properties)
    //  and larger memory occupancy.
    const double Ebinary = utilities::G*binstar->getstar(0)->getp_0(Mass::ID)*binstar->getstar(1)->getp_0(Mass::ID)/(2*a);

    return Ebinary>Ekin_extern_average;
}

int HardeningFastCluster::evolve(Binstar *binstar) {

    //Now if the binary is hard (|Ebinary|>Ekin_average) integrate the hardening, otherwise just skip
    //TODO Should we put this check in the binary evolution to disable any hardening? I prefer to take this here because
    //  1- We don-t pollute the binary evolution method
    //  2- In principle other hardening options can deal with this check in a different way (or even not have the check at all)
    if (is_hard_binary(binstar)){
        const double &dt = binstar->getp(BTimestep::ID)*utilities::Myr_to_yr; //Dt in yr
        double sigma =  binstar->get_svpar_num("hard_sigma")*utilities::kms_to_RSunyr;
        const double &a = binstar->getp(Semimajor::ID);
        double rhoc = binstar->get_svpar_num("hard_rhoc")/(utilities::parsec_to_Rsun*utilities::parsec_to_Rsun*utilities::parsec_to_Rsun);

        double Ga_sigma = utilities::G*rhoc*a/sigma;

        double DA = -2.0*M_PI*binstar->get_svpar_num("hard_xi")*Ga_sigma*a*dt; //Eq. 11 Mapelli+21
        double DE = 2*M_PI*binstar->get_svpar_num("hard_xi")*binstar->get_svpar_num("hard_kappa")*Ga_sigma*dt; //Eq. 13 Mapelli+21

        set(Semimajor::ID,DA);
        set(Eccentricity::ID,DE);
    }
    else{
        set(Semimajor::ID,0.);
        set(Eccentricity::ID,0.);
    }



    return EXIT_SUCCESS;
}


