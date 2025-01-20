//
// Created by Lorenzo Valentini on 03/03/2023.
//

#include <Windaccretion.h>
#include <binstar.h>
#include <star.h>
#include <lookup_and_phases.h>

Windaccretion* Windaccretion::Instance(_UNUSED IO* _io) {
    auto it = GetStaticMap().find(_io->winds_mode);
    if (it != GetStaticMap().end()) {
        GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    }
    else {
        std::string available_options;
        for (auto const& option : GetStaticMap()) available_options += option.first + ", ";
        svlog.critical("Option " + _io->winds_mode + " not available for process " + name() + ". The available options are: \n" + available_options,
                       __FILE__, __LINE__, sevnstd::notimplemented_error(" "));
    }
    return  (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

//----------------  Hurley implementation  -----------------------------
double WindaccretionHurley::estimate_accreted_mass(double DM, Star *donor, Star *accretor, Binstar *binstar) const {
    double G = utilities::G;//grav const to be in Rsun, yr, Msun
    double e = binstar->getp(Eccentricity::ID);
    double a = binstar->getp(Semimajor::ID);
    double Mtot = binstar->getstar(0)->getp_0(Mass::ID) + binstar->getstar(1)->getp_0(Mass::ID);
    double vorb2 = G * Mtot / a;
    double _rsqrt_1_m_e2 = 1.0 / std::sqrt(1.0 - e * e);

    //Estimate the maximum accretion efficiency
    double f_MT_wind = estimate_maximum_accretion_efficiency(donor,accretor,binstar);

    //Properties at the beginning of the timestep
    double Mdonor = donor->getp_0(Mass::ID);
    double Radx = std::min(donor->getp_0(Radius::ID), binstar->Radx(donor->get_ID())); //If RLO is active (R>RL), so the RL (or the Rcore if larger) will be used as an effective radius
    //this property is at the end of the timestep (the secondary accretes mass after it has lost its own)
    double Maccretor = accretor->getp(Mass::ID);

    ///Estimate Maccreted
    double Maccreted;
    //TODO Let pureHe accrete?
    //Naked Helium or Naked CO does not accrete
    if (accretor->aminakedhelium() or accretor->aminakedco())
        Maccreted = 0;
        //Check if a is 0
    else if (a <= 0.) {
        svlog.warning("Semimajor axis is <0, winds accretion diverges. Wind accretion set to its maximum,", __FILE__, __LINE__);
        Maccreted = f_MT_wind * fabs(DM);
    }
    else {
        double vw2 = (2.0 * betaw * G * Mdonor) / Radx; //eq 9 of hurley 2002,
        double _1_p_vnorm2 = vorb2 / vw2 + 1.0;

        /*begin eq. 6 of Hurley 2002 */
        Maccreted = -_rsqrt_1_m_e2;
        Maccreted *= G * G * (Maccretor * Maccretor) / (vw2 * vw2);
        Maccreted *= alphaw / (2.0 * a * a);
        Maccreted *= 1.0 / std::sqrt(_1_p_vnorm2 * _1_p_vnorm2 * _1_p_vnorm2);
        Maccreted *= DM;
        Maccreted = std::min(Maccreted, f_MT_wind * fabs(DM)); //limit wind accretion efficiency at f_MT_wind%


        //Limit Mass accretion on degenerate accretors (WD,BH,NS) objects due to Eddington limit
        if (accretor->amIdegenerate_0()) {
            const double dt = binstar->getp(BTimestep::ID) * utilities::Myr_to_yr; //Elapsed time in this timestep in yr
            const double Max_Mass_eddington = dt * utilities::dMdt_Eddington_accretion(donor, accretor,
                                                                                       binstar->get_svpar_num("eddington_factor")); //Max Mass that can be accreted in the current timestpe due to the Eddington limit
            Maccreted = std::min(Maccreted, Max_Mass_eddington);

        }

    }

    return Maccreted;
}

int WindaccretionHurley::accrete_mass(Binstar* binstar) {

    double dt = binstar->getp(BTimestep::ID);

    donor = binstar->getstar(0);
    accretor = binstar->getstar(1);

    //Evolve the single stellar properties
    for (size_t i = 0; i < 2; i++) {

        ///Estimate Maccreted
        double DM = donor->getp(dMdt::ID) * dt; //mass loss of primary due to winds during timestep
        double Maccreted = estimate_accreted_mass(DM,donor,accretor,binstar);
        double Radx = std::min(donor->getp_0(Radius::ID), binstar->Radx(donor->get_ID())); //If RLO is active (R>RL), so the RL (or the Rcore if larger) will be used as an effective radius

        if (Maccreted < 0.0)
            svlog.critical("Maccreted cannot be [" + std::to_string(Maccreted) + "]: it's less than 0!", __FILE__, __LINE__);

        set_var(accretor->get_ID(), Mass::ID, Maccreted); //add the accreted mass to the accretor....
        set_var(accretor->get_ID(), dMcumul_binary::ID, Maccreted); //add the accreted mass to the accretor....
        set_var(accretor->get_ID(), dMcumulacc_wind::ID, Maccreted); //add the accreted mass to the accretor....


        //Add angular momentum due to the accreted mass
        double Laccreted = 0.6666666666666666 * Radx * Radx * Maccreted * donor->getp_0(OmegaSpin::ID);
        set_var(accretor->get_ID(), AngMomSpin::ID, muw * Laccreted);

        //Correct AngmomLost by the donor if Radx!=Radius
        if (Radx != donor->getp_0(Radius::ID)) {
            //1-Recover the Mass lost through winds
            double Mlost = donor->getp(Mass::ID) - donor->getp_0(Mass::ID); //Notice, it should be Mlost<=0
            //2-Recover the  Llost through sse (the one obtained assuming R=Radius)
            double dLlost_sse = 0.6666666666666666 * donor->getp_0(Radius::ID) * donor->getp_0(Radius::ID) * Mlost * donor->getp_0(OmegaSpin::ID);
            //3-Estimate the effective L angular momentum lost assuming Radius=radx
            double dLlost_effective = 0.6666666666666666 * Radx * Radx * Mlost * donor->getp_0(OmegaSpin::ID);
            //4-Set the change removing the dLlost_sse and adding dLlost_effective
            set_var(donor->get_ID(), AngMomSpin::ID, dLlost_effective - dLlost_sse);
            //Notice Mlost<0, therefore we are actually adding dLlost_sse (to remove the Llost in SSE) and
            //subtracting dLlost_effective.
        }

        ///Uncomment this if you want to take into account SpinUp due to mass acrrtion on NSs
        //Estimate change of Bman and OmegaRem if this a Neutron star
        //if (accretor->amiNS() and Maccreted!=0){
        //    handle_NS_massaccretion(binstar,accretor,Maccreted);
        //}

        utilities::swap_stars(donor, accretor);

    }

    return 0;

}

int WindaccretionHurley::evolve(Binstar* binstar) {

    //Evolve stellar properties
    accrete_mass(binstar);

    //Evolve the binary properties (Note: stars are swaped inside DA, DE)
    set(Semimajor::ID, DA(binstar, ID));
    set(Eccentricity::ID, DE(binstar, ID));

    svlog.debug("TOT DA " + utilities::n2s(get(Semimajor::ID), __FILE__, __LINE__), __FILE__, __LINE__);
    svlog.debug("TOT DE " + utilities::n2s(get(Eccentricity::ID), __FILE__, __LINE__), __FILE__, __LINE__);
    //utilities::wait(__FILE__,__LINE__);

    return 0;

}

double WindaccretionHurley::DA(Binstar* b, int procID) {
    double a = b->getp(Semimajor::ID);
    double e = b->getp(Eccentricity::ID);
    double dt = b->getp(BTimestep::ID);
    Star* donor = b->getstar(0);
    Star* accretor = b->getstar(1);
    double Mb = donor->getp_0(Mass::ID) + accretor->getp_0(Mass::ID);

    double Maccretor;
    double Mdonated;
    double Maccreted;
    double _DA = 0;

    //evolve 'a' for winds (from both stars)
    for (size_t i = 0; i < 2; i++) {

        Maccretor = accretor->getp_0(Mass::ID);  //Mass before the single stellar evolution
        Mdonated = donor->getp(dMdt::ID) * dt,
                Maccreted = b->getprocess(procID)->get_var(accretor->get_ID(), Mass::ID);

        _DA += a * (-Mdonated / Mb - ((2.0 - e * e) / Maccretor + (1.0 + e * e) / Mb) * Maccreted / (1.0 - e * e));
        utilities::swap_stars(donor, accretor);
    }

    return _DA;
}

double WindaccretionHurley::DE(Binstar* b, int procID) {

    double e = b->getp(Eccentricity::ID);
    Star* donor = b->getstar(0);
    Star* accretor = b->getstar(1);
    double Mb = donor->getp_0(Mass::ID) + accretor->getp_0(Mass::ID);

    double Maccretor;
    double Maccreted;
    double _DE = 0;

    //evolve 'a' for winds (from both stars)
    for (size_t i = 0; i < 2; i++) {

        Maccretor = accretor->getp_0(Mass::ID); //Mass before the single stellar evolution
        Maccreted = b->getprocess(procID)->get_var(accretor->get_ID(), Mass::ID);

        _DE += -e * Maccreted * (1.0 / Mb + 0.5 / Maccretor);
        utilities::swap_stars(donor, accretor);

    }

    return _DE;
}

double WindaccretionHurleyTT24mod::estimate_maximum_accretion_efficiency(Star *donor,
                                                                  Star *accretor,
                                                                  _UNUSED Binstar *binstar) const {

    double macc = accretor->getp_0(Mass::ID);
    double Mtot = macc + donor->getp_0(Mass::ID);

    return macc*macc/(Mtot*Mtot);
}

//------------------------------------------------------------------------------------------

int WindaccretionDisabled::evolve(_UNUSED Binstar* binstar) {
    return 0.0;
}