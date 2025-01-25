//
// Created by spera on 13/02/19.
//

#include <Processes.h>
#include <Orbit.h>
#include <binstar.h>
#include <supernova.h>
#include <utilities.h>
#include <remnant.h>

//TODO We should refactored all the Process section, eliminating the Orbit.h using a method similar
//to the one we used for the SN for example. We have already

//TODO NOVA should be a process not something that happens inside RLO because it should take into account all the mass transfered

Process::~Process() {
        delete orb_change;
        orb_change=nullptr;
}


double MaccretionProcess::NS_DBmag_accretion(_UNUSED Binstar *b, Star *s, double DM) const {
    ///Checks
    if (!s->amiNS())
        svlog.critical("Called NS_DBmag_accretion but the input star is not a NS",__FILE__,__LINE__,sevnstd::sanity_error());
    else if (DM<0)
        svlog.critical("Called NS_DBmag_accretion but the input DM is negative (it should be positive, it is supposed to be an accretion event)",__FILE__,__LINE__,sevnstd::sanity_error());

    //Recast the remnant pointer to access info about NS
    const NSrem *ns = static_cast<const NSrem*>(s->get_staremnant());

    const double &B0 = s->getp_0(Bmag::ID); //Our B0 is the one at the beginning of the step in Gauss
    const double &Bmin = ns->get_Bmin(); //Bmin in Gauss
    //DMNS in the Equation is the input parameter DM in Msun
    const double DMd = s->get_svpar_num("ns_magnetic_mscale"); //Magnetic field decay mass-scale in Msun

    double Bnew; //CE Estimate this value using Eq. 6
    Bnew = (B0-Bmin)*std::exp(-DM/DMd) + Bmin;
    double DBmag = Bnew-B0;//We have to return a difference a DBmag.


    return DBmag;
}

double MaccretionProcess::NS_DOmegaRem_accretion(Binstar *b, Star *s, double DM) const {
    ///Checks
    if (!s->amiNS())
        svlog.critical("Called NS_DBmag_accretion but the input star is not a NS",__FILE__,__LINE__,sevnstd::sanity_error());
    else if (DM<0)
        svlog.critical("Called NS_DBmag_accretion but the input DM is negative (it should be positive, it is supposed to be an accretion event)",__FILE__,__LINE__,sevnstd::sanity_error());

    //useful constants
    const double yr = 3.1557600e7; //yr in s
    const double rsun = 6.95700e10; //rsun in cm
    const double msun = 1.98892e33; //msun in g
    const double G_cgs = utilities::G*rsun*rsun*rsun/(msun*yr*yr); //cm^3 s^-2 g^-1
    //Recast the remnant pointer to access info about NS
    //const NSrem *ns = static_cast<const NSrem*>(s->get_staremnant());
    const double &Omega0 = s->getp_0(OmegaRem::ID); //Our Omega0 is the one at the beginning of the step in 1/s
    const double Inertia = s->getp_0(Inertia::ID)*msun*rsun*rsun; //Inertia in g cm^2

    //const double &B0 = s->getp_0(Bmag::ID); //B0 at the beginning of the step
    //const double &dt = b->getp(BTimestep::ID)*1e6*yr; //Timestep in s
    //const double R = s->getp_0(Radius::ID)*rsun; //Radius of the NS in cm
    const double M = s->getp_0(Mass::ID)*msun; //Mass of the NS in g
    double DM_g = DM*msun; //Mass difference in g


    //double factor1 = std::pow(1.0/(G_cgs*8), 1.0/7.0);
    //double factor2 = std::pow(pow(R, 6.0)*dt/(DM_g*pow(M, 0.5)), 2.0/7.0);
    //double factor3 = std::pow(B0, 4.0/7.0);
    //double RA=0.5*factor1*factor2*factor3; //Magnetic radius in cm
    //Use the new functon in utilites
    double RA = 0.5*utilities::R_Alfven(s,DM/b->getp(BTimestep::ID),true)*utilities::Rsun_cgs; //Magnetic radius in cm




    double Vdiff = std::pow(G_cgs*M/(RA*RA*RA), 0.5)-Omega0; //in s^-1
    double DOmega=Vdiff*RA*RA*DM_g/Inertia; //Variation in the angular frequency of the NS in s^-1

    //std::cout<<" "<<DOmega<<" "<<Omega0+DOmega<<" "<<factor3<<__FILE__<<" "<<__LINE__<<std::endl;
    //utilities::hardwait(" CHECk omega",DM,RA,Vdiff,DOmega,Omega0,__FILE__,__LINE__);


    //MNS is the input parameter DM in MSUN

    //NOTICE: Eq. 6 in Chattopadhyay relate dOmega / dt to dM / dt, however SEVN works with discrete timestep and what we are interested to is just DOmega,
    //therefore integrating we obtain just
    //DOmega  propto DM that is the input parameter of this function

    return DOmega;
}

void MaccretionProcess::handle_NS_massaccretion(Binstar *b, Star *s, double DM) {

    if (b->get_svpar_bool("ns_detailed_accretion")){
        double DBmag = NS_DBmag_accretion(b,s,DM);
        double DOmegaRem = NS_DOmegaRem_accretion(b,s,DM);
        set_var(s->get_ID(),Bmag::ID,DBmag);
        set_var(s->get_ID(),OmegaRem::ID,DOmegaRem);
    }
}

/*******************************************
********     SN KICKS    *************
*******************************************/
int SNKicks::special_evolve(Binstar *binstar) {

    //NB return 0 if the system is not destroyed, 1 if the system is destroyed

    ///If disabled just exit
    if (orb_change->name()==Lookup::snkmap_name.at(Lookup::SNKickMode::_SNKickdisabled))
        return 0;

    ///Initiliase orbit changer (here DA and DE are filled with values)
    orb_change->init(binstar);
    auto *orb_change_casted = dynamic_cast<Orbital_change_SNKicks*>(orb_change);

    //Check if the systes has been destroyed by SNkick
    bool check_broken_da = ( binstar->getp(Semimajor::ID) + orb_change->DA(binstar,ID) )<=0; //a_fin<=0
    bool check_broken_de = ( binstar->getp(Eccentricity::ID) + orb_change->DE(binstar,ID) )>=1; //ecc_fin>=1


    double a_fin=0.,e_fin=0., cosnu=0., vcom=0.;
    std::string w;
    //If the system is broken just set the binary to the broken state
    if(check_broken_da or check_broken_de){
        binstar->broken = true;
        a_fin=std::nan("");
        e_fin=std::nan("");
        cosnu=std::nan("");
        vcom=std::nan("");
        svlog.pdebug("System destroyed by SNkick(s)",__FILE__,__LINE__);
        w = SNKicks::log_message(binstar,a_fin,e_fin,cosnu,vcom);
        binstar->print_to_log(w);
        set_event(Lookup::EventsList::SNBroken);
        return 1;
    }
    //If DA or DE !=0 set the DA, DE and disable the Timestep a and e check
    else if( orb_change->DA(binstar,ID)!=0 or orb_change->DE(binstar,ID)!=0){
        a_fin=binstar->getp(Semimajor::ID)+orb_change->DA(binstar, ID);
        e_fin=binstar->getp(Eccentricity::ID)+orb_change->DE(binstar, ID);
        cosnu=orb_change_casted->get_cos_nu();
        vcom=orb_change_casted->get_vcom();
        set(Semimajor::ID, orb_change->DA(binstar, ID));
        set(Eccentricity::ID, orb_change->DE(binstar, ID));
        //We don't want to repeat the timestep due to the a and e change after a sn explosion
        //since it is an impulsive variation
        binstar->disable_e_check=true;
        binstar->disable_a_check=true;
    }
    w = SNKicks::log_message(binstar,a_fin,e_fin,cosnu,vcom);
    binstar->print_to_log(w);



    return 0;


}

std::string SNKicks::log_message(Binstar *binstar,double a_fin,double e_fin, double cos_nu, double vcom) {

    std::string star_sn, star_other;

    if (binstar->getstar(0)->getp_0(Phase::ID)<Lookup::Remnant and binstar->getstar(0)->getp(Phase::ID)==Lookup::Remnant){
        star_sn = utilities::log_star_info(binstar->getstar(0),true);
        star_other= utilities::log_star_info(binstar->getstar(1),true);
    }
    else{
        star_sn = utilities::log_star_info(binstar->getstar(1),true);
        star_other= utilities::log_star_info(binstar->getstar(0),true);
    }

    std::string w = utilities::log_print("BSN", binstar,star_sn,star_other,
            binstar->getp(Semimajor::ID),binstar->getp(Eccentricity::ID),
            a_fin,e_fin,cos_nu,vcom);

    return w;

}


/*******************************************
********     GW RADIATION    *************
*******************************************/

int GWrad::evolve(Binstar *binstar){

    //When the condition are_not_merged is false, it means that the semimajor is smaller than the larger Rs, so the two
    //stars are merging. In that case, for now we stop the GWrad effect.
    ///Check condition to stop GWrad
    Star* s1 = binstar->getstar(0);
    Star* s2 = binstar->getstar(1);
    //Set the maximum GWtime (in Hubble time) for which the GWrad are enavled
    double gw_time_check = binstar->get_svpar_num("gw_tshold");
    //The last stable circular orbit (LSCO) is 3*Rs, we compare a with the sum of the two LSCO
    double LSCO = 3*(s1->getp(Rs::ID)+s2->getp(Rs::ID));
    bool are_merged = binstar->getp(Semimajor::ID)<=LSCO; //If a<3Rs flags the merger
    bool check_tgw_scale = binstar->getp(GWtime::ID)<gw_time_check*utilities::tH; //If tscale_GW<factor*THubble
    bool is_not_disabled=orb_change->name()!=Lookup::gwmap_name.at(Lookup::GWMode::_GWdisabled); //Check if they are disabled



    ///If are_mergerd id true a<=3Rs, flags the mix and exit.
    if (are_merged){
        binstar->mix=true;
        set_event((double)Lookup::EventsList::GW_Merger);
        utilities::wait("Merging from GW",binstar->mix,binstar->getstar(0)->getp(Mass::ID),__FILE__,__LINE__);
        return 0;
    }

    bool check_bco = !binstar->get_svpar_bool("gw_onlybco") or (s1->amiremnant() and s2->amiremnant());

    ///DO the evolution if: i) GW are enabled, ii) GWtime<f*THubble iii) if gw_onlybco is true check if this is a bco (binary compact object).
    if (is_not_disabled and check_tgw_scale and check_bco){
        //Evolve the binary properties
        set(Semimajor::ID, orb_change->DA(binstar, ID));
        set(Eccentricity::ID, orb_change->DE(binstar, ID));
        if (first_time_GW){
            set_event((double)Lookup::EventsList::GWBegin);
            first_time_GW=false;
        }
    }


    //svlog.debug("TOT DA "+utilities::n2s(get(Semimajor::ID),__FILE__,__LINE__),__FILE__,__LINE__);
    //svlog.debug("TOT DE "+utilities::n2s(get(Eccentricity::ID),__FILE__,__LINE__),__FILE__,__LINE__);



    return 0 ;


}



/*******************************************
********         TIDES       *************
*******************************************/

int Tides::evolve(Binstar *binstar) {

	//Evolve the binary properties (Note: stars are swaped inside DA, DE)

	/// Init stuff for tides
    orb_change->init(binstar);

	set(Semimajor::ID,orb_change->DA(binstar, ID));
	set(Eccentricity::ID,orb_change->DE(binstar, ID));

    set_var(binstar->getstar(0)->get_ID(),AngMomSpin::ID,orb_change->DAngMomSpin(binstar,ID,0));
    set_var(binstar->getstar(1)->get_ID(),AngMomSpin::ID,orb_change->DAngMomSpin(binstar,ID,1));

    //std::cout<<" Spin1 " <<orb_change->DOSpin(binstar,ID,0)<<" "<<binstar->getstar(0)->getp(OmegaSpin::ID)<<std::endl;
    //std::cout<<" Spin2 " <<orb_change->DOSpin(binstar,ID,1)<<" "<<binstar->getstar(1)->getp(OmegaSpin::ID)<<std::endl;

	//svlog.debug("DA "+utilities::n2s(orb_change->DA(binstar, ID),__FILE__,__LINE__),__FILE__,__LINE__);
	//svlog.debug("DE "+utilities::n2s(orb_change->DE(binstar, ID),__FILE__,__LINE__),__FILE__,__LINE__);
	//utilities::wait(__FILE__,__LINE__);
	return 0;

}


/*******************************************
******** ROCHE LOBE OVERFLOW *************
*******************************************/

int RocheLobe::evolve(Binstar *binstar) {

    if (orb_change->name()==Lookup::rlmap_name.at(Lookup::RLMode::_RLdisabled))
        return EXIT_SUCCESS;



    orb_change->init( binstar);
    auto *orb_change_casted = dynamic_cast<Orbital_change_RL*>(orb_change);

    //If the stars are not filling the Roche Lobe just return
    if (orb_change->is_process_ongoing()){

        ///precheck for swapped cases
        //This pre-check is needed because sometime in the same phase a RLO ends and a new one starts
        //his can happens when, for example, the donor becomes a nakedHe so that R<RL in the same phase
        //in which the other starts trigger a R>RL.
        //Here we check if the  swapping condition is verified and we set the various flags
        check_and_handle_swapping_condition(binstar);

        if(!RLO_last_step){
            //Reset dMcumul because we just start a new RLO
            reset_dMcumul_RLO_in_stars(binstar);
            //Reset dMcumul
            reset_DM_global();
            std::string w=log_message_start(binstar,orb_change_casted->get_mt_value(),orb_change_casted->get_mt_tshold());
            binstar->print_to_log(w);
            set_event((double)Lookup::EventsList::RLOBegin);
        }


        Star* s     = binstar->getstar(0); //Star with the first star as donor
        Star* s_other  = binstar->getstar(1); //Start with the second star as accretor




        //If just one (or none) of the two stars overflow the Roche Lobe
        for (size_t i=0; i<2; i++){

            double DM=orb_change->DM(binstar, ID, s->get_ID());

            ///Check peculiar accretion before to start
            //If the star accreting mass is a nakedHe  and it is accreting Helium, update all the MHE
            if (DM>0 and s->aminakedhelium() and s_other->whatamidonating_0()==Material::He) {
                set_var(s->get_ID(), MHE::ID, DM);
            }
            //If the star accreting mass is a nakedCO or a nakedHe and the accreted mass is not He
            //just raise a not implemented error. We do not allow this atm.
            else if (DM>0 and (s->aminakedhelium() or s->aminakedco())){
                svlog.critical("Accretion on nakedCo or non He accretion on a nakedHe is not "
                               "implemented",__FILE__,__LINE__,sevnstd::notimplemented_error());
            }

            set_var(s->get_ID(), Mass::ID, DM);
            set_var(s->get_ID(), dMcumul_binary::ID, DM);
            set_var(s->get_ID(), dMcumul_RLO::ID, DM);
            set_var(s->get_ID(),AngMomSpin::ID, orb_change->DAngMomSpin(binstar,ID,s->get_ID()));

            //Estimate change of Bman and OmegaRem if this an accreting Neutron star
            if (s->amiNS() and DM>0){
                handle_NS_massaccretion(binstar,s,DM);
            }

            if (s->get_ID()==0) DM_global_0+=DM;
            else DM_global_1+=DM;

            //If the star losing mass is a Naked Helium or a Naked CO, we should remove the mass also from the HE and CO core
            //TODO Myabe this should be handeld directly in star
            ///Donor Check
            //Remember if it is nakedco it is also nakedhe, the reverse is not true, so
            //we check nakedco before nakedhelium
            if (DM<0 and s->aminakedco()){
                set_var(s->get_ID(), MHE::ID, DM);
                set_var(s->get_ID(), MCO::ID, DM);
            }
            else if(DM<0 and s->aminakedhelium()){
                set_var(s->get_ID(), MHE::ID, DM);
            }


            utilities::swap_stars(s, s_other); //Swap the stars and check the RLO for the other star.
        }




        ///ORBITAL CHANGES
        //DA
        set(Semimajor::ID, orb_change->DA(binstar, ID));
        set(Eccentricity::ID, orb_change->DE(binstar, ID));

        RLO_last_step=true;
    }
    else if (orb_change_casted->is_colliding()){
        std::string w=log_message_end(binstar);
        binstar->print_to_log(w);
        set_event((double)Lookup::EventsList::Collision);
    }
    else {
        if(RLO_last_step){
            std::string w=log_message_end(binstar);
            binstar->print_to_log(w);
            //Set the event only if it has a lower priority. All the events with higher priority ends the RLO,
            //so in this case this event is useless
            if (binstar->getp(BEvent::ID)<Lookup::EventsList::RLOEnd) set_event((double)Lookup::EventsList::RLOEnd);
        }
        RLO_last_step=false;
        id_donor=-1; //reset id_donor
        RLOswapped=false;
    }

    return EXIT_SUCCESS;

}

//Roche Lobe Special evolve for shallowd stars
int RocheLobe::special_evolve(Binstar *binstar){

    if (orb_change->name()==Lookup::rlmap_name.at(Lookup::RLMode::_RLdisabled))
        return EXIT_SUCCESS;

    std::string w;


    Star* donor;
    Star* accretor;
    //Select donor, accretor
    if (binstar->is_swallowed[0]){
        donor = binstar->getstar(0);
        accretor = binstar->getstar(1);
    }
    else if (binstar->is_swallowed[1]){
        donor = binstar->getstar(1);
        accretor = binstar->getstar(0);
    }
    else{
        svlog.critical("Called RocheLobe special evolve (dynamic swallowing), but neither Stars have been flagged as swallowed.",
                       __FILE__,__LINE__,sevnstd::rl_error(""));
    }

    //Check1: the swallowing should be called just for donor without an already developed core
    if (donor->aminakedhelium() and donor->getp(MCO::ID)>0 and !donor->aminakedco()){
        svlog.critical("Called RocheLobe special evolve (dynamic swallowing), but the donor is a nakedHelium with a CO core, this is not expected",
                       __FILE__,__LINE__,sevnstd::rl_error(""));
    }
    //TODO
    //17/02/2023: GI changed this line:
    //else if (donor->getp(MHE::ID)>0 and !donor->aminakedco() and !donor->aminakedhelium()){
    //with the new one below.
    //At this point the idea is that just stars without a core could trigger a swallowing,
    //In the current implementation this happens fro stars with BSE phase <=1. However
    //in the SEVN-BSE conversion we consider phase 1 also stars we a tiny amount of He in the core at the beginning of the HE growth
    //(to avoid to have a common envelope with a star with a very tiny core (<1e-5))
    //Myabe we have just to thing a more robust implementation of the RLO-CE, so that if a CE is triggered
    //and the core is  tiny we just merge everything (as in the case of swallowing).
    //08/03/2023: GI changed this line:
    //else if ( donor->getp(Phase::ID)>1 and !donor->aminakedco() and !donor->aminakedhelium() and !donor->amiremnant()){
    //with the following to take into account the situations in which the He core is tiny and we continue to consider this star as a MS
    else if ( donor->getp(Phase::ID)>1 and donor->getp(PhaseBSE::ID)>1 and !donor->aminakedco() and !donor->aminakedhelium() and !donor->amiremnant()){
        svlog.critical("Called RocheLobe special evolve (dynamic swallowing), but the donor is a H-star with a (developed) He core, this is not expected",
                       __FILE__,__LINE__,sevnstd::rl_error(""));
    }
    //Ok for remnant as donor, they never have a core

    //Call the special evolve to set DM_accreted
    orb_change->speciale_evolve(binstar);
    double DM_accreted  = orb_change->DM(binstar,ID,accretor->get_ID());
    double DM_swallowed = orb_change->DM(binstar,ID,donor->get_ID());


    //Check-2, check that DM_swallowed=Mass_donor
    if(!utilities::areEqual(fabs(DM_swallowed),donor->getp(Mass::ID))){
        svlog.critical("In the RocheLobe special evolve (dynamic swallowing), the donated mass is not equal to the total mass of the donor, this is not expected",
                       __FILE__,__LINE__,sevnstd::rl_error(""));
    }

    utilities::sn_explosion fate = binstar->check_accretion_on_compact(donor->get_ID(), accretor->get_ID(), DM_accreted);
    if (fate==utilities::SNIA_EXPLODE){
        w=RocheLobe::log_message_swallowed(binstar,DM_accreted);
        //Binstar set to broken and accretor set to empy inside check_accretion_on_compact, so just donor to empy
        donor->set_empty();
    }
    //TODO ATM, we do not allow a nakedCO to accrete anything during a swallowing,
    //it is not clear what happens here, it is has to be considered as a compact object.
    //But even in this case if it accretes mass should it jump to a new H or pureHe tracks?
    //For now, for the sake of simplicity, we just assume that any accreted mass will be just removed quickly due to wind mass loss
    else if(accretor->aminakedco()){
        DM_accreted=0;
        w=RocheLobe::log_message_swallowed(binstar,DM_accreted);
        //Do not accret anything, just set stuff and raise a Warning
        donor->set_empty();
        binstar->set_broken();
        binstar->set_onesurvived();
        svlog.warning("In the RocheLobe special evolve (dynamic swallowing), the accretor is nakedCO and the accreted mass is >0. ATM SEVN is not handling this situation, the accreted mass will be set to 0.",
                      __FILE__,__LINE__);
    }
    else{

        //Save  the log message before changing the binary properties
        w=RocheLobe::log_message_swallowed(binstar,DM_accreted);


        if (DM_accreted>0){
            //Add all or part of the swallowed material to the accretor
            accretor->update_from_binary(Mass::ID, DM_accreted);
            accretor->update_from_binary(dMcumul_binary::ID, DM_accreted);

            //If the donor is a naked Helium, the accretor should have an Helium core or be a remnant to add the new material
            if (donor->aminakedhelium() and !accretor->amiremnant() and accretor->getp(MHE::ID)>0){
                accretor->update_from_binary(MHE::ID, DM_accreted); //Add the new material in the core
            } else if(donor->aminakedhelium() and !accretor->amiremnant()){
                svlog.critical("In the RocheLobe special evolve (dynamic swallowing), the swallowed star is a nakedHe, but the accretor has not an He-core and it is not a remnant, this is not expected",
                               __FILE__,__LINE__,sevnstd::rl_error(""));
            }
            //If it is a remnant it is not needed to update the Hecore, just the Mass is stored

            //If the donor is nakedCO, the accretor should have a CO core or be a remnant to add the new material
            if (donor->aminakedco() and !accretor->amiremnant() and accretor->getp(MCO::ID)>0){
                accretor->update_from_binary(MCO::ID, DM_accreted); //Add the new material in the CO core
                accretor->update_from_binary(MHE::ID, DM_accreted); //Add the new material in the HE core
            }else if(donor->aminakedco() and !accretor->amiremnant()){
                svlog.critical("In the RocheLobe special evolve (dynamic swallowing), the swallowed star is a nakedco, but the accretor has not an CO-core and it is not a remnant, this is not expected",
                               __FILE__,__LINE__,sevnstd::rl_error(""));
            }
            //If it is a remnant it is not needed to update the COcore and the Hecore, just the Mass is stored


            //If Hydrogen is added to a pureHe, we have to jump
            if (accretor->aminakedhelium() and !donor->aminakedhelium()){
                accretor->jump_to_normal_tracks();
            }
                //If the core is changed (H-stars that accreted from a nakedHe) jump to new tracks as it is a merger
            else if(!accretor->aminakedhelium() and donor->aminakedhelium()){
                accretor->find_new_track_after_merger();
            }
        }



        //Set donor to empty and binstar to broken
        donor->set_empty();
        binstar->set_broken();
        binstar->set_onesurvived();
        //Print log message
        binstar->print_to_log(w);
    }



    //RLO can set also RLOBegin, so if the dynamic swallowing happens at the onset of RLO, the event RLOB is already set
    //to avoid to overwrite it if it is already set we set the composite event RLOB+Swallowing
    if (get_event()==Lookup::EventsList::RLOBegin){
        set_event((double)Lookup::EventsList::RLOB_Swallowed);
    } else{
        set_event((double)Lookup::EventsList::Swallowed);
    }



    ///LAST STUFF
    //In any case reset swallowed
    binstar->is_swallowed[0]=binstar->is_swallowed[1]=false; //Reset swallowed attribute
    return EXIT_SUCCESS;

}

std::string RocheLobe::_log_message(Binstar *binstar, bool swap_donor_accretor) {

    std::string star_donor,star_accretor;

    if (binstar->getstar(0)->getp_0(Radius::ID)>=binstar->getp(RL0::ID)){
        star_donor=utilities::log_star_info(binstar->getstar(0),true) + ":"
                + utilities::n2s(binstar->getstar(0)->getp_0(Radius::ID),__FILE__,__LINE__) + ":"
                + utilities::n2s(binstar->getp(RL0::ID),__FILE__,__LINE__);
        star_accretor=utilities::log_star_info(binstar->getstar(1),true) + ":"
                   + utilities::n2s(binstar->getstar(1)->getp_0(Radius::ID),__FILE__,__LINE__) + ":"
                   + utilities::n2s(binstar->getp(RL1::ID),__FILE__,__LINE__);
    } else{
        star_accretor=utilities::log_star_info(binstar->getstar(0),true) + ":"
                   + utilities::n2s(binstar->getstar(0)->getp_0(Radius::ID),__FILE__,__LINE__) + ":"
                   + utilities::n2s(binstar->getp(RL0::ID),__FILE__,__LINE__);
        star_donor=utilities::log_star_info(binstar->getstar(1),true) + ":"
                      + utilities::n2s(binstar->getstar(1)->getp_0(Radius::ID),__FILE__,__LINE__) + ":"
                      + utilities::n2s(binstar->getp(RL1::ID),__FILE__,__LINE__);
    }

    if (swap_donor_accretor){
        return star_accretor + ":" +  star_donor;
    }

    return star_donor + ":" +  star_accretor;
}


std::string RocheLobe::log_message_start(Binstar *binstar,double q,double qcrit,bool swap_donor_accretor) {

    std::string other_msg = _log_message(binstar,swap_donor_accretor);

    std::string w = utilities::log_print("RLO_BEGIN", binstar,other_msg, q,  qcrit, binstar->getp(Semimajor::ID),binstar->getp(Eccentricity::ID));

    return w;
}

std::string RocheLobe::log_message_end(Binstar *binstar, bool swap_donor_accretor) {

    std::string other_msg = _log_message(binstar,swap_donor_accretor);

    std::string w = utilities::log_print("RLO_END", binstar,other_msg,std::min(DM_global_0,DM_global_1),std::max(DM_global_0,DM_global_1),binstar->getp(Semimajor::ID),binstar->getp(Eccentricity::ID));

    return w;
}

std::string RocheLobe::log_message_swallowed(Binstar *binstar, double DM_accreted) {

    std::string swallowed_stars;
    std::string other;

    Star *swallowed_star = binstar->is_swallowed[0] ? binstar->getstar(0) : binstar->getstar(1);
    Star *other_star = binstar->is_swallowed[0] ? binstar->getstar(1) : binstar->getstar(0);

    return log_message_swallowed(binstar, swallowed_star, other_star,DM_accreted);
}

std::string RocheLobe::log_message_swallowed(Binstar *binstar, Star *swallowed, Star *other, double DM_accreted) {

    std::string swallowed_star = utilities::log_star_info(swallowed);
    std::string other_star = utilities::log_star_info(other);
    std::string w = utilities::log_print("SWALLOWED", binstar,swallowed_star,other_star,DM_accreted);

    return w;
}

bool RocheLobe::check_and_handle_swapping_condition(Binstar *binstar) {

    Star* s     = binstar->getstar(0); //Star with the first star as donor
    Star* s_other  = binstar->getstar(1); //Start with the second star as accretor

    //This check is needed because sometime in the same phase a RLO ends and a new one starts
    //his can happens when, for example, the donor becomes a nakedHe so that R<RL in the same phase
    //in which the other starts trigger a R>RL.
    //Here we check if the  swapping condition is verified and we set the various flags
    for (size_t i=0; i<2; i++){

        double DM=orb_change->DM(binstar, ID, s->get_ID());
        const int id_star = int(s->get_ID());

        //Set the current id_donor if id_donor is not set (=-1)
        if (id_donor==-1 and DM<0){
            id_donor=id_star;
            break; //We can exit
        }
            //check swapped
        else if(DM<0 and  id_donor!=-1 and id_donor!=id_star){
            id_donor=id_star;
            RLOswapped=true;
        }

        if (RLOswapped){
            //Message: the old RLO ends
            //NOTICE, here we swap the donor and accretor, because
            //the end refers to the previous RLO phase before the donor/accretor swapping
            //so the swap is needed to be consistent with the previous RLOB begins message
            std::string w=log_message_end(binstar,true);
            binstar->print_to_log(w);
            RLO_last_step=false;
            //reset_DM_global(); //Reset DM
            //Message: a new RLO starts
            //binstar->print_to_log(w);
            RLOswapped=false;
            //Reset dMcumul because we just start a new RLO
            reset_dMcumul_RLO_in_stars(binstar);
            break;
        }

        utilities::swap_stars(s, s_other); //Swap the stars and check the RLO for the other star.
    }

    return RLOswapped;
}

int RocheLobe::reset_dMcumul_RLO_in_stars(Binstar *binstar) {

    for (int i=0; i<2; i++){
        Star* s=binstar->getstar(i);
        s->update_from_binary(int(dMcumul_RLO::ID),-s->getp(dMcumul_RLO::ID),binstar);
    }

    return  EXIT_SUCCESS;
}

/*******************************************
******** MIX *************
*******************************************/

int Mix::special_evolve(Binstar *binstar){

    if (binstar->mix){
        utilities::wait("Hey I have to mix",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);

        Star *donor = binstar->getstar(0); //Donor is the star that will set to empty
        Star *accretor = binstar->getstar(1); //Accretor is the star that will remain as results of the mix

        ///Choose the star that remains and the one that is set to empty
        //First handle the general case: the accretor is the more evolved star,
        // or the more compact (massive) remnant if both are remnant (handled in get_id_more_evolved_star).
        if(binstar->get_id_more_evolved_star()==0)
            //Now the star 0 becomes the accretor and the star 1 the donor
            utilities::swap_stars(donor,accretor);
        //Now handle a special case: If the accretor is a naked helium and the donor not and they have the same innermost core
        //we swap so that the accretor will be the normal star.  This is done because jumping to a normal
        //track from a pureHe is more dangerous since if there are problems on finding a good match the star
        //does not jump and the code raises an error.  In case of a normal track, if we don't find a match
        //we can just continue following the original track.
        //GI 18/02/21: bug fix, to really swap the donor have to have at least a Helium core (the CO core is handled inside).
        //We check >1E-3 instead of >0 to avoid to select the accretor star as the one that is just starting to grow its HE core.
        //GI 02/09/22: We further extend the bug fix to disable the donor swab if the H-star is growing its core from 0, i.e. it is the TMS phase
        //
        //if (accretor->aminakedhelium() and (!donor->aminakedhelium() and donor->getp(MHE::ID)>1E-3)){
        if (accretor->aminakedhelium() and (!donor->aminakedhelium() and
            donor->getp(Phase::ID)>Lookup::TerminalMainSequence and
                donor->getp(Phase::ID)!=Lookup::TerminalCoreHeBurning)){
            unsigned int id_inner_core_accretor = accretor->getp(MCO::ID)!=0 ? MCO::ID : MHE::ID;
            unsigned int id_inner_core_donor =  donor->getp(MCO::ID)!=0 ? MCO::ID : MHE::ID;
            if (id_inner_core_accretor==id_inner_core_donor)
                utilities::swap_stars(donor,accretor);
        }



        ///Let's mix
        //Case 1 - Mix between not remnant stars
        if (!accretor->amiremnant() and !donor->amiremnant()){

            //Set new maximumCO and minmum HE
            accretor->set_MCO_max_aftermerge(donor);
            accretor->set_MHE_min_aftermerge(donor);

            //Mass from the donor
            double DM_new     = donor->getp(Mass::ID);
            double DMHE_new   = donor->getp(MHE::ID);
            double DMCO_new   = donor->getp(MCO::ID);
            double MCO_old    = accretor->getp(MCO::ID);
            double MHE_old    = accretor->getp(MHE::ID);

            //Update masses of the accretor
            accretor->update_from_binary(Mass::ID, DM_new);
            accretor->update_from_binary(dMcumul_binary::ID, DM_new);
            accretor->update_from_binary(MHE::ID, DMHE_new);
            accretor->update_from_binary(MCO::ID, DMCO_new);

            ///Jump to a new track
            //Sanity check on MCO
            if (MCO_old==0 and DMCO_new>0){
                svlog.critical("This is an embarrassing situation, in Mix the donor has a CO core and the accretor not",
                               __FILE__,__LINE__,sevnstd::ce_error());
            }
            // Sanity check on MHE
            else if (MHE_old==0 and DMHE_new>0){
                svlog.critical("This is an embarrassing situation, in Mix the donor has a HE core and the accretor not",
                               __FILE__,__LINE__,sevnstd::ce_error());
            }
            // The accretor is a nakedHelium and the donor not, so we have to trigger a jump to a normal track
            // There is no possibility that the donor has a CO core and the naked helium not because we already check this
            // at the beginning
            else if (accretor->aminakedhelium() and !donor->aminakedhelium()){
                accretor->jump_to_normal_tracks();
            }
            else{
                accretor->find_new_track_after_merger();
            }
            binstar->set_onesurvived(); //Onesurveived is true
        }
        //Case 2 - Mix between WDs
        //TODO Notice that here mixing with a WD is treated as mixing with a NS/BH, in Hurley we have a more complicate outcomes
        else if (accretor->amiWD() and donor->amiWD()) {

            //Check if we have two HeWD
            bool check_double_HeWD = donor->getp(RemnantType::ID)==Lookup::Remnants::HeWD && accretor->getp(RemnantType::ID)==Lookup::Remnants::HeWD;

            //Case 2a - SNI explosion
            if (check_double_HeWD){
                accretor->explode_as_SNI();
                if (binstar->onesurvived) binstar->set_onesurvived(); //Why this condition? I forgot, we have to check
            }
            else{
                //Update the Mass
                accretor->update_from_binary(Mass::ID, donor->getp(Mass::ID));
                binstar->set_onesurvived();
            }

        }
        //Case 3 -  NS, BH or WD mixing with a NS/BH
        else if (accretor->amiCompact() and donor->amiremnant()){
            //Update the Mass
            accretor->update_from_binary(Mass::ID, donor->getp(Mass::ID));
            binstar->set_onesurvived();
        }
        //Case 4 - A not remnant donor over a NS, BH, WD
        //-> Do not accreate nothing, just set the donor to empty, see below
        //From commonenvelope::main_stellar_collision in common_envelope.cpp in SEVN1


        ///LOGPRINT AND EVENT SET
        std::string w =Mix::log_message(binstar,accretor,donor);
        binstar->print_to_log(w);
        //utilities::hardwait(" Event before ",get_event());
        set_event((double)Lookup::EventsList::Merger);
        //utilities::hardwait(" Event After ",get_event());

        ///LAST STUFF
        //In any case set the binary to broken and put the donor to empty
        //donor->set_empty_in_bse(binstar);
        donor->set_empty();
        binstar->mix=false; //Reset mix attribute
        binstar->set_broken(); //Broken is true,


    }


    return EXIT_SUCCESS;

}

/**
 * Log message print
 * @param binstar  Pointer to binstar
 * @param accretor  Pointer the accretor (the star surviving the merger)
 * @param donor Pointer to the donor (the star set to empty after the merger)
 * @return A string message from utilities::log_print containing
 * IDaccretor:Maccretor:MHEaccretor:MCOaccretor:Phaseaccretor:RemTypeaccretor:IDdonor:Mdono:MHEdonor:MCOdonor:Phasedonor:RemTypedonor
 * @Notice: This function assume that whe it is called, the properties of the accretor have been already changed
 * while the donor has not been set to empty yet.
 */
std::string Mix::log_message(Binstar *binstar, Star *accretor, Star *donor) {


    std::string info_accretor_old = utilities::log_star_info(accretor,true);
    std::string info_donor = utilities::log_star_info(donor,false);

    std::string w = utilities::log_print("MERGER",binstar,
                                         info_accretor_old,
                                         accretor->getp(Radius::ID),
                                         info_donor,
                                         donor->getp(Radius::ID),
                                         accretor->getp(Mass::ID),
                                         binstar->getp(Semimajor::ID),
                                         binstar->getp(Eccentricity::ID));

    return w;
}

/*******************************************
******** CE *************
*******************************************/
//TODO: Note as we wrote it the special evolve should be after the update of the single stellar property and before the update of the binary property
int CommonEnvelope::special_evolve(Binstar *binstar){


    //Check common envelope
    if(!binstar->comenv)
        return EXIT_FAILURE;

    //Check if disabled
    if (orb_change->name()==Lookup::cemap_name.at(Lookup::CEMode::_CEdisabled)){
        binstar->comenv=false; //Reset comenv attribute
        return EXIT_FAILURE;
    }



    utilities::wait("CE started",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);
    //Initialise
    init(binstar);


    std::string logmess_first = CommonEnvelope::log_mess(binstar,primary,secondary);

    //Check
    if (!isinisialised)
        svlog.critical("CE process has not been initialised",__FILE__,__LINE__,sevnstd::ce_error());

    double afin_initial = a_fin;

    //Decide the final fate of the system
    double a_ini = binstar->getp(Semimajor::ID);
    double Mass_primary_ini = primary->getp(Mass::ID);
    double Mass_secondary_ini = secondary->getp(Mass::ID);
    double Radius_primary_ini = primary->getp(Radius::ID);
    double Radius_secondary_ini = secondary->getp(Radius::ID);
    int Phase_old_primary = primary->getp(Phase::ID);
    int Phase_old_secondary = secondary->getp(Phase::ID);
    int Phase_bse_old_secondary = secondary->get_bse_phase();
    double Lambda_secondary_old=secondary->getp(Lambda::ID);
    double Lambda_primary_old=primary->getp(Lambda::ID);
    double Ebind_old_primary=primary->getp(Ebind::ID);
    double Ebind_old_secondary=secondary->getp(Ebind::ID);


    bool coalesce_done=false;
    std::string merger_mess;
    //If the stars are still filling the Roche Lobe start a coalescence
    if (R_core_primary>RL_primary_final || R_core_secondary>RL_secondary_final){
        utilities::wait("Coalesce",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);
        merger_mess=main_coalesce(binstar);
        coalesce_done=true;
    }
    //otherwise just loose the envelope
    else {
        utilities::wait("Lose the envelope",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);
        lose_the_envelope(binstar);
        //Update Orbital parameters
        double DA=a_fin-binstar->getp(Semimajor::ID);
        double DE=ecc_fin-binstar->getp(Eccentricity::ID);
        set(Semimajor::ID, DA);
        set(Eccentricity::ID, DE);

        //TODO This part below takes into account the mass accretion on NS
        //during the CE. This is just a prototype and it is now commented by default.
        //In oder to make it a standard feature we should:
        // 1- Enable this mass accretion for all the compact remnants (maybe also stars)
        // 2- Create runtime parameters to tune this mass accretion and a boolean flag to enable/disabe this option
        // 3- Atm we handle the mass accretion changing the spin and the Bmag but we do not really accrete the mass on the NS
        /*
        double DM_accrete_compact_remnant_min{0.04};
        double DM_accrete_compact_remnant_max{0.1};
        for (unsigned int i=0; i<2; i++){
            if (binstar->getstar(i)->amiNS()){
                Star *ns = binstar->getstar(i);
                double DM_on_NS=(std::uniform_real_distribution<>(DM_accrete_compact_remnant_min,
                                                                  DM_accrete_compact_remnant_max))(utilities::mtrand);
                //When we switch to a complete mass accretion (including change of Omega) remove the comment here
                //and comment the next two lines
                //handle_NS_massaccretion(binstar,ns,DM_on_NS);
                double DBmag = NS_DBmag_accretion(binstar,ns,DM_on_NS);
                set_var(ns->get_ID(),Bmag::ID,DBmag);
            }
        } */


        //TODO here we need something better before to star to pollute the class of this disable etc etc
        //for example we can move this check directly the BinaryProperty
        //anche check if with check_disabled that return true or false and automatically reset to false if it was true
        binstar->disable_a_check=true;
        binstar->disable_e_check=true;
        binstar->disable_stellar_rotation_check=true;
    }


    ///LOG PRINT
    //TODO Here we have to use afin instead of getp because we call special evolve we don't evolve the Binary, maybe we have to rethink this
    logmess_first += ":" + utilities::n2s(a_ini,__FILE__,__LINE__) +":"+ utilities::n2s(a_fin,__FILE__,__LINE__); //Add semimajor initial and final
    logmess_first += ":"+utilities::n2s(coalesce_done,__FILE__,__LINE__);//Add Coalesce done
    binstar->print_to_log(logmess_first);

    if (coalesce_done) binstar->print_to_log(merger_mess);

    //Set Event
    const Lookup::EventsList& event = coalesce_done ? Lookup::EventsList::CE_Merger : Lookup::EventsList::CE;
    set_event((double)event);






    svlog.pdebug("Binary ID:", binstar->get_ID(), "Name:", binstar->get_name(),
                 "started CE at Worldtime",binstar->getp(BWorldtime::ID),
                 "\n ID primary", primary->get_ID(),
                 "\n ain", binstar->getp(Semimajor::ID), "afin",afin_initial, "afin_effective", a_fin,
                 "\n ein", binstar->getp(Eccentricity::ID), "efin", ecc_fin,
                 "\n Mass primary before CE", Mass_primary_ini,
                 "\n Mass secondary before CE", Mass_secondary_ini,
                 "\n Mass primary after CE", primary->getp(Mass::ID),
                 "\n Mass secondary after CE", secondary->getp(Mass::ID),
                 "\n Core primary", M_core_primary,
                 "\n Core secondary", M_core_secondary,
                 "\n Envelope primary", M_env_primary,
                 "\n Envelope secondary", M_env_secondary,
                 "\n Radius primary",  Radius_primary_ini,
                 "\n Radius secondary", Radius_secondary_ini,
                 "\n Core Radius primary", R_core_primary,
                 "\n Core Radius secondary", R_core_secondary,
                 "\n Roche Lobe primary  initial", binstar->getp(RL0::ID),
                 "\n Roche Lobe secondary final", binstar->getp(RL1::ID),
                 "\n Roche Lobe primary  final", RL_primary_final,
                 "\n Roche Lobe secondary final", RL_secondary_final,
                 "\n Ebind ini", Ebind_ini,
                 "\n Ebind old primary", Ebind_old_primary,
                 "\n Ebind old secondary", Ebind_old_secondary,
                 "\n Lambda primary",   primary->getp(Lambda::ID),
                 "\n Lambda primary old",   Lambda_primary_old,
                 "\n Lambda secondary", secondary->getp(Lambda::ID),
                 "\n Lambda secondary old", Lambda_secondary_old,
                 "\n Phase primary",   primary->getp(Phase::ID),
                 "\n Phase primary BSE",   primary->get_bse_phase(),
                 "\n Phase primary old",   Phase_old_primary,
                 "\n Naked Helium",   primary->aminakedhelium(), primary->getp(Mass::ID)-primary->getp(MHE::ID),
                 "\n Phase secondary",   secondary->getp(Phase::ID),
                 "\n Phase secondary BSE",   secondary->get_bse_phase(),
                 "\n Phase secondary old",   Phase_old_secondary,
                 "\n Phase secondary BSE old",   Phase_bse_old_secondary,
                 "\n Naked Helium",   secondary->aminakedhelium(),
                 "\n FILE",__FILE__,
                 "\n LINE",__LINE__);
    utilities::wait(__FILE__,__LINE__);

    binstar->comenv=false; //Reset comenv attribute


    return EXIT_SUCCESS;

}

int CommonEnvelope::init(Binstar *binstar) {

    double a = binstar->getp(Semimajor::ID);
    double ecc = binstar->getp(Eccentricity::ID);

    //Initialise primary-secondary
    whoisprimary(binstar);

    //NOTICE: all energies are divided by G
    //TODO Check if we really need M_env
    //Parameters
    //lambda = binstar->get_svpar_num("ce_lambda");
    alpha = binstar->get_svpar_num("ce_alpha");
    //Mass and radii
    M_env_primary    = primary->aminakedhelium() ? primary->getp(Mass::ID) - primary->getp(MCO::ID) : primary->getp(Mass::ID) - primary->getp(MHE::ID);
    M_core_primary   = primary->getp(Mass::ID) - M_env_primary;
    R_core_primary   = primary->aminakedhelium() ? primary->getp(RCO::ID) : primary->getp(RHE::ID);

    //if s2 does not have a REAL core or is a remnant (M_env_secondary==secondary->getp(Mass::ID)),
    //the entire star is considered as a core and the envelope is 0
    //s2 could have mcore = 0.0. In this case we consider mcore = M_real (the entire star s2 is a core).
    if (secondary->haveienvelope()){
        M_env_secondary  = secondary->aminakedhelium() ? secondary->getp(Mass::ID) - secondary->getp(MCO::ID) : secondary->getp(Mass::ID) - secondary->getp(MHE::ID);
        R_core_secondary = secondary->aminakedhelium() ? secondary->getp(RCO::ID) : secondary->getp(RHE::ID);
    }
    else{
        M_env_secondary  = 0;
        R_core_secondary = secondary->getp(Radius::ID);
    }
    M_core_secondary = secondary->getp(Mass::ID) - M_env_secondary;


    //Initial binding energies of envelope. This is the reservoir of energy that can be used in this phase
    Ebind_ini = (primary->getp(Ebind::ID) + secondary->getp(Ebind::ID)); //Initial binding energy of the envelopes (Eq. 69 Hurley02)
    //Orbital energy of the two bodies (cores)
    Eorb_ini = - (M_core_primary * M_core_secondary)/(2.0*a); //Eq. 70 Hurley02
    //Final Orbital Energy: the initial one plus a fraction of the initial reservoir of binding energy
    Eorb_fin = Eorb_ini + Ebind_ini/alpha;

    //Semimajor axis evolution
    a_fin = - (M_core_primary * M_core_secondary)/(2.0*Eorb_fin); //New semimajor axes

    //Eccentricity evolution (following BSE/MOBSE)
    //Check if any eccentricity remains in the orbit by first using
    //energy to circularise the orbit before removing angular momentum.
    //note this should not be done in case of CE SN ...
    //TODO See above, what about SN?
    //E of a circular orbit at same a
    double  Ecirc = Eorb_ini/(1.0 - ecc*ecc);
    //See if there is room for any residual eccentricity
    if(Eorb_fin > Ecirc) //Remember E is negative, so Eorb_fin>Ecirc means |Eorb_fin|<|Ecirc| therefore Eorb_fin/Ecirc<1
        ecc_fin = sqrt(1 - Eorb_fin/Ecirc);
    else
        ecc_fin = 0;
    //TODO I am not sure about this assumption, Ecirc is the Energy you have to loose to have  a circular orbit at fixed a, here a is evolving too.

    //Estimate new Roche Lobe Radius
    RL_primary_final   = utilities::roche_lobe_Eg(M_core_primary,  M_core_secondary, a_fin);
    RL_secondary_final = utilities::roche_lobe_Eg(M_core_secondary, M_core_primary, a_fin);

    isinisialised = true;

    return EXIT_SUCCESS;
}

int CommonEnvelope::whoisprimary(Binstar *binstar){


    Star *star1 = binstar->getstar(0);
    Star *star2 = binstar->getstar(1);




    //Note that this special evolve happens outside the normal cycle of the Processe

    double mcore1 = star1->aminakedhelium() ? star1->getp(MCO::ID) : star1->getp(MHE::ID);
    double mcore2 = star2->aminakedhelium() ? star2->getp(MCO::ID) : star2->getp(MHE::ID);

    //Notice, We check Radius of stars at time t0 because the process that can trigger the CE are estimated using the T0 property
    bool check_core1 =mcore1 !=0.0;
    bool check_core2 =mcore2 !=0.0;


    //GI 13/03: Before this if-elseif-else part included also checks on the R>RL.
    //This was a logic check needed because before the COllision could be triggered only during the RLO
    ///Set primary and secondary
    //If both have core: the primary is the more evolved star
    if (check_core1 and check_core2){

        unsigned int _id_more_evolved = binstar->get_id_more_evolved_star();
        if (_id_more_evolved==0){
            primary = star1; secondary = star2;
        }
        else{
            primary = star2; secondary = star1;
        }
    }
    //Now check if star1 has a core and if filling the RL, it is the primary
    else if (check_core1){
        primary   = star1; secondary = star2;
    }
    //Now check if star2 has a core and if filling the RL, it is the primary
    else if (check_core2){
        primary   = star2; secondary = star1;
    }
    //this system should not have started the Common Envelope.
    else{
        svlog.critical("Entered in the Common Envelope process, but neither of the two stars has a core"
                       "\n M0: "+ utilities::n2s(star1->getp(Mass::ID),__FILE__,__LINE__) + " MHE0: " + utilities::n2s(star1->getp(MHE::ID),__FILE__,__LINE__) +
                       " MCO0: " + utilities::n2s(star1->getp(MCO::ID),__FILE__,__LINE__) +
                       " R0: "+ utilities::n2s(star1->getp(Radius::ID),__FILE__,__LINE__) + " RHE0: "+ utilities::n2s(star1->getp(RHE::ID),__FILE__,__LINE__) +
                       " RCO0: "+ utilities::n2s(star1->getp(RCO::ID),__FILE__,__LINE__) +
                       " RLO: " + utilities::n2s(binstar->getp(RL0::ID),__FILE__,__LINE__) +
                       "\n M1: "+ utilities::n2s(star2->getp(Mass::ID),__FILE__,__LINE__)  + " MHE1: " + utilities::n2s(star2->getp(MHE::ID),__FILE__,__LINE__) +
                       " MCO1: " + utilities::n2s(star2->getp(MCO::ID),__FILE__,__LINE__) +
                       " R1: "+ utilities::n2s(star2->getp(Radius::ID),__FILE__,__LINE__) + " Rc1: "+ utilities::n2s(star2->getp(RHE::ID),__FILE__,__LINE__) +
                       " RCO1: "+ utilities::n2s(star2->getp(RCO::ID),__FILE__,__LINE__) +
                       " RL1: " + utilities::n2s(binstar->getp(RL1::ID),__FILE__,__LINE__),
                       __FILE__, __LINE__, sevnstd::ce_error());
    }


    return EXIT_SUCCESS;


}

std::string CommonEnvelope::main_coalesce(Binstar *binstar){
    //Check
    if (!isinisialised)
        svlog.critical("CE process has not been initialised",__FILE__,__LINE__,sevnstd::ce_error());

    ///START
    //In this situation, the two cores merge before to reach the a_fin needed to expel all the envelope
    //therefore the outcome will be a star with the core that is the mix of the two cores and an envelope that is
    //the material that has not been expelled yet.

    std::string logmess="";

    ///Initial check
    //There are two cases in which we can immediately exit from the function:
    //1-Safety check, the primary should have a core+envelope
    if(  (primary->getp(MHE::ID)==0)  or (primary->aminakedhelium() and primary->getp(MCO::ID)==0)  )
        svlog.critical("The primary star in binary " + binstar->get_name()
                       + " (ID: " + utilities::n2s(binstar->get_ID(),__FILE__,__LINE__) + ")"
                       + " has not envelope to loose"
                       + " ( M: " + utilities::n2s(primary->getp(Mass::ID),__FILE__,__LINE__) + " Msun,"
                       + " MHE: " + utilities::n2s(primary->getp(MHE::ID),__FILE__,__LINE__) + " Msun,"
                       + " MCO: " + utilities::n2s(primary->getp(MHE::ID),__FILE__,__LINE__) + " Msun)."
                       + "\nThis is akward since this star started the CE",
                       __FILE__,__LINE__,sevnstd::ce_error());
    //2-If the secondary is a BH or a NS the result is a unstable Thorne object.
    //that leaves only the core (i.e. the secondary)
    //At this moment also a secondary WD is considered unstable and leaves only the core.
    //TODO The inclusion of the WD type here is just a placehorder, we have to found a better handling
    //Select NS, BH and WD.
    else if (secondary->amiremnant() and !primary->amiremnant()){
        //primary->set_empty_in_bse(binstar); //Primary is destroyed, secondary is not changed
        //std::string logmess = Mix::log_message(binstar,primary,secondary);
        logmess=log_message_swallowed(binstar); //Log Swallowed
        primary->set_empty(); //Primary is destroyed, secondary is not changed
    }
    //Note the primary cannot ne a NS or BH or WD, This is always the case because these objects have not envelopes
    //We put just a a check to be sure
    else if(secondary->amiremnant() and primary->amiremnant()){
        svlog.critical("CE colascence between two remnants. This is impossible: a system with two remnants "
                       "cannot trigger the CE",__FILE__,__LINE__,sevnstd::ce_error());
    }
    else{
        //Update maximum co and minimum mhe before coalescence
        primary->set_MCO_max_aftermerge(secondary);
        primary->set_MHE_min_aftermerge(secondary);

        //3A-Coalesce with the old SEVN formalism using binding energy
        if(binstar->get_svpar_num("ce_kce")==-1 && binstar->get_svpar_num("ce_knce")==-1)
            logmess=coalesce_with_bindingEnergy(binstar);
        //3B- Coalesce wit ht new SEVN2 formalism
        else
            logmess=coalesce(binstar);
    }



    //In any case after the coalescence there is anymore a binary system but only a new single star.
    binstar->set_broken(); //System is broken
    binstar->set_onesurvived(); //Onesurveived is true, one star (the secondary) survived

    return logmess;

}

std::string CommonEnvelope::coalesce_with_bindingEnergy(Binstar *binstar) {

    //Check
    if (!isinisialised)
        svlog.critical("CE process has not been initialised",__FILE__,__LINE__,sevnstd::ce_error());


    ///ESTIMATE THE FINAL ENVELOPE ENERGY
    //1-Estimate the effective a_fin, it is the one for wich RL=rc
    //   - We know RL = f  * a_fin
    //   - We want Rc = f * a_fin_true so:
    //   - Rc = RL/a_fin * a_fin_true and
    //   - a_fin_true = Rc/RL * a_fin
    if(R_core_primary/RL_primary_final > R_core_secondary/RL_secondary_final)
        a_fin *= R_core_primary / RL_primary_final;
    else
        a_fin *= R_core_secondary / RL_secondary_final;
    //2-Estimate the final orbital energy
    Eorb_fin = -0.5*M_core_primary*M_core_secondary/a_fin;
    //3-Estimate new envelope binding energy
    // We have that Delta Ebin = - alpha*Delta Eorb (What we gain in orbital energy, we lose from envelope binding)
    // E_bin_f - Ebind_i = - alpha*(Eorb_f - Eorb_i)
    // therefore E_bin_f = Ebind_i - alpha*(Eorb_f - Eorb_i).
    // However i SEVN1 it is writte as E_bin_f = lambda*(Ebind_i - alpha*(Eorb_f - Eorb_i))
    // This is justified as:
    //      doing that I can find the best match in track using just the new core and then
    //      present mass and radius eq. 74, Hurley et al. 2002, MNRAS 329, 897

    binstar->CE_E_tomatch =  (Ebind_ini - alpha*(Eorb_fin-Eorb_ini));
    svlog.pdebug("The final stars will an Envelope binding energy",binstar->CE_E_tomatch,lambda*Ebind_ini,M_core_primary,secondary->getp(Mass::ID),__FILE__,__LINE__);


    ///ESTIMATE THE CORE OF THE NEW STAR

    //Test Minimum mass afger the coalescence
    //double Min_mass = secondary->getp(Phase::ID)==Lookup::Phases::MainSequence? primary->getp(MHE::ID) + secondary->getp(Mass::ID) : primary->getp(MHE::ID)+secondary->getp(MHE::ID);
    double Min_mass = 0;//TODO This is not used inside find_track_after_CE_Ebinding
    double Max_mass = primary->getp(Mass::ID) + secondary->getp(Mass::ID);

    //In this situation the cores of the two stars are merging, so first we set the
    //mass of the new core, the mass of the envelope is found  in change track
    //Notice that even if a pureHe main sequence is treated as a core, we let all the helium to be accreted following Hurley, 2002
    primary->update_from_binary(MCO::ID, secondary->getp(MCO::ID));
    if (!primary->aminakedhelium() or !secondary->aminakedhelium()){
        primary->update_from_binary(MHE::ID, secondary->getp(MHE::ID)); //MHE is always HE+CO
        //Find new track,
        primary->find_track_after_CE_Ebinding(binstar->CE_E_tomatch,Min_mass,Max_mass);
    }
    else {
        //Find new track, ofr pureHE
        primary->find_track_after_CE_Ebinding(binstar->CE_E_tomatch,Min_mass,Max_mass,true);
    }

    std::string logmess = Mix::log_message(binstar,primary,donor);

    //After the coalescence the secondary does not exist anymore
    //secondary->set_empty_in_bse(binstar);
    secondary->set_empty();

    return logmess;
}

std::string CommonEnvelope::coalesce(Binstar *binstar) {
    //Check
    if (!isinisialised)
        svlog.critical("CE process has not been initialised",__FILE__,__LINE__,sevnstd::ce_error());

    ///Estimate the effective a_fin, it is the one for wich RL=rc
    //   - We know RL = f  * a_fin
    //   - We want Rc = f * a_fin_true so:
    //   - Rc = RL/a_fin * a_fin_true and
    //   - a_fin_true = Rc/RL * a_fin
    if(R_core_primary/RL_primary_final > R_core_secondary/RL_secondary_final)
        a_fin *= R_core_primary / RL_primary_final;
    else
        a_fin *= R_core_secondary / RL_secondary_final;


    //ESTIMATE THE CORE OF THE NEW STAR
    //In this situation the cores of the two stars are merging, so first we set the
    //mass of the new core, the mass of the envelope is found  in change track outside this function (maybe)
    //Notice that even with a pureHe main sequence, we let all the helium to be accreted following Hurley, 2002
    double DeltaMHE_final =  secondary->getp(MHE::ID); //MHE finale - MHE initial
    double DeltaMCO_final =  secondary->getp(MCO::ID); //MCO finale - MCO initial
    //ESTIMATE THE CORE OF THE NEW STAR
    double Mfinal  =  final_mass_after_coalescence();
    //Check
    if (Mfinal>(primary->getp(Mass::ID)+secondary->getp(Mass::ID))){
        svlog.pwarning("Binstar",binstar->get_name(),"ID(",binstar->get_ID(),")","In coalesce the final mass is larger than the initial total mass. It will be lowered"
                      " to the initial total mass.",__FILE__,__LINE__);
        Mfinal=primary->getp(Mass::ID)+secondary->getp(Mass::ID);
    }


    double DeltaMf =  Mfinal - primary->getp(Mass::ID); //MF final - Mf initial
    //Update variation
    double MCO_old = primary->getp(MCO::ID);
    double MHE_old = primary->getp(MHE::ID);

    // Other checks
    if (MCO_old==0 and DeltaMCO_final>0)
        svlog.critical("This is an embarrassing situation, in CE the secondary has a CO core and the primary not",
                       __FILE__,__LINE__,sevnstd::ce_error());
    else if (MHE_old==0 and DeltaMHE_final>0)
        svlog.critical("This is an embarrassing situation, in CE the secondary has a HE core and the primary not",
                       __FILE__,__LINE__,sevnstd::ce_error());
        //In this case we are changing the mass of the inner core (e.g. accreting He with a CO core, or accreting just total mass with a He core)
    else if ((DeltaMCO_final>0 or (MCO_old==0 and DeltaMHE_final>0))
            or ( (MCO_old!=0 and DeltaMHE_final>0) or (DeltaMCO_final==0 and DeltaMHE_final==0)))
        ;// This is ok
    else {
        svlog.pdebug(DeltaMCO_final,MCO_old,DeltaMHE_final,primary->getp(MHE::ID),primary->getp(MCO::ID),secondary->getp(MHE::ID),secondary->getp(MCO::ID));
        svlog.critical("Unknow situation in coalescence after CE",__FILE__,__LINE__,sevnstd::ce_error());
    }

    // Now handle the merger creating the new star
    primary->update_from_binary(Mass::ID, DeltaMf);
    primary->update_from_binary(MCO::ID, DeltaMCO_final);
    primary->update_from_binary(dMcumul_binary::ID, DeltaMf);
    secondary->update_from_binary(dMcumul_binary::ID, -DeltaMf);


    //If both are naked helium, we update the Helium Mass (we already updated the total mass)
    //In the following force jump is always set to true
    if (primary->aminakedhelium() and secondary->aminakedhelium()){
        primary->update_from_binary(MHE::ID, primary->getp(Mass::ID)-primary->getp(MHE::ID));
        primary->find_new_track_after_merger();
    }
    else{
        primary->update_from_binary(MHE::ID, DeltaMHE_final);
        if(primary->aminakedhelium())
            primary->jump_to_normal_tracks(); //do not need to  set is_merging=true, since this is always true when using jump_to_normal_tracks
        else
            primary->find_new_track_after_merger();
    }

    std::string logmess = Mix::log_message(binstar,primary,secondary);

    //After the coalescence the secondary does not exist anymore
    secondary->set_empty();

    return logmess;
}

int CommonEnvelope::lose_the_envelope(Binstar *binstar) {

    //Check
    if (!isinisialised)
        svlog.critical("CE process has not been initialised",__FILE__,__LINE__,sevnstd::ce_error());



    //The return of lose_the_envelope is EXIT_FAILURE if the envelope
    //has not been removed, e.g. because the star has not a core-envelope segregation
    //if a CE is happening the primary NEED to have a core-envelope structure, therefore
    //if the return is EXIT_FAILURE for the primary return an error.
    if (primary->lose_the_envelope()==EXIT_FAILURE)
        svlog.critical("The primary star in binary " + binstar->get_name()
                       + " (ID: " + utilities::n2s(binstar->get_ID(),__FILE__,__LINE__) + ")"
                       + " has not envelope to loose"
                       + " ( M: " + utilities::n2s(primary->getp(Mass::ID),__FILE__,__LINE__) + " Msun,"
                       + " MHE: " + utilities::n2s(primary->getp(MHE::ID),__FILE__,__LINE__) + " Msun,"
                       + " MCO: " + utilities::n2s(primary->getp(MHE::ID),__FILE__,__LINE__) + " Msun)."
                       + "\nThis is akward since this star started the CE",
                       __FILE__,__LINE__,sevnstd::ce_error());
    secondary->lose_the_envelope();

    utilities::wait("Envelope loosed",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);


    return EXIT_SUCCESS;
}

//This is a standalone function, do not need common variabile initialisation/
double CommonEnvelope::final_mass_after_coalescence(){

    //Check
    if (!isinisialised)
        svlog.critical("CE process has not been initialised",__FILE__,__LINE__,sevnstd::ce_error());

    //utilities::hardwait("COALESCE",__FILE__,__LINE__);

    double kce = primary->get_svpar_num("ce_kce");
    double knce = primary->get_svpar_num("ce_knce");


    bool check_ms =  secondary->getp(Phase::ID)==Lookup::Phases::MainSequence or secondary->getp(Phase::ID)==Lookup::Phases::CoreHeBurning;
    bool check_both_HE   = primary->aminakedhelium() and secondary->aminakedhelium();
    double Mf=0, Mc=0, Mnce=0, Mce=0;

    //TODO NEED TO BE TESTED
    if (check_both_HE and check_ms){
        Mc = primary->getp(MCO::ID);
        Mnce = secondary->getp(Mass::ID);
        Mce = primary->getp(Mass::ID)-primary->getp(MCO::ID);
    }
    else if (check_both_HE){
        Mc  = primary->getp(MCO::ID) + secondary->getp(MCO::ID);
        Mce = (primary->getp(Mass::ID)-primary->getp(MCO::ID)) + (secondary->getp(Mass::ID)-secondary->getp(MCO::ID));
    }
    else{
        //Primary
        if(primary->aminakedhelium()){
            Mc+=primary->getp(Mass::ID);
        }
        else{
            Mc+=primary->getp(MHE::ID);
            Mce+=primary->getp(Mass::ID)-primary->getp(MHE::ID);
        }

        //Secondary
        if(secondary->aminakedhelium())
            Mc+=secondary->getp(Mass::ID);
        else if (check_ms)
            Mnce+=secondary->getp(Mass::ID);
        else{
            Mc+=secondary->getp(MHE::ID);
            Mce+=secondary->getp(Mass::ID)-secondary->getp(MHE::ID);
        }
    }



    //Options
    //Use the Hurley, 2002 equation 77 without mass limit
    //The mass limit actually are implicit in the hurley_final_mass between the mass of the cores and the total masses
    if (knce==-1 && kce>=0){
        double afin_eff;
        if(R_core_primary/RL_primary_final > R_core_secondary/RL_secondary_final)
            afin_eff = a_fin* R_core_primary / RL_primary_final;
        else
            afin_eff = a_fin* R_core_secondary / RL_secondary_final;

        Eorb_fin = -0.5*M_core_primary*M_core_secondary/afin_eff;
        double Ebind_final = (Ebind_ini - alpha*(Eorb_fin-Eorb_ini));
        Mf = hurley_final_mass(Ebind_final);

        //TODO what happens, when we have Helium? it should go to the core even if is a evolved pureHe star,check!
    }
    //Use the Hurley, 2002 eq 77, but only for the CE part rescaling for the minimum and maximum mass.
    //Use the factor knce for the mass non participating to the CE (e.g. MS secondary)
    else if (kce==-1 && knce>=0){
        double afin_eff;
        if(R_core_primary/RL_primary_final > R_core_secondary/RL_secondary_final)
            afin_eff = a_fin* R_core_primary / RL_primary_final;
        else
            afin_eff = a_fin* R_core_secondary / RL_secondary_final;

        Eorb_fin = -0.5*M_core_primary*M_core_secondary/afin_eff;
        double Ebind_final = (Ebind_ini - alpha*(Eorb_fin-Eorb_ini));
        Mf = hurley_final_mass(Ebind_final);
        //rescale
        double MF_BSE_min = primary->getp(MHE::ID)+secondary->getp(MHE::ID);
        double MF_BSE_max = primary->getp(Mass::ID)+secondary->getp(Mass::ID);
        double KCE =  (Mf - MF_BSE_min) / (MF_BSE_max - MF_BSE_min);

        Mf = Mc + knce * Mnce + KCE*Mce;

    }
    //New formalism
    else if(kce>=0 && knce>=0){
        Mf = Mc + knce * Mnce + kce * Mce;
    }
    //Raise error
    else{
        svlog.critical("ce_kce and ce_knce are both negative, cannot find a final mass after coalesecence in CE",
                __FILE__,__LINE__,sevnstd::sevnerr());
    }

    return Mf;
}

std::string CommonEnvelope::log_mess(Binstar *binstar, Star *primary, Star *secondary) {

    std::string info_primary_old = utilities::log_star_info(primary,true);
    std::string info_secondary= utilities::log_star_info(secondary,false);


    std::string w = utilities::log_print("CE", binstar,
                                         info_primary_old,
                                         info_secondary);

    return w;

}

std::string CommonEnvelope::log_message_swallowed(Binstar *binstar) {
    //If called here, it means that the whole donor has been removed without mergning so the DM accreted is 0
    return RocheLobe::log_message_swallowed(binstar,primary,secondary,0.0);
}

