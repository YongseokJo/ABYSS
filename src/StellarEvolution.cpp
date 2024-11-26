#include "global.h"
#include "defs.h"

void setBHspin(Particle* ptcl) {
    std::random_device rd; // Obtain a random number from hardware
    std::mt19937 mt(rd()); // Seed the generator
    std::uniform_real_distribution<> distr(0.0, 1.0); // Define the range (0 to 1)
    REAL phi = 2 * M_PI * distr(mt);
    REAL theta = M_PI * distr(mt);

    ptcl->a_spin[0] = ptcl->star->getp(Xspin::ID) * sin(theta)*cos(phi);
    ptcl->a_spin[1] = ptcl->star->getp(Xspin::ID) * sin(theta)*sin(phi);
    ptcl->a_spin[2] = ptcl->star->getp(Xspin::ID) * cos(theta);
}

void UpdateEvolution(Particle* ptcl) {

    if (!ptcl->star->amiremnant()) {
        ptcl->dm += ptcl->Mass - ptcl->star->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->star->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->star->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
    }
    else if (ptcl->star->amiWD()) {
        fprintf(SEVNout, "WD. PID: %d, Mass: %e Msol, Radius: %e Rsol, Time: %e Myr\n", ptcl->PID, ptcl->star->getp(Mass::ID), ptcl->star->getp(Radius::ID), ptcl->EvolutionTime);
        ptcl->dm += ptcl->Mass - ptcl->star->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->star->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->star->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        ptcl->EvolutionTime = NUMERIC_FLOAT_MAX;
        if (ptcl->star->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->star->vkick[0], ptcl->star->vkick[1], ptcl->star->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->star->vkick[i]/(velocity_unit/yr*pc/1e5);
        }
    }
    else if (ptcl->star->amiNS()) {
        fprintf(SEVNout, "NS. PID: %d, Mass: %e Msol, Radius: %e Rsol, Time: %e Myr\n", ptcl->PID, ptcl->star->getp(Mass::ID), ptcl->star->getp(Radius::ID), ptcl->EvolutionTime);
        ptcl->dm += ptcl->Mass - ptcl->star->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->star->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->star->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        ptcl->EvolutionTime = NUMERIC_FLOAT_MAX;
        if (ptcl->star->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->star->vkick[0], ptcl->star->vkick[1], ptcl->star->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->star->vkick[i]/(velocity_unit/yr*pc/1e5);
        }
    }
    else if (ptcl->star->amiBH()) {
        fprintf(SEVNout, "BH. PID: %d, Mass: %e Msol, Radius: %e Rsol, Time: %e Myr\n", ptcl->PID, ptcl->star->getp(Mass::ID), ptcl->star->getp(Radius::ID), ptcl->EvolutionTime);
        setBHspin(ptcl);
        fprintf(SEVNout, "\tDimless spin. mag: %e, (%e, %e, %e)\n", ptcl->star->getp(Xspin::ID), ptcl->a_spin[0], ptcl->a_spin[1], ptcl->a_spin[2]);
        ptcl->dm += ptcl->Mass - ptcl->star->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->star->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->star->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        ptcl->EvolutionTime = NUMERIC_FLOAT_MAX;
        ptcl->ParticleType = Blackhole+SingleParticle;
        if (ptcl->star->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->star->vkick[0], ptcl->star->vkick[1], ptcl->star->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->star->vkick[i]/(velocity_unit/yr*pc/1e5);
        }
    }
    else if (ptcl->star->amiempty()) {
        fprintf(SEVNout, "Empty. PID: %d, Mass: %e Msol, Time: %e Myr\n", ptcl->PID, ptcl->star->getp(Mass::ID), ptcl->EvolutionTime);
        ptcl->dm += ptcl->Mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = 0;
        ptcl->isErase = true;
        MasslessList.push_back(ptcl);
    }
    fflush(SEVNout);

}

// Reference: int Mix::special_evolve(Binstar *binstar) in Processes.cpp of SEVN
void Mix(Star* star1, Star* star2) {

    // utilities::wait("Hey I have to mix",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);

    Star *donor = star1; //Donor is the star that will set to empty
    Star *accretor = star2; //Accretor is the star that will remain as results of the mix

    ///Choose the star that remains and the one that is set to empty
    //First handle the general case: the accretor is the more evolved star,
    // or the more compact (massive) remnant if both are remnant (handled in get_id_more_evolved_star).

    //If one of the star is empty return the other one
    int get_id_more_evolved_star;

    //In the stars have the same phase (not remnant handled before), return the more evolved is the one with the larger plife
    if (donor->getp(Phase::ID)==accretor->getp(Phase::ID))
        get_id_more_evolved_star = donor->plife()>=accretor->plife() ? 0 : 1;
    else
        get_id_more_evolved_star = donor->getp(Phase::ID)>=accretor->getp(Phase::ID) ? 0 : 1;
    if(get_id_more_evolved_star==0)
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
            throw std::runtime_error("This is an embarrassing situation, in Mix the donor has a CO core and the accretor not");
            // svlog.critical("This is an embarrassing situation, in Mix the donor has a CO core and the accretor not",
            //                 __FILE__,__LINE__,sevnstd::ce_error());
        }
        // Sanity check on MHE
        else if (MHE_old==0 and DMHE_new>0){
            throw std::runtime_error("This is an embarrassing situation, in Mix the donor has a HE core and the accretor not");
            // svlog.critical("This is an embarrassing situation, in Mix the donor has a HE core and the accretor not",
            //                 __FILE__,__LINE__,sevnstd::ce_error());
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
    }
    //Case 2 - Mix between WDs
    //TODO Notice that here mixing with a WD is treated as mixing with a NS/BH, in Hurley we have a more complicate outcomes
    else if (accretor->amiWD() and donor->amiWD()) {

        //Check if we have two HeWD
        bool check_double_HeWD = donor->getp(RemnantType::ID)==Lookup::Remnants::HeWD && accretor->getp(RemnantType::ID)==Lookup::Remnants::HeWD;

        //Case 2a - SNI explosion
        if (check_double_HeWD){
            accretor->explode_as_SNI(); // Eunwoo: this should be treated carefully!
            // if (binstar->onesurvived) binstar->set_onesurvived(); //Why this condition? I forgot, we have to check
        }
        else{
            //Update the Mass
            accretor->update_from_binary(Mass::ID, donor->getp(Mass::ID));
            // binstar->set_onesurvived();
        }

    }
    //Case 3 -  NS, BH or WD mixing with a NS/BH
    else if (accretor->amiCompact() and donor->amiremnant()){
        //Update the Mass
        accretor->update_from_binary(Mass::ID, donor->getp(Mass::ID));
    }
    //Case 4 - A not remnant donor over a NS, BH, WD
    //-> Do not accreate nothing, just set the donor to empty, see below
    //From commonenvelope::main_stellar_collision in common_envelope.cpp in SEVN1


    ///LOGPRINT AND EVENT SET
    // std::string w =Mix::log_message(binstar,accretor,donor);
    // binstar->print_to_log(w);
    //utilities::hardwait(" Event before ",get_event());
    // set_event((double)Lookup::EventsList::Merger);
    //utilities::hardwait(" Event After ",get_event());

    ///LAST STUFF
    //In any case set the binary to broken and put the donor to empty
    //donor->set_empty_in_bse(binstar);
    donor->set_empty();

}

// Staremnant* remnant;

// if (ptcl->star->amiWD())
// 	remnant = new WDrem(ptcl->star, ptcl->star->getp(Mass::ID), ptcl->EvolutionTime);
// else if (ptcl->star->amiNS())
// 	remnant = new NSrem(ptcl->star, ptcl->star->getp(Mass::ID), ptcl->EvolutionTime);
// else if (ptcl->star->amiBH())
// 	remnant = new BHrem(ptcl->star, ptcl->star->getp(Mass::ID), ptcl->EvolutionTime);

// ptcl->star->set_staremnant(remnant);

// fprintf(SEVNout, "Breaktrigger!(Remnant phase: %d), PID: %d, Mass: %e Msol, Time: %e Myr\n", ptcl->star->get_staremnant()->get_remnant_type(), ptcl->PID, ptcl->star->getp(Mass::ID), ptcl->EvolutionTime);
// fprintf(SEVNout, "Breaktrigger!(Remnant phase: %d), PID: %d, Mass: %e Msol, Time: %e Myr\n", ptcl->star->getp(RemnantType::ID), ptcl->PID, ptcl->star->getp(Mass::ID), ptcl->EvolutionTime); // not working well