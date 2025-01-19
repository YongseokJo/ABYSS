//
// Created by spera on 30/12/18.
//



#include <supernova.h>
#include <pairinstability.h>
#include <star.h>
#include <remnant.h>

supernova::supernova(Star *s){

    fallback_frac            = -1.0; //it means "not-set, yet"
    Average_remnant          = -1.0; //it means "not-set, yet"
    Average_ejected          = -1.0; //it means "not-set, yet"
    Mremnant_preCorrection   = -1.0; //it means "not-set, yet"
    pisn_correction          = 1;
    remnant_type             = Lookup::Remnants::NotARemnant; //not a remnant
    sn_type                  = Lookup::SNExplosionType::Unknown;

    if(s != nullptr) { //This means that I am calling the SN constructor from the specific instance of the SN type, and not just for registring

        //Set the kick model
        std::string kicksmodel = s->get_svpar_str("sn_kicks");
        kick = Kicks::Instance(kicksmodel);
        if (kick == nullptr)
            svlog.critical("Unknown SN kick model: [" + kicksmodel + "]", __FILE__, __LINE__);

        //Set the pair instability model
        std::string pairinstabilitymodel = s->get_svpar_str("sn_pairinstability");
        pairinstability = PairInstability::Instance(pairinstabilitymodel);

        //Set the neutrino mass loss model
        neutrinomassloss = NeutrinoMassLoss::Instance(s->get_svpar_str("sn_neutrinomloss"));

    }


}

supernova::~supernova(){
    delete kick;
    kick = nullptr;
    delete pairinstability;
    pairinstability= nullptr;
    delete neutrinomassloss;
    neutrinomassloss = nullptr;
    delete Mass_preSN;
    Mass_preSN= nullptr;
}

void supernova::initialise_remnant(Star *s, double Mass_remnant, Lookup::Remnants Remnant_type){

    ///Check input
    if (Mass_remnant>0 and Remnant_type==Lookup::Remnants::Empty)
        svlog.critical("It is not possible to use a Mass (" + utilities::n2s(Mass_remnant,__FILE__,__LINE__) +
        ") larger than 0 with Remnant_type=Empty",__FILE__,__LINE__,sevnstd::sevnio_error(""));
    else if (Mass_remnant<=0 and Remnant_type!=Lookup::Remnants::Empty)
        svlog.critical("It is not possible to use a Mass<=0 (" + utilities::n2s(Mass_remnant,__FILE__,__LINE__) +
        ") if the Remnant_type is not Empty",__FILE__,__LINE__,sevnstd::sevnio_error(""));
    else if (Mass_remnant>s->get_svpar_num("sn_Mchandra") and (Remnant_type==Lookup::Remnants::HeWD or Remnant_type==Lookup::Remnants::COWD or Remnant_type==Lookup::Remnants::ONeWD))
        svlog.critical("It is not possible to use a Mass (" + utilities::n2s(Mass_remnant,__FILE__,__LINE__) +
        ") larger than sn_Mchandra with Remnant_type=WD",__FILE__,__LINE__,sevnstd::sevnio_error(""));
    else if (Mass_remnant>s->get_svpar_num("sn_max_ns_mass") and (Remnant_type==Lookup::Remnants::NS_CCSN or Remnant_type==Lookup::Remnants::NS_ECSN) )
        svlog.critical("It is not possible to use a Mass (" + utilities::n2s(Mass_remnant,__FILE__,__LINE__) +
        ") larger than sn_max_ns_mass with Remnant_type=NS",__FILE__,__LINE__,sevnstd::sevnio_error(""));


    ///
    if (Mass_remnant>0.0 ) {
        Mremnant = Mass_remnant;
        remnant_type =  Remnant_type;
        set_staremnant(s);
        s->set_remnant();
    }
    else {
        Mremnant = 0.0;
        remnant_type = Lookup::Remnants::Empty;
        s->set_empty();
    }

    Mejected = 0;


    return;

}

supernova* supernova::Instance(std::string const &name, Star *s){

    auto it = GetStaticMap().find(name);


    //If the name is in the list update the counter
    if (it != GetStaticMap().end()){
        GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    }
    //otherwise return an error
    else{
        SevnLogging svlog;
        std::string available_options;
        for (auto const& option : GetStaticMap()) available_options+=option.first+", ";
        svlog.critical("Option "+name+" not available for Process Supernova. The available options are: \n"+available_options,
                       __FILE__,__LINE__,sevnstd::notimplemented_error(" "));
    }


    return (it->second)->instance(s); //when I call the instance create a new object with the non trivial constructor

}


double supernova::neutrino_mass_loss(const double mass, Star *s) const {
    return neutrinomassloss->apply(mass,s);
}

void supernova::main(Star *s) {

    double Mtot_before=s->getp(Mass::ID);
    double MHE_before=s->getp(MHE::ID);
    double MCO_before=s->getp(MCO::ID);

    remnant_properties(s);
    s->set_remnant();
    std::string w="";

    if (remnant_type== Lookup::Remnants::HeWD or remnant_type== Lookup::Remnants::COWD or remnant_type== Lookup::Remnants::ONeWD){
        w = utilities::log_print("WD",s,Mtot_before,MHE_before,MCO_before,Mremnant,remnant_type);
    }
    else { //Apply kicks
        kick->apply(s);
        w = utilities::log_print("SN",s,Mtot_before,MHE_before,MCO_before,Mremnant,remnant_type,sn_type,
                kick->get_random_kick(),s->vkick[3],s->vkick[0],s->vkick[1],s->vkick[2]);
    }



    s->print_to_log(w);

}

void supernova::WDformation(_UNUSED Star *s) {

    double preSNMCO = get_preSN_Mass(MCO::ID);
    double preSNMHE = get_preSN_Mass(MHE::ID);

    //GI: 26/04/2023
    //TODO: Is this the best implementation we can have? It is true that a stripped low mass star never ignite Helium?
    //Following Hurley, a HeWD is formed when a low-mass star (a star with MZAMS<MHEf, mass helium flash) loses the envelope when it is
    //has not yet started the core Helium burning (see setting kw = 10 in hrdiag.f).
    //So here we just assume that a HeWD is formed if the stars arrive at this point with a MCO=0. ATM, it is not really possible
    //to arrive here through normal SSE with MCO=0, but this method can be triggered externally by some event in binary.
    //This is what we actually implemented in the method used to transform the star in a pureHe (Star::jump_to_pureHE_tracks), if it has not a CO core yet and Mzams< approx 2,
    //the turn_into_remnant method is called and we enter in the first condition forming a HeWD.
    if (preSNMCO<=utilities::TINY){
        Mremnant = preSNMHE; //WD mass taken as the mass of the HE core
        remnant_type = Lookup::Remnants::HeWD;
    }
    else {
        Mremnant = preSNMCO; //WD mass taken as the mass of the CO core
        if (preSNMHE < 1.6) //TODO 1.6 limit is taken from Hurley (M_{c,BAGB}), Section 6 of https://academic.oup.com/mnras/article/315/3/543/972062
            remnant_type = Lookup::Remnants::COWD;
        else
            remnant_type = Lookup::Remnants::ONeWD;
    }

    //set the sn_type
    sn_type = Lookup::SNExplosionType::WDformation;
}


void supernova::ECSN(Star *s) {
    double preSNMCO = get_preSN_Mass(MCO::ID);
    Mremnant = Mass_corrections_after_explosion(preSNMCO, s); //TODO s->getp(MCO::ID) or CO_LOWER_ECSN ?????
    remnant_type = Lookup::Remnants::NS_ECSN;
    //fallback_frac = 0.; //Unset
    //set the sn_type
    sn_type = Lookup::SNExplosionType::ElectronCapture;
}

void supernova::CCexplosion(Star *s){
    //std::cout<<" Before "<< Mremnant << " " << remnant_type <<" "<<sn_type<<std::endl;
    explosion(s); //let SN explode and set Mremnant
    //std::cout<<" After explosion "<< Mremnant << " " << remnant_type <<" "<<sn_type<<std::endl;
    Mremnant = Mass_corrections_after_explosion(Mremnant, s); //Correct Mremnant for PISN + neutrino mass loss
    //std::cout<<" After correction "<< Mremnant << " " << remnant_type <<" "<<sn_type<<std::endl;
    set_remnant_type_after_explosion(s,Mremnant,sn_type); //Finally set the remnant type
    //std::cout<<" After rem type "<< Mremnant << " " << remnant_type <<" "<<sn_type<<std::endl;
}

void supernova::set_staremnant(Star *s) {

    if(remnant_type == Lookup::Remnants::NS_CCSN) {
        s->set_staremnant(new NSCCrem(s,Mremnant));
    }
    else if (remnant_type == Lookup::Remnants::NS_ECSN){
        s->set_staremnant(new NSECrem(s,Mremnant));
    }
    else if (remnant_type == Lookup::Remnants::HeWD){
        s->set_staremnant(new HeWDrem(s,Mremnant));
    }
    else if(remnant_type == Lookup::Remnants::COWD ){

        s->set_staremnant(new COWDrem(s,Mremnant));
    }
    else if(remnant_type == Lookup::Remnants::ONeWD){
        s->set_staremnant(new ONeWDrem(s,Mremnant));
    }
    else if(remnant_type == Lookup::Remnants::BH){
        s->set_staremnant(new BHrem(s,Mremnant));
    }
}

void supernova::remnant_properties(Star *s) {

    double ecsn_tshold = s->aminakedhelium() and s->get_svpar_num("sn_co_lower_ecsn_pureHe")>0 ?
            s->get_svpar_num("sn_co_lower_ecsn_pureHe")  : s->get_svpar_num("sn_co_lower_ecsn");
    double sn_tshold = s->aminakedhelium() and s->get_svpar_num("sn_co_lower_sn_pureHe")>0 ?
                         s->get_svpar_num("sn_co_lower_sn_pureHe")  : s->get_svpar_num("sn_co_lower_sn");

    //0-Initialse preSN Mass
    set_preSN_Masses(s);
    double preSNMCO= get_preSN_Mass(MCO::ID);
    double preSNMass= get_preSN_Mass(Mass::ID);

    ///1-Set Mremnant e Mejected
    if(preSNMCO < ecsn_tshold) //WD formation.. no SN
        WDformation(s);
    else if(preSNMCO < sn_tshold) //ECSN.. NSs
        ECSN(s);
    else{ //SN explosion (NSs or BHs)
        CCexplosion(s);
    }

    Mejected = s->amiempty() ? 0 :  preSNMass - Mremnant; //Total ejected mass including neutrinos

    if (Mejected<0){
        svlog.critical("Mejected after a SN in negative: total Mass preSN="+
        utilities::n2s(preSNMass ,__FILE__,__LINE__)+", Remnant mass="+
        utilities::n2s(Mremnant ,__FILE__,__LINE__),__FILE__,__LINE__);
    }


    ///2-Set the staremant object
    set_staremnant(s);

}

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
Lookup::Remnants  supernova::set_remnant_type_after_explosion(Star *s, double mremnant, Lookup::SNExplosionType sntype){

    //Sanity check
    if (mremnant!=Mremnant or sntype!=sn_type){
        svlog.critical("The input of the method supernova::set_remnant_type_after_explosion (" +
                       utilities::n2s(mremnant,__FILE__,__LINE__)+", "+ utilities::n2s(sntype,__FILE__,__LINE__) +
                       "is not consistent with the supernova member Mremnant="+ utilities::n2s(Mremnant,__FILE__,__LINE__)
                        + " and sn_type="+ utilities::n2s(sn_type,__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::sanity_error(""));
    }

    ///Set Remnant type
    if (mremnant<1e-10){
        //If Mremants is 0 after correction the SN destroyed the stars (here we check a small number rather directly 0)
        remnant_type = Lookup::Remnants::Empty;
        s->set_empty();
    }
    else if(Mremnant >= s->get_svpar_num("sn_max_ns_mass")) //If the mass is larger than the max for a NS, this is a BH
        remnant_type = Lookup::Remnants::BH;
    else if(sntype==Lookup::SNExplosionType::ElectronCapture) //If EC this is a NS_ECSN
        remnant_type = Lookup::Remnants::NS_ECSN;
    else //Otherwise this is a Core Collapse NS
        remnant_type = Lookup::Remnants::NS_CCSN;


    return remnant_type;
}



void supernova::explosion_SNI(Star *s) {
    s->set_empty();
    remnant_type = Lookup::Remnants::Empty;
    //set the sn_type
    sn_type = Lookup::SNExplosionType::Ia;
}

void supernova::explosion_SNI(Star *s, Binstar *b) {
    s->set_empty_in_bse(b);
    remnant_type = Lookup::Remnants::Empty;
    //set the sn_type
    sn_type = Lookup::SNExplosionType::Ia;
}

void supernova::set_Average_for_Unified(Star *s, double default_Average_Mremnant, double default_Average_Mejected) {

    Average_remnant = s->get_svpar_num("sn_Mremnant_average")!=-1.? s->get_svpar_num("sn_Mremnant_average") : default_Average_Mremnant;
    Average_ejected = s->get_svpar_num("sn_Mejected_average")!=-1.? s->get_svpar_num("sn_Mejected_average") : default_Average_Mejected;

}

double supernova::get_preSN_Mass(size_t massid) const {

    if (Mass_preSN== nullptr){
        svlog.critical("get_preSN_Mass is returning a 0, this means that Mass_preSN is not set yet");
    }

    if (massid==Mass::ID) return Mass_preSN->Mass;
    else if (massid==MHE::ID)  return Mass_preSN->MHE;
    else if (massid==MCO::ID)  return Mass_preSN->MCO;
    else svlog.critical("get_preSN_Mass has been called with an id that is not Mass::ID, MHE::ID or MCO::ID",
                        __FILE__,__LINE__,sevnstd::sanity_error(""));


    return 0.;


}

//Set preSN mass
int supernova::_set_preSN_Masses(Star *s, bool allow_pisn) {

    //If it is not a nullptr, we have to delete the old object to avoid memory leak
    if (Mass_preSN != nullptr){
        delete Mass_preSN;
        Mass_preSN = nullptr;
        svlog.warning("preSN Mass set more than  once",__FILE__,__LINE__);
    }

    if (allow_pisn){
        Mass_preSN = new  utilities::MassContainer(pairinstability->apply_beforeSN(s));
    } else{
        Mass_preSN = new utilities::MassContainer{s->getp(Mass::ID),s->getp(MHE::ID),s->getp(MCO::ID)};
    }

    return EXIT_SUCCESS;
}

/**
 * Estimate PISN and neutrino mass loss correction
 * @param mass Mremnant mass to correct
 * @param s Pointer to the exploding star
 * @return the final remnant mass after correction for pisn and neutrino mass loss
 */
double supernova::Mass_corrections_after_explosion(const double mass, Star *s) {
    Mremnant_preCorrection = mass; //Store the mass of the remnant before the correction
    double mremnant = mass;
    mremnant = pisn(mremnant, s);
    M_neutrinos = mremnant - neutrino_mass_loss(mremnant, s);
    mremnant = mremnant - M_neutrinos;
    return mremnant;
}

///PISN options
int PisnON::set_preSN_Masses(Star *s) {
    _set_preSN_Masses(s,true);
    return EXIT_SUCCESS;
}

/**
 * Apply the PPSIN correction to the remnant mass
 * @param mass  mass of the remnant
 * @param s Pointer to the star that generated the remnant
 * @return The new mass of the remnant corrected for the PISN  (it can also be 0), it also set the SN type
 */
double PisnON::pisn(const double mass, Star *s) {
    auto pireturn = get_pairinstability()->apply_afterSN(s, mass);
    pisn_correction=pireturn.Mremnant_after_pisn/mass; //set pisn_correction
    sn_type = pireturn.SNType;
    return pireturn.Mremnant_after_pisn;
}


int PisnOFF::set_preSN_Masses(Star *s) {
    _set_preSN_Masses(s,false);
    return EXIT_SUCCESS;
}

/**
 * DO NOT Apply the PPSIN correction to the remnant mass.
 * Jut set SNtype to CoreCOllapse and return the input mass
 * @param mass  mass of the remnant
 * @param s Pointer to the star that generated the remnant
 * @return The new mass of the remnant corrected for the PISN  (it can also be 0), it also set the SN type
 */
double PisnOFF::pisn(const double mass, _UNUSED Star *s) {
    pisn_correction=1.0; //set pisn_correction
    sn_type = Lookup::SNExplosionType::CoreCollapse;
    return mass;
}

//NSfromGau
double NSfromGau::get_NS_mass(Star *s){
   double Mrem;
   Mrem=generate_random_gau(s->get_svpar_num("sn_Mremnant_average_NS"),s->get_svpar_num("sn_Mremnant_std_NS"));
   //Limit the mass to the maximum allowed NS mass
   Mrem=std::max(1.1,std::min(Mrem,s->get_svpar_num("sn_max_ns_mass")));
   return Mrem;
}

void NSfromGau::ECSN(Star *s) {
    Mremnant = get_NS_mass(s);
    //Limit to the total mass of the star
    //Mremnant = std::min(Mremnant, get_preSN_Mass(Mass::ID));
    Mremnant = std::min(Mremnant, get_preSN_Mass(MCO::ID));
    remnant_type = Lookup::Remnants::NS_ECSN;
    sn_type = Lookup::SNExplosionType::ElectronCapture;
}

void NSfromGau::CCexplosion(Star *s){
    supernova::CCexplosion(s);
    if (remnant_type==Lookup::Remnants::NS_ECSN or remnant_type==Lookup::Remnants::NS_CCSN){
        Mremnant = get_NS_mass(s);
        Mremnant = std::min(Mremnant, get_preSN_Mass(Mass::ID));
    }
}