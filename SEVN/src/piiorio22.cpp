//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <piiorio22.h>
#include <star.h>


///Iorio22 model

//Ctor
PIIorio22::PIIorio22(bool reg) : PairInstability(false) {
    if (reg){
        Register(this);
    }
}

utilities::MassContainer PIIorio22::apply_beforeSN(Star *s) const {
    double MCO=s->getp(MCO::ID);
    double MHE=s->getp(MHE::ID);
    double Mass=s->getp(Mass::ID);
    double Mass_preSN{Mass}, MHE_preSN{MHE}, MCO_preSN{MCO};

    //Pulsation, remove the envelope Hydrogen and Helium envelope
    if (MCO>=38 and MHE<=135){
        Mass_preSN=MHE_preSN=MCO_preSN;
    }

    return utilities::MassContainer{Mass_preSN,MHE_preSN,MCO_preSN};
}

PISNreturn PIIorio22::apply_afterSN(Star *s, double mremnant) const {
    double MCO=s->getp(MCO::ID);
    double MHE=s->getp(MHE::ID);

    double Mremnant_final;
    auto sntype = Lookup::SNExplosionType::Unknown;

    //CCore collapse
    if (MCO<38 or MHE>135){
        Mremnant_final = mremnant;
        sntype = Lookup::SNExplosionType::CoreCollapse;
    }
    //PPISN
    else if (MCO>=38 and MCO<=60){
        Mremnant_final = ppisn_mass_limitation(mremnant,s);
        sntype = Lookup::SNExplosionType::PPISN;
    }
    //PISN
    else if (MCO>60){
        Mremnant_final = 0.;
        sntype = Lookup::SNExplosionType::PISN;
    }
    else{
        svlog.critical("Unknown case in PIIorio22::apply_afterSN, MCO="
                                +utilities::n2s(MCO,__FILE__,__LINE__)+
                                ", MHE="+utilities::n2s(MCO,__FILE__,__LINE__)+
                               ", MCO="+utilities::n2s(MCO,__FILE__,__LINE__),__FILE__,__LINE__,
                               sevnstd::sanity_error(""));
    }

    return PISNreturn{Mremnant_final, sntype};
}


PIIorio22Limited::PIIorio22Limited(bool reg) : PIIorio22(reg) {
    if (reg){
        Register(this);
    }
}

/**
 * Apply a mass limit to the BH mass after ppisn
 * Limit the BH mass to the MBH_max obtained with Eq. 2 in Farmer+19:
 * MBH_max=35.1 -3.9 log10(Z) -0.31 log10(Z)*log10(Z), where Z is the metallicity
 * @param Mremnant Remnant mass before the application of the mass limit
 * @param s Pointer to the exploding star
 * @return New value of the remnant mass after the application of the limits
 */
double PIIorio22Limited::ppisn_mass_limitation(double Mremnant, Star *s) const {

    const double coeffb1{35.1}, coeffb2{-3.9}, coeffb3{-0.32};
    const double logZ=std::log10(s->get_Z());

    double MBH_max = coeffb1 + coeffb2*logZ + coeffb3*logZ*logZ; //Eq. 2 in Farmer+19

    return std::min(Mremnant,MBH_max);
}



