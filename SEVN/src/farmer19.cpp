//
// Created by Giuliano Iorio on 02/08/2022.
//

#include <pifarmer19.h>
#include <star.h>

//Ctor
PIFarmer19::PIFarmer19(bool reg) : PairInstability(false) {
    if (reg){
        Register(this);
    }
}

PISNreturn PIFarmer19::apply_afterSN(Star *s, double mremnant) const {
    double MCO=s->getp(MCO::ID);
    double MHE=s->getp(MHE::ID);

    double Mremnant_final=mremnant;
    auto sntype = Lookup::SNExplosionType::Unknown;

    //CCore collapse
    if (MCO<38 or MHE>135){
        Mremnant_final = mremnant;
        sntype = Lookup::SNExplosionType::CoreCollapse;
    }
    //PISN
    else if (MCO>60){
        Mremnant_final = 0.;
        sntype = Lookup::SNExplosionType::PISN;
    }
    //PISN
    else if (MCO>=38 and MCO<=60){
        sntype = Lookup::SNExplosionType::PPISN;
        //GI 03/02/2023: I added a Z limit since this is the lowest Z simulated by Farmer+19, using a lower value of Z
        //we are extrapolating from their function, so we decide to limit the value of Z
        constexpr double Zlimit_farmer=1E-5;
        const double Z = std::max(s->get_Z(),Zlimit_farmer);
        //Constant for equation  A1 in Farmer+19
        const double a1{-0.096}, a2{8.564}, a3{-2.07}, a4{-152.97};
        //Minimum BH mass shown in the Farmer+19 plots
        const double Minmass{11.0};

        //Equation A1 in Farmer+19
        double MBHfarmer = a1*MCO*MCO + a2*MCO + a3*std::log10(Z) + a4;
        //Limit to the Minimum BH mass shown in the Farmer+19 plots
        MBHfarmer = std::max(Minmass,MBHfarmer);

        //Do not allow BH mass larger than the one coming from the SN model
        Mremnant_final = std::min(Mremnant_final,MBHfarmer);
    }
    else{
        svlog.critical("Unknown case in PIFarmer19::apply_afterSN, MCO="
                       +utilities::n2s(MCO,__FILE__,__LINE__)+
                       ", MHE="+utilities::n2s(MCO,__FILE__,__LINE__)+
                       ", MCO="+utilities::n2s(MCO,__FILE__,__LINE__),__FILE__,__LINE__,
                       sevnstd::sanity_error(""));
    }

    return PISNreturn{Mremnant_final, sntype};
}
