//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <pimapelli20.h>
#include <star.h>


///Mapelli20model

PIMapelli20::PIMapelli20(bool reg) : PairInstability(false) {
    if (reg){
        Register(this);
    }
}


PISNreturn PIMapelli20::apply_afterSN(Star *s, double mremnant) const {

    Lookup::SNExplosionType sntype = Lookup::SNExplosionType::Unknown; //initialise to Unknown
    double mremnant_after_pisn=pisn_correction(mremnant,s,sntype);
    return     PISNreturn{mremnant_after_pisn, sntype};
}

double PIMapelli20::pisn_correction(const double mass, Star *s, Lookup::SNExplosionType& sntype) {

    //Appendix A, Mapelli+20 (https://ui.adsabs.harvard.edu/abs/2020ApJ...888...76M/abstract)
    // (based on fit to the Woosley17 tables, https://ui.adsabs.harvard.edu/abs/2017ApJ...836..244W/abstract)

    //default: no correction, i.e. pisn_correction factor = 1.0, SN type core collapse;
    double pisn_correction = 1.0;
    sntype=Lookup::SNExplosionType::CoreCollapse;

    //Helium fraction
    double he_fin = s->getp(MHE::ID);
    double m_fin = s->getp(Mass::ID);
    double he_frac = he_fin/m_fin;

    //k_param (pisn_correction dependence upon helium fraction)
    double k_param = 0.67*he_frac + 0.1;

    //pulsation pair-instability supernova: pulses
    if(he_fin > 32.0 && he_fin < 64.0){
        sntype = Lookup::SNExplosionType::PPISN; //pulsation pair-instability supernova: pulses
        //use Table 2 fit from Woosley 2017
        if(he_frac < 0.9){
            if(he_fin <= 37.0)
                pisn_correction = (k_param - 1.0)/5.0*he_fin + (37.0 - 32.0*k_param)/5.0;
            else if (he_fin > 37.0 && he_fin <= 60.0)
                pisn_correction = k_param;
            else
                pisn_correction = -(k_param/4.0)*he_fin + 16.0*k_param;
        }
            //use WR table 1 fit from Woosley 2017
        else{
            if(he_fin <= 37.0)
                pisn_correction = (0.5226*he_frac-0.52974)*(he_fin-32.0) + 1.0;
            else if (he_fin > 37.0 && he_fin <= 56.0){
                double val = (0.5226*he_frac-0.52974)*5.0 + 1.0;
                if(val < 0.82916)
                    pisn_correction = val;
                else
                    pisn_correction = (-0.1381*he_frac+0.1309)*(he_fin-56.0) + 0.82916;
            }
            else
                pisn_correction = -0.103645*he_fin + 6.63328;
        }

    }

        //pair-instability supernova: disintegrated
    else if (he_fin >= 64.0 && he_fin < 135.0) {
        sntype = Lookup::SNExplosionType::PISN; //pair-instability supernova: disintegrated
        pisn_correction = 0.0;
    }
        //end of PISN... again standard direct collapse
    else{
        sntype = Lookup::SNExplosionType::CoreCollapse; //core-collapse (direct-collapse)
        pisn_correction = 1.0;
    }

    double corrected_mass=pisn_correction*mass;

    //Disable the creation of

    //The minimum remnant in the Woosley17 tables is 4.5 Msun, we assume
    //that if the PISN correction (pisn_correction!=1) create a final mass
    //below this number the system is actually destroyed by a PISN
    if (corrected_mass<4.5 and pisn_correction!=1){
        sntype = Lookup::SNExplosionType::PISN; //pair-instability supernova: disintegrated
        corrected_mass = 0.0;
    }
    //TODO The alternativa idea below is more similar to what we use now in Farmer (but in Farmer this problem never happens)
    //Alternative idea, now if the corrected mass is below 4.5 Msun, 4.5 Msun is used instead
    //corrected_mass = std::max(corrected_mass,4.5);


    return corrected_mass;
}

