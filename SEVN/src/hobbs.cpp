//
// Created by Giuliano Iorio on 26/07/2022.
//
#include <hobbs.h>
#include <star.h>
#include <supernova.h>

/********* Hobbs Pure **************/

void HobbsPure::_apply(Star *s) {
    double hobbs_std = s -> get_svpar_num("sn_kick_velocity_stdev");

    //maxwellian velocities with sigma = get_svpar_num("sn_kick_velocity_sigma"); usual: 265 km/s, as in Hobbs et al., 2005, MNRAS 360 974
    s->vkick[0] = draw_from_gaussian(hobbs_std);
    s->vkick[1] = draw_from_gaussian(hobbs_std);
    s->vkick[2] = draw_from_gaussian(hobbs_std);
    s->vkick[3] = random_velocity_kick = sqrt( s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);

}


/********* Hobbs  **************/

void Hobbs::_apply(Star *s) {

    double hobbs_std = s -> get_svpar_num("sn_kick_velocity_stdev");

    //maxwellian velocities with sigma = get_svpar_num("sn_kick_velocity_sigma"); usual: 265 km/s, as in Hobbs et al., 2005, MNRAS 360 974
    //all kicks are rescaled accordingly to the fallback fraction (direct collapse = no kicks)
    s->vkick[0] = draw_from_gaussian(hobbs_std);
    s->vkick[1] = draw_from_gaussian(hobbs_std);
    s->vkick[2] = draw_from_gaussian(hobbs_std);
    s->vkick[3] = random_velocity_kick = sqrt( s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);

    //TODO This kick model is based on the fallback fraction that is properly defined only
    //in the rapid and delayed model. So either we disable this kick model for the others sn model
    //or we can use the approximate method below to estimate the fallback frack
    double f_fb_correction = estimate_fallback_correction(s);

    for (auto& vkick : s->vkick)
        vkick*=f_fb_correction;


}

double Hobbs::estimate_fallback_correction(Star *s) {

    double f_fb_correction;

    /** New proptotype for fb estimate ***/
    //Atm if the SN model do not se the fallback frac is it is equal to -1, below there is a simple solution
    //to estimate the fallback frack assuming always a proto compact object mass of 1.1 Msun.
    //Correct new
    if (s->get_supernova()->get_fallback_frac()<0){
        //NOTICE: in the classical Rapid and Delayed formalism, the fallback is estimated
        //considering the compact remnant mass before the correction for the neutrino mass loss
        //and pair instability, therefore here we emulate this formalism considering the Mass of the remnant
        //not corrected.
        const double& Mrem   = s->get_supernova()->get_Mremnant_preCorrection();
        const double& MpreSN = s->get_supernova()->get_preSN_Mass(Mass::ID);
        const double Mproto  = 1.1;

        double fb = std::max((Mrem - Mproto),0.) / std::max(MpreSN - Mproto,Mproto);
        f_fb_correction = 1 - fb;
    }  else{
        f_fb_correction = (1.0 - s->get_supernova()->get_fallback_frac());
    }


    /*** Old implementation ***/
    //f_fb_correction =  s->get_supernova()->get_fallback_frac() != -1 ? (1.0 - s->get_supernova()->get_fallback_frac()) : 1.0;


    if (f_fb_correction>1 or f_fb_correction<0){
        svlog.critical("Fallback correction ("+utilities::n2s(f_fb_correction,__FILE__,__LINE__)+
                       ") is outside the range [0,1]",__FILE__,__LINE__,sevnstd::sanity_error(""));
    }

    return  f_fb_correction;

}
