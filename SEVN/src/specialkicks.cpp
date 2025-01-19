//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <specialkicks.h>
#include <star.h>
#include <supernova.h>

/********* CC15 **************/

void CC15::_apply(Star *s) {
    //maxwellian velocities with sigma = 265 km/s, as in Hobbs et al., 2005, MNRAS 360 974
    //all kicks are rescaled accordingly to the fallback fraction (direct collapse = no kicks)

    s->vkick[0] = gaussian15(utilities::mtrand);
    s->vkick[1] = gaussian15(utilities::mtrand);
    s->vkick[2] = gaussian15(utilities::mtrand);
    s->vkick[3] = random_velocity_kick =sqrt(s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);

    if ( (s->get_supernova()->get_remnant_type()==Remnants::BH)){
        double f_fb_correction = estimate_fallback_correction(s);
        for (auto& v : s->vkick)
            v=v*f_fb_correction;
    }
}


/********* EC15CC265 **************/

void EC15CC265::_apply(Star *s) {

    //If ECSN explosion, use the ECSN Maxwellian, for CCSN use the CCSN Maxwellian
    //all kicks are rescaled accordingly to the fallback fraction (direct collapse = no kicks)

    if(s->get_supernova()->get_remnant_type()==Remnants::NS_ECSN){
        s->vkick[0] = gaussian_ecsn(utilities::mtrand);
        s->vkick[1] = gaussian_ecsn(utilities::mtrand);
        s->vkick[2] = gaussian_ecsn(utilities::mtrand);
        s->vkick[3] = random_velocity_kick =sqrt(s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);
    } else{
        s->vkick[0] = gaussian_ccsn(utilities::mtrand);
        s->vkick[1] = gaussian_ccsn(utilities::mtrand);
        s->vkick[2] = gaussian_ccsn(utilities::mtrand);
        s->vkick[3] = random_velocity_kick =sqrt(s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);
    }

    if ( (s->get_supernova()->get_remnant_type()==Remnants::BH)){
        double f_fb_correction = estimate_fallback_correction(s);
        for (auto& v : s->vkick)
            v=v*f_fb_correction;
    }

}


/********* ECUS30 **************/

void ECUS30::_apply(Star *s) {

    //A SN is considered an Ultra Stripped SN if the mass other than the CO core mass is smaller than 0.1
    //This is always true for nakedCO stars.
    //Notice here we can use getp since the stars propertie have not yet been changed when kicks apply is called
    //TODO Call Mass and MCO here can be dangerous, we have to have a clear standard on the developer guide when the properies will change after a SN
    bool isUSSN = (s->getp(Mass::ID)-s->getp(MCO::ID)<=0.1);

    //First case, ECSN or USSN, they receive a kick from a Maxwellian with sigma 30 km/s
    if (s->get_supernova()->get_remnant_type()==Remnants::NS_ECSN or isUSSN){
        s->vkick[0] = draw_from_gaussian(sigma_ecsn);
        s->vkick[1] = draw_from_gaussian(sigma_ecsn);
        s->vkick[2] = draw_from_gaussian(sigma_ecsn);
    }
        //Second case, CCSN (Notice we are not correcting for the BH fallback)
    else{
        double hobbs_std = s -> get_svpar_num("sn_kick_velocity_stdev");
        s->vkick[0] = draw_from_gaussian(hobbs_std);
        s->vkick[1] = draw_from_gaussian(hobbs_std);
        s->vkick[2] = draw_from_gaussian(hobbs_std);
    }


    //Vkick magnitude
    s->vkick[3] = random_velocity_kick =sqrt(s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);


}
