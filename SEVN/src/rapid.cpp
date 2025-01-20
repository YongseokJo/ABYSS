//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <supernova.h>
#include <star.h>
#include <rapid.h>

///Rapid
rapid::rapid(Star *s) : supernova(s) {

    if ( s == nullptr) {
        Register(this);
    }
    else {
        //Default Mremnant, Mejected for this SN model
        double _default_Mremnant_average=1.27;
        double _default_Mejected_average=10.9;
        set_Average_for_Unified(s,_default_Mremnant_average,_default_Mejected_average);
    }

}
void rapid::explosion(_UNUSED Star *s) {

    double final_CO = get_preSN_Mass(MCO::ID);
    double final_mass = get_preSN_Mass(Mass::ID);

    double mproto = 1.1;
    double alpha_R = 0.25 - 1.275/(final_mass - mproto);
    double beta_R = 1. - 11.*alpha_R;

    if(final_CO < 2.5){
        fallback_frac = 0.2/(final_mass-mproto);
    }
    else if  (final_CO >= 2.5 && final_CO < 6.) {
        fallback_frac = (0.286*final_CO - 0.514)/(final_mass - mproto);
    }
    else if (final_CO >= 6. and final_CO < 7.){
        fallback_frac = 1;
    }
    else if ( final_CO >= 7. and final_CO < 11. ){
        fallback_frac = alpha_R*final_CO + beta_R;
    }
    else if (final_CO >= 11.){
        fallback_frac = 1;
    }
    else{
        svlog.critical("Unexpected final_CO value in delayed::explosion",__FILE__,__LINE__,sevnstd::sn_error());
    }

    Mremnant = mproto + fallback_frac*(final_mass-mproto);
}


///Rapid with NS mass drawn from a Gaussian
rapid_gauNS::rapid_gauNS(Star *s) : supernova(s), rapid(s) {
    if(s == nullptr) {
        Register(this);
    }
    else {
        //Default Mremnant, Mejected for this SN model
        double _default_Mremnant_average=s->get_svpar_num("sn_Mremnant_average_NS"); //Input parameter
        double _default_Mejected_average=Average_ejected; //Same of the parent delayed class that has been already initialised
        set_Average_for_Unified(s,_default_Mremnant_average,_default_Mejected_average);
    }
}
