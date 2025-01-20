//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <supernova.h>
#include <star.h>
#include <delayed.h>

///Delayed
delayed::delayed(Star *s) : supernova(s) {
    if(s == nullptr) {
        Register(this);
    }
    else {
        //Default Mremnant, Mejected for this SN model
        double _default_Mremnant_average=1.36;
        double _default_Mejected_average=10.45;
        set_Average_for_Unified(s,_default_Mremnant_average,_default_Mejected_average);
    }
}

void delayed::explosion(_UNUSED Star *s) {

    double final_CO = get_preSN_Mass(MCO::ID);
    double final_mass = get_preSN_Mass(Mass::ID);


    double mproto, alpha_D, beta_D;

    if(final_CO < 2.5){
        mproto = 1.15;
        fallback_frac = 0.2/(final_mass - mproto);
    }
    else if(final_CO >= 2.5 && final_CO < 3.5){
        mproto = 1.15;
        fallback_frac = (0.5*final_CO - 1.05)/(final_mass - mproto);
    }
    else if(final_CO >= 3.5 && final_CO < 6){
        mproto = 1.2;
        alpha_D = 0.133 - 0.084/(final_mass - mproto);
        beta_D = 1. - 11.*alpha_D;
        fallback_frac = alpha_D*final_CO + beta_D;
    }
    else if(final_CO >= 6. && final_CO < 11.){
        mproto = 1.3;
        alpha_D = 0.133 - 0.084/(final_mass - mproto);
        beta_D = 1. - 11.*alpha_D;
        fallback_frac = alpha_D*final_CO + beta_D;
    }
    else if(final_CO >= 11.){
        mproto = 1.5;
        fallback_frac = 1.;
    }
    else{
        svlog.critical("Unexpected final_CO value in delayed::explosion",__FILE__,__LINE__,sevnstd::sn_error());
    }


    Mremnant = mproto + fallback_frac*(final_mass-mproto);


    return;
}


///Delayed with NS mass drawn from a Gaussian
delayed_gauNS::delayed_gauNS(Star *s) : supernova(s), delayed(s) {
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
