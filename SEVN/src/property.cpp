//
// Created by mario on 09/02/19.
//




#include <property.h>
#include <star.h>
#include <binstar.h>
#include <iomanip>
#include <supernova.h>
#include <remnant.h>
#include <lambda_nanjing.h>

void Property::evolve_fake(Star *s){

    ///Set v0=v
    set_0_fk(value); //save the previous value, before the evolution

    ///Move to new value in the interpolating track
    FOR4 {
        val_in[_i] = val_ref[_i] + s->pos[_i];
    }

    ///Find new values at times t for the interpolatring tracks using linear interpolation
    for(int i = 0; i < 4; i++) {
        //cout<<"cycle = "<<i<<endl;
        double *_val = val_in[i];
        double *_time = s->times_in[i];

        double slope = ((*(_val + 1)) - *_val) / ((*(_time +1)) - *_time);
        double intercept = *_val - slope * (*_time);

        interpolating_values[i] = slope * s->ttimes[i] + intercept; //ttimes come from function tracktimes()

    }

    ///Find the interpolating z value, M1*v1 + M2*v2
    double val_zlow = interpolating_values[0]*wM[0] + interpolating_values[1]*wM[1];
    double val_zhigh = interpolating_values[2]*wM[2] + interpolating_values[3]*wM[3];

    ///Find the final interpolated value, v=Z1*v1 + Z2*v2
    set_fk(wZ[0]*val_zlow + wZ[1]*val_zhigh);


    //if (s->getp(Worldtime::ID)>9.8947394131 and name()=="MCO" and s->get_ID()==0)
    //    utilities::wait("AAc",value0,value, V0, V,__FILE__,__LINE__);

}

/**Table properties**/
void TableProperty::set_refpointers(Star *s) {

    for(int i = 0; i < 4; i++) {
        val_ref[i] = s->get_table_pointer_atpos(TabID(), i, 0);
        if (val_ref[i]== nullptr){
            svlog.critical("One of the interpolating tracks did not found tables for the  property " + name(),
                           __FILE__,__LINE__,sevnstd::sanity_error());
        }
    }
    table_loaded=true;

}

void TableProperty::set_w(Star *s) {
    set_wM(s);
    set_wZ(s);
}

void TableProperty::set_wM(Star *s) {
    std::string wM_option = s->get_svpar_str("ev_setwM");

    if (wM_option=="linear"){
        set_wM_linear(s);
    } else if (wM_option=="rational"){
        set_wM_rational(s);
    } else if (wM_option=="log"){
        set_wM_log(s);
    } else{
        svlog.critical("Unknown set M weights options for property "+
                       name(),__FILE__,__LINE__,sevnstd::sanity_error(""));
    }

}

void TableProperty::set_wZ(Star *s) {
    //Set the metallicity weights for interpolation

    double *Ztrack = s->get_ztrack();
    double Z = s->get_Z();

    wZ[0] = (Ztrack[1] - Z) / (Ztrack[1] - Ztrack[0]);
    wZ[1] = (Z - Ztrack[0]) / (Ztrack[1] - Ztrack[0]);
}

void TableProperty::set_wM_linear(Star *s) {
    //Set linear mass weights for interpolation
    double *Mtrack = s->get_mtrack();
    double mzams = s->get_zams();
    //set also the default weights
    for(int i = 0; i < 4; i+=2) {
        wM[i] = (Mtrack[i + 1] - mzams) / (Mtrack[i + 1] - Mtrack[i]);
        wM[i+1] = (mzams - Mtrack[i]) / (Mtrack[i + 1] - Mtrack[i]);
    }
}

void TableProperty::set_wM_rational(Star *s) {
    //Set  mass weights for interpolation
    double *Mtrack = s->get_mtrack();
    double mzams = s->get_zams();
    //set also the default weights
    for(int i = 0; i < 4; i+=2) {
        wM[i] = (Mtrack[i]/mzams)*(Mtrack[i + 1] - mzams) / (Mtrack[i + 1] - Mtrack[i]);
        wM[i+1] = (Mtrack[i+1]/mzams)*(mzams - Mtrack[i]) / (Mtrack[i + 1] - Mtrack[i]);
    }
}

void TableProperty::set_wM_log(Star *s) {
    //Set log mass weights for interpolation
    double *Mtrack = s->get_mtrack();
    double lmzams = std::log10(s->get_zams());
    //set also the default weights
    for(int i = 0; i < 4; i+=2) {
        wM[i] = (std::log10(Mtrack[i + 1]) - lmzams) / (std::log10(Mtrack[i + 1]) - std::log10(Mtrack[i]));
        wM[i+1] = (lmzams - std::log10(Mtrack[i])) / (std::log10(Mtrack[i + 1]) - std::log10(Mtrack[i]));
    }

}



/**Optional table properties**/
void OptionalTableProperty::set_refpointers(Star *s) {

    table_loaded=true; //Preset table_loaded to true, if some of the table are not loaded set table_loaded to false below
    for(int i = 0; i < 4; i++) {

        val_ref[i] = s->get_table_pointer_atpos(TabID(), i, 0);

        if (val_ref[i]==nullptr){
            V=V0=value=value0=std::nan("");
            table_loaded=false;
            break;
        }

    }

}

void OptionalTableProperty::update_derived(Star *s) {
    //If acting as derived properties
    if (amiderived()){
        //Here we want to evolve the derived properties, without setting V0
        //1-Save the V0 value
        double _V0= V0;
        //2-Evolve (depending on the stellar type)
        //First check if it is nakedCO (nakedCO is not a remnant)
        if(s->aminakedco()) evolve_nakedco(s);
        //Second, check if it is not a remnant , in this case just evolve
        else if(!s->amiremnant()) evolve(s);
        //Finally, if it is not empty, we call evolve remnant
        else if(!s->amiempty()) evolve_remnant(s);
        //Notice: do nothing if the star is empty
        //3-Restore the V0 value
        set_0(_V0);
    }
    //else use the common update_derived
    else{
        Property::update_derived(s);
    }

}

void SurfaceAbundanceTable::set_refpointers(Star *s) {

    //If we don't want to use the tables just set them to false
    if (!s->get_svpar_bool("tabuse_Xsup"))
        table_loaded=false;
    else
        OptionalTableProperty::set_refpointers(s);

    if (!table_loaded and s->get_svpar_bool("tabuse_Xsup"))
        svlog.critical("The table "+name()+" has not be loaded, but the option tabuse_envconv is set to true."
                                           "This table is needed if the option is enabled.",__FILE__,__LINE__,sevnstd::sanity_error());
}

void ConvectiveTable::set_refpointers(Star *s) {

    //If we don't want to use the tables just set them to false
    if (!s->get_svpar_bool("tabuse_envconv"))
        table_loaded=false;
    else
        OptionalTableProperty::set_refpointers(s);

    if (!table_loaded and s->get_svpar_bool("tabuse_envconv"))
        svlog.critical("The table "+name()+" has not be loaded, but the option tabuse_envconv is set to true."
                                           "This table is needed if the option is enabled.",__FILE__,__LINE__,sevnstd::sanity_error());

}

/**Mass object**/
void Mass_obejct::update_from_binary(Star* s, const double &DV, _UNUSED  Binstar *b) {

    if (std::isnan(DV) or std::isinf(DV))
        svlog.critical("Update from Binary of property " +
                       name()+"  is nan or infinite (Star ID:"+
                       utilities::n2s(s->get_ID(),__FILE__,__LINE__)+").",
                       __FILE__,__LINE__);

    V+=DV;

    //TODO problem: the mass lost in the processes can be larger than the total mass (cumulative)
    // In this case we set the mass to a tiny number. This is true for all the mass property.
    // Anyway, usually this is due to a sudden very high mass loss that is then correct by the
    // adaptive time-step criterium. Still left as a TO DO to understand if something better can be done
    if (V<=0 and std::abs(DV)>0) {
        svlog.warning("Update from Binary for Mass property is using a DV larger than the current Mass. The "
                      " star will be set to a tiny value. Notice that this is safely taken into account"
                      " if the update is from binary processes.", __FILE__, __LINE__);
        //s->set_empty();
        V = utilities::TINY;
    } else{

    }



    /*
    // This is just a temporary patch, if V_DV<0, we have to put this star empty.
    double V_tmp = V+DV;

    if (DV>V or std::isnan(DV) or std::isinf(DV))
        svlog.critical("Update from Binary of property " +
        name()+" is larger than the property value or is nan or infite",__FILE__,__LINE__);



    //utilities::wait(DV,__FILE__,__LINE__);
    V = V_tmp>0 ? V_tmp : 1e-2*V;
    svlog.debug("update_from_binary AFTER  " + name() + " "+ std::to_string(V),__FILE__,__LINE__);

    //if (s->getp(Phase::ID)==7 and s->getp(Worldtime::ID)>464.30)
        //utilities::wait(V,DV,s->getp(dMcumul_binary::ID)-s->getp_0(dMcumul_binary::ID),__FILE__,__LINE__);
    */
}


/**Mass**/
void Mass::set_remnant(Star *s) {

    V0 = V;
    value0 = value;

    V = s->get_staremnant()->get_Mremnant_at_born();
    value = V;
}

void Mass::changed_track(_UNUSED Star* s,Star* s_jtrack){
    //copy the fake tracks value  from s_jtrack
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}

void Mass::correct_interpolation_errors(_UNUSED Star *s)  {

    //1-Mass cannot increase in a timestep
    if ( (value>value0) && (value0>0.0) ){
        svlog.warning("Interpolation Warning: total mass cannot increase in "
        "single stellar evolution.\n   value0 is "+utilities::n2s(value0,__FILE__,__LINE__)
        +   ", value is "+utilities::n2s(value,__FILE__,__LINE__)+
        ".\n   value set to value0 and V to V0.",__FILE__,__LINE__);
        set_fk(value0);
        set(V0);
        utilities::wait();
    }


}

void Mass::correct_interpolation_errors_real( Star *s)  {


        if ( (V>V0) && (V0>0.0) ){
            svlog.warning("Interpolation Warning: total mass cannot increase in "
            "single stellar evolution.\n    V0 is "+utilities::n2s(V0,__FILE__,__LINE__)
            +", V is "+utilities::n2s(V,__FILE__,__LINE__)+
            ".\n    V set to V0.",__FILE__,__LINE__);
            set(V0);
        }

        //If I am a naked helium, I cannot go below the minimum HE mass set the first time the star develope a CO core
        //unless this was already decreased below this minimum by binary proceses (this is checked with get_0>=get_MHE_min(), if it is false the reset of the mass is not done)
        if (!std::isnan(s->get_MHE_min()) and s->aminakedhelium() and get_0()>=s->get_MHE_min()){
            set(get() < s->get_MHE_min() ?  s->get_MHE_min() : get() );
        }

}

/**Inertia**/
void Inertia::changed_track(Star *s, Star *s_jtrack) {


    if (table_loaded){
        // When we change traclk we are using the old mass but the new raidius
        // We assume that I propto MR so the new inertia propto Mold Rnew^2
        // while the new stars as I propto Mnew Rnew^2, therefore we rescale for a factor Mold/Mnew

        //GI 23/03, I found that sometime when the new Inertia is small (close to 1) there could be
        //a problem in evolving the Inertia with a mismatch between V,VO and value,value0.
        //THis is due to the fact that internally we use a logInertia and when Inertia<1, logInertia<0
        //It could happens that sometime DValue is negative but also V0 is negative and as a result
        //the inertia will increase instead of decrease.
        //In order to avoid this problem, we set V to Vjtrack if Inertia is smaller than 1.5
        if (s_jtrack->getp(ID)<1.5){
            V=std::log10(s_jtrack->getp(ID));
        } else{
            double weight_factor=s->getp(Mass::ID)/s_jtrack->getp(Mass::ID);
            weight_factor = fabs(weight_factor-1) < 0.1 ? 1. : weight_factor; // If the masses are very similar do not use this factor
            V=std::log10(weight_factor*s_jtrack->getp(ID)); //Remember get is 10**V that is stored in log
        }
    } else {
        V=estimate_logInertia(s); //Here we have to estimate again the inertia on s using the updated value of Masses and Radius, already updated.
    }


    value=log10(s_jtrack->getp_fk(ID));
    value0=log10(s_jtrack->getp_fk0(ID));

}
void Inertia::set_remnant(Star *s) {

    V0 = V;
    value0 = value;

    V = log10(s->get_staremnant(Inertia::ID));
    value = log10(s->get_staremnant(Inertia::ID));

}
void Inertia::evolve_remnant(Star *s) {
    V0 = V;
    V = log10(s->get_staremnant(Inertia::ID));
}

void Inertia::evolve_nakedco(Star *s) {
    V0 = V;
    V = log10(estimate_Inertia_homogeneous_sphere(s));
}

void Inertia::set_refpointers(Star *s) {

    //If we don't want to use the tables just set them to false
    if (!s->get_svpar_bool("tabuse_inertia")){
        table_loaded=false;
        //Set the function to be used
        if (s->get_svpar_str("inertiamode")=="hsphere")
            inertia_func = &Inertia::estimate_Inertia_homogeneous_sphere;
        else if (s->get_svpar_str("inertiamode")=="hspherecore")
            inertia_func = &Inertia::estimate_Inertia_homogeneous_sphere_wcore;
        else if (s->get_svpar_str("inertiamode")=="Hurley")
            inertia_func = &Inertia::estimate_Inertia_Hurley;
        else if (s->get_svpar_str("inertiamode")=="DeMink"){
            inertia_func = &Inertia::estimate_Inertia_DeMink;
        }
        else
            svlog.critical("inertiamode "+s->get_svpar_str("inertiamode")+" not allowed.",__FILE__,__LINE__,sevnstd::params_error());
    }
    else
        OptionalTableProperty::set_refpointers(s);

    if (!table_loaded and s->get_svpar_bool("tabuse_inertia"))
        svlog.critical("The table "+name()+" has not be loaded, but the option tabuse_inertia is set to true."
                                           "This table is needed if the option is enabled.",__FILE__,__LINE__,sevnstd::sanity_error());
}

void Inertia::correct_interpolation_errors_real(_UNUSED Star *s){

    if (get()<=1e-6){
        V=-6; //Log10(0.1), set a minimum value
    }

}


void Inertia::evolve_without_table(Star *s) {
    set_0(std::log10(get()));
    set( estimate_logInertia(s) );
    //Set V,V0 but also v, v0. This is needed because even if the tables are loaded
    //the user can choice do not use them.
    synch_v_value_without_table();
}
double Inertia::estimate_logInertia(Star *s) {
    return std::log10((*this.*inertia_func)(s));
}

double Inertia::estimate_Inertia_homogeneous_sphere(double Mass,double Outer_radius, double Inner_radius){


    if (Outer_radius==0 and Inner_radius==0)
        return 0.0;
    else if (Inner_radius>=Outer_radius)
        svlog.critical("In estimate_Inertia_homogeneous_sphere Inner radius ("+
        utilities::n2s(Inner_radius,__FILE__,__LINE__)+") is equal or larger  than Outer radius ("+utilities::n2s(Outer_radius,__FILE__,__LINE__)+")");


    double den =  Outer_radius*Outer_radius*Outer_radius; //Rout^3
    double num =  den*Outer_radius*Outer_radius; //Rout^5

    if (Inner_radius>0){
        double Rin3=Inner_radius*Inner_radius*Inner_radius;
        num = num - Rin3*Inner_radius*Inner_radius; //Rout^5-Rin^5
        den = den - Rin3; //Rout^3-Rin^3
    }


    return 0.4*Mass*num/den; //0.4 is 2/5 for a homogeneous sphere
}
double Inertia::estimate_Inertia_homogeneous_sphere(Star *s){
    return estimate_Inertia_homogeneous_sphere(s->getp(Mass::ID),s->getp(Radius::ID));
}
double Inertia::estimate_Inertia_homogeneous_sphere_wcore(Star *s){

    double inertia_core     = estimate_Inertia_homogeneous_sphere(s->Mcore(),s->Rcore());
    double inertia_envelope = estimate_Inertia_homogeneous_sphere(s->Menvelope(),s->getp(Radius::ID),s->Rcore());
    return inertia_core + inertia_envelope;
}

//TODO This is taken directly from SEVN1, but in bse k2 depends on stellar phase
double Inertia::estimate_Inertia_Hurley(_UNUSED Star *s){

    double k2=0.1, k3=0.21;
    double ifac = 2.5; //Factor to eliminate the 2/5 factor used inside estimate_Inertia_homogeneous_sphere

    double inertia_core     = k3*ifac*estimate_Inertia_homogeneous_sphere(s->Mcore(),s->Rcore());
    double inertia_envelope = k2*ifac*estimate_Inertia_homogeneous_sphere(s->Menvelope(),s->getp(Radius::ID),s->Rcore());

    return inertia_core + inertia_envelope;

}
//TODO This is taken directly from SEVN1, but in DeMink13 it seems that this formalism is valid only for MS
double Inertia::estimate_Inertia_DeMink(Star *s){

    const double Radius = s->getp(Radius::ID);
    //TODO It is not clear if after a jump we should use mzams-rzams or mzams0-rzmas0
    const double mzams = s->get_zams();
    const double rzams = s->get_rzams();

    double logmz = std::log10(mzams);
    double logm = std::log10(s->getp(Mass::ID));
    //Estimate k0, Eq. A1
    double k = std::min(0.21, std::max(0.09-0.27*logmz, 0.037+0.033*logmz));
    //correct for logmz>1/3
    if(logmz > 1.3) k-= 0.055*(logmz -1.3)*(logmz -1.3);

    //Estimate exponent C, Eq. A3
    //TODO In the paper is not clear if in this part we should use logm or logmz, in the old SEVN verison we used logmz, but in a  meeting (04/06/21) we decided to use logm
    double C;
    if (logm<0){
        C=-2.5;
    }
    else if(logm > 0 and logm < 0.2){
        C=-2.5+5*logm;
    }
    else{
        C=-1.5;
    }

    //Estimate final k, Eq. A2
    k = (k - 0.025)*std::pow(Radius/rzams, C) +  0.025 * pow(Radius/rzams, -0.1);

    return k*s->getp(Mass::ID)*Radius*Radius;

}

void Inertia::set_wM(Star *s) {
    std::string wM_option = s->get_svpar_str("ev_setwM_log");

    if (wM_option=="linear"){
        set_wM_linear(s);
    } else if (wM_option=="rational"){
        set_wM_rational(s);
    } else if (wM_option=="log"){
        set_wM_log(s);
    } else{
        svlog.critical("Unknown set M weights options for property "+
                       name(),__FILE__,__LINE__,sevnstd::sanity_error(""));
    }
}


/**Luminosity**/
void Luminosity::changed_track(_UNUSED Star* s,Star* s_jtrack){

    //Use the raw V,value,value0 to avoid to lose precision thorugh the log10 conversion
    V=s_jtrack->getp_raw(ID);
    value=s_jtrack->getp_fk_raw(ID);
    value0=s_jtrack->getp_fk0_raw(ID);

    ///Copy the Luminosity from the jtrack star
    //V=log10(s_jtrack->getp(ID)); //Remember get is 10**V that is stored in log
    //V0=log10(s_jtrack->getp(ID));
    //value=log10(s_jtrack->getp_fk(ID));
    //value0=log10(s_jtrack->getp_fk0(ID));
    //svlog.pdebug("Radius changed track", get(), get_0(), value, value0,__FILE__,__LINE__);
    //utilities::wait();
}
void Luminosity::set_remnant(Star *s) {
        V0 = V;
        value0 = value;

        V = log10(s->get_staremnant(Luminosity::ID));
        value = log10(s->get_staremnant(Luminosity::ID));
}
void Luminosity::evolve_remnant(Star *s) {
    V0 = V;
    V = log10(s->get_staremnant(Luminosity::ID));
}

void Luminosity::set_wM(Star *s) {
    std::string wM_option = s->get_svpar_str("ev_setwM_log");

    if (wM_option=="linear"){
        set_wM_linear(s);
    } else if (wM_option=="rational"){
        set_wM_rational(s);
    } else if (wM_option=="log"){
        set_wM_log(s);
    } else{
        svlog.critical("Unknown set M weights options for property "+
                       name(),__FILE__,__LINE__,sevnstd::sanity_error(""));
    }
}

/**Radius**/
void Radius::evolve(Star *s){

    if (s->amifollowingQHE()){
        set_0(get());
    }
    else{
        Property::evolve(s);
    }


    return;
}

void Radius::evolve_fake(Star *s){

    //If ev_R_linear_interp normal evolve_fake using the values stored in log
    if (!s->get_svpar_bool("ev_R_linear_interp")){
        Property::evolve_fake(s);
    }
    //Otherwise, transform the store value to lin and perform the linear interpolation
    else{
        ///Set v0=v
        set_0_fk(value); //save the previous value, before the evolution

        ///Move to new value in the interpolating track
        FOR4 {
            val_in[_i] = val_ref[_i] + s->pos[_i];
        }

        ///Find new values at times t for the interpolating tracks using linear interpolation
        for(int i = 0; i < 4; i++) {
            double *_val = val_in[i];
            double *_time = s->times_in[i];

            //Key difference!! Transform the stored values to linear values
            double _val_lin[2]={std::pow(10,*_val),std::pow(10,*(_val+1))};

            double slope = (_val_lin[1] - _val_lin[0]) / ((*(_time +1)) - *_time);
            double intercept = _val_lin[0] - slope * (*_time);

            interpolating_values[i] = slope * s->ttimes[i] + intercept; //ttimes come from function tracktimes()
        }

        ///Find the interpolating z value, M1*v1 + M2*v2
        double val_zlow = interpolating_values[0]*wM[0] + interpolating_values[1]*wM[1];
        double val_zhigh = interpolating_values[2]*wM[2] + interpolating_values[3]*wM[3];

        ///Find the final interpolated value, v=Z1*v1 + Z2*v2
        //Key difference!! Transform back the interpolated value in log
        set_fk(std::log10(wZ[0]*val_zlow + wZ[1]*val_zhigh));

    }

    return;
}



void Radius::changed_track(_UNUSED Star* s,Star* s_jtrack){

    //If we are followng the QHE does not change Radius
    if (s->amifollowingQHE())
        return;

    //Use the raw V,value,value0 to avoid to lose precision thorugh the log10 conversion
    V=s_jtrack->getp_raw(ID);
    value=s_jtrack->getp_fk_raw(ID);
    value0=s_jtrack->getp_fk0_raw(ID);

    //V=log10(s_jtrack->getp(ID)); //Remember get is 10**V that is stored in log
    //V0=log10(s_jtrack->getp(ID));
    //value=log10(s_jtrack->getp_fk(ID));
    //value0=log10(s_jtrack->getp_fk0(ID));
    //svlog.pdebug("Radius changed track", get(), get_0(), value, value0,__FILE__,__LINE__);
    //utilities::wait();
}
void Radius::set_remnant(Star *s) {

    V0 = V;
    value0 = value;

    V = log10(s->get_staremnant(Radius::ID));
    value = log10(s->get_staremnant(Radius::ID));

}
void Radius::evolve_remnant(Star *s) {
    V0 = V;
    V = log10(s->get_staremnant(Radius::ID));
}

void Radius::set_wM(Star *s) {

    std::string wM_option;

    //If the ev_R_linear_interp is true we are using linear values to interpolate the radius
    //so set the weights considering the setwM option
    if (s->get_svpar_bool("ev_R_linear_interp")){
        wM_option = s->get_svpar_str("ev_setwM");
    }
        //while if ev_R_linear_interp is false, use the weights for log quantities
    else{
        wM_option = s->get_svpar_str("ev_setwM_log");
    }


    if (wM_option=="linear"){
        set_wM_linear(s);
    } else if (wM_option=="rational"){
        set_wM_rational(s);
    } else if (wM_option=="log"){
        set_wM_log(s);
    } else{
        svlog.critical("Unknown set M weights options for property "+
                       name(),__FILE__,__LINE__,sevnstd::sanity_error(""));
    }

}

void Radius::update_from_binary(_UNUSED Star *s,  const double &DV, _UNUSED Binstar *b) {

    if (std::isnan(DV) or std::isinf(DV))
        svlog.critical("Update from Binary of property " +
                       name()+" is nan or infinite",__FILE__,__LINE__);

    // The Radius is internally stored in log, so transform DV to log before updating it
    V = std::log10(pow(10,V) + DV);
}


/**Phase**/
//this property evolves outside the main loop
void Phase::special_evolve(Star *s){

    const double *tphase = s->get_tphase();

    //cout<<" OLD phase = "<<get()<<"  subphase = "<<get_fk()<<endl;
    set_0(get());
    set_0_fk(get_0_fk());

    set(-1);
    double localtime = s->getp(Localtime::ID);

//TODO for WR stars maybe we should refer to different phases (same for pureHE tracks)
//if(amiwr()) then use the tphaseWR vector coupled with the PhasesWR enum

    svlog.pdebug("Local time in Phase evolve ",localtime,V0,std::setprecision(20),__FILE__,__LINE__);
    svlog.pdebug(std::setprecision(20),"Tphase "+std::to_string(tphase[0]),std::setprecision(20),__FILE__,__LINE__);
    svlog.pdebug("Check "+std::to_string(localtime < tphase[0]),std::setprecision(20),__FILE__,__LINE__);
    svlog.pdebug(std::setprecision(20),"Local time in Phase evolve ",localtime,
    "Tphase "+std::to_string(tphase[0]),"Check "+std::to_string(localtime < tphase[0]),std::setprecision(20),__FILE__,__LINE__);

    if(localtime < tphase[0])
        svlog.critical("Initial time is likely negative (Star ID " + utilities::n2s(s->get_ID(),__FILE__,__LINE__) + ")",__FILE__,__LINE__);

    for(int i = 0; i < Nphases-1; i++) {
        svlog.pdebug(std::setprecision(20),"LTRRRRR",localtime,tphase[i],tphase[i+1],s->get_zams(),__FILE__,__LINE__);
        if (localtime < tphase[i + 1] && localtime >= tphase[i])
            set(i);
    }

    //TODO the next two line are for debug remove them
    //if (s->get_ID()==0 and s->getp(Worldtime::ID)!=0)
    //    std::cout<<std::setprecision(20)<<s->getp(Worldtime::ID)<<" "<<localtime<<" "<<tphase[Nphases-1]<<std::endl;
    if(localtime >= tphase[Nphases-1]){
        //set(Nphases-1);
        set(Lookup::Phases::Remnant);
    }


    if(get() == -1)
        svlog.critical("Cannot initialize the phase of the star", __FILE__, __LINE__);

    if(tphase[(int)get()] == -1)
        svlog.critical("Cannot set the phase of the star, tphase not set. Interpolating track with different phases? Check in the look-up tables",
                     __FILE__, __LINE__);

}
void Phase::changed_track(_UNUSED Star *s, Star *s_jtrack){
    V=s_jtrack->getp(Phase::ID);
}

void Phase::set_wM(Star *s) {
    std::string wM_option = s->get_svpar_str("ev_setwM_tphase");

    if (wM_option=="linear"){
        set_wM_linear(s);
    } else if (wM_option=="rational"){
        set_wM_rational(s);
    } else if (wM_option=="log"){
        set_wM_log(s);
    } else{
        svlog.critical("Unknown set M weights options for property "+
                       name(),__FILE__,__LINE__,sevnstd::sanity_error(""));
    }
}

//Notice here we don't do nothing because as for Localtime and  Worldtime this property has been already evolved in the special evolve
//void Phase::set_remnant(Star *s) {
    //return;
    //Notice here we don't do nothing because as for Localtime and  Worldtime this property has been already evolved in the special evolve

    //V0 = V;  Notice here we do not use
    //value0 = value;
    //V = Lookup::Phases::Remnant;
    //value = Lookup::Phases::Remnant;
//}

/**Localtime**/
void Localtime::special_evolve(Star *s){
    svlog.debug("Old Localtime T_0= "+utilities::n2s(get_0(),__FILE__,__LINE__),__FILE__,__LINE__);
    svlog.debug("Old Localtime T= "+utilities::n2s(get(),__FILE__,__LINE__),__FILE__,__LINE__);

    set_0(get());
    set(get() + s->getp(Timestep::ID));

    svlog.debug("New Localtime T_0= "+utilities::n2s(get_0(),__FILE__,__LINE__),__FILE__,__LINE__);
    svlog.debug("New Localtime T= "+utilities::n2s(get(),__FILE__,__LINE__),__FILE__,__LINE__);
}
void Localtime::changed_track(_UNUSED Star* s,Star* s_jtrack){
    V=V0=value=value0=s_jtrack->getp(ID);
}

/**Temperature**/
void Temperature::evolve(Star *s){
    set_0(get());

    double Teff4=s->getp(Luminosity::ID)/(4*M_PI*s->getp(Radius::ID)*s->getp(Radius::ID)*utilities::Sigma_StefBoltz);
    set(pow(Teff4,0.25));
}
void Temperature::changed_track(_UNUSED Star *s, Star *s_jtrack) {
    //V=s_jtrack->getp(ID);
    //V0=s_jtrack->getp(ID);
    V = s_jtrack->getp(Luminosity::ID)/(4*M_PI*s_jtrack->getp(Radius::ID)*s_jtrack->getp(Radius::ID)*utilities::Sigma_StefBoltz);
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}

/*Mass objects*/
double Mass_obejct::get_Vlast(const Star *s) const {

    double last_values [4]; //last values of the 4 interpolating track

    for(int i = 0; i < 4; i++){
        last_values[i] = s->get_last_fromtable((Lookup::Tables)TabID(),i);
    }

    double val_zlow = last_values[0]*wM[0] + last_values[1]*wM[1];
    double val_zhigh = last_values[2]*wM[2] + last_values[3]*wM[3];

    ///Find the final interpolated value, v=Z1*v1 + Z2*v2
   return wZ[0]*val_zlow + wZ[1]*val_zhigh;

}

/**MCO**/
void MCO::evolve(Star *s){

    if(s->getp(Phase::ID) >= TerminalCoreHeBurning)
        Property::evolve(s);
    else
        return;
}

void MCO::correct_interpolation_errors(Star *s) {
    double v =  get_fk();
    double v0 = get_0_fk();
    //1-In some rare case (so far only on a Mac laptop with GCC-11 and O3 optimization)
    //numeric errors on the interpolations return return values smaller than 0
    //In this case if the Star is the Main Sequence just set 0
    //otherwise an arbitrary small number
    if (get_fk()<0){
        if (s->getp(Phase::ID) >= TerminalCoreHeBurning) set_fk(1E-20);
        else set_fk(0);
    }
    //2-CO mass cannot decrease
    set_fk((v < v0) ? v0 : v);
}


void MCO::correct_interpolation_errors_real(_UNUSED Star *s)  {

    double v =  get();
    double v0 = get_0();


    ///1-CO cannot decrease
    set((v < v0) ? v0 : v); //CO mass cannot decrease

    //2-CO cannot grow more than maximum CO (if set)
    //unless some binary processes give to the star more MCO than the maximum.
    //This condition is checked looking at v0 (that is the last value after a complete bse step, is the star is a binary).
    //If this value is  larger tha maximum CO, the reset to the maximum CO is not done (see next if).
    if (!std::isnan(s->get_MCO_max()) and v0<=s->get_MCO_max()){
            set(v > s->get_MCO_max() ?  s->get_MCO_max() : v);
    }
    //If the MCO was already decreased above  the maximumjust set it to this value (the one after the last step of BSE)
    //do not allow to become larger due to SSE
    else if(!std::isnan(s->get_MCO_max())){
        set(v > s->get_MCO_max() ?  v0 : v);
    }


    ///3-MCO cannot be larger than MHE or total Mass
    //If MCO is larger than the Mass in the current timestep there was a point in which the decreasing Mass
    //and the increasing  MCO cross each other. Since the interpolation is linear, we can easily find this point
    //Then we set Mass=MHE=MCO=Mcross.
    //TODO In order to be totally consistent we shoud trigger a repeat setting the newtimestep as the one to reach the crossig point.
    if (get()>s->getp(Mass::ID)){
        double m_mass, m_mco, b_mass, b_mco;
        utilities::find_line(0,s->getp(Timestep::ID), s->getp_0(Mass::ID),s->getp(Mass::ID),m_mass,b_mass);
        utilities::find_line(0,s->getp(Timestep::ID), get_0(),get(),m_mco,b_mco);
        double t_match = (b_mco-b_mass)/(m_mass-m_mco);
        double M_match = m_mass*t_match+b_mass;
        set(M_match);
        s->copy_property_from(Mass::ID,MCO::ID);
        //The following is a trick to don't have exactly Mass=MHE (this will force the radius to be equal with a sequent likely repeat due to large Radius variation)
        M_match=M_match-std::min(1e-5,0.1*s->get_svpar_num("ev_naked_tshold"));
        set(M_match);
        s->copy_property_from(MHE::ID,MCO::ID);
        //The following is a trick to don't have exactly Mass=MHE (this will force the radius to be equal with a sequent likely repeat due to large Radius variation)
        M_match=M_match-std::min(1e-5,0.1*s->get_svpar_num("ev_naked_tshold"));
        set(M_match);
        //If it is a pureHe set also the new MHE mass
        if (s->amifollowingpureHetrack()){
            s->copy_property_from(MHE::ID,Mass::ID);
        } else{
            s->copy_property_from(MHE::ID,MCO::ID);
            //The following is a trick to don't have exactly MCO=MHE (this will force the radius to be equal with a sequent likely repeat due to large Radius variation)
            M_match=M_match-std::min(1e-5,0.1*s->get_svpar_num("ev_naked_tshold"));
            set(M_match);
        }
    }
    else if(get()>s->getp(MHE::ID)){
        double m_mhe, m_mco, b_mhe, b_mco;
        utilities::find_line(0,s->getp(Timestep::ID), s->getp_0(MHE::ID),s->getp(MHE::ID),m_mhe,b_mhe);
        utilities::find_line(0,s->getp(Timestep::ID), get_0(),get(),m_mco,b_mco);
        double t_match = (b_mco-b_mhe)/(m_mhe-m_mco);
        double M_match = m_mhe*t_match+b_mhe;
        set(M_match);
        s->copy_property_from(MHE::ID,MCO::ID);
        set(M_match-std::min(1e-5,0.1*s->get_svpar_num("ev_naked_tshold")));
    }

}



void MCO::changed_track(_UNUSED Star *s, Star *s_jtrack) {

    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
    //TODO THis is just a patch to avoid the problem in which we wrongly have an overmassive or undermassive core
    //if (s->getp(Phase::ID) == TerminalCoreHeBurning and V<1e-5)
    //    V=  s_jtrack->getp(ID)>0.0 ? s_jtrack->getp(ID) : 1e-7;


}


/**MHE**/
void MHE::evolve(Star *s){

    //Check if it is naked helium
    //PureHe tracks from Parsec have Mtot and R that are slighly different wrt MHE and RHE, this becasue there is
    //an Helium atmosphere. In order to simplify the pureHE star we force the evolution to just follow the total Mass evolution
    if (s->amifollowingpureHetrack()){
        set_0(get());
        set_0_fk(get_0_fk());
        set_fk(s->getp_fk(Mass::ID));
        set(s->getp(Mass::ID));
        return;
    }
    else if (s->getp(Phase::ID) >= TerminalMainSequence)
        Property::evolve(s);
    else
        return;

    //cout<<" NEW MHE = "<<get()<<endl;
}

void MHE::correct_interpolation_errors(Star *s) {
    //In some rare case (so far only on a Mac laptop with GCC-11 and O3 optimization)
    //numeric errors on the interpolations return return values smaller than 0
    //In this case if the Star is the Main Sequence just set 0
    //otherwise an arbitrary small number
    if (get_fk()<0){
        if (s->getp(Phase::ID)>=TerminalMainSequence) set_fk(1E-20);
        else set_fk(0);
    }
}


void MHE::correct_interpolation_errors_real(Star *s){


    //I cannot go below the minimum HE mass set the first time the star develope a CO core
    //unless this was already decreased below this minimum by binary proceses (this is checked with get_0>=get_MHE_min(), if it is false the reset of the mass is not done)
    //see next else if
    if (!std::isnan(s->get_MHE_min()) and get_0()>=s->get_MHE_min()){
        set(get() < s->get_MHE_min() ?  s->get_MHE_min() : get() );
    }
    //TODO Notice that using it we essentially disable the pureHE winds if the star that becomes a pureHE start with a MHE smaller than the MHE_min
    //If the MHE was already decreased below the minimum by the binary process, just set it to this value (the one after the last step of BSE)
    //do not allow to become even smaller
    else if (!std::isnan(s->get_MHE_min())){
        set(get() < s->get_MHE_min() ?  get_0() : get() );
    }

    ///Now compare with Mass
    //If MHE is larger than the Mass in the current timestep there was a point in which the decreasing Mass
    //and the increasing or decreasing MHE cross each other. Since the interpolation is linear, we can easily find this point
    //Then we set Mass=MHE=Mcross.
    //TODO In order to be totally consistent we shoud triggere a repeat setting the newtimestep as the one to reach the crossig point.
    if (get()>s->getp(Mass::ID)){
        double m_mass, m_mhe, b_mass, b_mhe;
        utilities::find_line(0,s->getp(Timestep::ID), s->getp_0(Mass::ID),s->getp(Mass::ID),m_mass,b_mass);
        utilities::find_line(0,s->getp(Timestep::ID), get_0(),get(),m_mhe,b_mhe);
        double t_match = (b_mhe-b_mass)/(m_mass-m_mhe);
        double M_match = m_mass*t_match+b_mass;
        set(M_match);
        s->copy_property_from(Mass::ID,MHE::ID);
        //Set real M_match
        //The following is a trick to don't have exactly Mass=MHE (this will force the radius to be equal with a sequent likely repeat due to large Radius variation)
        set(M_match-std::min(1e-5,0.1*s->get_svpar_num("ev_naked_tshold")));
    }


}
void MHE::changed_track(Star *s, Star *s_jtrack) {

    //If Naked Helium, the MHE has to follow the total mass
    if (s->amifollowingpureHetrack())
        V=s->getp(Mass::ID);
    //TODO THis is just a patch to avoid the problem in which we wrongly have an overmassive or undermassive core
    //else if (s->getp(Phase::ID) == TerminalMainSequence and V<1e-5)
    //    V=s_jtrack->getp(ID)>0.0 ? s_jtrack->getp(ID) : 1e-7;


    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);



}


/**RHE**/
void RHE::evolve(Star *s){

    //First think, check if it is naked helium. In this case RHE is just equal to the total Radius
    if (s->amifollowingpureHetrack()){
        set_0(get());
        set_0_fk(get_0_fk());
        set_fk(s->getp_fk(Radius::ID));
        set(s->getp(Radius::ID));
        return;
    }
    //second check if He core actually exists and table are loaded(Phase larger than Terminal Main Sequence)
    else if(s->getp(Phase::ID) >= TerminalMainSequence and table_loaded){
        Property::evolve(s);
    }
    //third  if He core actually exists and but table are not loaded evolve using functional form
    else if(s->getp(Phase::ID) >= TerminalMainSequence and !table_loaded){
        evolve_without_table(s);
    }
    //fourth in all the other cases do nothing
    else{
        return;
    }
    return;
}

void RHE::correct_interpolation_errors(Star *s) {
    //In some rare case (so far only on a Mac laptop with GCC-11 and O3 optimization)
    //numeric errors on the interpolations return return values smaller than 0
    //In this case if the Star is the Main Sequence just set 0
    //otherwise an arbitrary small number
    if (get_fk()<0){
        if (s->getp(Phase::ID)>=TerminalMainSequence) set_fk(1E-20);
        else set_fk(0);
    }

    //Check if RHE>Radius from interpolating tracks
    //THis can happen sometime because R is stored and interpolated in log, while RHE and RCO use the linear scale.
    if (s->getp(MHE::ID)<s->getp(Mass::ID) and s->getp_fk(Radius::ID)<=s->getp_fk(RHE::ID)){
        set_fk(0.9999*s->getp_fk(Radius::ID));
    }




}

void RHE::correct_interpolation_errors_real(Star *s){

    if (s->getp(MHE::ID)<s->getp(Mass::ID)){
        //RHE cannot be larger than R, for the real star
        set(get() > 0.9999*s->getp(Radius::ID) ? 0.9999*s->getp(Radius::ID) : get());
    }
    else if (s->getp(MHE::ID)==s->getp(Mass::ID) and !s->aminakedhelium()){
        //This is (will become) a naked Helium set the Radius to the RHE
        s->copy_property_from(Radius::ID,RHE::ID);
    }

}
void RHE::changed_track(Star *s, Star *s_jtrack) {

    //If Naked Helium, the RHE has to follow the Radius
    if (s->amifollowingpureHetrack())
        V=s->getp(Radius::ID);
    else if (table_loaded)
        V=s_jtrack->getp(ID);
    else
        V=estimate_Rcore(s); //Use the updated MHE mass

    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);

}
void RHE::evolve_without_table(Star *s){

    set_0(get());
    set( estimate_Rcore(s) );
    //Set V,V0 but also v, v0. This is needed because even if the tables are loaded
    //the user can choice to not use them.
    synch_v_value_without_table();

    return;
}
void RHE::set_refpointers(Star *s) {
    //If we don't want to use the tables just set them to false
    if (!s->get_svpar_bool("tabuse_rhe"))
        table_loaded=false;
    else
        OptionalTableProperty::set_refpointers(s);

    if (!table_loaded and s->get_svpar_bool("tabuse_rhe"))
        svlog.critical("The table "+name()+" has not be loaded, but the option tabuse_rhe is set to true."
                                           "This table is needed if the option is enabled.",__FILE__,__LINE__,sevnstd::sanity_error());
}
double RHE::estimate_Rcore(Star *s) {

    double best_scale_factor=0.10746657; //RHE at MHE=1 Msun from a fit done on the SEVNtracks_G tables  (see Wiki)
    double Rcore_from_func = CoreRadius::estimate_Rcore(s->getp(MHE::ID),best_scale_factor);

    //Do not allow Rcore larger than total Radius
    return Rcore_from_func >= s->getp(Radius::ID) ? 0.99*s->getp(Radius::ID) : Rcore_from_func;
}

/**RCO**/
void RCO::evolve(Star *s){

    if(s->getp(Phase::ID) >= TerminalCoreHeBurning and table_loaded)
        Property::evolve(s);
    else if(s->getp(Phase::ID) >= TerminalCoreHeBurning and !table_loaded)
        evolve_without_table(s);
    else
        return;

}

void RCO::correct_interpolation_errors(Star *s) {
    //In some rare case (so far only on a Mac laptop with GCC-11 and O3 optimization)
    //numeric errors on the interpolations return return values smaller than 0
    //In this case if the Star is the Main Sequence just set 0
    //otherwise an arbitrary small number
    if (get_fk()<0){
        if (s->getp(Phase::ID) >= TerminalCoreHeBurning) set_fk(1E-20);
        else set_fk(0);
    }

    //Check if RCO>RHE from interpolating tracks
    //THis can happen sometime because RHE can be modified because it is larger than the Radius
    // since R is stored and interpolated in log, while RHE and RCO use the linear scale.
    if (s->getp(MCO::ID)<s->getp(MHE::ID) and s->getp_fk(RHE::ID)<=s->getp_fk(RCO::ID)){
        set_fk(0.9999*s->getp_fk(RHE::ID));
    }


}

void RCO::correct_interpolation_errors_real(Star *s){

    //RCO cannot be larger than RHE, for the real star
    if (get()>0.9999*s->getp(Radius::ID))
        set(0.9999*s->getp(Radius::ID));
    else if (get()>s->getp(RHE::ID))
        set(0.9999*s->getp(RHE::ID));

    //TODO Notice that setting the radius to the RCO radius in this phase we prevent possible RLOs (the star is set to naked CO only at the end)
    //If MCO==MHE the radii have to be the same, we set the radius to the radius of the CO (if the star is a pureHE the Radius is set to RCO too)
    //In any case this star will become a nakedCO star
    if (s->getp(MCO::ID)==s->getp(MHE::ID)){
        s->copy_property_from(RHE::ID,RCO::ID);
        if (s->aminakedhelium())  s->copy_property_from(Radius::ID,RCO::ID);
    }



    //Sometimes During the initial stage of the TerminalHecore, CCO is 0, in that case we force it to be not zero so that it is consistent we
    //other section of the code
    //if (s->getp(Phase::ID) == TerminalCoreHeBurning and V==0.0)
    //    V=1e-7;

}
void RCO::changed_track(_UNUSED Star *s, Star *s_jtrack) {

    if (table_loaded)
        V=s_jtrack->getp(ID);
    else
        V=estimate_Rcore(s); //Use the updated MHE mass

    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}
void RCO::set_refpointers(Star *s) {
    //If we don't want to use the tables just set them to false
    if (!s->get_svpar_bool("tabuse_rco"))
        table_loaded=false;
    else
        OptionalTableProperty::set_refpointers(s);

    if (!table_loaded and s->get_svpar_bool("tabuse_rco"))
        svlog.critical("The table "+name()+" has not be loaded, but the option tabuse_rco is set to true."
                                           "This table is needed if the option is enabled.",__FILE__,__LINE__,sevnstd::sanity_error());
}
void RCO::evolve_without_table(Star *s){

    set_0(get());
    set( estimate_Rcore(s) );
    //Set V,V0 but also v, v0. This is needed because even if the tables are loaded
    //the user can choice to not use them.
    synch_v_value_without_table();

    return;
}
double RCO::estimate_Rcore(Star *s) {

    double best_scale_factor=0.04149729; //RCO at MCO=1 Msun from a fit done on the SEVNtracks_G tables  (see Wiki)
    double Rcore_from_func = CoreRadius::estimate_Rcore(s->getp(MCO::ID),best_scale_factor);

    //Do not allow Rcore larger than the RHE radius
    return Rcore_from_func >= s->getp(RHE::ID) ? 0.99*s->getp(RHE::ID) : Rcore_from_func;
}



/**Bmag**/
void Bmag::evolve(_UNUSED Star *s) {
    set_0(get());
}
void Bmag::set_remnant(Star *s) {
    V0 = V;
    value0 = value;

    V = s->get_staremnant(Bmag::ID);
    value = s->get_staremnant(Bmag::ID);
}
void Bmag::evolve_remnant(Star *s) {
    set_0(get());
    set(s->get_staremnant(Bmag::ID));
}

/**OmegaRem**/
void OmegaRem::set_remnant(Star *s) {
    V0 = V;
    value0 = value;

    V = s->get_staremnant(OmegaRem::ID);
    value = s->get_staremnant(OmegaRem::ID);
}
void OmegaRem::evolve_remnant(Star *s) {
    set_0(get());
    set(s->get_staremnant(OmegaRem::ID));
}




/**Rs**/
void Rs::evolve( Star *s){

    set_0(get());
    set(utilities::R_Schwarzschild(s->getp(Mass::ID)));
    //set(rs);
}

/**Spin**/
void Spin::evolve(Star *s) {
    set_0(get());
    set(Spin_from_OmegaSpin(s->getp(OmegaSpin::ID),s->getp(Mass::ID),s->getp(Radius::ID)));
}


double Spin::Spin_from_OmegaSpin(double OmegaSpin, double Mass, double Radius) {
    return OmegaSpin/utilities::omega_crit(Mass,Radius);
}


/**Worldtime**/
void Worldtime::special_evolve(Star *s){

    svlog.debug("Old Worldtime T_0= "+utilities::n2s(get_0(),__FILE__,__LINE__),__FILE__,__LINE__);
    svlog.debug("Old Worldtime T= "+utilities::n2s(get(),__FILE__,__LINE__),__FILE__,__LINE__);

    set_0(get());
    set(get() + s->getp(Timestep::ID));

    svlog.debug("New Worldtime T_0= "+utilities::n2s(get_0(),__FILE__,__LINE__),__FILE__,__LINE__);
    svlog.debug("New Worldtime T= "+utilities::n2s(get(),__FILE__,__LINE__),__FILE__,__LINE__);

}

/**NextOutput**/
void NextOutput::special_evolve(Star *s){
    set_0(get());

    //Condition to avoid that the get() + s remains always lower than the current worldtime.
    get() + s->get_dtout()>s->getp(Worldtime::ID) ? set(get() + s->get_dtout()) : set(s->getp(Worldtime::ID) + s->get_dtout());
}

/**Timestep**/
void Timestep::evolve(Star *s){

    set_0(get());
    double max_variation = s->get_svpar_num("ts_maximum_variation");

    //bool dprint=s->getp(Worldtime::ID)>1.686631 and s->get_ID()==0;
    bool  dprint=false;

    s->repeatstep = false; //Reset flag repeatstep

    //GI 09/02/2024 prototype
    // If the timestep is tiny two stuff:
    // - never repeat
    // - next timestep is always a fixed value to avoid that the very tiny timestep produce numerical instability on the
    //      predictions of the next timestep:
    if (get_0()<=utilities::TINY){
        double new_timestep = 1e-7;
        //Check and modify if new_timestep>max_dt or new_timestep<min_dt
        check_dt_limits(new_timestep, s);
        //Always accept the step and propose a new one larger of about 1 month.
        set(new_timestep);
        return;
    }


    ///Check and repeat
    //Perform check and repeat only if dt is not equal to the minimum timestep
    if ( get_0()!=s->get_svpar_num("ts_min_dt") ){

        //check if total mass has varied a lot, despite the timestep control
        if(check_repeat(s, s->getp_0(Mass::ID), s->getp(Mass::ID), s->getp(dMdt::ID), 0.0)){

            //if (s->get_ID()==0) std::cout<<" Mass repeat "<<__FILE__<<" "<<__LINE__<<std::endl;

            if (dprint)
                utilities::wait("Mass repeat",s->get_ID(),s->getp_0(Mass::ID), s->getp(Mass::ID), s->getp(Phase::ID), s->getp(Localtime::ID), s->getp_0(MHE::ID), s->getp(MHE::ID),s->hestart(),__FILE__,__LINE__);
            handle_star_check_repeat(s);
            return;
        }


        //check if helium has varied a lot, despiting the timestep control
        //If the star is a pureHE, disable the MHE check
        if (!s->aminakedhelium()){
            if(check_repeat(s, s->getp_0(MHE::ID), s->getp(MHE::ID), s->getp(dMHEdt::ID), s->hestart())){

                //if (s->get_ID()==0) std::cout<<" MHE repeat "<<__FILE__<<" "<<__LINE__<<std::endl;


                if (dprint)
                    utilities::wait("MHE repeat",s->get_ID(), s->getp(Phase::ID), s->getp(Localtime::ID), s->getp_0(MHE::ID), s->getp(MHE::ID),s->hestart(),__FILE__,__LINE__);
                handle_star_check_repeat(s);
                return;
            }
        }


        //check if carbon-oxygen has varied a lot, despiting the timestep control
        if(check_repeat(s, s->getp_0(MCO::ID), s->getp(MCO::ID), s->getp(dMCOdt::ID), s->costart())){

            //if (s->get_ID()==0) std::cout<<" MCO repeat "<<__FILE__<<" "<<__LINE__<<std::endl;


            if (dprint)
                utilities::wait("MCO repeat", s->getp(Timestep::ID),__FILE__,__LINE__);
            handle_star_check_repeat(s);
            return;
        }

        //check if Radius has varied a lot, despiting the timestep control
        if(check_repeat(s, s->getp_0(Radius::ID), s->getp(Radius::ID), s->getp(dRdt::ID), 0.0)){

            //if (s->get_ID()==0) std::cout<<" Radius repeat "<<__FILE__<<" "<<__LINE__<<std::endl;
            if (dprint){
                utilities::wait("R repeat",s->getp(Phase::ID),s->getp_0(Radius::ID),s->getp(Radius::ID),s->getp(dRdt::ID),__FILE__,__LINE__);
                utilities::wait("info",s->getp(Mass::ID),s->getp(MHE::ID),s->getp(MCO::ID),s->aminakedco(),s->amiremnant(),__FILE__,__LINE__);
            }

            handle_star_check_repeat(s);
            return;

        }

        //Check OmegaSpin
        double dSpindt = (s->getp(OmegaSpin::ID)-s->getp_0(OmegaSpin::ID))/s->getp(Timestep::ID);
        //Make this check only of if allowed by the runtime option (ts_check_spin) and the velocity is significant (> 1 km/s)
        //Also do not check it if it is a remnant, OmegaSpin does not make any sense for remnant, they have their own properties (e.g. Xspin or OmegaRem)
        bool dSpin_check =!s->amiremnant() and s->get_svpar_bool("ts_check_spin")  and  s->getp_0(OmegaSpin::ID)>1.;
        if( dSpin_check and check_repeat(s,s->getp_0(OmegaSpin::ID),s->getp(OmegaSpin::ID),dSpindt,0.0)){
            handle_star_check_repeat(s);
            return;
        }


    }


    //Still we want some control if something strange is happening (we don't expect to have difference larger than
    //the machine precision)

    ///Check on Wolrdtime
    //If Worldtime>=tf, we are just at the end of the simulation
    //This can happen because  of lose of machine precision, especially if we want
    //to end the simulation after the SN explosion. In this case the code will try to do a single
    //large timestep and the final Worldtime can be slightly larger than tf (order of 1e-16).
    //If this happens just return
    double time_to_end=s->get_tf() - s->getp(Worldtime::ID);
    if (utilities::areEqual(time_to_end,0.) and time_to_end<=0){
        return; //This is the last step not need to estimate a new timestep.
    }
    //Still we want some control if something strange is happening (we don't expect to have a difference larger than
    //the machine precision)
    else if (time_to_end<=0){
        svlog.critical("The Worldtime is unexpectedly  larger than the final time",
                       __FILE__,__LINE__,sevnstd::sanity_error(""));
    }



    double new_timestep = std::min(s->get_dtmax(), time_to_end);

    svlog.debug("Max timestep"+std::to_string(new_timestep));
    if(dprint)
        utilities::wait("before", s->getp_0(Mass::ID), s->getp(Mass::ID),s->get_ID(),new_timestep, s->getp(Worldtime::ID),s->get_dtmax(),s->get_tf() - s->getp(Worldtime::ID),s->getp(Mass::ID),__FILE__,__LINE__);

    //determine next dt by imposing a maximum relative variation of stellar mass and radius
    if(fabs(s->getp(dMdt::ID)) != 0.0) {
        new_timestep = std::min(new_timestep,
                                max_variation * (s->getp(Mass::ID) / fabs(s->getp(dMdt::ID))));
        //std::cout<<std::scientific<<s->getp(Mass::ID)<<std::endl;
        //std::cout<<std::scientific<<s->getp(MHE::ID)<<std::endl;
        //svlog.pdebug("Entrato",__FILE__,__LINE__);
        if(dprint)
            utilities::wait("dMdt", new_timestep, s->getp(Worldtime::ID),s->getp(dMdt::ID),max_variation * (s->getp(Mass::ID) / fabs(s->getp(dMdt::ID))),__FILE__,__LINE__);
    }

    svlog.debug("timestep dM/dt end"+std::to_string(new_timestep)+"  " + std::to_string(s->getp(Mass::ID))+ " " + std::to_string(s->getp(dMdt::ID)));

    if(fabs(s->getp(dRdt::ID)) != 0.0)
        new_timestep = std::min(new_timestep, max_variation*( s->getp(Radius::ID)/fabs(s->getp(dRdt::ID)) ) );

    if(dprint)
        utilities::wait("dRdt", new_timestep, s->getp(Radius::ID),s->getp_0(Radius::ID),s->getp(Worldtime::ID),s->getp(dRdt::ID), max_variation*( s->getp(Radius::ID)/fabs(s->getp(dRdt::ID)) ),__FILE__,__LINE__);

    svlog.debug("timestep dR/dt end "+std::to_string(new_timestep)+"  " + std::to_string(s->getp(dRdt::ID))+ " " + std::to_string(s->getp(dRdt::ID)));


    //Disable if the star is a Naked Helium
    if( (s->getp(MHE::ID) != 0.0 or s->getp_0(MHE::ID) != 0.0) and !s->aminakedhelium()) {
        if (s->getp_0(MHE::ID) > 4.0 && s->getp(MHE::ID) > 4.0) {
            new_timestep = std::min(new_timestep,
                                    max_variation * (s->getp(MHE::ID) / fabs(s->getp(dMHEdt::ID))));
        } else {//fixed maximum variation under 2.0 Msun (avoid small time steps just because Mhe or MCO are small)
            double new_dt = 0.2 / s->getp(Timestep::ID);
            new_timestep = std::min(new_timestep, new_dt);
        }
    }

    if(dprint)
        utilities::wait("dMHEdt", new_timestep, s->getp(Worldtime::ID),s->getp(dMHEdt::ID),max_variation * (s->getp(MHE::ID) / fabs(s->getp(dMHEdt::ID))),__FILE__,__LINE__);

    svlog.debug("timestep dMhe/dt end "+std::to_string(new_timestep)+"  " + std::to_string(s->getp(MHE::ID))+ " " + std::to_string(s->getp(dMHEdt::ID)));


    if(s->getp(MCO::ID) != 0.0 || s->getp_0(MCO::ID) != 0.0) {
        if (s->getp_0(MCO::ID) > 4.0 && s->getp(MCO::ID) > 4.0) {
            new_timestep = std::min(new_timestep,
                                    max_variation * (s->getp(MCO::ID) / fabs(s->getp(dMCOdt::ID))));
        } else {//fixed maximum variation under 2.0Msun (avoid small time steps just because Mhe or MCO are small)
            double new_dt = 0.2 / s->getp(Timestep::ID);
            new_timestep = std::min(new_timestep, new_dt);
        }
    }
    //svlog.pinfo("New timestep end",new_timestep,__FILE__,__LINE__);

    if(dprint)
        utilities::wait("dMCO", new_timestep, s->getp(Worldtime::ID),__FILE__,__LINE__);

    ///Predict timestep trying to avoid the Mass to become smaller than the MHE
    //Notice if the Masses are already within the ev_naked_tshold we skip this.
    //We need the skip to a void a sort of "Achille s paradox" in which we continue to infinitely reduce the timestep to avoid that the Mass becomes larger than MHE
    //In practie, here we reduce the timestep to slowly approach the moment in which the star will becomes pureHE or nakedCO
    //Here we start from:
    // M = M0 + dMdt * DT
    // m = m0 + dmdt *DT
    // where for example M=Mass and m=MHE
    // now we want a DT so that m=M, so
    //     m0 + dmdt *DT = M0 + dMdt * DT
    // DT (dmdt - dMdt ) = M0-m0
    // so DT = (M0-m0) / (dmdt - dMdt )
    // hence we set DT_new = 0.5* (M0-m0) / (dmdt - dMdt )
    bool check_almost_naked_Hstar= !s->aminakedhelium() and s->getp(Mass::ID)>s->getp(MHE::ID)+std::max(1e-6,s->get_svpar_num("ev_naked_tshold"));
    bool check_almost_naked_HEstar= s->aminakedhelium() and s->getp(Mass::ID)>s->getp(MCO::ID)+std::max(1e-6,s->get_svpar_num("ev_naked_tshold"));


    if (check_almost_naked_Hstar and s->getp(MHE::ID)>=0.9*s->getp(Mass::ID)){
        double new_dt= 0.5*( s->getp(Mass::ID)-s->getp(MHE::ID) ) / (s->getp(dMHEdt::ID)-s->getp(dMdt::ID));
        if (new_dt>0) new_timestep = std::min(new_timestep,new_dt);
    }
    if (check_almost_naked_HEstar and s->getp(MCO::ID)>=0.9*s->getp(Mass::ID)){
        double new_dt= 0.5*( s->getp(Mass::ID)-s->getp(MCO::ID) ) / (s->getp(dMCOdt::ID)-s->getp(dMdt::ID));
        if (new_dt>0) new_timestep = std::min(new_timestep,new_dt);
    }
    /*****************************/

    ///Spins
    double dOSpin = fabs(s->getp(OmegaSpin::ID) - s->getp_0(OmegaSpin::ID))/s->getp(Timestep::ID);
    //Make this check only of if allowed by the runtime option (ts_check_spin) and the velocity is significant (> 1 km/s)
    //Also do not check it if it is a remnant, OmegaSpin does not make any sense for remnant, they have their own properties (e.g. Xspin or OmegaRem)
    bool dSpin_check =!s->amiremnant() and s->get_svpar_bool("ts_check_spin")  and  s->getp_0(OmegaSpin::ID)>1. and dOSpin != 0.0;
    //Notice: Make this check just when OmegaSpin is significant (in order to avoi a large number of timesteps with OmegaSpin close to 0)
    if(dSpin_check){
        new_timestep = std::min(new_timestep, max_variation*( s->getp(OmegaSpin::ID)/dOSpin ) );
    }




    //reduce the time-step at every critical passage (e.g. formation of the He core, CO core, remnant, ...)
    //TODO maybe we should take into account also the transformation into a WR star??

    double referencetime = 0.0;

    if(s->getp(Localtime::ID)+new_timestep > s->get_next_tphase()) {
        referencetime = s->get_next_tphase();
    }


    svlog.pdebug("referencetime",std::setprecision(20),referencetime,new_timestep,s->getp(Localtime::ID),s->get_next_tphase());

    if(dprint)
        utilities::wait("referencetime", std::setprecision(20),referencetime,new_timestep,s->getp(Localtime::ID),s->get_next_tphase(),__FILE__,__LINE__);

    //TODO the factor 1e10 is needed to follows the last part of the evolution.
    //TODO It seems to work better than just using a percentange, but it is maybe worth checking at a certain point
    if(referencetime != 0.0) {
        if (s->changedphase) {
            //cout << "A timestep/localtime  = " << new_timestep  << "   " << referencetime << " " << s->getp(Localtime::ID) << endl;
            //utilities::hardwait();
            new_timestep = std::min(new_timestep,  (referencetime - s->getp(Localtime::ID)+1e-10));
            s->changedphase = false;
            svlog.pdebug("A New timestep end",new_timestep,referencetime,s->getp(Localtime::ID),s->getp(Mass::ID),__FILE__,__LINE__);
            if (dprint)
                utilities::wait("A New timestep end",new_timestep,referencetime,s->getp(Localtime::ID),s->getp(Mass::ID),__FILE__,__LINE__);

        } else {
            s->changedphase = true;
            //cout << " changed phase true" << endl;
            //cout << "P timestep/localtime = " << new_timestep  << "   " << referencetime << " " << s->getp(Localtime::ID) << endl;
            new_timestep = std::min(new_timestep,   (referencetime - s->getp(Localtime::ID)-1e-10));
            if (new_timestep<0) new_timestep= (referencetime - s->getp(Localtime::ID))*0.99;
            //cout << " timestep/localtime = " << new_timestep  << "   " << referencetime << " " << s->getp(Localtime::ID) << " " << s->getp(MCO::ID)  << endl;
            //utilities::wait();
            svlog.pdebug("B New timestep end",new_timestep,referencetime,s->getp(Localtime::ID),__FILE__,__LINE__);
            if (dprint)
                utilities::wait("B New timestep end",new_timestep,referencetime,s->getp(Localtime::ID),s->getp(Mass::ID),__FILE__,__LINE__);
        }
    }


    //Check and modify if new_timestep>max_dt or new_timestep<min_dt
    check_dt_limits(new_timestep, s);


    if (new_timestep==0)
        svlog.critical("Proposed timestep is equal to 0",__FILE__,__LINE__,sevnstd::sse_error());
    else if(new_timestep<0)
        svlog.critical("Proposed timestep is negative "+
                       utilities::n2s(new_timestep,__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::sse_error());
    else if (new_timestep<utilities::TINY)
        svlog.warning("New dt estimated in Timestep evolve is extremely small ("+
                      utilities::n2s(new_timestep,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                                     "if this message it is frequent during the evoluton of a system, it "
                                                                     "can be an hint that somethin is broken in the sse.",__FILE__,__LINE__);


    set(new_timestep);

}
//TODO All the hardcoded numbers should be input parameters
//TODO Write different check and repeat for MHE/MCO and for Radius, currently we use tstart to discriminate
inline bool Timestep::check_repeat(Star *s, const double &V0, const double &V, const double &derivative, const double tstart) {

    double new_dt = s->get_dtmax();
    bool shouldIrepeat = false;
    double max_variation =   s->get_svpar_num("ts_maximum_variation");


    if(tstart!=0.0){



        if(V <= 4.0 or V0 <= 4.0){
            if(fabs(V-V0) > 0.4){
                new_dt = max_variation*(0.4/fabs(derivative)); //TODO this should be weighted with the mass of the star!! 0.4 is now hardcoded!!
                shouldIrepeat = true;
           }
            else
            shouldIrepeat = false;
        }

         else if(fabs(V-V0)/V > 2.0 * max_variation && V != 0.0){ //repeat tge step if we got twice the variation!
            new_dt = max_variation*(V/fabs(derivative));
            shouldIrepeat = true;
        }
        else
            shouldIrepeat = false;

        //TODO Is something as below needed? It should assure that if already are in a new phase we do not came brack in the previous phase
        //new_dt = std::max(new_dt,tstart-s->getp_0(Localtime::ID)+new_dt);

    }
    else{
            if(fabs(V-V0)/V > 2.0 * max_variation) {
                shouldIrepeat = true;
                new_dt = get()*0.5;//max_variation*(m/fabs(derivative));

        }
    }
    //TODO it should be better to use (m-m0)/m0 as relative variation (relative wrt the starting point)


    //if (s->getp(Worldtime::ID)>9.8947394131 and s->get_ID()==0)
    //    utilities::wait("CHeck and r",V,V0,fabs(V-V0),shouldIrepeat,__FILE__,__LINE__);

    //Check and modify if new_dt>max_dt or new_dt<min_dt
    check_dt_limits(new_dt, s);


    if(new_dt >= 0.8*get()) return false;


    if(shouldIrepeat) {

        if (new_dt<utilities::TINY)
            svlog.warning("New dt estimated in check and repeat is extremely small ("+
            utilities::n2s(new_dt,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                     "if this message it is frequent during the evoluton of a system, it "
                                                     "can be an hint that somethin is broken in the sse.",__FILE__,__LINE__);

        set(new_dt);
        set_0(get());
        return true;
    }
    else
        return false;

}
inline bool Timestep::check_repeat_almost_naked(Star *s) {

    if (checked_almost_naked){
        checked_almost_naked=false;
        return false;
    }


    double new_dt = s->get_dtmax();
    bool shouldIrepeat = false;



    //Here we start from:
    // M = M0 + dMdt * DT
    // m = m0 + dmdt *DT
    // where for example M=Mass and m=MHE
    // now we want a DT so that m<M, so
    //     m0 + dmdt *DT < M0 + dMdt * DT
    // DT (dmdt - dMdt ) < M0-m0
    // so DT < (M0-m0) / (dmdt - dMdt )
    // hence we set DT_new = 0.9* (M0-m0) / (dmdt - dMdt )

    if (!s->aminakedhelium() and std::abs(s->getp(Mass::ID)-s->getp(MHE::ID))<=s->get_svpar_num("ev_naked_tshold")){
        new_dt = 0.5*( s->getp_0(Mass::ID)-s->getp_0(MHE::ID) ) / (s->getp(dMHEdt::ID)-s->getp(dMdt::ID));
        shouldIrepeat=true;
    }
    else if  (std::abs(s->getp(Mass::ID)-s->getp(MCO::ID))<=s->get_svpar_num("ev_naked_tshold")){
        new_dt = 0.5*( s->getp_0(Mass::ID)-s->getp_0(MCO::ID) ) / (s->getp(dMCOdt::ID)-s->getp(dMdt::ID));
        shouldIrepeat=true;
    }
    else if (s->getp(MHE::ID)!=0 and s->getp(MCO::ID)!=0 and std::abs(s->getp(MHE::ID)-s->getp(MCO::ID))<=s->get_svpar_num("ev_naked_tshold")){
        new_dt = 0.5*( s->getp_0(MHE::ID)-s->getp_0(MCO::ID) ) / (s->getp(dMCOdt::ID)-s->getp(dMHEdt::ID));
        shouldIrepeat=true;
    }

    if (shouldIrepeat){
        //Check and modify if new_dt>max_dt or new_dt<min_dt
        check_dt_limits(new_dt, s);


        if(new_dt >= 0.8*get()){
            return false;
        }

        if (new_dt<utilities::TINY)
            svlog.warning("New dt estimated in check and repeat is extremely small ("+
                          utilities::n2s(new_dt,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                                   "if this message it is frequent during the evoluton of a system, it "
                                                                   "can be an hint that somethin is broken in the sse.",__FILE__,__LINE__);

        set(new_dt);
        set_0(get());

        checked_almost_naked=true;

        return true;
    }

    return false;


}

void Timestep::handle_star_check_repeat(Star *s){
    s->repeatstep = true; //Set the repeatstep flag to true
    //TODO I remove the following command, because the reset of changephase is now handled dreictly by the method Star::restore, that is more genreal since it can bel called also after a repeat triggered by the Binary evolution
    //Myabe we can remove this function and use s->repeatstep = true directly in Timestep::evolve
    //s->changedphase = false; //In case the previous step start the reduction of the timestep for the next change of phase, reset it since we are repeating
    return;
}
void Timestep::resynch(Star *s) {
    set(s->getp(Timestep::ID));
}

double Timestep::timestep_remnant(Star *s) {
    double new_timestep;

    if (s->get_COnaked() and s->getp(Phase::ID) == Remnant)
        svlog.critical(
                "Inside Timestep::set_remnant the star flags remnant and isconaked are both true, this is not allowed",
                __FILE__, __LINE__, sevnstd::sse_error());
    else if (s->get_COnaked()) {
        new_timestep = 1.001 * (s->get_tphase()[Remnant] - s->getp(Localtime::ID));
    } else if (s->getp(Phase::ID) == Remnant) {
        new_timestep = s->get_tf() - s->getp(Worldtime::ID);
        if (new_timestep == 0)
            new_timestep = get_0(); //If we reach tf, just put a fake timestep to avoid to raise the warning
    } else if (!s->get_COnaked() and s->getp(Phase::ID) != Remnant) {
        svlog.critical(
                "Inside Timestep::set_remnant when isconaked is false and the Phase is not Remnant, this is not allowed",
                __FILE__, __LINE__, sevnstd::sse_error());
    } else
        svlog.critical(
                "Inside Timestep::set_remnant the star flags  remnant and isconaked are both false, this is not allowed",
                __FILE__, __LINE__, sevnstd::sse_error());

    ///Check on Wolrdtime
    //If Worldtime>tf, we are just at the end of the simulation
    //This can happen because  of lose of machine precision, especially if we want
    //to end the simulation after the SN explosion. In this case the code will try to do a single
    //large timestep and the final Worldtime can be slightly larger than tf (order of 1e-16).
    //If this happens just put a fake small timestep.
    if (utilities::areEqual(s->getp(Worldtime::ID), s->get_tf()))
        return utilities::smallestSignificativeStep(s->getp(Worldtime::ID));
    //Still we want some control if something strange is happening (we don't expect to have a difference larger than
    //the machine precision)
    else if (s->getp(Worldtime::ID)>s->get_tf())
        svlog.critical("The Worldtime is unexpectedly  larger than the final time",
                       __FILE__,__LINE__,sevnstd::sanity_error(""));

    new_timestep = std::min(new_timestep,s->get_dtmax()); //Use the mimimum

    return new_timestep;
}
void Timestep::set_remnant(Star *s) {

    set_0(get());

    double new_timestep = timestep_remnant(s);

    if (new_timestep==0)
        svlog.critical("Proposed timestep is equal to 0",__FILE__,__LINE__,sevnstd::sse_error());
    if(new_timestep<0)
        svlog.critical("Proposed timestep is negative " +
                       utilities::n2s(new_timestep, __FILE__, __LINE__), __FILE__, __LINE__, sevnstd::sse_error());
    else if (new_timestep<utilities::TINY)
        svlog.warning("New dt estimated in Timestep evolve is extremely small ("+
                      utilities::n2s(new_timestep,__FILE__,__LINE__)+"). If this happens just few times, it can be ignored, but"
                                                                     "if this message it is frequent during the evoluton of a system, it "
                                                                     "can be an hint that somethin is broken in the sse.",__FILE__,__LINE__);

    set(new_timestep);

}
void Timestep::set_remnant_in_bse(Star *s, _UNUSED Binstar *b){
    double old_dt0 = get_0();
    set_remnant(s);
    set_0(old_dt0); //Restore t0
}
void Timestep::set_empty_in_bse(Star *s, Binstar *b){
    set_remnant_in_bse(s,b);
}

void Timestep::check_dt_limits(double &dt, Star *s){

    if (s->get_svpar_num("ts_max_dt")>0)
        dt=std::min(dt, s->get_svpar_num("ts_max_dt"));

    if (s->get_svpar_num("ts_min_dt")>0)
        dt=std::max(dt, s->get_svpar_num("ts_min_dt"));
}
void Timestep::changed_track(_UNUSED Star* s,Star* s_jtrack){
    //Here we set the next time step rescaling the dt proposed for the ration between the DT_phase of the new track wrt to the old one.
    //double rescale = (s_jtrack->get_next_tphase()-s_jtrack->get_current_tphase()) / (s->get_next_tphase()-s->get_current_tphase());
    //V=V*rescale;
    //GI 24/05/22, new implementation just get the minimum between the two timestep, this is more conservative
    //and solve a problem with the print timestep
    V=std::min(V,s_jtrack->getp(Timestep::ID));
}



/**dRdt**/
void dRdt::evolve(Star *s){
    set_0(get());
    set((s->getp(Radius::ID) - s->getp_0(Radius::ID))/s->getp(Timestep::ID));
}

/**dMCOdt**/
void dMCOdt::evolve(Star *s){
    set_0(get());
set((s->getp(MCO::ID) - s->getp_0(MCO::ID))/s->getp(Timestep::ID));
}

/**dMHEdt**/
void dMHEdt::evolve(Star *s){
    set_0(get());
set((s->getp(MHE::ID) - s->getp_0(MHE::ID))/s->getp(Timestep::ID));
}

/**dMdt**/
void dMdt::evolve(Star *s){
    set_0(get());
    svlog.debug("Inside dM/dt "+ std::to_string(s->getp(Timestep::ID)));
    //utilities::wait();
    set((s->getp(Mass::ID) - s->getp_0(Mass::ID))/s->getp(Timestep::ID));
}



/**Chemical composition on the surface**/
void Hsup::changed_track(_UNUSED Star *s, Star *s_jtrack) {
    //copy the fake tracks value  from s_jtrack
    V=s_jtrack->getp(ID);
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}
void HEsup::changed_track(_UNUSED Star *s, Star *s_jtrack) {
    //copy the fake tracks value  from s_jtrack
    V=s_jtrack->getp(ID);
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}
void Csup::changed_track(_UNUSED Star *s, Star *s_jtrack) {
    //copy the fake tracks value  from s_jtrack
    V=s_jtrack->getp(ID);
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}
void Nsup::changed_track(_UNUSED Star *s, Star *s_jtrack) {
    //copy the fake tracks value  from s_jtrack
    V=s_jtrack->getp(ID);
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}
void Osup::changed_track(_UNUSED Star *s, Star *s_jtrack) {
    //copy the fake tracks value  from s_jtrack
    V=s_jtrack->getp(ID);
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}

/**Convective envelope properties**/
void Qconv::changed_track(_UNUSED Star *s, Star *s_jtrack) {
    if (table_loaded){
        V=s_jtrack->getp(ID);
    } else {
        V=estimate_Qconv(s_jtrack);
    }
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);

}
void Qconv::evolve_without_table(Star *s) {

    set_0(get());
    set(estimate_Qconv(s));
    //Set V,V0 but also v, v0. This is needed because even if the tables are loaded
    //the user can choice do not use them.
    synch_v_value_without_table();

}
double Qconv::estimate_Qconv(Star *s) {


    //We follow Hurley02 with some minor change:
    //-We assume all pureHE stars are radiative (see below)
    //-We use the SEVN phase instead of the bse phases (we need this to track the lifetime in each phase)
    //-We assume all the stars <1.25 are completely convective. In BSE stars <0.3 are completely convective,
    //while stars between 0.3 and 1.25 are partially convective. However to set the convective Radius we should
    //have an estimate of the zams Radius for a 0.3 Msun, that depends on Models and metallcities. Given that
    //this is in any case a rough estiamte of the Dconv, we decided to use 1.25 as limit


    double Mcnv=0;

    //We assume all the pureHE are radiative.
    //We check with parsec tracks that this is totally consistent:
    //all the pureHE are completely radiative or have a very light convective envelope (Qconv<1e-3) toward the end
    //of their life. Metal rich pureHE stars (Z>0.02) are purely radiate up to the end.
    if (s->aminakedhelium()){
        return 0.0;
    }
    //From Hurley+00
    else if (s->getp(Phase::ID)==Lookup::MainSequence){
        if (s->get_zams()>1.25){
            return 0.0;
        }
        else{
            Mcnv = s->getp(Mass::ID);
            Mcnv*=pow((1-s->plife()),0.25); //plife is t/tMS as in Hurley00, Sec. 7.2
        }
    }
    else if (s->getp(Phase::ID)==Lookup::TerminalMainSequence){
        Mcnv=s->plife()*(s->getp(Mass::ID)-s->getp(MHE::ID)); //plife = t-tMS/(tGB_start - tMS) as in Hurley00, Sec. 7.2
    }
    //All the Giants
    else{
        Mcnv=s->getp(Mass::ID)-s->getp(MHE::ID);
    }

    return std::min(Mcnv/s->getp(Mass::ID),1.0); //Depthconv is the Mass of the envelope layer normalised over the total mass

}
void Tconv::changed_track(_UNUSED Star *s, Star *s_jtrack) {
    ///Copy the Tconv from the jtrack star
    V=s_jtrack->getp(ID);
    //V0=log10(s_jtrack->getp(ID));
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);
}
void Depthconv::changed_track(_UNUSED Star *s, Star *s_jtrack) {

    if (table_loaded){
        V=s_jtrack->getp(ID);
    } else {
        V=estimate_Dconv(s_jtrack);
    }
    value=s_jtrack->getp_fk(ID);
    value0=s_jtrack->getp_fk0(ID);

}
void Depthconv::evolve_without_table(Star *s) {
    set_0(get());
    set(estimate_Dconv(s));
    //Set V,V0 but also v, v0. This is needed because even if the tables are loaded
    //the user can choice do not use them.
    synch_v_value_without_table();
}
double Depthconv::estimate_Dconv(Star *s) {

    //We follow Hurley02 with some minor change:
    //-We assume all pureHE stars are radiative (see below)
    //-We use the SEVN phase instead of the bse phases (we need this to track the lifetime in each phase)
    //-We assume all the stars <1.25 are completely convective. In BSE stars <0.3 are completely convective,
    //while stars between 0.3 and 1.25 are partially convective. However to set the convective Radius we should
    //have an estimate of the zams Radius for a 0.3 Msun, that depends on Models and metallcities. Given that
    //this is in any case a rough estiamte of the Dconv, we decided to use 1.25 as limit

    double Dcnv=0.0;

    //We assume all the pureHE are radiative.
    //We check with parsec tracks that this is totally consistent:
    //all the pureHE are completely radiative or have a very light convective envelope (Qconv<1e-3) toward the end
    //of their life. Metal rich pureHE stars (Z>0.02) are purely radiate up to the end.
    if (s->aminakedhelium()){
        return 0.0;
    }
    //From Hurley+02
    else if (s->getp(Phase::ID)==Lookup::MainSequence){
        if (s->get_zams()>1.25){
            return 0.0;
        }
        else{
            Dcnv = s->getp(Radius::ID); //RHE is 0 becasue MainSequence
            Dcnv*=pow((1-s->plife()),0.25); //plife is t/tMS, from Hurley02 Eq. 37
        }
    }
    else if (s->getp(Phase::ID)==Lookup::TerminalMainSequence){
        Dcnv=std::sqrt(s->plife())*(s->getp(Radius::ID)-s->getp(RHE::ID)); //plife = t-tMS/(tGB_start - tMS)  from Hurley02 Eq. 39.
    }
    //All the Giants
    else{
        Dcnv=s->getp(Radius::ID)-s->getp(RHE::ID);
    }


    return std::min(Dcnv/s->getp(Radius::ID),1.0); //Depthconv is normalised over the total radius
}

/**Spin**/

void AngMomSpin::evolve(_UNUSED Star *s){

    set_0(get());
    //If 0 just return, it cannot become larger than 0 because Ltot=0
    if (get_0()==0)
        return;

    double new_spin = evolve_angmom(s);

    set(new_spin);
}

double AngMomSpin::evolve_angmom(_UNUSED Star *s){
    double dM_wind = s->getp(Mass::ID)-s->getp_0(Mass::ID);
    const double & R0 = s->getp_0(Radius::ID);
    double OmegaSpin0  = s->getp(OmegaSpin::ID); //OmegaSpin has to be evolved


    ///1- Stellar angular momentum loss through winds
    // In order to estimate the angular momentum loss by winds, we assume that the wind removed
    // a thin shell of material with Radius R0 and mass DM_wind, so dL = 2/3 * R0^2 * R0^2 * dM * OmegaSpin0
    //0.6666666666666666=2/3, notice dM_wind is negative, so dLwind is negative
    double dLwind = 0.6666666666666666*R0*R0*dM_wind*OmegaSpin0;
    ///2- Stellar angular momentum loss trough magnetic braking
    //Here we follow Hurley+02, Eq. 50 and MOBSE, in MOBSE the comment is:
    // Include magnetic braking for stars that have appreciable convective
    // envelopes. This includes MS stars with M < 1.25, HG stars near the GB
    // and giants. MB is not allowed for fully convective MS stars.
    //Notice that even if stated there is not a proper check on the amount of convective envelope
    //TODO Use the amoun of convective envelope mass instead of the whole convective envelope
    double dLmbraking=0;
    if (s->getp_0(Mass::ID)>=0.35 and s->get_bse_phase_0()<10){
        dLmbraking = -(5.83e-16*s->Menvelope(true)*R0*R0*R0*OmegaSpin0*OmegaSpin0*OmegaSpin0/s->getp_0(Mass::ID))*s->getp(Timestep::ID)*1e-6; //Timestep is in Myr, Eq. is in yr
    }
    ///2- Estimate new angmom
    //V is the current L
    double Lnew  = V + dLwind + dLmbraking; //Both dLwind and dLmbraking are <0
    //Limit to 0 or to maximum Lcrit allowe by critical rotational velocity
    double Lcrit = s->getp(Inertia::ID)*utilities::omega_crit(s->getp(Mass::ID),s->getp(Radius::ID));
    Lnew = std::min(std::max(Lnew,0.0),Lcrit);

    return Lnew;
}

void AngMomSpin::changed_track(Star *s, _UNUSED Star *s_jtrack) {
    //We consider that angular momentum is conserved, however we have to check
    //if the new angmom is below the critical L, otherwise set is to Lcrit
    //Notice Inertia, Mass and Radius alread changed
    //if (s->get_ID()==0) utilities::hardwait("Jc", V, s->getp(Inertia::ID),utilities::omega_crit(s->getp(Mass::ID),s->getp(Radius::ID)));
    double Lcrit = s->getp(Inertia::ID)*utilities::omega_crit(s->getp(Mass::ID),s->getp(Radius::ID));
    V = std::min(V,Lcrit);

}


void OmegaSpin::evolve(_UNUSED Star *s){

    set_0(get());
    set(s->getp(AngMomSpin::ID)/s->getp(Inertia::ID));

}

/**Xspin (black hole spin) properties**/
void Xspin::set_remnant(Star *s){
    V0 = V;
    value0 = value;
    V = s->get_staremnant(Xspin::ID);
    value = s->get_staremnant(Xspin::ID);
}

void Xspin::evolve_remnant(Star *s){
    set_0(get());

    set(s->get_staremnant(Xspin::ID));
}


/**RemnantType properties**/
void RemnantType::set_remnant(Star *s) {

    V = value = (double)s->get_staremnant()->get_remnant_type();
}

void RemnantType::set_empty(_UNUSED Star *s){
    V0 = V = Lookup::Remnants::Empty;
    value0 = value = Lookup::Remnants::Empty;
}

void RemnantType::evolve(_UNUSED Star *s) {
    //Do not evolve if it is not a remnant
    //V = V0 = value = value0 = s->get_supernova()->get_remnant_type();
}

/**JIT properties**/
//Lambda for envelope binding energy
double Lambda::get(const Star *s) {


    //If evolve_number>last_evolve_number, we have to obtain a new estimate for lambda
    if(new_estimate_needed()){
        //When we don't have a clear core-envelope separation return nan
        if (!s->haveienvelope())
            set(std::nan(""));
        else
            set(estimate_lambda(s)); //Set the new V and set evolve_number=last_evolve_number
    }


    return V;
}

double Lambda::estimate_lambda (const Star *star){
    const double &lambda_option = star->get_svpar_num("star_lambda");

    if (lambda_option==-1)
        return estimate_lambda_BSE(star, false,true);
    else if (lambda_option==-11)
        return estimate_lambda_BSE(star, false,false);
    else if (lambda_option==-12)
        return estimate_lambda_BSE(star, true,true);
    else if (lambda_option==-13)
        return estimate_lambda_BSE(star, true,false);
    else if(lambda_option==-2)
        return estimate_lambda_Claeys14(star,false);
    else if(lambda_option==-21)
        return estimate_lambda_Claeys14(star,true);
    else if(lambda_option==-3)
        return estimate_lambda_Izzard04(star,false);
    else if(lambda_option==-31)
        return estimate_lambda_Izzard04(star,true);
    else if(lambda_option==-4)
        return estimate_lambda_Klencki21(star,true);
    else if(lambda_option==-41)
        return estimate_lambda_Klencki21(star,false);
    else if(lambda_option==-5)
        return  estimate_lambda_Nanjing(star,true);  //As in Compas but interpolating between XU&Li10 and Dominik+12 tracks
    else if(lambda_option==-51)
        return  estimate_lambda_Nanjing(star,false); //Exactly as in Compas
    else if (lambda_option>0){
        int bse_phase = star->get_bse_phase(); //NakedHelium or WR
        if (bse_phase==7 or bse_phase==8 or bse_phase==9) return star->get_svpar_num("star_lambda_pureHe");
        else return star->get_svpar_num("star_lambda");
    }
    else
        svlog.critical("Option " + utilities::n2s(lambda_option,__FILE__,__LINE__)+
        " not implemented",__FILE__,__LINE__,sevnstd::sevnio_error());

    return std::nan("");
}

double Lambda::estimate_lambda_Izzard04(const Star *star,bool whole_cenv) {

    //Lambda estimate following Claeys et al. 2014 (Appendix A)

    ///Initial checks
    //Extra check: we should not call this fuction for Main Sequence Star and Remnants
    if (!star->haveienvelope())
        svlog.critical("Trying to estimate the lambda for envelope binding energy, "
                       "for a star without a clear Envelope/Core separation:"
                       "\nStar ID: " + utilities::n2s(star->get_ID(),__FILE__,__LINE__) +
                       "\nStar Phase " + utilities::n2s((int)star->getp(Phase::ID),__FILE__,__LINE__),
                       __FILE__,__LINE__,sevnstd::ce_error());


    //To be consistent with MOBSE, use the BSE phases
    int bse_phase = star->get_bse_phase();



    //If it is a Naked He, Clayes has not fitting formula, she just assumed 0.5
    if (bse_phase==7 or bse_phase==8 or bse_phase==9)
        return star->get_svpar_num("star_lambda_pureHe"); //No fit just return 0.5 as done in MOBSE/BSE/Claeys



    ///Start Estimate
    double lambda_ce; //Final results

    /**
     * The final lambda depends on the Mass of the convective envelope.
     * If we are using convective tables just use the table value, otherwise estimate the Mcenv
     * using analytic approximation from Hurley (see Qconv::estimate_without_table)
     * Mcenv_env=Qconv*Mass,
     * If whole_cenv is true consider the whole envelope as convective
     */
    double Mcenv= whole_cenv ? star->getp(Mass::ID)-star->getp(MHE::ID) : star->getp(Mass::ID)*star->getp(Qconv::ID);//Mass of the convective envelope

    //Mcenv = whole_cenv ? 1 : star->getp(Mass::ID)*star->getp(Qconv::ID)/(star->getp(Mass::ID)-star->getp(MHE::ID));

    //Other variables
    double logL = std::log10(star->getp(Luminosity::ID));
    double logM = std::log10(star->getp(Mass::ID));


    double lam1=0, lam2=0; //lambda1 e lambda2 in Claeys et al. 2014 (Appendix A)



    ///Estimate lam2 (needed only for Mcenv<1)
    if (Mcenv<1){
        //Calculate Rzams
        double Rzams = get_Rzams(star);
        lam2 = 0.42*pow(Rzams/star->getp(Radius::ID),0.4);
    }

    ///Estimate lam1 (needed only for Mcenv>0)
    if (Mcenv>0){

        //HG and RGB (eq. A3)
        if (bse_phase<=3){
            double proposed_lam1 = 3/(2.4 + pow(star->getp(Mass::ID),-3./2.));
            proposed_lam1 -= 0.15*logL;
            lam1 = std::min(0.8, proposed_lam1);
        }
            //CHEB, E-AGB, TP-AGB (eq. A4-A5)
        else if (bse_phase<=6){
            double lam3 = std::min(0.9, 0.58 + 0.75*logM) - 0.08*logL; //Eq. A4 (typo in the paper: -0.9)
            //CHEB, E-AGB
            if(bse_phase<6) lam1 = std::min(0.8, std::min(1.25-0.15*logL,lam3)); //Eq. A5a
                //TP-AGB
            else lam1= std::min(std::max(-3.5 -0.75*logM+logL, lam3),1.0); //Note maximum of 1 is written in the text.



        }
        else
            svlog.critical("Star with BSE Phase "+utilities::n2s(bse_phase,__FILE__,__LINE__)+
                           "should enter here.",__FILE__,__LINE__,sevnstd::ce_error());



    }

    ///Multiply by 2
    //Notice, in the Claeys paper it seems that this multiplication by 2 has to be made only at the end
    //after the correction for ionization energy, however in bse and in the original izzard phd (http://personal.ph.surrey.ac.uk/~ri0005/doc/thesis/thesis.pdf, Appendix E)
    //is the opposite, lam1 is estimated multplying by 2 and then it is used in the ion machinary.
    //Therefore here we multiply both lam1 and lam2 by 2
    lam2=2*lam2;
    lam1=2*lam1;

    ///Estimate ionization contribution //Eq. A6
    double f_ion=star->get_svpar_num("star_lambda_fth"); //This is the equivalent of lambda_ion in Claeys et al. 2014 (Appendix A)
    if (f_ion>0 and Mcenv>0) {
        double lam4, lam5; //Eq. A6, A8
        //Parameter to estimate lam5
        double a_ion, b_ion, c_ion, d_ion;
        b_ion = std::max(1.5, 3 - 5 * logM);
        c_ion = std::max(3.7 + 1.6 * logM, 3.3 + 2.1 * logM);
        if (bse_phase <= 3) { //HG,RGB
            a_ion = std::min(-0.5, 1.2 * (logM - 0.25) * (logM - 0.25) - 0.7); // Eq. A9a
            d_ion = std::max(0., std::min(0.15, 0.15 - 0.25 * logM)); //Eq. 12a
        } else if (bse_phase <= 6) { //CHEB,EAGB.TAGB
            a_ion = std::max(-0.5, -0.2 - logM); // Eq. A9b
            d_ion = 0.; //Eq. 12b
        } else {
            svlog.critical("Unexpected BSE type " + utilities::n2s(bse_phase, __FILE__, __LINE__), __FILE__, __LINE__,
                           sevnstd::sanity_error());
        }

        lam5 = 1 / ( a_ion + std::atan(b_ion*(c_ion-logL)) + d_ion*(logL-2.) ); //Eq. A8
        lam4 = std::max( std::min(lam5,100.),lam1 ); //Eq. A7


        //Eq. A6
        if (f_ion==1) lam1 = lam4;
        else lam1 = lam1 + f_ion*(lam4-lam1); //Eq. A6


    }

    ///Final lambda_ce estimate //Eq. A1 without 2, see comment above
    if (Mcenv==0) lambda_ce=lam2;
    else if (Mcenv>=1) lambda_ce =lam1;
    else lambda_ce=(lam2 + std::sqrt(Mcenv)*(lam1-lam2));


    return lambda_ce;
}

double Lambda::estimate_lambda_Claeys14(const Star *star,bool whole_cenv) {

    //Lambda estimate following Claeys et al. 2014 (Appendix A)

    ///Initial checks
    //Extra check: we should not call this fuction for Main Sequence Star and Remnants
    if (!star->haveienvelope())
        svlog.critical("Trying to estimate the lambda for envelope binding energy, "
                       "for a star without a clear Envelope/Core separation:"
                       "\nStar ID: " + utilities::n2s(star->get_ID(),__FILE__,__LINE__) +
                       "\nStar Phase " + utilities::n2s((int)star->getp(Phase::ID),__FILE__,__LINE__),
                       __FILE__,__LINE__,sevnstd::ce_error());


    //To be consistent with MOBSE, use the BSE phases
    int bse_phase = star->get_bse_phase();



    //If it is a Naked He, Clayes has not fitting formula, she just assumed 0.5
    //TODO we can do it with PARSEC
    if (bse_phase==7 or bse_phase==8 or bse_phase==9) //NakedHelium or WR
        return star->get_svpar_num("star_lambda_pureHe"); //No fit just return 0.5 as done in MOBSE/BSE/Claeys



    ///Start Estimate
    double lambda_ce; //Final results

    /**
     * The final lambda depends on the Mass of the convective envelope.
     * If we are using convective tables just use the table value, otherwise estimate the Mcenv
     * using analytic approximation from Hurley (see Qconv::estimate_without_table)
     * Mcenv_env=Qconv*Mass,
     * If whole_cenv is true consider the whole envelope as convective
     */
    double Mcenv= whole_cenv ? star->getp(Mass::ID)-star->getp(MHE::ID) : star->getp(Mass::ID)*star->getp(Qconv::ID);//Mass of the convective envelope

    //Mcenv = whole_cenv ? 1 : star->getp(Mass::ID)*star->getp(Qconv::ID)/(star->getp(Mass::ID)-star->getp(MHE::ID));

    //Other variables
    double logL = std::log10(star->getp(Luminosity::ID));
    double logM = std::log10(star->getp(Mass::ID));


    double lam1=0, lam2=0; //lambda1 e lambda2 in Claeys et al. 2014 (Appendix A)



    ///Estimate lam2 (needed only for Mcenv<1)
    if (Mcenv<1){

        //Calculate Rzams
        double Rzams = get_Rzams(star);
        lam2 = 0.42*pow(Rzams/star->getp(Radius::ID),0.4);
    }

    ///Estimate lam1 (needed only for Mcenv>0)
    if (Mcenv>0){

        //HG and RGB (eq. A3)
        if (bse_phase<=3){
            double proposed_lam1 = 3/(2.4 + pow(star->getp(Mass::ID),-3./2.));
            proposed_lam1 -= 0.15*logL;
            lam1 = std::min(0.8, proposed_lam1);
        }
            //CHEB, E-AGB, TP-AGB (eq. A4-A5)
        else if (bse_phase<=6){
            double lam3 = std::min(0.9, 0.58 + 0.75*logM) - 0.08*logL; //Eq. A4 (typo in the paper: -0.9)
            //CHEB, E-AGB
            if(bse_phase<6) lam1 = std::min(0.8, std::min(1.25-0.15*logL,lam3)); //Eq. A5a
                //TP-AGB
            else lam1= std::min(std::max(-3.5 -0.75*logM+logL, lam3),1.0); //Note maximum of 1 is written in the text.



        }
        else
            svlog.critical("Star with BSE Phase "+utilities::n2s(bse_phase,__FILE__,__LINE__)+
                           "should enter here.",__FILE__,__LINE__,sevnstd::ce_error());



    }



    ///Estimate ionization contribution //Eq. A6
    double f_ion=star->get_svpar_num("star_lambda_fth"); //This is the equivalent of lambda_ion in Claeys et al. 2014 (Appendix A)
    if (f_ion>0 and Mcenv>0) {
        double lam4, lam5; //Eq. A6, A8
        //Parameter to estimate lam5
        double a_ion, b_ion, c_ion, d_ion;
        b_ion = std::max(1.5, 3 - 5 * logM);
        c_ion = std::max(3.7 + 1.6 * logM, 3.3 + 2.1 * logM);
        if (bse_phase <= 3) { //HG,RGB
            a_ion = std::min(-0.5, 1.2 * (logM - 0.25) * (logM - 0.25) - 0.7); // Eq. A9a
            d_ion = std::max(0., std::min(0.15, 0.15 - 0.25 * logM)); //Eq. 12a
        } else if (bse_phase <= 6) { //CHEB,EAGB.TAGB
            a_ion = std::max(-0.5, -0.2 - logM); // Eq. A9b
            d_ion = 0.; //Eq. 12b
        } else {
            svlog.critical("Unexpected BSE type " + utilities::n2s(bse_phase, __FILE__, __LINE__), __FILE__, __LINE__,
                           sevnstd::sanity_error());
        }

        lam5 = 1 / ( a_ion + std::atan(b_ion*(c_ion-logL)) + d_ion*(logL-2.) ); //Eq. A8
        lam4 = std::max( std::min(lam5,100.),lam1 ); //Eq. A7


        //Eq. A6
        if (f_ion==1) lam1 = lam4;
        else lam1 = lam1 + f_ion*(lam4-lam1); //Eq. A6


    }

    ///Final lambda_ce estimate //Eq. A1 without 2, see comment above
    if (Mcenv==0) lambda_ce=2*lam2;
    else if (Mcenv>=1) lambda_ce =2*lam1;
    else lambda_ce=2*(lam2 + std::sqrt(Mcenv)*(lam1-lam2));


    return lambda_ce;
}

double Lambda::estimate_lambda_Parsec(const Star *star) {
    return 0.1*estimate_lambda_BSE(star,false);
}

double Lambda::estimate_lambda_BSE(const Star *star, bool whole_cenv, bool m0_as_hurley) {

    //Lambda estimate following exactly the implementation of the function CELAMF in zfuncs.f  in MOBSE ()
    //NOTICE: this should come from the Appendix A in Claeys+14, but there is a typo fixed (it was already fixed
    //also in estimate_lambda_Claeys14) and several other differences.

    ///Initial checks
    //Extra check: we should not call this fuction for Main Sequence Star and Remnants
    if (!star->haveienvelope())
        svlog.critical("Trying to estimate the lambda for envelope binding energy, "
                       "for a star without a clear Envelope/Core separation:"
                       "\nStar ID: " + utilities::n2s(star->get_ID(),__FILE__,__LINE__) +
                       "\nStar Phase " + utilities::n2s((int)star->getp(Phase::ID),__FILE__,__LINE__),
                       __FILE__,__LINE__,sevnstd::ce_error());


    //To be consistent with MOBSE, use the BSE phases
    int bse_phase = star->get_bse_phase();

    //If it is a Naked He, Clayes has not fitting formula, she just assumed 0.5
    if (bse_phase==7 or bse_phase==8 or bse_phase==9) //NakedHelium or WR
        return star->get_svpar_num("star_lambda_pureHe"); //No fit just return 0.5 as done in MOBSE/BSE/Claeys




    ///Start Estimate
    double lambda_ce; //Final results

    /**
     * The final lambda depends on the  fraction  of the mass convective envelope over the total envelope mass.
     * If we are using convective tables just use the table value, otherwise estimate the Mcenv
     * using analytic approximation from Hurley (see Qconv::estimate_without_table)
     * fcenv=Qconv,
     * If whole_cenv is true consider the whole envelope as convective fcenv=1
     */
    double fcenv = whole_cenv ? 1 : star->getp(Mass::ID)*star->getp(Qconv::ID) / (star->getp(Mass::ID) - star->getp(MHE::ID));//Mass of the convective envelope



    //Other variables
    //Notice, in Clayes+2014, the mass seems the current mass while in BSE/MOBSE it is M0
    //M0 is the initial effective mass (see Sec. 7.1 in Hurley+00), during the MS, HG and nakedHelium MS
    //M0=M, but for later stages M0 is not changed except when the star starts the ZAHB (M0=M) or become a pureHE
    //This was explained in Sec. 7.1  because the age depends mainly on the core and in these phases the evolution of the envelope
    //and of the core is decoupled.
    //However (GI opinion, for the binding energy of the envelope the current mass of the envelope DOES matter)
    double Mass = star->getp(Mass::ID);

    //Exploring MOBSE/BSE we realise that during MS and nakedHe MS M0=M, in all the other cases
    //M0 is freezed to the value at the end of the MS and of the nakedHE MS
    //THerefore, the first time we enter here in a phase outside MS (we do not have to consider Cheb since the lambda for them is set to 0.5)
    //we can set a value M0 that can be used later, in this way we will be in agreement with MOBSE/BSE
    if (star->getp(Phase::ID)>1 and m0_as_hurley){
        Mass = get_M0_BSE(star);
    }

    double lam1=0, lam2=0; //lambda1 e lambda2 in Claeys et al. 2014 (Appendix A)

    if (fcenv>0){

        double logL = std::log10(star->getp(Luminosity::ID));
        double logM = std::log10(Mass);

        //Comment from MOBSE: Formulae for giant-like stars; also used for HG and CHeB stars close
        // to the Hayashi track.
        if (bse_phase<=5){
            double m1 =  Mass;
            if (bse_phase>3) m1 = 100.0;
            lam1 = 3./(2.4 + 1./std::pow(m1,3./2.)) - 0.15*logL;
            lam1 = std::min(lam1,0.8);
        }
        else{
            lam1 = -3.5 - 0.75*logM + logL;
        }

        if (bse_phase>3){
            //Comment from MOBSE: Eq A.4  lam2 = lam3
            lam2 = std::min(0.9,0.58 + 0.75*logM) - 0.08*logL;

            if (bse_phase<6){
                lam1 = std::min(lam2,lam1);
            }else{
                lam1 = std::max(lam2,lam1);
                lam1 = std::min(lam1,1.);
            }
        }

        lam1 = 2.*lam1;

        double f_ion=star->get_svpar_num("star_lambda_fth"); //This is the equivalent of lambda_ion in Claeys et al. 2014 (Appendix A)
        if (f_ion>0){
        //Comment from MOBSE:     Use a fraction FAC of the ionization energy in the energy balance.
        // Note that from eq. A.6 to eq. A.12 lam5 = 1/lam2
        // Eq. A.9
            double aa=0, bb=0, cc=0, dd=0;
            if (bse_phase<=3){
                double logM_m_025= (logM - 0.25);
                aa = std::min(1.2*logM_m_025*logM_m_025 - 0.7,-0.5);
            } else{
                aa = std::max(-0.2 - logM,-0.5);
            }

            //Comment from MOBSE: Eq. A.10
            bb = std::max(3. - 5.*logM,1.5);
            //Comment from MOBSE: Eq. A.11
            cc = std::max(3.7 + 1.6*logM,3.3 + 2.1*logM);
            lam2 = aa + std::atan(bb*(cc - logL));
            //Comment from MOBSE: Eq. A.12
            if(bse_phase<=3){
                dd = std::max(0.,std::min(0.15,0.15 - 0.25*logM));
                lam2 = lam2 + dd*(logL - 2.);
            }

            lam2 = std::max(lam2,0.01);
            lam2 = std::max(1./lam2,lam1);

            if(f_ion>=1.0){
                lam1 = lam2;
            } else{
                //Comment from MOBSE: Eq. A.6
                lam1 = lam1 + f_ion*(lam2 - lam1);
            }
        }
    }

    if (fcenv<1){
        double Rzams = get_Rzams(star);
        lam2 = 0.42*std::pow(Rzams/star->getp(Radius::ID),0.4);
        lam2 = 2*lam2;
    }

    if (fcenv<=0.0){
        lambda_ce = lam2;
    }
    else if (fcenv>=1){
        lambda_ce = lam1;
    } else{
        //Comment from MOBSE: Interpolate between HG and GB values depending on conv. envelope mass.
        lambda_ce = lam2 + std::sqrt(fcenv) * (lam1-lam2);
    }



    return lambda_ce;
}



double Lambda::get_M0_BSE(const Star *s) {
    if (Mzams_cachedM!=s->get_zams() or Zmet_cachedM!=s->get_Z()){
        Star s_auxiliary(s, 0, s->get_zams(), s->get_Z(), "tms"); //Notice no need to seet pureHe false of true since pureHE cannot enter here
        M0_cached = s_auxiliary.getp(Mass::ID);
        Mzams_cachedM=s->get_zams();
        Zmet_cachedM=s->get_Z();
    }

    return M0_cached;
}


double Lambda::get_Rzams(const Star *s) {
    if (Mzams_cachedR!=s->get_zams() or Zmet_cachedR!=s->get_Z()){
        Star s_auxiliary(s, 0, s->get_zams(), s->get_Z(), "zams"); //Notice no need to seet pureHe false of true since pureHE cannot enter here
        Rzams_cached = s_auxiliary.getp(Radius::ID);
        Mzams_cachedR=s->get_zams();
        Zmet_cachedR=s->get_Z();
    }

    return Rzams_cached;
}

double Lambda::estimate_lambda_Klencki21(const Star *star, bool interpolate) {

    if (star->aminakedhelium()){ //In this case we consider really only the pureHe stars since we have a fit for WR direcly from the tracks
        return star->get_svpar_num("star_lambda_pureHe"); //No fit just return 0.5 as done in MOBSE/BSE/Claeys
    }

    //If it is the first time initialise the unique_ptr with the Lambda Klencki class
    //and reset the flag
    if (first_call){
        if (interpolate)  lambda_base_ptr = utilities::make_unique<Lambda_Klencki_interpolator>(star);
        else lambda_base_ptr = utilities::make_unique<Lambda_Klencki>(star);
        first_call=false;
    }

    //Return lambda from Klencki
    return lambda_base_ptr->operator()(star);
}

double Lambda::estimate_lambda_Nanjing(const Star *star, bool interpolate) {
    //If it is the first time initialise the unique_ptr with the Lambda Nanjing class
    //and reset the flag
    if (first_call){
        if (interpolate) lambda_base_ptr = utilities::make_unique<Lambda_Nanjing_interpolator>(star);
        else lambda_base_ptr = utilities::make_unique<Lambda_Nanjing>(star);
        first_call=false;
    }

    //Return lambda from Nanjing
    return lambda_base_ptr->operator()(star);
}


//Envelope binding energy
double Ebind::get(const Star *s) {

    //If evolve_number>last_evolve_number, we have to obtain a new estimate for lambda
    if(new_estimate_needed()){
        //When we don't have a clear core-envelope separation return nan
        if (s->haveienvelope()){
            double Menv  = s->aminakedhelium() ? s->getp(Mass::ID) - s->getp(MCO::ID) : s->getp(Mass::ID) - s->getp(MHE::ID);
            double E = - (1/s->getp(Lambda::ID))*s->getp(Mass::ID)*Menv/(s->getp(Radius::ID));
            set(E); //Set the new V and set evolve_number=last_evolve_number
        }
        else
            set(0.); //Set the new V and set evolve_number=last_evolve_number
    }
    return V;
}


//Phase BSE
double PhaseBSE::get(const Star *s){
    //If evolve_number>last_evolve_number, we have to obtain a new estimate for lambda
    if(new_estimate_needed()){
        set(s->get_bse_phase());
    }

    return V;
}

//Zams
double Zams::get(const Star *s){
    //If evolve_number>last_evolve_number, we have to obtain a new estimate for lambda
    if(new_estimate_needed()){
        set(s->get_zams());
    }

    return V;
}

double Zmet::get(const Star* s){

    if(new_estimate_needed()){
        set(s->get_Z());
    }

    return V;
}

void dMcumul_RLO::update_from_binary(_UNUSED Star *s, const double &DV, _UNUSED  Binstar *b){

        if (std::isnan(DV) or std::isinf(DV))
            svlog.critical("Update from Binary of property " +
                           name()+" is nan or infinite",__FILE__,__LINE__);


        ///GI 21/12/22, I removed this check because, it was not really robust
        ///for example at a certain point we could add that over-critical rotating stars do not accrete
        ///during a RLO even if the RLO is still ongoing, therefore the check that DV is 0 is not robust
        ///For this reason we decided to add the Binary as a optional parameter of update_form_binaries (make sense)
        ///so that we can directly use the information from the binary. In this case we can check if the RLO is ongoing
        //if it is not we set V to 0

        //Check if we have just stopped a RLO phase
        //e.g. when DV from RLO is 0 and dMcumu_RLO is not 0
        //If this is the case reset the property V (i.e. V=0).
        //We do not use directly the method  reset()
        //because we don't want to set V0 to 0.
        //A non 0 zero value of V0 is needed if we have to restore the property after a repetitions
        //if (DV==0. and V!=0.)
        //    V=value=0.0;
        //else
        //    V = V + DV;
        //Call without binary, not allowed
        if (b== nullptr){
            svlog.critical("Called update from binary for dMcumul_RLO but without including the binary object, this is not allowed",
                           __FILE__,__LINE__,sevnstd::sanity_error(""));
        }
        //Call with binary, but RLO is not ongoing and DV is not 0
        else if (!utilities::areEqual(DV,0.0) and !b->getprocess(RocheLobe::ID)->is_process_ongoing()){
            svlog.critical("Called update from binary for dMcumul_RLO but RLO is not ongoing, this is not allowed",
                           __FILE__,__LINE__,sevnstd::sanity_error(""));
        }
        //RLO is not ongoing, just set V and value to 0.0
        else if (!b->getprocess(RocheLobe::ID)->is_process_ongoing()){
            V=value=0.0;
        }
        //in the other cases, update it
        else{
            V = V + DV;
        }

    }


double dMRLOdt::get(const Star *s) {

    //If evolve_number>last_evolve_number, we have to obtain a new estimate for lambda
    if(new_estimate_needed()) {

        //If dMcumul_RLO is currently 0, it means that no RLO is ongoing
        if (s->getp(dMcumul_RLO::ID)==0.0) {
            set(0.0);
        }
        //Deal with cases of a RLO swap (when in a RLO ends and in the next step a RLO stars)
        //for the RLO process the RLO was always ongoing and the donor and accretor swaps
        //So if dMcumu_RLO before times dMcumu_RLO after is negative or 0, we just start a RLO
        //or we in a case of donor/accretor swap. IN such cases the rate is just the current value dMcumul/timestep
        else if (s->getp(dMcumul_RLO::ID)*s->getp_0(dMcumul_RLO::ID)<=0){
            const auto dt = s->get_last_Timestep(); //Timestep in Myr
            set(s->getp(dMcumul_RLO::ID)/dt);
        }
        else{
            const auto dMRLO = s->getp(dMcumul_RLO::ID) - s->getp_0(dMcumul_RLO::ID);  //dM in Myr in the current timestep, it could positive or negative
            //Here we use get_last_Timestep, because this is  JIT properties therefore sometime it could be called before
            //the Timestep evolution (and so dt should be get(Timestep::ID)) and some time after (and so dt should be get_0(Timestep::ID))
            //get_last_Timestep takes into account these differences
            const auto dt = s->get_last_Timestep(); //Timestep in Myr
            set(dMRLO/dt);
        }

    }

    return V;
}

double dMaccwinddt::get(const Star *s) {
    //If evolve_number>last_evolve_number, we have to obtain a new estimate for lambda
    if(new_estimate_needed()) {
        double dMwind = s->getp(dMcumulacc_wind::ID) -  s->getp_0(dMcumulacc_wind::ID);
        //Here we use get_last_Timestep, because this is  JIT properties therefore sometime it could be called before
        //the Timestep evolution (and so dt should be get(Timestep::ID)) and some time after (and so dt should be get_0(Timestep::ID))
        //get_last_Timestep takes into account these differences
        const auto dt = s->get_last_Timestep(); //Timestep in Myr
        set(std::max(0.0,dMwind/dt));
    }

    return V;
}




double Event::get(const Star *s) {

    //If evolve_number>last_evolve_number, we have to obtain a new estimate
    if(new_estimate_needed()){

        if (!is_qhe_set and s->amifollowingQHE()){
            is_qhe_set=true; //To avoid to enter here again
            set(Lookup::EventsList::QHE);
        }
        else if (s->getp_0(Phase::ID)!=s->getp(Phase::ID) and s->getp(Worldtime::ID)>0){
            set(Lookup::EventsList::ChangePhase);
        }
        else if (s->getp_0(RemnantType::ID)!=s->getp(RemnantType::ID) and s->getp(Worldtime::ID)>0){
            set(Lookup::EventsList::ChangeRemnant);
        }
        else{
            set(Lookup::EventsList::NoEvent);
        }
    }

    return V;
}

void NSsalpha::set_remnant(Star *s) {
    //First check if the star is a NS
    if (s->amiNS()){
        //2-We are sure the Star is a NS, so downcast from Staremant* to NS* to use the get_salpha method
        const NSrem *n = static_cast<const NSrem*> (s->get_staremnant()); //Cast to NSrem pointer
        sinalpha = n->get_salpha();
    }
}


double Plife::get(const Star *s) {
    //If evolve_number>last_evolve_number, we have to obtain a new estimate for lambda
    if(new_estimate_needed()) {
        if (s->getp(Phase::ID)<Lookup::Phases::Remnant) set(s->plife());
    }
    return V;
}

void Derived_Property::update_derived(Star *s) {
    //Here we want to evolve the derived properties, without setting V0
    //1-Save the V0 value
    double _V0=V0;
    //2-Evolve (depending on the stellar type)
    if(!s->amiremnant()) evolve(s);
    else if(s->aminakedco()) evolve_nakedco(s);
    else if(!s->amiempty()) evolve_remnant(s);
    //Notice: do nothing if the star is empty
    //3-Restore the V0 value
    set_0(_V0);
}


