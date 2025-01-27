//
// Created by iorio on 7/5/21.
//

#include <remnant.h>
#include <star.h>


Staremnant::Staremnant(Star *s, double Mremnant) : born_time(s->getp(Worldtime::ID)),Mremnant_at_born(Mremnant) {}

double Staremnant::get(Star *s, size_t ID) const {

    if (ID==Mass::ID)
        return Mass(s);
    else if(ID==Radius::ID)
        return Radius(s);
    else if(ID== Luminosity::ID)
        return Luminosity(s);
    else if(ID== OmegaRem::ID)
        return OmegaRem(s);
    else if(ID== Inertia::ID)
        return Inertia(s);
    else if(ID== Bmag::ID)
        return Bmag(s);
    else if(ID == Xspin::ID)
        return Xspin(s);
    else
        svlog.critical("Property with ID="+ utilities::n2s(ID,__FILE__,__LINE__)+" not implemented in Staremnant",__FILE__,__LINE__,sevnstd::notimplemented_error());

    return 0.;
}

double Staremnant::InertiaSphere(Star *s) const {
    return 0.4*s->getp(Mass::ID)*s->getp(Radius::ID)*s->getp(Radius::ID);
}

double Staremnant::age(Star *s) const {
    return s->getp(Worldtime::ID) - born_time;
}

double Staremnant::Mass(Star *s) const {
    return s->getp(Mass::ID);
}

/****BH****/
double BHrem::Radius(Star *s) const {
    return utilities::R_Schwarzschild(s->getp(Mass::ID));
}

double BHrem::Xspin(Star *s) const {

    if (s -> get_svpar_str("xspinmode")=="accretion"){
        return XspinAccretion(s);
    }

    return xspin;
}

double BHrem::estimate_Xspin(Star *s) const {

    const auto &spin_mode = s -> get_svpar_str("xspinmode");


    if (spin_mode == "disabled")
            return std::nan("");
    
    const double &z0 = s -> get_Z();
    const double &mco = s -> getp_0(MCO::ID); //Last value before the SN explosion, the SN explosion has already reset MCO at this point

    if (spin_mode == "geneva")
        return XspinGeneva(z0, mco);
    else if (spin_mode == "mesa")
        return XspinMESA(z0, mco);
    else if (spin_mode == "fuller")
        return XspinFuller();
    else if (spin_mode == "maxwellian"){
        const double &sigma_xspin = s -> get_svpar_num("xspin_sigma_maxwell");
        return XspinMaxwellian(sigma_xspin);
    }
    else if (spin_mode == "zeros"){
        return XspinZeros();
    }
    else if (spin_mode == "accretion"){
        return 0.0;
    }
    else{
        svlog.critical("xspinmode "+spin_mode+" not allowed",__FILE__,__LINE__,sevnstd::params_error(""));
    }

    return std::nan("");
}

void BHrem::default_initialiser(Star *s) {
    remnant_type=Lookup::Remnants::BH;
    set_Xspin(estimate_Xspin(s));
}

int BHrem::apply_Bavera_correction_to_Xspin(double period, double mass_wr) {
    double alpha, beta, period_days,_xspin;
    period_days = period * 365.25;


    if (period_days > 1.0)
        _xspin=xspin;
    else {
        /* These values are used if period and mass_wr are taken at C depletion. (Here we use the instant previous to the explotion.) */
        alpha = 0.029928 + exp(-0.282998 * mass_wr);
        alpha = -0.051237 / alpha;
        beta = 0.010905 + exp(-0.422213 * mass_wr);
        beta = -0.027090 / beta;
        _xspin = log10(period_days) * (alpha * log10(period_days) + beta);
        _xspin = std::min(std::max(0.0,_xspin),1.0); //Force xspin between 0 and 1
    }

    //Set the new value of _xspin
    set_Xspin(_xspin);

    return EXIT_SUCCESS;
}

double BHrem::XspinAccretion(Star *s) const {

    double Mratio = get_Mremnant_at_born() / s->getp(Mass::ID);
    double _xspin=1;
    if (Mratio>0.999999){
        _xspin=xspin;
    }
    else if (1/Mratio<=std::sqrt(6)){
            _xspin = std::sqrt(2./3.) * Mratio * (4. - std::sqrt(18*Mratio*Mratio-2));
    }

    return _xspin;
}


/****NS****/
double NSrem::Luminosity(Star *s) const {
    //cooling curve  LNS=0.02*M^(2/3)/(max(t,0.1)^2) from Eq. 93 using Hurley, 2000
    double Dt = std::max(age(s),0.1); //Time from the remnant creation
    return 0.02*std::pow(s->getp(Mass::ID),2./3.)/(Dt*Dt);
}

double NSrem::Inertia(Star *s) const {
    //Simple the ienrtia of a sphere 2/5 M R^2
    return InertiaSphere(s);
}

double NSrem::OmegaRem(_UNUSED Star *s) const {
    //CE: Fill the function, the return of this function is the Omega (in s) after DT from the NS birth
    //Useful constants
    const double yr = 3.1557600e7; //yr in seconds
    const double rsun = 6.95700e10; //rsun in cm
    const double msun = 1.98892e33; //msun in g

    //The first time it is called it is just Omega0
    if(age(s)==0){
        return Omega0;
    }
    //return s->getp(OmegaRem::ID);

    //TODO Radius and Inertia should be the one at the beginning of the evolution not their last values
    //In any case it does not matter a lot since Radius and Inertia will be constant during SSE
    const double R = s->getp(Radius::ID)*rsun; //Stellar radius at the beginning of the evolution
    const double I = s->getp(Inertia::ID)*rsun*rsun*msun; //Inertia  at the beginning of the evolution

    const double &Bi =  s->getp_0(Bmag::ID); //Initial B in Gauss in this timestep, Notice B should have been already evolved
    const double &Bf =  s->getp(Bmag::ID);//Final B in Gauss, Notice B should have been already evolved, we can also use Bmag(s)
    const double &Omegai = s->getp(OmegaRem::ID); //Initial value for Omega, Notice we use getp because we are evolving it now, so the  value of the previous step is still in getp
    const double & dt = s->getp(Timestep::ID)*1e6*yr; //Timestep
    const double taud = tau_magnetic*1e6*yr; //Magnetic field decay time scale in Myr;
    const double c_cm_s = utilities::c*rsun/yr; //Speed of light in cm/s

    double factor1 = 2*std::pow(R, 6)*sinalpha*sinalpha/(3*std::pow(c_cm_s, 3)*I);
    //I think it should be 8*M_PI*std::pow(R, 6)*sinalpha*sinalpha/(3*std::pow(c_cm_s, 3)*I);

    double factor2 = Bmin*Bmin*dt-taud*Bmin*(Bf-Bi)-taud/2*(Bf*Bf-Bi*Bi);

    double Omega = std::pow(pow(Omegai, -2)+2*factor1*factor2, -0.5);

    return Omega;
}

double NSrem::Bmag(_UNUSED Star *s) const {

    //CE: Fill the function, the return of this function is the Bmag (in Gauss) after DT from the NS birth
    //The first time it is called it is just B0
    if(age(s)==0){
        return B0;
    }
    //return s->getp(Bmag::ID);


    //Estimate vaule for this timestep
    const double& dt = s->getp(Timestep::ID);
    const double& Bold = s->getp(Bmag::ID); //Old value of the magnetic field
    const double taud = tau_magnetic; //Magnetic field decay time scale in Myr;

    double B = (Bold-Bmin)*std::exp(-dt/taud)+Bmin;

    return B;
}


void NSrem::default_initialiser(Star *s) {
    remnant_type=Lookup::Remnants::NS_CCSN;

    //Drawn a random alpha (we store the sinalpha)
    //Assuming an isotropic magnetic axis direction we have to uniformly sample in cos(alpha) between -1 and 1.
    double cosalpha = generate_uniform(-1,1);
    sinalpha = std::sqrt(1-cosalpha*cosalpha);

    //TODO Different options for the B0 and Omega0 initilisation are now hardcoded.
    //we should add a new set of parameters to do it at runtime or a more sophisticated method to add alternative NS models

    //Initiliase B0
    //Generate uniform (see Sgalletta+23)
    double B0min = 1e10; //CE: this is the minimum value for the uniform distribution of B0 in Gauss
    double B0max = 1e13; //CE: this is the maximum value for the uniform distribution of B0 Gauss
    Bmin = 1e8; //Bmin in Gauss
    B0 = generate_uniform(B0min,B0max);
    //Generate flat-in-log (see Sgalletta+23)
    //double logB0 = generate_uniform(10.0, 13.0);
    //B0 = std::pow(10.0, logB0);
    //Generate Faucher-Giguere: log normal distribution
    //double logB0 =  (std::normal_distribution<>{12.65, 0.55})(utilities::mtrand);
    //B0 = std::pow(10.0, logB0);

    //Initialise Omega0
    double Period0min = 0.01; //CE: this is the minimum value for the uniform distribution of Omega0 in s  (10 ms)
    double Period0max = 0.1; //CE: this is the maximum value for the uniform distribution of Omega0 in s (100 ms)
    double Period0 = generate_uniform(Period0min,Period0max);
    //Generate Faucher-Giguere: normal distrution with mean 0.3 s and variance 0.15 s (see Sgalletta+23)
    //double Period0 = (std::normal_distribution<>{0.3, 0.15})(utilities::mtrand);

    Omega0 = 2*M_PI/Period0;

    //Initialise
    tau_magnetic = s->get_svpar_num("ns_magnetic_tscale"); //Magnetic field decay time scale in Myr
}

void NSrem::print_log_message(Star *s) {
    std::string w = utilities::log_print("NS",s,remnant_type,get_Mremnant_at_born(),B0,Omega0,sinalpha);
    s->print_to_log(w);
}


/****WD****/
double WDrem::Radius(Star *s) const {
    //Eq.91 from Hurley+00
    const double RNS = NSrem::Rns; //Neutron star radius 10 km in Rsun, from Hurley, 2000
    const double &Mch = s->get_svpar_num("sn_Mchandra"); //Mchandra
    const double &Mass = s->getp(Mass::ID);

    double RWD = 0.0115*std::sqrt(pow(Mch/Mass,0.6666666667) -  pow(Mass/Mch,0.6666666667));

    return std::max(RNS,RWD);
}

double WDrem::Luminosity(Star *s) const {
    //Eq. 90 in Hurley+00
    double num = 635*s->getp(Mass::ID)*pow(s->get_Z(),0.4);
    double den = A_luminosity* pow(age(s) + 0.1,1.4);

    return num/den;
}


void Zombierem::initialise(Star *s) {
    remnant_type=Lookup::Remnants::Zombie;
    luminosity = s->getp(Luminosity::ID);
    radius = s->getp(Radius::ID);
}

