//
// Created by spera on 21/02/19.
//

#include <Orbit.h>
#include <Processes.h>
#include <binstar.h>
#include <errhand.h>
#include <supernova.h>
#include <utilities.h>


/*******************************************
********         TIDES       *************
*******************************************/

Orbital_change_Tides* Orbital_change_Tides::Instance(const std::string  &name)
{
	auto it = GetStaticMap().find(name);

	if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;

	return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

double Orbital_change_Tides::Mcnv(Star *star_wtides) const{
    //Just use property, they can automatically use tabled or fitting formulas from Hurley depending on the option tabuse_envconv
    return star_wtides->getp_0(Qconv::ID)*star_wtides->getp_0(Mass::ID);
}

double Orbital_change_Tides::Dcnv(Star *star_wtides) const{
    //Just use property, they can automatically use tabled or fitting formulas from Hurley depending on the option tabuse_envconv
    return star_wtides->getp_0(Depthconv::ID)*star_wtides->getp_0(Radius::ID);
}


double Orbital_change_Tides::compute_tconv_rasio96(Star *star_wtides, _UNUSED Binstar *b) {
    // Rasio+ 1996 Eq 4
    /// Should be depth of convective envelope
    double Denv_cnv = Dcnv(star_wtides);
    double Radx = std::min(star_wtides->getp_0(Radius::ID),b->Radx(star_wtides->get_ID()));
    double Rc = star_wtides->aminakedhelium() ? star_wtides->getp_0(RCO::ID) : star_wtides->getp_0(RHE::ID);
    Denv_cnv = std::max(std::min(Denv_cnv,Radx-Rc),1E-10); //Same as row 2339-2340 in MOBSE
    /// Should be mass of convective envelope
    double Menv_cnv = Mcnv(star_wtides);

    //cout << "Radius " << utilities::n2s(star_wtides->getp_0(Radius::ID), __FILE__, __LINE__) <<" "<< star_wtides->getp_0(RHE::ID)  << "; Mass " << utilities::n2s(star_wtides->getp_0(Mass::ID), __FILE__, __LINE__) << endl;
    //cout << "Denv_cnv " << utilities::n2s(Denv_cnv, __FILE__, __LINE__) << "; Menv_cnv " << utilities::n2s(Menv_cnv, __FILE__, __LINE__) << endl;
    //cout << "Mass " << utilities::n2s(star_wtides->getp_0(Mass::ID), __FILE__, __LINE__) << "; MHE " << utilities::n2s(star_wtides->getp_0(MHE::ID), __FILE__, __LINE__) << endl;
    //cout << "Mass N " << utilities::n2s(star_wtides->getp(Mass::ID), __FILE__, __LINE__) << "; MHE " << utilities::n2s(star_wtides->getp(MHE::ID), __FILE__, __LINE__) << endl;
    //cout << "Mzams " << star_wtides->get_zams()<< star_wtides->getp(Worldtime::ID)<<endl;

    /// Modification from Hurley+ 2002 Eq 31
    double t_cnv = Menv_cnv * Denv_cnv * (star_wtides->getp_0(Radius::ID) - 0.5 * Denv_cnv);
    t_cnv /= 3 * star_wtides->getp_0(Luminosity::ID) * utilities::LSun_to_Solar;  ///< l_t from LSun is converted to MSun RSun^2 yr^-3
    t_cnv = pow(t_cnv, 0.333333333333333333);

    //t_cnv *= 0.411; ///< random constant found in Hurley 2002


    if (t_cnv == 0 or std::isinf(t_cnv) or std::isnan(t_cnv)) {
        svlog.critical("Tconv has bad value: " + utilities::n2s(t_cnv, __FILE__, __LINE__), __FILE__, __LINE__);
    }
    return t_cnv;
}

double Orbital_change_Tides::compute_tconv_pwheel18(Star *star_wtides, _UNUSED Binstar *b) {
    // Price-Whelan+ 2018 Eq 9
    /// Should be mass of convective envelope
    double Menv_cnv = Mcnv(star_wtides);

    double t_cnv = 0.5 * pow(Menv_cnv, 0.333333333333333333)
                   * pow(star_wtides->getp_0(Temperature::ID)*2e-4, -1.333333333333333333); ///< In years
    if (t_cnv == 0 or std::isinf(t_cnv) or std::isnan(t_cnv)) {
        svlog.critical("Tconv has bad value: " + utilities::n2s(t_cnv, __FILE__, __LINE__), __FILE__, __LINE__);
    }
    return t_cnv;
}

double Orbital_change_Tides::compute_kt_zahn_conv(double Mass, double Menv_cnv, double t_cnv, double tide_freq){

    /*
    double conv_freq = 2*M_PI/t_cnv;
    //TODO This is an altenative implementation put here by Ale, talk to him and create a new tide option
    //GI 24/05/2022, I restored the old MOBSE-like implementaiton, I already check that it seems that the differences are negligible
    /// Vidal+ 2020
    double fratio = tide_freq/conv_freq;
    double fredfac;
    if (fratio < 1) {
        fredfac = 1.0;
    } else if (fratio < 5) {
        fredfac = fratio;
    } else {
        fredfac = fratio * fratio * 5;
    }
    */
    double fredfac=std::min(1.0,0.5*0.5*tide_freq*tide_freq/(t_cnv*t_cnv));


    double kt = 2 * fredfac * Menv_cnv / (21*t_cnv*Mass);
    return kt; //1/Myr

}

double Orbital_change_Tides::compute_kt_zahn_conv(Star *star_wtides, _UNUSED Binstar *binstar) {



    /// Should be mass of convective envelope
    double Menv_cnv = Mcnv(star_wtides);

    double t_cnv;
    if (star_wtides->table_loaded(Tconv::ID) and star_wtides->get_svpar_bool("tabuse_envconv")){
        t_cnv = star_wtides->getp_0(Tconv::ID);
    }else{
        t_cnv = get_tconv(star_wtides,binstar);
    }

    double omegaorb = 2*M_PI/binstar->getp(Period::ID); //yr^-1
    double  omegaspin = star_wtides->getp_0(OmegaSpin::ID); //yr^-1
    double  Ptid = 2*M_PI/fabs(omegaorb-omegaspin); //Equivalent to Eq. 33 in Hurley+02

    double fconv=std::min(1.0,0.25*Ptid*Ptid/(t_cnv*t_cnv)); //Eq. 32 in Hurley+02 (0.25=1/2^2)


    //return compute_kt_zahn_conv(star_wtides->getp_0(Mass::ID),Menv_cnv,t_cnv,Ptid);

    return  2 * fconv * Menv_cnv / (21*t_cnv*star_wtides->getp_0(Mass::ID)); //Eq. 30 in Hurley+02 yr^-1
}




double Orbital_change_Tides::compute_kt_zahn_rad(Star *star_wtides, Star *star_pert, Binstar *binstar)
{

    double E2;
    E2 = 1.58313e-9 * pow(star_wtides->getp_0(Mass::ID), 2.84);

    const double &mwtides = star_wtides->getp_0(Mass::ID);
    double Radx = std::min(star_wtides->getp_0(Radius::ID),binstar->Radx(star_wtides->get_ID()));
    //Use Radx instead of Radius same as row 2334 in MOBSE

    const double &a = binstar->getp(Semimajor::ID);

    double q = star_pert->getp_0(Mass::ID) / mwtides;


    double kt = sqrt(utilities::G * mwtides*Radx*Radx / (a*a*a*a*a)) *pow(1.+q, 0.833333333333333333)*E2;
    return kt;
}


double Orbital_change_Tides::dadt_hut(Star* star_wtides, Star *star_pert, Binstar *binstar,
                       double kt, double spin) {

    //If kt is 0 just return 0
    if (kt==0) return 0;

    double mwtides = star_wtides->getp_0(Mass::ID);
    double mpert = star_pert->getp_0(Mass::ID);
    double Radxwtides = std::min(star_wtides->getp_0(Radius::ID),binstar->Radx(star_wtides->get_ID()));
    double e = binstar->getp(Eccentricity::ID);
    double a = binstar->getp(Semimajor::ID);

    double q = mpert / mwtides;
    double ecc2 = e * e;
    double ome2 = 1.0 - ecc2;
    double r_over_a = Radxwtides/binstar->getp(Semimajor::ID);
    double omegaorb = 2*M_PI/binstar->getp(Period::ID);

    double da = -6 *kt*q*(1+q) * pow(r_over_a, 8) * a / pow(ome2, 7.5) * (f1(ecc2) - pow(ome2, 1.5) * f2(ecc2) * spin/omegaorb);

    return da;
}

double Orbital_change_Tides::dedt_hut(Star *star_wtides, Star *star_pert, Binstar *binstar, double kt, double spin)
{
    //If kt is 0 just return 0
    if (kt==0) return 0;

    const double &mwtides = star_wtides->getp_0(Mass::ID);
    const double &mpert = star_pert->getp_0(Mass::ID);
    double Radxwtides = std::min(star_wtides->getp_0(Radius::ID),binstar->Radx(star_wtides->get_ID()));
    const double &e = binstar->getp(Eccentricity::ID);

    double q = mpert / mwtides;
    double ecc2 = e * e;
    double ome2 = 1.0 - ecc2;
    double r_over_a = Radxwtides/binstar->getp(Semimajor::ID);
    double omegaorb = 2*M_PI/binstar->getp(Period::ID);

    double de = -27 *kt*q*(1+q) * pow(r_over_a, 8) * e/pow(ome2, 6.5) * (f3(ecc2) - 11./18.*pow(ome2, 1.5) * f4(ecc2) *spin/omegaorb);
    return de;
}

double Orbital_change_Tides::dspindt_hut(Star *star_wtides, Star *star_pert, Binstar *binstar, double kt, double spin)
{
    const double &mwtides = star_wtides->getp_0(Mass::ID);
    const double &mpert = star_pert->getp_0(Mass::ID);
    const double &e = binstar->getp(Eccentricity::ID);
    double Radxwtides = std::min(star_wtides->getp_0(Radius::ID),binstar->Radx(star_wtides->get_ID()));

    double q = mpert / mwtides;
    double ecc2 = e * e;
    double ome2 = 1.0 - ecc2;
    double r_over_a = Radxwtides/binstar->getp(Semimajor::ID);
    double omegaorb = 2*M_PI/binstar->getp(Period::ID);
    //Notice, in BSE/MOBSE we have that the Inertia is:
    // Inerta=k2*(M-Mc)*R*R + k3*Mc*Rc*rc and
    // rg2= k2 for MS stars, CHeburning (phase 4) and pureHE Ms (phase 7) (why for phase4 too?)
    // rg2= k2*(M-Mc)/M (why?) for star with phase <9 and not 1,4,7
    // rg2= k3 in the other cases.
    // Here in SEVN we just assume the real definition rg = I/(M*R*R), we use R instaed of radx because inertia is defined with R
    double rg2 = star_wtides->getp_0(Inertia::ID)/(mwtides * star_wtides->getp_0(Radius::ID) * star_wtides->getp_0(Radius::ID));


    double dospin = 3*kt*q*q/rg2 * pow(r_over_a, 6) * omegaorb/pow(ome2, 6) * (f2(ecc2) - pow(ome2, 1.5)*f5(ecc2) * spin/omegaorb);
    return dospin;
}


void Tides_simple::init(Binstar *binstar) {
    //Reset values
    _DA=0.;
    _DE=0.;
    _DLspin[0]=_DLspin[1]=0.;

    Star *s1 = binstar->getstar(0);
    Star *s2 = binstar->getstar(1);

    const double &e = binstar->getp(Eccentricity::ID);
    double omegaorb = 2*M_PI/binstar->getp(Period::ID);
    double dt = binstar->getp(BTimestep::ID)*1e6; //Myr to yr

    /// Here we ignore the effect of spin
    syncspin = pseudosync_spin(e*e, omegaorb); //Equilibrium angular momentum


    for(size_t i = 0; i < 2; i++) {

        /// Here we check if there are tides and, in case,
        /// Compute kt for each star
        if (s1->amiremnant() or s1->aminakedco()){
            kt[i] = 0;
        }
        else if (Mcnv(s1)==0){
            kt[i] = compute_kt_zahn_rad(s1, s2, binstar);
        }
        else{
            kt[i] = compute_kt_zahn_conv(s1, binstar);
        }

        if (kt[i]!=0){
            double dadt = dadt_hut(s1,s2,binstar,kt[i],s1->getp_0(OmegaSpin::ID));
            double dedt = dedt_hut(s1,s2,binstar,kt[i],s1->getp_0(OmegaSpin::ID));
            double dOSdt = dspindt_hut(s1,s2,binstar,kt[i],s1->getp_0(OmegaSpin::ID));



            //Spinlimiter
            //Notice that the  equilibrium tides assume that there is a Equilibrium rotational
            //velocity (Eq. 34 in Hurley+02) for which no angular momentum can be transferred (syncspin value).
            //In order to take into account this limit we use a method already used in SEVN1:
            //1-We estimate the time needed to reach the equilibrium rotation
            //considering the current rate of spin change (dOSdt)
            //2-If this timestep is smaller than the current timestep (BTimestep) we use it
            //as an effective timestep. In pratictice we consider that the tides acts until the spin
            //is equal to the equilibrium spin.
            //TODO see below
            //Notice: given our synchronous evolution, in principle  we should wait the effect of other processes
            //and correct a-posteriori the effect of tides. However this is maybe too complicated consider
            //the effect. Moreover, when the tides are very strong they usually dominated the choice of the Timestep
            //considering the change in the semimajor axis, eccentricity and spins.
            double limited_dt = (syncspin-s1->getp_0(OmegaSpin::ID))/dOSdt;
            //Notice, Eq. 26 in Hurley is always positive if syncspin>s1->getp_0(OmegaSpin::ID)
            //otherwise always negative, therefore here we check if limited_dt>0 otherwise we throw an error
            if (limited_dt < 0) svlog.critical("Change of angular velocity due to tides has a wrong sign",
                                               __FILE__,__LINE__);

            double dt_effecttive =  std::min(limited_dt,dt);


            _DA+=dadt*dt_effecttive;
            _DE+=dedt*dt_effecttive;


            //Rescale inertia to take into account Radx
            double Inertiax = s1->getp_0(Inertia::ID);
            double Radx = binstar->Radx(s1->get_ID());
            //Since Inertia propto R^2, rescale it by a factor Radx^2/Radius^2
            //if Radius>Radx
            if (s1->getp_0(Radius::ID)>Radx){
                Inertiax = Inertiax * (Radx*Radx)/(s1->getp_0(Radius::ID)*s1->getp_0(Radius::ID));
            }
            _DLspin[i]+=dOSdt*dt_effecttive*Inertiax;
            //if (s1->get_ID()==1) utilities::hardwait(dt,dt_effecttive,dOSdt*dt_effecttive,syncspin,s1->getp_0(OmegaSpin::ID));
        }

        utilities::swap_stars(s1, s2);
    }
}

double Tides_simple::DA(_UNUSED Binstar *binstar, _UNUSED int procID){

	return _DA;
}

double Tides_simple::DE(_UNUSED Binstar *binstar, _UNUSED int procID){
	return _DE;
}

double Tides_simple::DAngMomSpin(_UNUSED Binstar *b, _UNUSED int procID, int starID) {
    return _DLspin[starID];
}

double Tides_simple_notab::Mcnv(Star *star_wtides) const{
    if (star_wtides->getp_0(Phase::ID)<=Lookup::TerminalMainSequence or star_wtides->amiremnant())
        return 0.; //pure radiative envelope
    else
        return star_wtides->getp_0(Mass::ID) - star_wtides->getp_0(MHE::ID); //Whole envelope is convective

}

double Tides_simple_notab::Dcnv(Star *star_wtides) const{


    if (star_wtides->getp_0(Phase::ID)<=Lookup::TerminalMainSequence or star_wtides->amiremnant())
        return 0.; //pure radiative envelope
    else
        return star_wtides->getp_0(Radius::ID) - star_wtides->getp_0(RHE::ID); //Whole envelope is convective
    //Notice pureHE stars will  not enter in the first if, but they will return 0 in any case since (Radius=RHE), this
    //is ok we want the pureHE stars to have a pure radiative envelope
}



/*******************************************/

/*******************************************
********   GW    *************
*******************************************/
Orbital_change_GW* Orbital_change_GW::Instance(const std::string  &name)
{
    auto it = GetStaticMap().find(name);

    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;

    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

double Peters_gw::DE(Binstar *b, _UNUSED int procID){

    ///Eq.5.7 in Peters64

    //Using the value before the SSE
    double M1 = b->getstar(0)->getp_0(Mass::ID);
    double M2 = b->getstar(1)->getp_0(Mass::ID);
    double a  = b->getp(Semimajor::ID);
    double e  = b->getp(Eccentricity::ID);
    double dt = b->getp(BTimestep::ID) * utilities::Myr_to_yr; //G and c that we are using are in units of yr
    double edot_e,_DE;

    edot_e=Peters_GW_de(M1,M2,a,e); //Eq. 5.7 in Peters64
    _DE = edot_e * e * dt;  // De= edot/e  * e * dt;

    return _DE;
}

double Peters_gw::DA(Binstar *b, _UNUSED int procID){

    //Using the value before the SSE
    double M1 = b->getstar(0)->getp_0(Mass::ID);
    double M2 = b->getstar(1)->getp_0(Mass::ID);
    double a  = b->getp(Semimajor::ID);
    double e  = b->getp(Eccentricity::ID);
    double dt = b->getp(BTimestep::ID) * utilities::Myr_to_yr; //G and c that we are using are in units of yr
    double adot_a,_DA;

    adot_a=Peters_GW_da(M1,M2,a,e); //Eq. 5.6 in Peters64
    _DA = adot_a * a * dt;  // Da= adot/a  * a * dt;

    return _DA;

}

/*******************************************/

/*******************************************
********   Roche Lobe    ******************
*******************************************/

/**** Orbital change ****/
Orbital_change_RL* Orbital_change_RL::Instance(const std::string &name) {

    auto it = GetStaticMap().find(name);

    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;

    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}




///q
double Orbital_change_RL::estimate_q(Star *donor, Star *accretor){
    return donor->getp_0(Mass::ID)/accretor->getp_0(Mass::ID);
}

double Orbital_change_RL::Hfrac(Star *s) {
    //TODO Now we should have this value in the table Hsup
    return utilities::Hfrac(s);
}

double Orbital_change_RL::dMdt_eddington(Star *donor, Star *accretor, Binstar *b) {
    return utilities::dMdt_Eddington_accretion(donor,accretor,b->get_svpar_num("eddington_factor"));
}

///Time scales
double Orbital_change_RL::kelvin_helmotz_tscale(Star *s) {

    double tkh;
    double Mass = s->getp_0(Mass::ID);
    double R = s->getp_0(Radius::ID);
    double L = s->getp_0(Luminosity::ID);
    //Menvelope mass:
    //if M=MCO (naked CO) Menvelope=Mass
    //if M=MHE Menvelope = M-MCO
    //else Menvelope=M-MHE (if MHE, Menvelope=Mass)
    //COnsistent with Eq. 61 in Hurley+02
    double Menvelope=s->Menvelope(true);

    tkh = 1.0e7 * Mass * Menvelope/(R*L);  //MM in yr, units hardcoded, eq 61 hurley+2002

    return tkh;
}

double Orbital_change_RL::dynamic_tscale(Star *s) {

    double tdyn;
    double Mass = s->getp_0(Mass::ID);
    double R = s->getp_0(Radius::ID);

    tdyn = 5.05e-5 * sqrt((R*R*R)/Mass); //MM in yr, units hardcoded, eq 63 hurley+ 2002!!!


    return tdyn;
}


///Mass transfer
int Orbital_change_RL::speciale_evolve(_UNUSED Binstar *b){
    return dynamic_swallowing(b);
}


/**** Hurley ****/
void Hurley_rl::init_common_variable(Star *donor_star, Star *accretor_star, Binstar *b){



    donor = donor_star;
    accretor = accretor_star;

    star_type1=donor->get_bse_phase_0(),
    star_type2=accretor->get_bse_phase_0(); //BSE phases


    //Radii
    frl1          = fRL_radius(donor,b);  // RL over R1
    rl1           = frl1  * donor->getp_0(Radius::ID); //Radius of the Roche Lobe for the donor
    is_RL_in_core = donor->aminakedhelium() ?  rl1<=donor->getp_0(RCO::ID) : rl1<=donor->getp_0(RHE::ID); //Is the RL inside the Core?
    ///mt stuff - NOTICE, the interpretation of these two values depend on the mt_stability option
    ///if the option is from the qcrit family, mt_value=q and mt_tshold=qcrit,
    ///if the option is from the zeta family, mt_value=zeta_adiabatic and mt_tshold=zeta_RL
    mt_value     = mt_stability->get(donor,accretor,b);//estimate_q(donor, accretor);
    mt_tshold    = mt_stability->get_tshold(donor,accretor,b);//qcrit(donor,accretor);


    ///Accretion efficency
    f_MT = b->get_svpar_num("rlo_f_mass_accreted");
    dmdt_edd = dMdt_eddington(donor, accretor, b); //Eddingtont rate of the accretion on the accretor

    ///Times
    dt    = b->getp(BTimestep::ID)*utilities::Myr_to_yr; //Equations are in yr
    tkh1  = kelvin_helmotz_tscale(donor); //KH timescale of the star filling the RL
    tdyn1 = dynamic_tscale(donor); //Dynamical timescale of the star filling the RL
    taum  = std::sqrt(tkh1*tdyn1); //Time scale for dynamic Mass Loss through RL
    tkh2  = kelvin_helmotz_tscale(accretor); //KH timescale of the accretor
}

void Hurley_rl::init(Binstar *b){

    //First be sure to reset the Mass accreted/Donated
    //Initialise to default some quantities
    unset_mix(); //set local mix to false
    unset_comenv(); //set local comenv to false
    reset_DM();  //Set DM to 0
    reset_swallowed(); //Set is_swallowed to 0;
    novae = super_eddington =false;
    _DA = _DE = 0; //reset DA, DE
    _DLspin[0]=_DLspin[1]=0.; //reset spin
    unset_is_RLO_happening(); //Reset is_RLO_happening to false
    unset_is_colliding();

    ///NOW check for Mix or Common envelope due to both stars  filling the Roche Lobe
    //Notice the  check_doubleRLO updates the relevant local quantities (mix, ce)
    if (check_doubleRLO(b)){
        set_is_RLO_happening(); //RLO happening to true
    }
    ///NOW analyse the RLO
    else{
        Star* donor_star=b->getstar(0); //Star with the first star as donor
        Star* accretor_star=b->getstar(1); //Start with the second star as accretor



        //If just one (or none) of the two stars overflow the Roche Lobe
        for (size_t i=0; i<2; i++){

            //This is enefficient because we calculate a lot of properties that can be maybe not needed if the RLO is not happening
            //GI mobed this after the fRL_Radius check
            //init_common_variable(donor_star, accretor_star, b); //Init all the common variables, and then init again given that we have swaped the stars



            if (fRL_radius(donor_star,b)<=1){

                init_common_variable(donor_star, accretor_star, b); //Init all the common variables, and then init again given that we have swaped the stars


                //utilities::hardwait("W","R",donor->getp(Radius::ID),"R0",donor->getp_0(Radius::ID), "ID", donor->get_ID(),
                //                    "RLO",b->getp(RL0::ID),"RLO_0",b->getp_0(RL0::ID),
                //                    "RL1",b->getp(RL1::ID),"RL1_0",b->getp_0(RL1::ID),__FILE__,__LINE__);

                //Notice Collision is not checked here since thre is now a specific process Kollision
                //Check if stars collide at periastron, in case set the proper flag and break the cycle (all the set are made inside the function)
                //if(b->get_svpar_bool("rlo_enable_collision") and check_collision_at_periastron(donor, accretor, b)){
                //    set_is_colliding();
                //    break;
                //}

                RLO(b);
                set_orbital(b); //This is done only if fRL_radius<=1, there is no danger DA is repeated two times

                set_is_RLO_happening(); //RLO happening to true
            }
            utilities::swap_stars(donor_star, accretor_star); //Swap the stars and check the RLO for the other star.
        }
    }

    ///Now SET the binary  mix or comenv outcome
    //Notice, that RLO use its own local mix and comenv variable, then the binary one are updated only if needed (so that we dont' overwrite)
    if (get_mix() && get_comenv())
        svlog.critical("RLO overlof (Hurley) set both mix and comenv to true, this is not allowed",
                __FILE__,__LINE__, sevnstd::rl_error());
    //Now change only if mix or comenv needs to be set to true to avoid to overwrite mix or comenv that are already true
    else if (get_mix())
        b->mix = get_mix();
    else if (get_comenv())
        b->comenv = get_comenv();

    //Reset _DA _DE


}


///Aux
bool Hurley_rl::check_doubleRLO(Binstar *b){

    Star* star1=b->getstar(0); //Star with the first star as donor
    Star* star2=b->getstar(1); //Start with the second star as accretor

    //Note fRL_radius is Rstar/R_Roche_Lobe
    bool check_double_RLO = (fRL_radius(star1,b)<=1 && fRL_radius(star2,b)<=1); //Check if both stars are overfilling the Roche Lobe.

    //If both stars overfill the Roche Lobe  we obtain a Common Envelope or  a Mix depending on the bse phase.
    if (check_double_RLO){
        svlog.pdebug("Both stars are filling the Roche Lobe",__FILE__,__LINE__);
        outcome_double_RLO(b); //Set mix or cmnenv to true
        if (get_mix())
            svlog.pdebug("Star are mixing",__FILE__,__LINE__);
        else if (get_comenv())
            svlog.pdebug("Common envelope begins",__FILE__,__LINE__);
        else
            svlog.critical("One between mix or common envelope has to be true",__FILE__,__LINE__,sevnstd::rl_error());
        //utilities::wait("Starting a double Roche Overflow at time "+utilities::n2s(donor->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
        return true;
    }

    return false;
}

bool Hurley_rl::check_collision_at_periastron(Star *donor, Star *accretor, Binstar *b){


    //Note fRL_radius is Rstar/R_Roche_Lobe
    double peri = b->getp(Semimajor::ID)*(1-b->getp(Eccentricity::ID)); //Pericentric radius of the binary.
    bool check_collision = (donor->getp_0(Radius::ID) + accretor->getp_0(Radius::ID))>peri;  //Check if the stars collide at the pericenter.



    //If stars collide we obtain a Common Envelope or  a Mix depending on the bse phase.
    if (check_collision){
        svlog.pdebug("Stars collide",__FILE__,__LINE__);

        double rlo_donor, rlo_accretor;
        if (donor->get_ID()==0){
            rlo_donor = b->getp(RL0::ID);
            rlo_accretor = b->getp(RL1::ID);
        } else{
            rlo_donor = b->getp(RL1::ID);
            rlo_accretor = b->getp(RL0::ID);
        }

        std::string w = utilities::log_print("COLLISION",b,
                donor->get_ID(),donor->getp_0(Mass::ID),donor->getp_0(Radius::ID),donor->getp_0(Phase::ID),
                accretor->get_ID(),accretor->getp_0(Mass::ID),accretor->getp_0(Radius::ID),accretor->getp_0(Phase::ID),
                b->getp(Semimajor::ID),b->getp(Eccentricity::ID),rlo_donor,rlo_accretor);
        b->print_to_log(w);

        outcome_collision(b);
        if (get_mix())
            svlog.pdebug("Star are mixing",__FILE__,__LINE__);
        else if (get_comenv())
            svlog.pdebug("Common envelope begins",__FILE__,__LINE__);
        else
            svlog.critical("One between mix or common envelope has to be true",__FILE__,__LINE__,sevnstd::rl_error());
        utilities::wait("Collision at time "+utilities::n2s(donor->getp(Worldtime::ID),__FILE__,__LINE__),__FILE__,__LINE__);
        return true;
    }

    return false;

}

int Hurley_rl::outcome_double_RLO(Binstar *b) {

    //In order to mantain the formalism for which the primary is the most massive star
    //NB this is not needed if the condizion on the primary and secondary are mirrored (as here)
    Star * primary    = b->getstar(0)->getp(Mass::ID) >= b->getstar(1)->getp(Mass::ID) ? b->getstar(0) : b->getstar(1);
    Star * secondary = b->getstar(0)->getp(Mass::ID) >= b->getstar(1)->getp(Mass::ID) ? b->getstar(1) : b->getstar(0);

    //Donor and accretor bse phase
    int pbse = primary->get_bse_phase_0();
    int sbse = secondary->get_bse_phase_0();


    if ( ( (pbse>=3) and (pbse<=9) and (pbse!=7) ) or ( (sbse>=3) and (sbse<=9) and (sbse!=7) ) )
        set_comenv();
    else
        set_mix();


    return EXIT_SUCCESS;
}

int Hurley_rl::outcome_collision(Binstar *b) {

    Star * primary    = b->getstar(0)->getp(Mass::ID) >= b->getstar(1)->getp(Mass::ID) ? b->getstar(0) : b->getstar(1);
    Star * secondary = b->getstar(0)->getp(Mass::ID) >= b->getstar(1)->getp(Mass::ID) ? b->getstar(1) : b->getstar(0);


    //Donor and accretor bse phase
    int pbse = primary->get_bse_phase_0();
    int sbse = secondary->get_bse_phase_0();

    if ( ( (pbse>=3) and (pbse<=9) and (pbse!=7) ) or ( (sbse>=2) and (sbse<=9) and (sbse!=7) ) )
        set_comenv();
    else
        set_mix();

    return EXIT_SUCCESS;
}

int Hurley_rl::outcome_collision(Star *donor, Star *accretor, _UNUSED Binstar *b) {


    //Donor and accretor bse phase
    int pbse = donor->get_bse_phase_0();
    int sbse = accretor->get_bse_phase_0();

    if ( ( (pbse>=3) && (pbse<=9) && (pbse!=7) ) | ( (sbse>=2) && (sbse<=9) && (sbse!=7) ) )
        set_comenv();
    else
        set_mix();

    return EXIT_SUCCESS;
}

//TODO in SEVN1 we have a sort of repeat when R1 is too large wrt RL, but this check should be done in Timestep evolve
void  Hurley_rl::RLO(Binstar *b){

    double DM_donor=0, DM_accreted=0;

    /*********************************************************************************
    **MM: Sec 2.6.4 hurley+ 2002: if a star is <0.7 Msun it is completely convective
    ** if companion is also low mass, mass transfer is always dynamical till entire star is swallowed.
     * This will be dynamical mass transfer of a similar nature to common-envelope evolution.
    *********************************************************************************/
    //bool critical = is_critical();
    bool critical = mt_stability->mt_unstable(donor,accretor,b);

    //bool critical2 = mt_stability->mt_unstable(donor,accretor,b);
    //std::cout<<std::setprecision(20)<<donor->getp(Worldtime::ID)<< " " <<donor->get_bse_phase_0() <<" " << " " << mt_stability->get_tshold(donor,accretor,b)<< " "<< critical << std::endl;


    //TODO At the moment we don't have  an equivalent bse phase 0 star
    if( (star_type1==0) and (critical) ){

        //Entire star is swallowed
        set_swallowed(donor->get_ID());
        b->is_swallowed[donor->get_ID()]=true;
        utilities::wait("Dynamic unstable accretion from full convective star",mt_value,mt_tshold,__FILE__,__LINE__);

    }
    else if(!b->get_svpar_bool("rlo_mtstable_ms") and (star_type1==1 or star_type1==7) and critical){
        //Entire star is swallowed
        set_swallowed(donor->get_ID());
        b->is_swallowed[donor->get_ID()]=true;
        utilities::wait("Dynamic unstable accretion from MS-like star",mt_value,mt_tshold,__FILE__,__LINE__);
    }
        /*********************************************************************************/
        /**MM: UNSTABLE ACCRETION ON DYNAMICAL TIMESCALE: CALL COMMON ENVELOPE OR MERGE**/
        /*********************************************************************************/
    else if ((star_type1==2 && critical) and b->get_svpar_bool("optimistic_scenario_hg")){
        set_comenv();
    }
    else if( (star_type1==2 && critical) || (star_type1<=2 && star_type2<=2 && critical) || (!is_giant_like(star_type1) && frl1>10) ){
        //Case for MS or HG stars when qcrit or case for non giant like star where the radius is much larger than the RL1:
        // Note, here I group here three conditions that were separated in SEVN1 (the first in row 456 (and 469, repeat) the other two in row 660 (branch neutron star))
        set_mix();//b->mix=true;
    }
    else if(   (((star_type1==3) || (star_type1==5) || (star_type1==6) || (star_type1==8) || (star_type1==9)) && ((critical) || (is_RL_in_core))) || ((star_type1==4) && (critical))){
        //Case of evolved type
        set_comenv();//b->comenv=true;
    }
    else if (star_type1==13 || star_type1==14){
        //Case for NS and BH
        check_compact_accretor(star_type1, star_type2); //Check if the accretor is "more compact" than the donor
        set_mix();//b->mix = true;
    } else if(star_type1>=10 and star_type1<=12 and critical){
        ///Case for WDs
        //svlog.pdebug("Type SEVN", donor->getp(Phase::ID), "Type BSE", donor->get_bse_phase(),__FILE__,__LINE__);
        check_compact_accretor(star_type1, star_type2); //Check if the accretor is "more compact" than the donor

        if (star_type1==10 and star_type2==10) {
            /// Assume the energy released by ignition of triple=alpha is enough to destroy the star(s)
            accretor->tobe_exploded_as_SNI();
            //donor->set_empty_in_bse(b);
            //donor->set_empty();
            //The star is now empty, the binary system becomes broken or empty
            //b->set_broken();
            //b->set_onesurvived(); //If onesurvived is already true, the system becomes empty

            utilities::wait("SN1a triggered by tripleHa",__FILE__,__LINE__);
        }
        //GI 02/09/2022, I have disabled this since it is an heritage of BSE that it seems not useful anymore
        //else if ((b->get_svpar_num("eddington_factor")>10.0)){ //Hardcoded number to de-facto disable the eddington accretion and let the stars mix
        //    set_mix();//b->mix = true;
        //}
        else{
            set_swallowed(donor->get_ID());
            b->is_swallowed[donor->get_ID()]=true;
        }


    }
    /******************************************************************************/
    /** STABLE MASS TRANSFER ON THERMAL TIME SCALE OR NUCLEAR (take the minimum) (Hurley+2002, Sec. 2.6.2 and Sec. 2.6.3)**/
    /******************************************************************************/
    else{
        double dm_nuclear = nuclear_dmdt(donor, b) * dt; //Nuclear mass transfer
        if(is_giant_like(star_type1)){
            //Case of Giant and Giant like stars, Upper limit dm from nuclear mass trasnfer
            double dm_thermal = thermal_dmdt(donor, b) * dt; //Thermal mass transfer
            DM_donor = std::min(dm_nuclear,dm_thermal); //Thermal mass transfer cannot be more than nuclear
            //utilities::wait("Nuclear or thermal kGiant", b->mix, b->comenv, dm_nuclear, dm_thermal,  DM_donor, DM_accreted, __FILE__,__LINE__);

        }
        else{
            double dm_dynamic = dynamic_dmdt(donor, b) * dt; //Dynamic time scale (Eq. 60 Hurley+2002)
            DM_donor = std::min(dm_nuclear,dm_dynamic); //Thermal mass transfer cannot be more than dynamic
            //utilities::wait("Nuclear or thermal No Giant", b->mix, b->comenv, DM_donor, DM_accreted, __FILE__,__LINE__);

        }
        DM_accreted = thermal_nuclear_accreted_mass(DM_donor,donor,accretor,b);
    }


    //Apply the special correction (if any)
    DM_accreted = correct_accreted_mass(DM_accreted, DM_donor, donor, accretor, b);
    //If mix or comenv we don't need to set the DM
    if (!get_mix() && !get_comenv()){
        //Be sure that we are not exceeding the total mass.
        //This is just a pre-check to avoid to have negative masses, however
        // the method check_nakedHe_or_nakedCO_after_binary_evolution will take into account this properly.
        //Moreover, sometime the adaptive time step can just force a repeat step with a lower timestep
        if (DM_donor>=donor->getp_0(Mass::ID)){
            DM_donor=0.999999*donor->getp_0(Mass::ID); //Just to avoid to have exactly the Mass of the donor
            svlog.warning("The donated mass (" + utilities::n2s(DM_donor,__FILE__,__LINE__) +
            " Msun) is larger than the initial stellar mass ("+utilities::n2s(donor->getp_0(Mass::ID),__FILE__,__LINE__)+")"
            ,__FILE__,__LINE__);
            //svlog.critical("The donated mass (" + utilities::n2s(DM_donor,__FILE__,__LINE__) +
            //" Msun) is larger than the initial stellar mass ("+utilities::n2s(donor->getp_0(Mass::ID),__FILE__,__LINE__)+")"
            //,__FILE__,__LINE__,sevnstd::rl_error());
        }

        //DM_donor    = min(DM_donor, donor->getp_0(Mass::ID));
        DM_accreted = std::min(DM_accreted, DM_donor); //Safe condition, limit the accretion to the actual mass loss

        set_DM(-DM_donor,donor->get_ID());
        set_DM(DM_accreted,accretor->get_ID());
    }


    //Summary
    svlog.pdebug("********** \nEND RLO \n**********",
            "\n*WorldTime", donor->getp(Worldtime::ID),
            "\n*ID Donor", donor->get_ID(),
            "\n*ID Accretor", accretor->get_ID(),
            "\n*Mix?", get_mix(),
            "\n*CE?", get_comenv(),
            "\n*DM Donor", get_DM(donor->get_ID()),
            "\n*DM Accretor", get_DM(accretor->get_ID()),
            "\n*Tdyn", dynamic_tscale(donor),
            "\n*Tkh", kelvin_helmotz_tscale(donor),
            "\n*Phase Donor", "BSE: ",star_type1,donor->get_bse_phase_0(),"SEVN:", int(donor->getp(Phase::ID)),
            "\n*Phase Accretor", "BSE: ",star_type2,"SEVN:", int(accretor->getp(Phase::ID)),
            "\n*Mass Donor", donor->getp_0(Mass::ID),
            "\n*Mass Accretor", accretor->getp_0(Mass::ID),
            "\n*Mass HE Donor", donor->getp_0(MHE::ID),
            "\n*Mass HE Accretor", accretor->getp_0(MHE::ID),
            "\n*Mass CO Donor", donor->getp_0(MCO::ID),
            "\n*Mass CO Accretor", accretor->getp_0(MCO::ID),
            "\n*q", estimate_q(donor, accretor),
            "\n*qcrit", mt_tshold,
            "\n*R donor", donor->getp(Radius::ID),
            "\n*R accretor", accretor->getp(Radius::ID),
            "\n*RL donor", fRL_radius(donor, b)*donor->getp(Radius::ID),
            "\n*RL accretor", fRL_radius(accretor, b)*accretor->getp(Radius::ID),
            "\n",__FILE__,__LINE__);


    if (get_mix() || get_comenv() || is_swallowed(0) || is_swallowed(1)){
        utilities::wait("Dynamic unstable accretion",mt_value,mt_tshold,is_RL_in_core,"Mix?",get_mix(), "CE?",
                        get_comenv(), "Swallowed 1?", is_swallowed(0), "Swallowed 2?", is_swallowed(1), "Time",
                        b->getp(BWorldtime::ID),__FILE__,__LINE__);
    }

    return;
}

void  Hurley_rl::set_orbital(_UNUSED Binstar *b){

    ///Load current Mass transfer
    double M1 = std::abs(get_DM(0)), M2 = std::abs(get_DM(1));
    //Mdonated

    ///First check if actual mass transfer is happening
    if (M1>0 or M2>0){

        /************* DA ******************************/
        double gamma=b->get_svpar_num("rlo_gamma_angmom"); //parameter to set the orbital change behaviour

        double DM = std::abs(M1-M2); //Mass lost from the system during RLO
        double DM_donor    = donor->get_ID()==0 ? M1 : M2; //Absolute value of Mass lost from the primary
        double DM_accretor = accretor->get_ID()==0 ? M1 : M2; //Absolute value of Mass accreted on the secondary

        double Oorb = 2.0 * M_PI/b->getp(Period::ID); //Omega Orb
        double a = b->getp(Semimajor::ID); //Semimajor axis
        double ecc2 = b->getp(Eccentricity::ID)*b->getp(Eccentricity::ID);

        double Mdonor    = donor->getp_0(Mass::ID);
        double Maccretor = accretor->getp_0(Mass::ID);
        double Mtot = Mdonor + Maccretor;

        double djorb; //dJ
        double djspin_donor=0.; //Stellar Angmom lost by the donor
        double djspin_accretor=0.; //Stellar Angmom gained by the accretor
        double Angmom = b->getp(AngMom::ID); //Old AngMom


        ///Angmom lost from the system
        /*    When mass is lost from the system during RLOF there are now
        *    three choices as to how the orbital angular momentum is
        *    affected: a) the lost material carries with it a fraction
        *    gamma of the orbital angular momentum, i.e.
        *    dJorb = gamma*dm*a^2*omega_orb*; b) the material carries with it
        *    the specific angular momentum of the primary, i.e.
        *    dJorb = dm*a_1^2*omega_orb; or c) assume the material is lost
        *    from the system as if a wind from the secondary, i.e.
        *    dJorb = dm*a_2^2*omega_orb.
        *    The parameter gamma is an input option.
        *    Total angular momentum is conserved in this model.
        *  NOTE GI: this a comment from BSE, we actually consider also (1-e*e), but
        *  e is instantaneously put to 0.
        */

        djorb =  DM;
        if(super_eddington  or novae or gamma<-1.5) //ANGMOM lost from the secondary
            djorb *= (Mdonor*Mdonor)/(Mtot*Mtot); //(M1/Mtot)^2, because a2=Mdonor/Mtot*a
        else if(gamma>=0.0) //fraction of ANG MOM lost from the system
            djorb *= gamma;
        else //ANGMOM lost from the primary
            djorb *= (Maccretor*Maccretor)/(Mtot*Mtot); //(M2/Mtot)^2, because a1=Maccretor/Mtot*a

        djorb *= a*a*sqrt(1.0-ecc2)*Oorb;


        ///Angmom from donor Orbit-spin coupling
        double radx1=std::min(donor->getp_0(Radius::ID),b->Radx(int(donor->get_ID())));
        djspin_donor = DM_donor*radx1*radx1*donor->getp_0(OmegaSpin::ID);//AngMom spin lost by the donor
        djorb = djorb - djspin_donor; //Losing from start and giving to orbit (total J constant) !!djorb is positive but it will be subctracted!!



        ////Angmom from accretor Orbit-spin coupling
        /*
        * Determine whether the transferred material forms an accretion
        * disk around the secondary or hits the secondary in a direct
        * stream, by using eq.(1) of Ulrich & Burger (1976, ApJ, 206, 509)
        * fitted to the calculations of Lubow & Shu (1974, ApJ, 198, 383).
        * If disc forms:
        * Alter spin of the degenerate secondary by assuming that material
        * falls onto the star from the inner edge of a Keplerian accretion
        * disk and that the system is in a steady state.
        * Else:
        * No accretion disk.
        * Calculate the angular momentum of the transferred material by
        * using the radius of the disk (see Ulrich & Burger) that would
        * have formed if allowed.
        */
        double q2=1/estimate_q(donor, accretor);
        double rdisc = 0.0425*a*pow((q2*(1.0+q2)),(1.0/4.0));
        double r_accretion = rdisc>accretor->getp_0(Radius::ID) ? accretor->getp_0(Radius::ID) :  1.7*rdisc;
        djspin_accretor = DM_accretor * sqrt(utilities::G*Maccretor*r_accretion); //AngMom accreted by the


        //djspint2 = djspint2 - djt;
        djorb = djorb + djspin_accretor; //Losing from orbit and giving to star (total J constant) !!djorb is positive but it will be subctracted!!


        ////Update Angmom
        Angmom -= djorb; //Variation of ANGMOM




        ///Now estimate _DA as anew - a, where anew estimated from the  updated AngMom
        //NB In order to be consistent the estimate of DA is made considering the new masses after the rlo (without considering other processes at this stage).
        double Mdonor_after_rlo = Mdonor - DM_donor;
        double Maccretor_after_rlo = Maccretor + DM_accretor;
        double Mtot_after_rlo = Mdonor_after_rlo + Maccretor_after_rlo;
        _DA = Angmom*Angmom * Mtot_after_rlo /  ( utilities::G * (1-ecc2) * Mdonor_after_rlo*Mdonor_after_rlo*Maccretor_after_rlo*Maccretor_after_rlo ) - a;
        /*******************************************/

        /************* DE  ******************************/
        //DE is put suddenly to 0
        //_DE = -1e30;
        //DE is simply not modified
        //b->disable_e_check=true; //The eccentricity here can drop from some value to 0, so disable the adaptive time step check.
        /*******************************************/

        /*************** DJSpin  ******************************/
        _DLspin[donor->get_ID()]    =  -djspin_donor; //djspin_donor is in absolute value, but it is lost from the donor
        _DLspin[accretor->get_ID()] =  djspin_accretor; //djspin_accretor is in absolute value

    } else{
        _DA = _DE = _DLspin[0] = _DLspin[1] = 0.;  //reset orbital
    }



}

double Hurley_rl::DA(_UNUSED Binstar *b, _UNUSED int procID){return _DA;}

/**
 * Variation of Eccentricity. In BSE derived from Hurley+20, the
 * eccentricity is set to 0 at onset of RLO. In order to take this
 * into account we set the DE to a large negative value so that e<0 and
 * in the subsequent checks it is forced to 0.
 * @param b Pointer to binary
 * @param procID ID of process calling DE.
 * @return
 */
double Hurley_rl::DE(_UNUSED Binstar *b, _UNUSED int procID){ return _DE;}

double Hurley_rl::DAngMomSpin(_UNUSED Binstar *b, _UNUSED int procID, int starID) {
    return _DLspin[starID];
}


//Mass transfer
int Hurley_rl::dynamic_swallowing(_UNUSED Binstar *b) {

    //SAFE CHECK: If the dynamic unstable mass transfer already trigger mix or ce skip the dynamic
    if (get_mix() or get_comenv()){
        return EXIT_SUCCESS;
    }


    Star *_donor_star = b->getstar(0);
    Star *_accretor_star = b->getstar(1);

    for (size_t i=0; i<2; i++){


        init_common_variable(_donor_star, _accretor_star, b); //Init all the common variables, and then  in the next cycle step init again given that we have swaped the stars

        if (is_swallowed(donor->get_ID())){


            //NB we use getp and not getp_0 because this function has to be called as special evolve, after normal binary evolution
            //THis mass takes already into account the Mass loss (or accreted) from Winds
            double DM_donor=donor->getp(Mass::ID); //The star has been swallowed. Transfer the current mass.
            double DM_accreted=0;

            //If it is a naked helium or a nakedco or a remnant do not accreate anything (we assume the excess mass is rapidly lost by winds or removed through jets)
            if (accretor->aminakedhelium() or accretor->aminakedco() or accretor->amiremnant()){
                DM_accreted=0;
            }
            //If it is not a degenerate object accrete just a fraction of the RLO mass from the donor.
            else if (star_type2<10 ){
                DM_accreted = f_MT*DM_donor;
            }
            //else { //From Steady Accretion on to degenerate object (Hurley+20)
            //    DM_accreted = std::min(dmdt_edd * taum, DM_donor);
            //    //Check if we are accreting at super eddington rate
            //    if (  (b->get_svpar_num("eddington_factor")>1)) super_eddington =true;
            //}

            set_DM(-DM_donor,donor->get_ID());
            set_DM(DM_accreted,accretor->get_ID());
            reset_swallowed();
            break; //Not needed to check the other star, Double RLO is taken into account elsewhere
        }
        else
            utilities::swap_stars(_donor_star, _accretor_star); //Swap the stars and check the RLO for the other star.

    }

    return EXIT_SUCCESS;


}

double Hurley_rl::dynamic_dmdt(Star *s, _UNUSED Binstar *b){
    //Eq. 62 in Hurley+2002
    double tdyn = dynamic_tscale(s);
    double M   = s->getp_0(Mass::ID);

    double dm_dynamic = M/tdyn;

    return dm_dynamic;
}

double Hurley_rl::nuclear_dmdt( Star *s, Binstar *b){
    //Note dm1 is from Hurley+2002 , Eq. 58, 59, but in the original version it is min(M,5)**2 instead of M**2
    //Howver in SEVN1 we used just M/
    double star_type = s->get_bse_phase_0();
    double frl1  = fRL_radius(s,b);  // RL over R1
    double lnR = std::log(1/frl1); //ln(R1/RL1)
    double M   =   std::min(s->getp_0(Mass::ID),s->get_svpar_num("rlo_max_nuclearmt"));
    double dm_nuclear = 3.0e-6*lnR*lnR*lnR*M*M; //Nuclear mass transfer

    if (star_type==2){
        //It seems not be present in Hurley+2002 paper, but it is in BSE and SEVN1
        double Mcore = s->getp_0(MHE::ID);
        double Mew = (s->getp_0(Mass::ID)-Mcore)/s->getp_0(Mass::ID);
        dm_nuclear *= std::max(Mew,0.01);
    }
    else if(star_type>=10){
        //The version in the paper (Hurley+2002, Sec. 2.6.2) is just 1.0e3/max(r,1.0e-4) without mass,
        //but both in BSE and SEVN1 we have 1.0e3*M/max(r,1.0e-4)
        dm_nuclear *= 1.0e3*M/std::max(s->getp_0(Radius::ID),1.0e-4);
    }

    return dm_nuclear;

}

double Hurley_rl::thermal_dmdt( Star *s, _UNUSED Binstar *b){
    //Eq. 60 in Hurley+2002
    double tkh = kelvin_helmotz_tscale(s);
    double M   = s->getp_0(Mass::ID);

    double dm_thermal = M/tkh;

    return dm_thermal;

}

double Hurley_rl::thermal_nuclear_accreted_mass(double DM_donor, Star *donor, Star *accretor, Binstar *b){

    //int star_type1 = donor->get_bse_phase_0();
    int star_type2 = accretor->get_bse_phase_0();
    double dt      = b->getp(BTimestep::ID)*utilities::Myr_to_yr; //Equations are in yr

    double epsnova = b->get_svpar_num("rlo_eps_nova");
    double f_MT = b->get_svpar_num("rlo_f_mass_accreted");
    double dmdt_edd = dMdt_eddington(donor, accretor, b); //Eddington rate of the accretion on the accretor

    ///Mass accretion

    //PureHe stars do not accrete. This because they will then loose this extra mass soo.
    if (accretor->aminakedhelium() or accretor->aminakedco())
        return 0.0;

    //Simplified accretion, Accrete just a fraction instead of the more complicated rules of Hurley+2002
    //MM 01/05/2020
    double DM_accretor=f_MT*DM_donor;

    //Accretion on a compact object limited by the Eddington accretion
    if (star_type2>=10){
        DM_accretor = std::min(dmdt_edd*dt, DM_accretor);

        //Check if we are accreting at super eddington rate
        if (  (b->get_svpar_num("eddington_factor")>1)) super_eddington =true;
    }

    //Accretion on a WD
    if (star_type2>=10 && star_type2<=12){

        //If the material comes from a non-WD stars,  the material is H-reach and the accretion creates a nova
        if (donor->whatamidonating()==Material::H){
        //if (star_type1<=6){
            DM_accretor *= epsnova;
            //This is a nova!
            novae = true;
        }

     }


    return DM_accretor;
 }

 double Hurley_rl::fRL_radius(Star *star, Binstar *b){

    //Remember binary proprety T are equivalent to  stellar property T_0 during the binary evolution (because sse have been already updated)
     if (star->get_ID()==0)
         return b->getp(RL0::ID)/star->getp_0(Radius::ID); //Star property previously of the SSE evolution
     else if (star->get_ID()==1)
         return b->getp(RL1::ID)/star->getp_0(Radius::ID);
     else
         svlog.critical("Got a star in binary with ID="+utilities::n2s(star->get_ID(),__FILE__,__LINE__)+
         ", only 0 and 1 allowed in a binary",__FILE__,__LINE__,sevnstd::bse_error());

     return -1.0;
 }


/**** Hurley  BSE ****/
double Hurley_rl_bse::thermal_nuclear_accreted_mass(double DM_donor, Star *donor, Star *accretor, Binstar *b){

     //int star_type1 = donor->get_bse_phase_0();
     int star_type2 = accretor->get_bse_phase_0();
     double dt      = b->getp(BTimestep::ID)*utilities::Myr_to_yr; //Equations are in yr

     double epsnova = b->get_svpar_num("rlo_eps_nova");
     double dmdt_edd = dMdt_eddington(donor, accretor, b); //Eddington rate of the accretion on the accretor

     double DM_accretor = 0.0;

     ///Mass accretion
    //Here we follow the MOBSE and old SEVN implementation, i.e. for MS,HG and CHEB star the accretion is limited by the
    //timescale of the secondary, we allow also accretion on nakedHe stars but only if the accreted mass is He
    if (star_type2<=2 or star_type2==4
    or (star_type2>=7 and star_type2<=9 and donor->whatamidonating_0()==Material::He)){
        double taum  = accretor->getp_0(Mass::ID)/DM_donor*dt; //Eq. 65 in Hurley+02
        double tkh2  = kelvin_helmotz_tscale(accretor);
        DM_accretor  = std::min(1.0,10.0*taum/tkh2)*DM_donor; //MOBSE row 1901
    }
    //PureHe stars do not accrete. This because they will then lose this extra mass soon.
    else if(accretor->aminakedhelium() or accretor->aminakedco()){
        return 0.0;
    }
    //For all the other stars add all the material to the giant envelope
    else{
        //From mobse: We have a giant whose envelope can aborb any transferred material.
        DM_accretor=DM_donor; //MOBSE row 2100
    }


     //Accretion on a compact object limited by the Eddington accretion
     if (star_type2>=10){
         DM_accretor = std::min(dmdt_edd*dt, DM_accretor);

         //Check if we are accreting at super eddington rate
         if (  (b->get_svpar_num("eddington_factor")>1)) super_eddington =true;
     }

     //Accretion on a WD
     if (star_type2>=10 && star_type2<=12){

         //If the material comes from a non-WD stars,  the material is H-reach and the accretion creates a nova
         if (donor->whatamidonating()==Material::H){
             //if (star_type1<=6){
             DM_accretor *= epsnova;
             //This is a nova!
             novae = true;
         }

     }

     return DM_accretor;
 }

int Hurley_rl_bse::dynamic_swallowing(Binstar *b) {
    //SAFE CHECK: If the dynamic unstable mass transfer already trigger mix or ce skip the dynamic
    if (get_mix() or get_comenv()){
        return EXIT_SUCCESS;
    }

    Star *_donor_star = b->getstar(0);
    Star *_accretor_star = b->getstar(1);

    for (size_t i=0; i<2; i++){

        //donor, accretor, taum, tkh2, dmdt_edd are all initiliased here
        init_common_variable(_donor_star, _accretor_star, b); //Init all the common variables, and then  in the next cycle step init again given that we have swaped the stars

        if (is_swallowed(donor->get_ID())){
            //NB we use getp and not getp_0 because this function has to be called as special evolve, after normal binary evolution
            //THis mass takes already into account the Mass loss (or accreted) from Winds
            double DM_donor=donor->getp(Mass::ID); //The star has been swallowed. Transfer the current mass.
            double DM_accreted=0;

            if (star_type2<=1){
                DM_accreted  = taum/tkh2*DM_donor;
            }
            //PureHe stars do not accrete. This because they will then loose this extra mass soon.
            else if(accretor->aminakedhelium() or accretor->aminakedco()){
                DM_accreted=0.0;
            }
            //From Steady Accretion on to degenerate object (Hurley+20)
            else if(star_type2>10){
                DM_accreted = std::min(dmdt_edd * taum, DM_donor);
                //Check if we are accreting at super eddington rate
                if (  (b->get_svpar_num("eddington_factor")>1)) super_eddington =true;
            }
            //For all the other stars add the material to the giant envelope
            else{
                DM_accreted=DM_donor;
            }

            set_DM(-DM_donor,donor->get_ID());
            set_DM(DM_accreted,accretor->get_ID());
            reset_swallowed();
            break; //Not needed to check the other star, Double RLO is taken into account elsewhere
        }
        else
            utilities::swap_stars(_donor_star, _accretor_star); //Swap the stars and check the RLO for the other star.

    }

    return EXIT_SUCCESS;

}



/**** Hurley  mod ****/
//TODO I think we don t need this anymore
void  Hurley_mod_rl::set_orbital(_UNUSED Binstar *b){

    ///Load current Mass transfer
    Star *donor = get_donor();
    Star *accretor = get_accretor();

    double M1 = std::abs(get_DM(0)), M2 = std::abs(get_DM(1));
    //Mdonated

    ///First check if actual mass transfer is happening
    if (M1>0 || M2>0){

        /************* DA ******************************
         *
         * da/a =  2 dJ/J + 2dMd/Md -2dMa/Ma  - (Md-Ma)/Mtot
         *
         * where dMd is the mass lost from the donor (absolute value)
         * and dMa the mass accreted on the accretor.
         */

        double DM = std::abs(M1-M2); //Mass lost from the system during RLO
        double DM_donor    = donor->get_ID()==0 ? M1 : M2; //dMd
        double DM_accretor = accretor->get_ID()==0 ? M1 : M2; //dMa
        double Mdonor    = donor->getp_0(Mass::ID); //Md
        double Maccretor = accretor->getp_0(Mass::ID); //Ma
        double Mtot = Mdonor + Maccretor; //Mtot
        double da_a = 0;  //da/a



        ///dJ
        double gamma=b->get_svpar_num("rlo_gamma_angmom"); //parameter to set the orbital change behaviour
        double Oorb = 2.0 * M_PI/b->getp(Period::ID); //Omega Orb
        double a = b->getp(Semimajor::ID); //Semimajor axis
        double ecc2 = b->getp(Eccentricity::ID)*b->getp(Eccentricity::ID);


        double djorb; //dJ
        double djt; //Variation of dj due to the stellar spin
        //double Angmom = b->getp(AngMom::ID); //Old AngMom


        ///Angmom lost from the system
        /*    When mass is lost from the system during RLOF there are now
        *    three choices as to how the orbital angular momentum is
        *    affected: a) the lost material carries with it a fraction
        *    gamma of the orbital angular momentum, i.e.
        *    dJorb = gamma*dm*a^2*omega_orb*; b) the material carries with it
        *    the specific angular momentum of the primary, i.e.
        *    dJorb = dm*a_1^2*omega_orb; or c) assume the material is lost
        *    from the system as if a wind from the secondary, i.e.
        *    dJorb = dm*a_2^2*omega_orb.
        *    The parameter gamma is an input option.
        *    Total angular momentum is conserved in this model.
        *  NOTE GI: this a comment from BSE, we actually consider also (1-e*e), but
        *  e is instantaneously put to 0.
        */

        djorb =  DM;
        if(super_eddington  || novae || gamma<-1.5) //ANGMOM lost from the secondary
            djorb *= (Mdonor*Mdonor)/(Mtot*Mtot); //(M1/Mtot)^2
        else if(gamma>=0.0) //fraction of ANG MOM lost from the system
            djorb *= gamma;
        else //ANGMOM lost from the primary
            djorb *= (Maccretor*Maccretor)/(Mtot*Mtot); //(M2/Mtot)^2

        djorb *= a*a*sqrt(1.0-ecc2)*Oorb;

        ///Angmom from donor Orbit-spin coupling
        double radx1=std::max(donor->getp_0(RHE::ID),get_rldonor());
        djt = DM_donor*radx1*radx1*donor->getp_0(Spin::ID)/donor->getp_0(Inertia::ID);
        //djspin= djspint -djt
        djorb = djorb - djt; //Losing from start and giving to orbit (total J constant) !!djorb is positive but it will be subctracted!!

        ////Angmom from accretor Orbit-spin coupling
        /*
        * Determine whether the transferred material forms an accretion
        * disk around the secondary or hits the secondary in a direct
        * stream, by using eq.(1) of Ulrich & Burger (1976, ApJ, 206, 509)
        * fitted to the calculations of Lubow & Shu (1974, ApJ, 198, 383).
        * If disc forms:
        * Alter spin of the degenerate secondary by assuming that material
        * falls onto the star from the inner edge of a Keplerian accretion
        * disk and that the system is in a steady state.
        * Else:
        * No accretion disk.
        * Calculate the angular momentum of the transferred material by
        * using the radius of the disk (see Ulrich & Burger) that would
        * have formed if allowed.
        */
        //double q2=estimate_q(accretor, donor);
        //double rdisc = 0.0425*a*pow((q2*(1.0+q2)),(1.0/4.0));
        //double r_accretion = rdisc>accretor->getp_0(Radius::ID) ? accretor->getp_0(Radius::ID) :  1.7*rdisc;
        djt =  DM_accretor * sqrt(utilities::G*Maccretor*accretor->getp_0(Radius::ID));
        //djspint2 = djspint2 - djt;
        djorb = djorb + djt; //Losing from orbit and giving to star (total J constant) !!djorb is positive but it will be subctracted!!


        ////Update a
        //NB djorb is assumed positive in the estimate, but it means the J lost, so its real value is negative.
        da_a += -2*djorb/b->getp(AngMom::ID);
        //DM donated
        da_a += 2*DM_donor/Mdonor;
        //DM accreted
        da_a -= 2*DM_accretor/Maccretor;
        //DM lost
        da_a -= DM/Mtot;


        ///Now estimate _DA as anew - a, where anew estimated from the  updated AngMom
        //NB actually it is not consistent to estimate anew in that way, Mdonor and Maccretor are the old one, but they have been updated from the binary evolution.
        set_DA(a*da_a);
        /*******************************************/

        /************* DE  ******************************/
        //DE is put suddenly to 0 following Hurley+02 and BSE -> syncronistation at the onset of RLO.
        set_DE(-1e30);
        b->disable_e_check=true; //The eccentricity here can drop from some value to 0, so disable the adaptive time step check.
        /*******************************************/

    } else{
        //reset orbital
        set_DA(0.0);
        set_DE(0.0);
    }



}

/*******************************************/


/*******************************************
********   Mix    *************
*******************************************/
Orbital_change_Mix* Orbital_change_Mix::Instance(const std::string  &name)
{
    auto it = GetStaticMap().find(name);

    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;

    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

/*******************************************/

/*******************************************
********   SN Kicks    *************
*******************************************/
Orbital_change_SNKicks* Orbital_change_SNKicks::Instance(const std::string &name) {
    auto it = GetStaticMap().find(name);

    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;

    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

//kicks
void Hurley_SNKicks::init(_UNUSED Binstar *b){

    double DA=0, DE=0;
    double a_in = b->getp(Semimajor::ID), ecc_in = b->getp(Eccentricity::ID);
    double a_fin=0, ecc_fin=0;
    Star *sn = b->getstar(0);
    Star *other = b->getstar(1);

    for (size_t i=0; i<2; i++){
        a_fin=ecc_fin=0;
        kick_star(sn, other, b, a_fin, ecc_fin);
        DA += a_fin - a_in;
        DE += ecc_fin - ecc_in;
        utilities::swap_stars(sn,other);//Swap the stars and check if the other star explodes
    }

    set_DA(DA);
    set_DE(DE);

}

double Hurley_SNKicks::DA(_UNUSED Binstar *b, _UNUSED int procID){return _DA;}
double Hurley_SNKicks::DE(_UNUSED Binstar *b, _UNUSED int procID){return _DE;}

int Hurley_SNKicks::kick_star(Star *sn, Star *other, Binstar *binstar, double& a_fin, double& ecc_fin) {

    //No explosion just return
    if(sn->getp(Phase::ID) < Lookup::Remnant or sn->kicked()){
        a_fin=binstar->getp(Semimajor::ID);
        ecc_fin=binstar->getp(Eccentricity::ID);

        return 0;
    }



    ///*************** Variables **********************///
    //Orbital parameters
    double ecc = binstar->getp(Eccentricity::ID);
    double a = binstar->getp(Semimajor::ID);
    double vorb, vfin; // relative velocity at beginning, relative speed at the end, velocity
    double r; //  distance of the two stars at the moment of the sn explosion

    // define some angles
    double  sin_beta, cos_beta; //Angle between the separation vector and the relative velocity of the stars (orbit in a plane)
    double  _UNUSED sin_phi, cos_phi; //Zenithal angle of the Kick velocity wrt to the orbital plane (Fig. A1 in Hurley+02)
    double sin_omega, cos_omega; //Azimuthal angle of the Kick velocity wrt to the orbital plane (Fig. A1 in Hurley+02)
    double _UNUSED cos_nu; //Angle between the old and the new orbital plane


    //Variables used for the cross product of the separation vector and the new relative velocity (specific ang mom, hn = r x vn) vn=(vi,vj,vk) but vj is not needed since r is aligned along r so hj=0.
    double hn2; //magnitude of the new ang mom (after kick) to the power of 2
    double vx=0.,vx2=0., vy=0., _UNUSED vy2=0, vz=0., vz2=0; //vi, vk components to the power of 2 of the new velocity  (after kick, Eq. A12 in Hurley 02).
    double vk, vx_kick=0., vy_kick=0., vz_kick=0.; //variable to store the kick velocity in new units
    double _UNUSED vcom,vcomx,vcomy,vcomz; //New centre-of-mass velocity
    double Mbin_new, Mbin_old;//new and old binary mass
    Mbin_new = sn->getp(Mass::ID) + other->getp(Mass::ID);
    Mbin_old = sn->getp_0(Mass::ID) + other->getp_0(Mass::ID);
    ///***********************************************///




    //************ find the initial separation randomly choosing a mean anomaly *******/
    double X1 = _uniform_real(utilities::mtrand); //between 0.0 and 1.0
    double mm = X1*2.0*M_PI;       // mean anomaly
    double em = utilities::kepler(ecc, mm);   // eccentric anomaly
    r = a*(1.0 - ecc*cos(em));  // separation
    ///***********************************************///

    //******** find the initial relative velocity vector (eq A2 & A3 in Hurley 2002) *******/
    vorb = sqrt( (G* Mbin_old)* (2.0/r - 1.0/a));
    ///***********************************************///

    //******** Angles between r and the initial v *******/
    sin_beta = sqrt ( (a*a*(1.0 - ecc*ecc)) / (r*(2.0*a - r)) ); //Eq. A2 Hurley+02
    cos_beta = (-1.0*ecc*sin(em)) / sqrt( 1.0 - ecc*ecc*cos(em)*cos(em) ); //Eq. A3 Hurley+02
    //Sanity check //TODO DEBUG to be removed
    if (sin_beta*sin_beta+cos_beta*cos_beta-1>1e-10)
        svlog.critical("Sanity check failed in kicks: sin beta^2 +cos beta^2 !=1",__FILE__,__LINE__,sevnstd::sanity_error());
    //utilities::wait("Check ecc",sin_beta*sin_beta+cos_beta*cos_beta,__FILE__,__LINE__);
    ///***********************************************///


    //******** convert kick units from km/s to Rsun/yr ******/
    if (std::isnan(sn->vkick[3])){
        vk=0;
        //Just force the components to be 0
        sn->vkick[0] = sn->vkick[1] = sn->vkick[2] = sn->vkick[3] = 0.0; //No kick
        vx_kick = vy_kick = vz_kick = vk = 0.0; //Redundant but put here for clarity
    } else{
        //Transform the kicks units in Rsun/yr to be consistent with the rest of the code
        vx_kick = sn->vkick[0]*utilities::kms_to_RSunyr;
        vy_kick = sn->vkick[1]*utilities::kms_to_RSunyr;
        vz_kick = sn->vkick[2]*utilities::kms_to_RSunyr;
        vk = sn->vkick[3]*utilities::kms_to_RSunyr; // modules
    }
    //utilities::hardwait("SNKicks",sn->vkick[0],sn->vkick[3],__FILE__,__LINE__);
    ///***********************************************///

    //Separate between 0 and not 0 kicks;
    if (vk==0){
        vfin = vorb; //No kick, Vfin is still the initial Vorb
        vx = -vorb*sin_beta; //i(or x) component  of the new velocity after kick (Eq. A12 in Hurley 02 or Fig. A1)
        vx2 = vx*vx; //i(or x) component  of the new velocity after kick (Eq. A12 in Hurley 02 or Fig. A1)
        vy = -vorb*cos_beta; //j(or y) component  of the new velocity after kick  (Eq. A12 in Hurley 02 or Fig. A1)
        vy2 = vy*vy; //j(or y) component  of the new velocity after kick  (Eq. A12 in Hurley 02 or Fig. A1)
        vz = vz2 = 0.; //k( or z) component  of the new velocity after kick   (Eq. A12 in Hurley 02 or Fig. A1)
    }
    else{
        double sqrt_xy = sqrt(vx_kick*vx_kick + vy_kick*vy_kick);
        cos_phi = sqrt_xy / vk;
        sin_phi = vz_kick / vk;
        cos_omega = vx_kick / sqrt_xy;
        sin_omega = vy_kick  / sqrt_xy;
        vfin = sqrt(vk*vk + vorb*vorb - 2.0*vk*vorb*(cos_omega*cos_phi*sin_beta + sin_omega*cos_phi*cos_beta)); //Module of the new velocity Eq. A8 in Hurley+2002

        ///New velocity
        //i-x component
        vx  = (vx_kick - vorb*sin_beta); //vx kick (vk cos_omega cos_phi) is just sn->vkick[0]
        vx2 = vx*vx; //i(or x) component  of the new velocity after kick (Eq. A12 in Hurley 02 or Fig. A1)
        //j-y component
        vy  = (vy_kick - vorb*cos_beta); //vy kick (vk sin_omega cos_phi) is just sn->vkick[0]
        vy2 =  vy*vy; //j(or y) component  of the new velocity after kick (Eq. A12 in Hurley 02 or Fig. A1)
        //z-k component
        vz = vz_kick; //vk*sin_phi (see definiton of sin_phi), k( or z) component of the new ang mom  (Eq. A12 in Hurley 02 or Fig. A1)
        vz2 = vz*vz; //vk*vk*sin_phi*sin_phi (see definiton of sin_phi), k( or z) component of the new ang mom  (Eq. A12 in Hurley 02 or Fig. A1)

    }
    // the magnitude of the new specific angular momentum vector
    hn2 = r*r*(vz2 + vx2); //Magnitude of the new specific ang mom |h|=|r||v| Notice vj is not considered  because j is aligned along r by construction so hj=0 (see Fig. A1 in Hurley+02)
    //cos of the angle between the new and the old orbital plane
    //cos nu is estiamted as cos nu = Lz/|L|= Lz/sqrt(Lx^2 + Lz^2 + Ly^2)
    //with Lx= yVz -z*Vx, Ly=-x*Vz+Vx*z, Lz=x*Vy-Vx*y
    //but considering the assumed frame of reference  x=0, z=0, y=r (Fig. A1 in Hurley+02):
    //so Lx = rVz, Ly=0, Lz=-rVx
    // so cos nu = - rVx/(r*sqrt(Vz^2+Vx^2)) = - Vx/(Vz^2+Vx^2)
    // and sin nu = Lx/|L| = Vz/(Vz^2+Vx^2)
    //TODO Add sin nu as output in the log (nu alone is not enough to obtain the sign)
    cos_nu = -vx/std::sqrt((vz2+vx2)); //Eq. A13 in Hurley02
    //New centre of mass velocity
    double vkick_ratio = sn->getp(Mass::ID)/Mbin_new;
    //TODO Here sn->get_supernova()->get_Mejected() is equal to getp_0(Mass) - getp_0(Mass), but in principle other processes like wind can change this difference,
    //We want really just the Mejected or take into account anything
    double vold_ratio = (sn->get_supernova()->get_Mejected()*other->getp_0(Mass::ID))/(Mbin_new*Mbin_old);
    //COM velocity Eq. A14 in Hurley+02
    vcomx = vkick_ratio*vx_kick - vold_ratio*vorb*sin_beta;
    vcomy = vkick_ratio*vy_kick - vold_ratio*vorb*cos_beta;
    vcomz = vkick_ratio*vz_kick;
    //TODO This equation make clear that the vcom change also without a kick, just beacuse of mass loss, therefore vcom should be a binary property that evolve after each
    //mass loss process (also winds).
    vcom  = std::sqrt(vcomx*vcomx + vcomy*vcomy + vcomz*vcomz);
    set_cos_nu(cos_nu);
    set_vcom(vcom/utilities::kms_to_RSunyr); //Transform back in km/s
    ///***********************************************///

    ///********* Estimate new orbital parameters *********/
    double G_Mn = G*Mbin_new; //G*Mb_n where Mb_n is the new total mass of the system (after the SN explosion).


    //Estimate new a,
    //remember v^2 = G M (2/r -1/a), so the new  a is
    //a = G M r / (2*G*M - v*v*r), where M is the new total mass of the binary (after sn explosion) and v is the new velocity after kick
    double denom = ( 2*G_Mn - vfin*vfin*r );
    a_fin = denom > 0 ? G_Mn *r / denom : - 99999.0;  // When denom<=0, we have vfin^2 >= 2*G*M/r that is the escape velocity, so the system is broken, set a negative a.
    //a_fin = 1.0 / (2.0/r - vfin*vfin/(G*(sn->getp(Mass::ID) + other->getp(Mass::ID))));
    //utilities::wait("afin inside kick", a_fin, "v_fin",vfin,"vk",vk,"vorb",vorb,__FILE__,__LINE__);
    //double Blauw_kick=std::sqrt(G/r*(sn->getp_0(Mass::ID) + other->getp_0(Mass::ID))) - std::sqrt(1/r*G_Mn);
    //std::cout<<binstar->get_ID()<<","<<sn->get_zams()<<","<<a<<","<<a_fin<<","<<vk/utilities::kms_to_RSunyr<<","<<Blauw_kick/utilities::kms_to_RSunyr<<","<<sn->getp_0(Mass::ID)-sn->getp(Mass::ID)<<","<<sn->getp_0(MCO::ID)<<","<<sn->getp(Mass::ID)<<std::endl;

    //Estimate new e,
    //remember  |h| = sqrt ( G*M*a*(1-e^2)), so
    // e = sqrt( 1-|h|^2 / G*M*a  ) where M is the new total mass of the binary (after sn explosion)
    // a is the new semimajor axis after the explosion  and h is the new specific ang mom   after the kick
    double e2 = 1 - hn2/ ( G_Mn*a_fin );
    ecc_fin =  e2 > 0 ? sqrt(e2) : 0.0; //If e2<=0 the system is not anymore an ellipse, the sytem is broken.
    //ecc_fin = sqrt(max(1.0 - hn2/(G*a_fin*(sn->getp(Mass::ID) + other->getp(Mass::ID))),0.0));
    //utilities::wait("ecc fin inside kick", ecc_fin, ecc,1.0 - hn2/(G*a*(sn->getp(Mass::ID) + other->getp(Mass::ID))),__FILE__,__LINE__);

    ///***********************************************///

    //********* Set SN as KICKED *********/
    sn->set_kicked();
    //***********************************/

    return EXIT_SUCCESS;
}
/*******************************************/


/*******************************************
********   CE    *************
*******************************************/

Orbital_change_CE*   Orbital_change_CE::Instance(const std::string &name){
    auto it = GetStaticMap().find(name);

    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;

    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}


/*******************************************/
double Hurley_rl_propeller::correct_accreted_mass(double DM_accreted,  _UNUSED double DM_donor, _UNUSED Star *donor, Star *accretor,
                                                  Binstar *b) {
    //Disable accretion on NS for propeller mechanism (Campana+18, https://www.aanda.org/articles/aa/pdf/2018/02/aa30769-17.pdf)
    if (accretor->amiNS()) {

        //Magnetic Radius
        double Rm = 0.5*utilities::R_Alfven(accretor,DM_accreted/b->getp(BTimestep::ID),true); //Magnetic radius in Rsun
        double Omega = accretor->getp_0(OmegaRem::ID)*utilities::yr_cgs;//OmegaRem in 1/yr
        double Rcorotation = std::pow(utilities::G*accretor->getp_0(Mass::ID)/(Omega*Omega),1./3.);
        //Do not accrete nothing given that we are in the propeller regime
        if (Rm>Rcorotation){
            DM_accreted=0;
        }
    }

    return DM_accreted;
}
