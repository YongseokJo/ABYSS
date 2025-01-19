//
// Created by mario on 13/02/19.
//
#include <utilities.h>
#include <IO.h>
#include <binstar.h>
#include <star.h>

thread_local std::mt19937_64 utilities::mtrand;

unsigned long utilities::gen_rseed() {
    thread_local std::random_device rd; //thread_local means that the lifetime of the variable is mantained until the end of the thread, (it is not re-defined each time the function is called)
    thread_local std::uniform_int_distribution<unsigned long> distribution(0,1000000000000000);

    return distribution(rd);
}


std::vector<std::string> utilities::split(const std::string& s, char delimiter){
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)){
        tokens.push_back(token);
    }
    return tokens;
}


double utilities::maxwellian_cdf(double x, double sigma){
    double a = x/sigma;
    return (erf(0.7071067811865475*a) - 0.7978845608028654*a*exp(-0.5*a*a));
}

double utilities::maxwellian_pdf(double x, double sigma){
    double a = x/sigma;
    return (0.7978845608028654*a*a/sigma*exp(-0.5*a*a));
}


double utilities::roche_lobe_Eg(double Mass_primary, double Mass_secondary, double a) {

    double q=Mass_primary/Mass_secondary;

    double qpow,A,B;

    qpow=pow(q,2./3.);
    A=0.49*qpow;
    B=0.6*qpow + log(1+pow(q,1./3.));

    return a * (A/B);
}


double utilities::R_Alfven(Star *s, double dMdt, bool get0)  {

    //Lambda func to use get or get_0 values
    const double& R = get0 ? s->getp_0(Radius::ID)*Rsun_cgs : s->getp(Radius::ID)*Rsun_cgs;
    const double& M = get0 ? s->getp_0(Mass::ID)*Msun_cgs : s->getp(Mass::ID)*Msun_cgs;
    const double& B = get0 ? s->getp_0(Bmag::ID) : s->getp(Bmag::ID);
    dMdt *= Msun_cgs/(1e6*yr_cgs);

    double factor1 = std::pow(1.0/(G_cgs*8), 1.0/7.0);
    double factor2 = std::pow(pow(R, 6.0)/(dMdt*pow(M, 0.5)), 2.0/7.0);
    double factor3 = std::pow(B, 4.0/7.0);


    double RA=factor1*factor2*factor3; //Alfven radius in Rcm

    return RA/Rsun_cgs; //Alfven radius in  cm
}




std::string utilities::get_name(Star* s){
    return s->get_name();
}
long utilities::get_ID(Star* s){
    return s->get_ID();
}
std::string utilities::get_name(Binstar* b){
    return b->get_name();
}
long utilities::get_ID(Binstar* b){
    return b->get_ID();
}
double utilities::get_current_time(Star* s){
    return s->getp(Worldtime::ID);
}
double utilities::get_current_time(Binstar* b){
    return b->getp_0(BWorldtime::ID);
}


std::string utilities::log_star_info(Star *s,bool oldtstep) {


    std::string fstring;

    fstring+=utilities::n2s(s->get_ID(),__FILE__,__LINE__)+":";
    if (oldtstep){
        fstring+=utilities::n2s(s->getp_0(Mass::ID),__FILE__,__LINE__,3)+":";
        fstring+=utilities::n2s(s->getp_0(MHE::ID),__FILE__,__LINE__,3)+":";
        fstring+=utilities::n2s(s->getp_0(MCO::ID),__FILE__,__LINE__,3)+":";
        fstring+=utilities::n2s((int)s->getp_0(Phase::ID),__FILE__,__LINE__)+":";
        fstring+=utilities::n2s((int)s->getp_0(RemnantType::ID),__FILE__,__LINE__);
    }
    else{
        fstring+=utilities::n2s(s->getp(Mass::ID),__FILE__,__LINE__,3)+":";
        fstring+=utilities::n2s(s->getp(MHE::ID),__FILE__,__LINE__,3)+":";
        fstring+=utilities::n2s(s->getp(MCO::ID),__FILE__,__LINE__,3)+":";
        fstring+=utilities::n2s((int)s->getp(Phase::ID),__FILE__,__LINE__)+":";
        fstring+=utilities::n2s((int)s->getp(RemnantType::ID),__FILE__,__LINE__);
    }

    return fstring;
}

double utilities::Hfrac(Star *s) {
    //TODO Now we should have this value in the table Hsup
    return 0.760 -3.0*s->get_Z();
}

double utilities::dMdt_Eddington_accretion(Star *donor, Star *accretor, double eddfact){

    //Check eddington factor
    if (eddfact<0){
        sevnstd::SevnLogging svlog;
        svlog.critical("Eddington factor cannot  negative (it is"+
        utilities::n2s(eddfact,__FILE__,__LINE__), __FILE__,__LINE__,sevnstd::sanity_error(""));
    }

    double Hf = Hfrac(donor);
    double R2 = accretor->getp_0(Radius::ID); //Remember we use always the single stellar property at time T-1 (in order to apply SSE and BSE processes in synch)
    double dMdt = 2.08e-03*(1.0/(1.0 + Hf))*R2; //Eq. 67 in Hurley+02.

    return  eddfact*dMdt;
}

bool utilities::areEqual(double x, double y) {

    //From http://realtimecollisiondetection.net/blog/?p=89
    //Absolute tolerance works well for values<1, relative tolerance for larger values
    double maxXYOne = std::max( 1.0, std::max(std::fabs(x) , std::fabs(y) )) ;

    return std::fabs(x - y) <= std::numeric_limits<double>::epsilon()*maxXYOne;
}

double utilities::smallestSignificativeStep(double x) {
    double maxXYOne = std::max( 1.0, fabs(x));
    return 1.01*std::numeric_limits<double>::epsilon()*maxXYOne;
}




std::unique_ptr<utilities::ListGenerator> utilities::ListGenerator::make_unique(double _vstep, double _vstep_max, double _vstep_min) {
    return std::unique_ptr<ListGenerator>(new ListGenerator(_vstep,_vstep_max,_vstep_min));
}

std::unique_ptr<utilities::ListGenerator> utilities::ListGenerator::make_unique(std::vector<double> _vlist) {
    return std::unique_ptr<ListGenerator>(new ListGenerator(_vlist));
}
