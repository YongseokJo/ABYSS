//
// Created by Giuliano Iorio on 10/08/2021.
//

#ifndef SEVN_BSEINTEGRATOR_H
#define SEVN_BSEINTEGRATOR_H

constexpr double cZsun = 0.02;

/****************** BSE COEFFICIENTS ********************/

/**
 * Parent class to store the coefficients from Appendix A in Hurley+00
 * it receives as input the metallicity and the alpha,beta,gamma,eta parameters and it returns overloading the () operator
 * the value of the corresponde coefficient
 */
class BSE_Coefficient{

public:
    BSE_Coefficient(double Z, double alpha, double beta=0, double gamma=0, double eta=0, double mu=0) : alpha{alpha}, beta{beta}, gamma{gamma}, eta{eta}, mu{mu}, Z{Z}{
        set_coeff();
    }

    inline double operator() () const { return coeff;}
    inline double getZ(){return Z;}
    inline void change_Z(double newZ){
        Z=newZ;
        set_coeff();
    }

    const double alpha;
    const double beta;
    const double gamma;
    const double eta;
    const double mu;

protected:
    const double Zsun=cZsun;
    double Z;
    double coeff;

    inline void set_coeff(){
        double Zlog = std::log10(Z/Zsun);
        //Appendix A Hurley+00
        coeff = alpha + Zlog*(beta + Zlog*(gamma + Zlog*(eta + Zlog*mu))); // alpha+beta*Zlog+gamma*Zlog^2+eta*Zlog^3+mu*Zlog^4
    }
};

//TODO it is betteer to create a single class containing a vector with all the coefficients
/**
 * The following classes derived from BSE Coefficient store automatically the values of the coefficients in Appendix A.
 */
class BSE_Coefficient_a1 : public BSE_Coefficient{
public:
    BSE_Coefficient_a1(double Zinput) : BSE_Coefficient{Zinput,1593.890,2053.038,1231.1226,232.7785}{

    }
};
class BSE_Coefficient_a2 : public BSE_Coefficient{
public:
    BSE_Coefficient_a2(double Zinput) : BSE_Coefficient{Zinput,2.706708e3,1.483131e3,5.772723e3,7.411230e1}{}
};
class BSE_Coefficient_a3 : public BSE_Coefficient{
public:
    BSE_Coefficient_a3(double Zinput) : BSE_Coefficient{Zinput,1.466143e2,-1.048442e2,-6.795374e1,-1.391127e1}{}
};
class BSE_Coefficient_a4 : public BSE_Coefficient{
public:
    BSE_Coefficient_a4(double Zinput) : BSE_Coefficient{Zinput,4.141960e-2,4.564888e-2,2.958542e-2,5.571483e-3}{}
};
class BSE_Coefficient_a5 : public BSE_Coefficient{
public:
    BSE_Coefficient_a5(double Zinput) : BSE_Coefficient{Zinput,3.426349e-1}{}
};
class BSE_Coefficient_a6 : public BSE_Coefficient{
public:
    BSE_Coefficient_a6(double Zinput) : BSE_Coefficient{Zinput,1.949814e1,1.758178,-6.008212,-4.470533}{}
};
class BSE_Coefficient_a7 : public BSE_Coefficient{
public:
    BSE_Coefficient_a7(double Zinput) : BSE_Coefficient{Zinput,4.903830}{}
};
class BSE_Coefficient_a8 : public BSE_Coefficient{
public:
    BSE_Coefficient_a8(double Zinput) : BSE_Coefficient{Zinput,5.212154e-2,3.166411e-2,-2.750074e-3,-2.271549e-3}{}
};
class BSE_Coefficient_a9 : public BSE_Coefficient{
public:
    BSE_Coefficient_a9(double Zinput) : BSE_Coefficient{Zinput,1.312179,-3.294936e-1,9.231860e-2,2.610989e-2}{}
};
class BSE_Coefficient_a10 : public BSE_Coefficient{
public:
    BSE_Coefficient_a10(double Zinput) : BSE_Coefficient{Zinput,8.073972e-1}{}
};

/********************************************************/


/******************BSE PROPERTIES**************************/

/**
 * Basic Abstract class to initialise a BSE property and evolve it.
 * When created a Class property gets all the needed coefficient, then the value can
 * be called with the overloaded operator (Mass) or (Mass,time).
 */
class BSE_Property{
public:
    BSE_Property(double Zinput): Z{Zinput}{}
    virtual  double operator() (double M) = 0;
    //TODO add also the operator with M and time
    //virtual  double operator() (double M, double t) = 0;
    virtual std::string units(){return "";}

    //Just a wrapper to be called easily by  pointers
    inline double eval(double M){return operator()(M);}

protected:
    double chaced_M=0;
    double chaced_logP=0;
    double Z;

    std::vector<double> coeff;

};

///TIMES PROPERTIES

/**
 * Stellar times properties. The default unit is Myr
 */
class BSEtimes : public BSE_Property{
public:
    BSEtimes(double Zinput) : BSE_Property(Zinput){}
    std::string units() override{return "Myr";}
};

/**
 * Time to reach the base of the giant branch.
 * Eq. 4 in Hurley+00
 */
class Tbgb :  public  BSEtimes{
public:

    Tbgb(double Zinput): BSEtimes(Zinput), a1{BSE_Coefficient_a1(Zinput)()}, a2{BSE_Coefficient_a2(Zinput)()}, a3{BSE_Coefficient_a3(Zinput)()}, a4{BSE_Coefficient_a4(Zinput)()}, a5{BSE_Coefficient_a5(Zinput)()}{

    }

    double operator() (double M) override{

        double M3=M*M*M;
        double M7=M3*M*M*M*M;
        double num = a1 + a2*M3 + a3*std::pow(M,5.5) + M7;
        double den = a4*M*M + a5*M7;

        return num/den;
    }

private:

    const double a1;
    const double a2;
    const double a3;
    const double a4;
    const double a5;

};

/**
 * Time to reach the hook of the MS
 * Eq. 5 in Hurley+00
 */
class Thook : public Tbgb{
public:
    Thook(double Zinput): Tbgb(Zinput), a6{BSE_Coefficient_a6(Zinput)()}, a7{BSE_Coefficient_a7(Zinput)()},a8{BSE_Coefficient_a8(Zinput)()}, a9{BSE_Coefficient_a9(Zinput)()}, a10{BSE_Coefficient_a10(Zinput)()}{

    }

    double operator() (double M) override{
        double tbgb = Tbgb::operator()(M);
        return  mu(M)*tbgb;
    }

    //Eq. 7 in Hurley+00
    double mu(double M){
        double mut = 1.0-0.01*std::max(a6/std::pow(M,a7),a8+a9/std::pow(M,a10));
        return std::max(0.5,mut);
    }

private:

    const double a6;
    const double a7;
    const double a8;
    const double a9;
    const double a10;
};

/**
 * MS lifetime
 * Eq. 5 in Hurley+00
 */
class Tms : public Thook{

public:
    Tms(double Zinput) : Thook(Zinput){
        double Zeta = std::log10(Zinput/cZsun);
        //Eq. 6 in Hurley+00
        x=std::max(0.95, std::min(0.95-0.03*(Zeta+0.30103),0.99));
    }

    double operator() (double M) override{
        double tbgb = Tbgb::operator()(M);
        double thook = Thook::operator()(M);

        return std::max(thook,x*tbgb);
    }

protected:

    double x;


};

#endif //SEVN_BSEINTEGRATOR_H
