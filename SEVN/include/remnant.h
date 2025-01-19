//
// Created by iorio on 7/5/21.
//

//Notice, the Mass property in staremnant returns just the mass of the star, there are no variables
//to internally track the current mass of the star inside the class. The member Mremnant_at_born just track
//the mass of the systems when it is created.


#ifndef SEVN_REMNANT_H
#define SEVN_REMNANT_H

#include <iostream>
#include <sevnlog.h>
#include <utilities.h>
#include <lookup_and_phases.h>
#include <random>

class Star;



class Staremnant {
public:

    Staremnant(_UNUSED Star *s, double Mremnant, double time) : born_time(time), Mremnant_at_born(Mremnant){}
    Staremnant(Star *s, double Mremnant);
    SevnLogging svlog;

    virtual ~Staremnant()=default;

    inline double get_born_time() const {return born_time;}
    inline double get_Mremnant_at_born() const {return Mremnant_at_born;}

    virtual double Mass(_UNUSED Star *s) const;
    virtual double Radius(_UNUSED Star *s) const = 0; //Pure virtual
    virtual double Luminosity(_UNUSED Star *s) const = 0; //Pure virtual
    virtual double OmegaRem(_UNUSED Star *s) const = 0; //Pure virtual
    virtual double Inertia(_UNUSED Star *s) const = 0; //Pure virtual
    virtual double Bmag(_UNUSED Star *s) const  = 0; //Pure virtual
    virtual double Xspin(_UNUSED Star *s) const = 0; //Pure virtual

    /**
     * Get the property with given ID. This a wrapper and a dispatcher to get the results of the various class methods.
     * @param s Pointer to the star
     * @param ID  ID of the Property
     * @return The value of the property as estimated for the remnant or throw a not_implemented_error if the property is not available in the class
     */
    double get(Star* s, size_t ID) const;

    inline Lookup::Remnants get_remnant_type() const {return remnant_type;}

    /**
     * Estimate the elapsed time from the remannt creation. It is estimated simply as Worldtime-creation time
     * @param s Pointe to the star
     * @return  Age of the remnant in Myr.
     */
    double age(Star *s) const;

    double InertiaSphere(Star *s) const;

protected:
    Lookup::Remnants remnant_type;

private:

    double born_time; /*!< Time of remnant creation */
    double Mremnant_at_born;  /*!< Mass of the remnant at the moment of the creation */

};




class BHrem :  public Staremnant{

public:

    BHrem(_UNUSED Star *s, double Mremnant, double time) : Staremnant(s, Mremnant, time) {
        default_initialiser(s);
    }
    BHrem(Star *s, double Mremnant) : Staremnant(s,Mremnant){
        default_initialiser(s);
    };

    /** Schwarzschild radius 2GM/c^2 **/
    double Radius(_UNUSED Star *s) const override;
    /** BH Luminosity from Eq. 96 in Hurley+00**/
    double Luminosity(_UNUSED Star *s) const override { return 1e-10;}
    double OmegaRem(_UNUSED Star *s) const override {return 0.;}
    double Inertia(_UNUSED Star *s) const override {return InertiaSphere(s);}
    double Bmag(_UNUSED Star *s) const override {return 0.;}
    double Xspin(_UNUSED Star *s) const override;

    /** Sets **/
    int apply_Bavera_correction_to_Xspin(double period, double mass_wr);

protected:

    void default_initialiser(Star *s);

    /**
     * Set the value of Xspin
     * @param value number between 0 and 1
     * @return EXIT SUCCESS, thrown an sanity_error if xspin <0 or >1
     */
    inline int set_Xspin(double value){
        if (value<0 or value>1){
            svlog.critical("Xspin in Bhrem has to be between 0 and 1, current value is "
                           + utilities::n2s(value,__FILE__,__LINE__),__FILE__,__LINE__,sevnstd::sanity_error());
        }
        xspin=value;
        return EXIT_SUCCESS;
    }




    /**
     * Estimate Xspin given the input option
     * @param s Pointer to the star
     * @return Value of Xspin
     */
    double estimate_Xspin(_UNUSED Star *s) const;

    /**
     * Xspin following Equation 3 in https://www.aanda.org/articles/aa/pdf/2020/04/aa36528-19.pdf
     * @param z0 Metallicity
     * @param mco CO mass in Msun
     * @return adimensional BH spin Xspin
     */
    inline double XspinGeneva(const double z0, const double mco) const {

        double a, b, m1, m2, alow;

        a = -0.088;
        if (z0 > 0.01)
            b = 2.258, m1 = 16.0, alow = 0.13; //m2 = 24.2
        else if (z0 > 0.004) //   0.004 < z0 <= 0.01
            b = 3.578, m1 = 31.0, alow = 0.25; //m2 = 37.8
        else if (z0 > 0.0012) //  0.0012 < z0 <= 0.004
            b = 2.434, m1 = 18.0, alow = 0.0; //m2 = 27.7
        else // z0 <= 0.0012
            b = 3.666, m1 = 32.0, alow = 0.25; //m2 = 38.8

        m2 = (alow - b) / a;
        if (mco <= m1)
            return 0.85;
        else if (mco < m2)
            return a * mco + b;
        else
            return alow;
    }

    /**
     * Xspin following Equation 4 in https://www.aanda.org/articles/aa/pdf/2020/04/aa36528-19.pdf
     *  Fig. 2 of Belczynski1+20 (https://www.aanda.org/articles/aa/pdf/2020/04/aa36528-19.pdf) shows
     *   that the fitting equations are valid only up to approx a MCO mass of 70 Msun, but Xspin seems
     *  to drop to 0 for large MCO, therefore we just extrapoalte putting a lower limit at 0 Msun
     * @param z0 Metallicity
     * @param mco CO mass in Msun
     * @return adimensional BH spin Xspin
     */
    inline double XspinMESA(const double z0, const double mco) const {
        //TODO Check with Gaston the metallicity bounadry used in this function

        double a1{0.}, b1{0.};

        if (z0 > 0.0012 && z0 <= 0.004){
            if (mco <= 12.09)
                a1 = 0.0076, b1 = 0.05;
            else
                a1 = -0.0019, b1 = 0.165;
        }
        else if (z0 > 0.01)
            a1 = -0.0016, b1 = 0.115;
        else if (z0 > 0.004) //  0.004 < z0 <= 0.01
            a1 = -0.0006, b1 = 0.105;
        else // z0 <= 0.0012
            a1 = -0.0010, b1 = 0.125;

        //Fig. 2 of Belczynski1+20 (https://www.aanda.org/articles/aa/pdf/2020/04/aa36528-19.pdf) shows
        //that the fitting equations are valid only up to approx a MCO mass of 70 Msun, but Xspin seems
        //to drop to 0 for large MCO, therefore we just extraploate putting a lower limit at 0 Msun
        return std::max(a1 * mco + b1,0.0);
    }

    /**
     * Xspin following Equation 5 in https://www.aanda.org/articles/aa/pdf/2020/04/aa36528-19.pdf
     * @param z0 Metallicity
     * @param mco CO mass in Msun
     * @return adimensional BH spin Xspin
     */
    inline double XspinFuller() const {
        return 0.01;
    }

    /**
     * Xspin drawn from a Maxwellian distribution
     * @param sigma_xspin Maxwellian 1D sigma
     * @return adimensional BH spin Xspin
     */
    inline double XspinMaxwellian(const double sigma_xspin) const {

        double x1, x;

        std::normal_distribution<double> gaussian_xspin{0.0, sigma_xspin};

        x1 = gaussian_xspin(utilities::mtrand);
        x1 *= x1;
        x = gaussian_xspin(utilities::mtrand);
        x *= x;
        x += x1;
        x1 = gaussian_xspin(utilities::mtrand);
        x1 *= x1;
        x += x1;

        x = sqrt(x);
        if (x > 1.0) /* Just in case the distribution throws x > 1 (possible although not too probable?)*/
            x = 1.0;

        return x;
    }

    /**
     * Just return 0 for the spin
     * @return 0.
     */
    inline double XspinZeros() const {return 0.0;}

    //TODO This option is a test, we have to decide if keep it or no
    //TODO In case we keep it, myabe we have to think to a better implementation (using processes since here the anchges are due to mass acretion)
    /**
     * Xspin following  Zevin&Bavera22 https://arxiv.org/pdf/2203.02515.pdf
     * All the BH are born with an Xspin=0, then it is increased just through mass accretion
     * Eq. 4 in Zevin&Bavera22.
     * Notice:: in Zevin&Bavera22, they consider only MT in RLO, in our implementation we consider
     * also accretion from winds.
     * @param s Pointer to the star
     * @return Xspin after accretion (or 0 if not accretion)
     */
    double XspinAccretion(Star *s) const;

private:
    double xspin=std::nan(""); // Value of the BH spin.
};



class NSrem :  public Staremnant{

public:

    NSrem(_UNUSED Star *s, double Mremnant, double time) : Staremnant(s, Mremnant, time), root_distribution(0.,1.0) {
        default_initialiser(s);


    }
    NSrem(Star *s, double Mremnant) : Staremnant(s,Mremnant), root_distribution(0.,1.0) {
        default_initialiser(s);
    };


    double Radius(_UNUSED Star *s) const override {return Rns;}
    /**
     * NS Luminosity, Eq. 93 Hurley 2000,
     * LNS=0.02*M^(2/3)/(max(t,0.1)^2) from Eq. 93 using Hurley, 2000
     * @param s Pointer to the star
     * @return Neutrpm star luminosity
     */
    double Luminosity(_UNUSED Star *s) const override;
    double OmegaRem(_UNUSED Star *s) const  override;
    double Inertia(_UNUSED Star *s) const override;
    double Bmag(_UNUSED Star *s) const override;
    double Xspin(_UNUSED Star *s) const override {return std::nan("");}

    /**
     * Returh the alpha angle, i.e. the angle between the rotation axis and the magnetic axis
     * @return alpha in radians
     */
    inline double get_alpha() const{return std::asin(sinalpha);}

    /**
     * Returh the sin of the alpha angle, i.e. the angle between the rotation axis and the magnetic axis
     * @return sin alpha
     */
    inline double get_salpha() const{return sinalpha;}

    /**
     * Return the value of Bmin
     * @return Bmin in Gauss
     */
    inline double get_Bmin() const{return Bmin;}

protected:

    std::uniform_real_distribution<double> root_distribution;

    inline double generate_uniform(double a, double b){
        return root_distribution(utilities::mtrand)*(b-a) + a;
    }

    void default_initialiser(_UNUSED Star *s);

    void print_log_message(_UNUSED Star *s);

public:
    constexpr static double Rns = 11*utilities::km_to_RSun; //Neutron star radius 11 km in Rsun, from Hurley, 2000

private:
    double sinalpha; //sin of the angle between the rotation axis and the magnetic ax
    double B0, Bmin; //Initial and minimum magnetic field in Gauss
    double Omega0; //Initial Spin
    double tau_magnetic; //Magnetic field decay time scale in Myr

};

class NSCCrem : public NSrem {

public:
    NSCCrem(_UNUSED Star *s, double Mremnant, double time) : NSrem(s, Mremnant, time) {
        print_log_message(s);
    }
    NSCCrem(Star *s, double Mremnant) : NSrem(s,Mremnant){
        print_log_message(s);
    };
};

class NSECrem : public NSrem {

public:
    NSECrem(_UNUSED Star *s, double Mremnant, double time) : NSrem(s, Mremnant, time) {
        remnant_type=Lookup::Remnants::NS_ECSN;
        print_log_message(s);
    }
    NSECrem(Star *s, double Mremnant) : NSrem(s,Mremnant){
        remnant_type=Lookup::Remnants::NS_ECSN;
        print_log_message(s);
    };
};

class WDrem :  public Staremnant{

public:

    WDrem(_UNUSED Star *s, double Mremnant, double time) : Staremnant(s,Mremnant, time) {
        default_initialiser();
    }
    WDrem(Star *s, double Mremnant) : Staremnant(s, Mremnant) {
        default_initialiser();
    }

    /**
     * Radius from Eq. 91 in Hurley+00
     * @param s Pointer to star
     * @return Radius of the WD
     */
    double Radius(_UNUSED Star *s) const override;
    /**
     * Estimate the luminosity evolution of the WD followint Eq. 90 in Hurley+00
     * @param s Pointer to the star
     * @return Luminosity of the WD at a given age in Lsun.
     */
    double Luminosity(_UNUSED Star *s) const override;
    double OmegaRem(_UNUSED Star *s) const override {return 0.;}
    double Inertia(_UNUSED Star *s) const override {return InertiaSphere(s);};
    double Bmag(_UNUSED Star *s) const override {return 0.;}
    double Xspin(_UNUSED Star *s) const override {return 0.;}

protected:
    double A_luminosity; //To be used in the estimate of LWD (See Hurley+00 Sec. 6.2)

    inline  void default_initialiser(){}

};

class HeWDrem :  public WDrem {

public:

    HeWDrem(_UNUSED Star *s, double Mremnant, double time) : WDrem(s, Mremnant, time){
        default_initialiser();
    }

    HeWDrem(_UNUSED Star *s, double Mremnant) : WDrem(s, Mremnant){
        default_initialiser();
    }

protected:

    inline  void default_initialiser()  {
        remnant_type=Lookup::Remnants::HeWD;
        A_luminosity = 4; //From Hurley+00 Sec 6.2
    }
};

class COWDrem :  public WDrem {

public:

    COWDrem(_UNUSED Star *s, double Mremnant, double time) : WDrem(s, Mremnant, time){
        default_initialiser();
    }

    COWDrem(_UNUSED Star *s, double Mremnant) : WDrem(s, Mremnant){
        default_initialiser();
    }

protected:
    inline  void default_initialiser(){
        remnant_type=Lookup::Remnants::COWD;
        A_luminosity = 15; //From Hurley+00 Sec 6.2
    }
};

class ONeWDrem :  public WDrem {

public:

    ONeWDrem(_UNUSED Star *s, double Mremnant, double time) : WDrem(s, Mremnant, time){
        default_initialiser();
    }

    ONeWDrem(_UNUSED Star *s, double Mremnant) : WDrem(s, Mremnant){
        default_initialiser();
    }
protected:
    inline  void default_initialiser()  {
        remnant_type=Lookup::Remnants::ONeWD;
        A_luminosity = 17; //From Hurley+00 Sec 6.2
    }
};


/***
 * A remnant class to handle objects that cannot be currently taken into account.
 * For example stars that should jump in a track that is not available
 */
class Zombierem : public Staremnant{

public:

    Zombierem(Star *s, double Mremnant, double time) : Staremnant(s, Mremnant, time) {
        initialise(s);
    }

    Zombierem(Star *s, double Mremnant) : Staremnant(s,Mremnant){
        initialise(s);
    };

    double Radius(_UNUSED Star *s) const override { return radius;}
    double Luminosity(_UNUSED Star *s) const override { return luminosity;}
    double OmegaRem(_UNUSED Star *s) const override {return 0.;}
    double Inertia(_UNUSED Star *s) const override {return InertiaSphere(s);}
    double Bmag(_UNUSED Star *s) const override {return 0.;}
    double Xspin(_UNUSED Star *s) const override {return 0.;}

protected:

    void initialise(Star *s);

private:
    double luminosity;
    double radius;
};



#endif //SEVN_REMNANT_H
