//
// Created by iorio on 12/14/21.
//

#ifndef SEVN_LAMBDA_NANJING_H
#define SEVN_LAMBDA_NANJING_H

#include <vector>
#include <sevnlog.h>
#include <utilities.h>
#include <lambda_base.h>

class Star;
class IO;

/**
 * Basic Lambda model from the work of Xu&Li10 (https://iopscience.iop.org/article/10.1088/0004-637X/716/1/114) +
 * its errata corrige (https://iopscience.iop.org/article/10.1088/0004-637X/722/2/1985).
 * The implementation is taken directly from the COMPAS code https://github.com/TeamCOMPAS/COMPAS that is based
 * on the STARTRACK implementation (see Dominik+12, https://ui.adsabs.harvard.edu/abs/2012ApJ...759...52D/abstract).
 * Lambda is thus estimated as:
 *
 *      y = a + b1*x + b2*x^2 + b3+*x^3 + b4*x^4 + b5*x^5
 *
 *      where y is in general lambda (gravitational or considering also the thermal energy), sometime it is 1/lambda
 *      or sometime it is log lambda.
 *      x in general is the Radius, but sometime it is the mass of the envelope.
 *      The a and b coefficient are estimate for a given list of Mzams (1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,20.,50.,100.)
 *      and they are split in two broad metallicity bin (PopI and PopII), we define the PopI stars as the one with
 *      Z<0.085 that is half of the Zsun=0.017. They are also split in the three different evolutionary stages:
 *
 *      - Stage 1 begins at the center H exhaustion and ends when the star starts to shrink (i.e., near center He ignition).
 *      - Stage 2 follows and ends when the star starts to expand again.
 *      - Stage 3  Stage 3 begins after that and continues until the end of the evolution.
 *
 *      Following what has been done in COMPAS, we  identify the three stages with the SEVN phases:
 *
 *      - Stage 1: TerminalMainSequence, Hshellburning
 *      - Stage 2: coreHeburning
 *      - Stage 3: TerminalcoreHeburning, Heshellburning
 *
 *      NOTICE: We also use the lambda estimate for pureHE stars written in COMPAS in the Nanjing implementation.
 *      It has not published yet and  it is
 *          lambda = 0.3*R^(-0.8) with Radius Range between 0.25-120 Rsun, radii smaller or larger of that limit
 *          are set to the mimimum or maximum allowed value.
 *
 * The lambda can be retrieved using the overload operator (Star *s)
 * Using this class the lambda is estimated in a quantised way, i.e. for a Mzamsnot present in the table, the
 * closed track with Mtrack>Mzams is used. If Mzams>100, Mtracks=100.
 */
class Lambda_Nanjing : public Lambda_Base {

public:
    Lambda_Nanjing(){};
    Lambda_Nanjing(_UNUSED const Star *s){};
    Lambda_Nanjing(_UNUSED const IO *io){};
    virtual ~Lambda_Nanjing() = default;

    /**
     * Estimate lambda for a given star
     * @param s Pointer to the star
     * @return lambda
     */
    double operator() (const Star* s) override;

protected:


    /** Methods to estimate lambda in various phases **/
    virtual double lambda_Giant(const Star* s);
    virtual double lambda_Cheb(const Star* s);
    virtual double lambda_AGB(const Star* s);
    virtual double lambda_pureHe(const Star* s);

    //Properly estimated the two lambdas, they call set_coeff** to set the coefficients of the fittin equations
    virtual double lambda_Giant(_UNUSED double Mzams, _UNUSED  double Z, _UNUSED  double Radius, _UNUSED double Menv, _UNUSED double lambda_th);
    virtual double lambda_Cheb(_UNUSED double Mzams, _UNUSED  double Z, _UNUSED  double Radius, _UNUSED double Menv, _UNUSED double lambda_th);
    virtual double lambda_AGB(_UNUSED double Mzams, _UNUSED  double Z, _UNUSED  double Radius, _UNUSED double Menv, _UNUSED double lambda_th);
    virtual double lambda_pureHe(_UNUSED double Mzams, _UNUSED  double Z, _UNUSED  double Radius, _UNUSED double Menv, _UNUSED double lambda_th);

    //Set the coefficient of the lambda fitting equations
    typedef unsigned int coeff_status;
    virtual coeff_status set_coeff_Giant(double Mzams, double Z, double Radius);  //From HG::CalculateLambdaNanjing in HG.cpp in COMPAS
    virtual coeff_status set_coeff_Cheb(double Mzams, double Z, double Radius);   //From CHeB::CalculateLambdaNanjing in CHeb.cpp in COMPAS
    virtual coeff_status set_coeff_AGB(double Mzams, double Z, double Radius);    //From EAGB::CalculateLambdaNanjing in EAGB.cpp in COMPAS

    /**
     * Polynomial fitting equation for lambda (Eq. 13 in Xu&Li10)
     * @param x argument of the fitting equation
     * @return results (2) of the fifth order polynomial equation, the coefficients are taken from coeff_lambdab and coeff_lambdag
     */
    std::vector<double> fitting_equation(double x){
        double x2 = x * x;
        double x3 = x2 * x;
        double x4 = x2 * x2;
        double x5 = x3 * x2;

        double y1 = coeff_lambdab[0] + (coeff_lambdab[1] * x) + (coeff_lambdab[2] * x2) + (coeff_lambdab[3] * x3) + (coeff_lambdab[4] * x4) + (coeff_lambdab[5] * x5);
        double y2 = coeff_lambdag[0] + (coeff_lambdag[1] * x) + (coeff_lambdag[2] * x2) + (coeff_lambdag[3] * x3) + (coeff_lambdag[4] * x4) + (coeff_lambdag[5] * x5);

        std::vector<double> y_vec = {y1,y2};

        return y_vec;
    }

    sevnstd::SevnLogging svlog;
    std::vector<double> coeff_lambdab; /*!< coefficient to estiamate Lambda from Xi&Liu considering both the gravitational + internal energy. Set in set_coeff** */
    std::vector<double> coeff_lambdag; /*!< coefficient to estiamate Lambda from Xi&Liu considering  the gravitational  energy only. Set in set_coeff** */
    std::vector<double> maxBG = {};  /*!< array storing the maximum value allows for lambdab and lambdag. Set in set_coeff** */
    std::vector<double> lambdaBG = {}; /*!< array storing the values estimated for lambdab and lambdag. Set in set_coeff** */

    /**
     * Reset (clear) members (coeff_lambdab, coeff_lambdag, maxBG, lambdaBG)
     * @return EXIT_SUCCESS
     */
    inline int reset(){
         coeff_lambdab.clear();
         coeff_lambdag.clear();
         maxBG.clear();
         lambdaBG.clear();
         return EXIT_SUCCESS;
     }

private:

    const double ZpopI = 0.0085; /** Metallicity boundary between popI and popII stars **/
    constexpr static unsigned int COEF_SET = 1;
    constexpr static unsigned int COEF_NOT_SET = 0;

    /**
     * Wrapper for the simple equation lambdab*lambda_th + lambdag*(1-lambda_th)
     * @param lambda_th Fraction of the internal energy to use (0-only gravitational energy, 1-gravitational+all internal energy)
     * @return final lambda
     */
    inline double lambda_estimate(double lambda_th){
        return lambda_th*lambdaBG[0] + (1-lambda_th)*lambdaBG[1];
    }
};


/**
 * Class to estimate the Lambda using fitting equations in Xi&Liu10.
 * See the documentation of the base class for more information.
 *
 * With respect to the base class, the lambda is not estimated in a quantised form, bur rather using
 * an interpolation considering the Mzams. The interpolating Mzams are the one reported in the
 * Xu&Li errata (https://iopscience.iop.org/article/10.1088/0004-637X/722/2/1985): (1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,20.,50.,100.) Msun.
 *
 * Each time a lambda estimate is needed for a star with Mzams=Mzams*, two interpolating masses are defined so that
 * Mzams_track0<=Mzams*<=Mzams_track1, if Mzams*<1 Msun or Mzams*>100 Msun, the lambda is estimated using Mzams=1 Msun and Mzams1=100 Msun respectively
 * The lambda of the interpolating tracks are estimated considering the current Radius, envelope mass and metallicity of the current track
 */
class Lambda_Nanjing_interpolator : public Lambda_Nanjing{

public:
    Lambda_Nanjing_interpolator(){};
    Lambda_Nanjing_interpolator(_UNUSED const Star *s){};
    Lambda_Nanjing_interpolator(_UNUSED const IO *io){};

    double operator() (const Star* s) override;

protected:

    double operator() (unsigned int Phase, double Mzams, double Z, double Radius, double Menv, double lambda_th);


    static std::vector<double> Mzams_list;  /*!< list of all zams masses in the Xu&Li10 and Dominik+12  tables */
    double Mzams_cache{0}, Z_cache{0}; /*!< Last Mzams and Z used to estimate lambda */

    void find_interpolators(double Mzams);

private:

    double Mzams_interpolators [2] = {0.,0.}; /*!< Mzams of the two interpolator */
    double wM [2] = {0., 0.}; /*!< interpolator weights for the two Mzams */

};

#endif //SEVN_LAMBDA_NANJING_H
