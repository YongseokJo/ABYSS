//
// Created by Giuliano Iorio on 07/12/2021.
// This headed contains the class Lambda Klencki and its derived classes.
// They are used to use the Lambda fit from the paper Klencki+21 (Appendix A,  https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..54K/abstract).
//

#ifndef SEVN_LAMBDA_KLENCKI21_H
#define SEVN_LAMBDA_KLENCKI21_H

#include <cmath>
#include <iostream>
#include <vector>
#include <sevnlog.h>
#include <set>
#include <lambda_base.h>

class Star;
class IO;

/**
 * Basic Lambda_Klencki class.
 * The Lambda fit comes from Appendix A in Klencki+21 (https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..54K/abstract).
 * Lambda is thus estimates as:
 *      log10(lambda) = a*x^3 + b*x^2 + c*x + d
 *      where x=log10(R)
 * The value of a,b,c,d depends on the Mass (Mzams), Z (metallicity) and Radius.
 * All these values are tabulated for a set of given tracks, the table is in auxiliary_data/lambda_fit_Klencki21.dat
 * (original link  https: //ftp.science.ru.nl/astro/jklencki/).
 *
 * During the first class initialisation in each thread, the table is loaded as a static 2D vector
 * The following class initialisations will just use the already loaded table.
 *
 * The lambda can be retrieved using the overload operator (Star *s) or (double Mzams, double Z, double R).
 *
 * Notice: the lambda is estimated by Klenci+21 using MESA tracks between 10 and 80 Msun and Z/Zsun 0.01 - 0.1,
 * where Zsun=0.017 (not all the tracks contains all the ranges in Z).
 * Using this class the lambda is estimated in a quantised way, i.e. for a Mzams and a Z not present in the table, the
 * closed track with Mtrack>Mzams and Ztrack>Z. If Mzams<10, Mtracks=10 is used, for Z/Zsun<0.01, Z/Zsun=0.01 is used.
 */
class Lambda_Klencki  : public Lambda_Base {

public:
    Lambda_Klencki(const Star *s);
    Lambda_Klencki(const IO *io);
    virtual ~Lambda_Klencki() = default;

    /**
     * Find the row   of the Klencki table that best match Mzams and Z in input.
     * The best match is the one for which Mtrack>Mzams and Ztrack>Z.
     * If Mzams<10, Mtracks=10 is used, for Z/Zsun<0.01, Z/Zsun=0.01 is used.
     * @param Mzams  zams mass (Msun) of the star for which we want to find the the related track in the Klencki table
     * @param Z Z of the star for which we want to find the the related track in the Klencki table
     * @return return the quantities  of the matched  table row as a vector with dimension 17
     */
    std::vector<double> find_row(double Mzams,double Z);

    /**
     * Estimate lambda for a given star. Internally it calles the operator(double Mzams, double Z, double R).
     * @param s Star for which we want to estimate lambda
     * @return value of lambda.  Notice: the lambda is estimated by Klenci+21 using MESA tracks between 10 and 80 Msun and Z/Zsun 0.01 - 0.1,
     * where Zsun=0.017 (not all the tracks contains all the ranges in Z).
     * The lambda is estimated in a quantised way, i.e. for a Mzams and a Z not present in the table, the
     * closest track with Mtrack>Mzams and Ztrack>Z. If Mzams<10, Mtracks=10 is used, for Z/Zsun<0.01, Z/Zsun=0.01 is used.
     */
    double operator() (const Star* s) override;
    /**
     * Estimate the value of lambda as function of Mzams, Z and R.
     * Some part of the lambda estimate are cached. In particular the coefficient from Klencki+21 depends on the Mzams and Z,
     * therefore these two values are cached and the coefficients (see find_row) are loaded only  if Mzams or Z are different
     * with respect to the last call.
     * @param Mzams Zams mass in Msun
     * @param Z metallicity
     * @param R radius in Rsun
     * @return value of lambda.  Notice: the lambda is estimated by Klenci+21 using MESA tracks between 10 and 80 Msun and Z/Zsun 0.01 - 0.1,
     * where Zsun=0.017 (not all the tracks contains all the ranges in Z).
     * The lambda is estimated in a quantised way, i.e. for a Mzams and a Z not present in the table, the
     * closed track with Mtrack>Mzams and Ztrack>Z. If Mzams<10, Mtracks=10 is used, for Z/Zsun<0.01, Z/Zsun=0.01 is used.
     */
    virtual double operator() (double Mzams, double Z, double R);

protected:

     sevnstd::SevnLogging svlog;
     double Mzams_cache{0}, Z_cache{0}; /*!< Last Mzams and Z used to estimate lambda */
     std::vector<double> vector_cache; /*!< Last coefficient used to estimate lambda */
     const double Zsun=0.017; /*!< Value of Z for the sun as assumed in Klencki+21 */

     /**
      * Estimate lambda from Klencki+21 using Eq. A.1
      * @param R radius in Rsun
      * @param coefficients coefficients of the Eq. A.1 in Klencki + R12,R23 and Rmax, it has to be a vector with dimension 17 with the following columns:
      *     - 0: Mzams
      *     - 1: Z
      *     - 2: R12
      *     - 3: R23
      *     - 4: Rmax
      *     - 5-8: a1,b1,c1,d1
      *     - 9-12: a2,b2,c2,d2
      *     - 13-16: a3,b3,c3,d3
      * @return @return value of lambda
      */
     static double estimate_lambda(double R, const std::vector<double>& coefficients);

     /**
      * Here now we define two quantities that are used to load the external table
      * only once per class per thread. In this way we can save memory not loading one table
      * for each Lambda instance, this is possibile using static variables.
      * At the same time using thread_local we force to have one separate definition for these variable
      * for each thread avoiding data racing. In this case it is equivalent to defining the variable as
      * threadprivate (see https://stackoverflow.com/questions/60932116/does-thread-local-work-for-openmp-threads),
      * Variables that appear in threadprivate directives or variables with the _Thread_local (in C) or thread_local (in C++) storage-class specifier are threadprivate
      *
      * In the class constructors the table is loaded only if the variable already_loaded is set to false, after the table
      * loading alread_loaded is set to false.
      */
     static thread_local std::vector<std::vector<double>> table; /*!< 2D table  170 rows, 17 columns */
     static thread_local bool already_loaded; /*!< flag to check if it is needed to load table */
};



/**
 * Derived class from Lambda_Klencki to estimate the lambda using a Mzams-Z interpolation from the tracks tabulated
 * by Klencki+21.
 */
class Lambda_Klencki_interpolator : public  Lambda_Klencki{

public:
    Lambda_Klencki_interpolator(const Star *s);
    Lambda_Klencki_interpolator(const IO *io);
    /**
     * Estimate the value of lambda as function of Mzams,Z and R.
     * The lambda is estimated interpolating the lambdas obtained for the four interpolating tracks  in the Klencki tables.
     * The interpolating tracks  are so that  Mzams_track1<= Mzams <Mzams_track2 and  Z_track12 <= Z < Z_track34.
     * So the four interpolating tracks have:
     *   1 - Mzams_track1, Ztrack1 -> lambda1
     *   2-  Mzams_track1, Ztrack2 -> lambda2
     *   3-  Mzams_track2, Ztrack3 -> lambda3
     *   4-  Mzams_track2, Ztrack4 -> lambda4
     *  The final lambda is estimated as
     *    wM1*(wZ1*lambda1 + wZ2*lambda2) + wM2*(wZ3*lambda3 + wZ4*lambda4),
     *    whew wM1,wM2,wZ1,wZ2,wZ3,wZ4 are the interpolateing weights estimates as
     *    wM1 = (Mzams_track2 - Mzams) /  (Mzams_track2 - Mzams_track1)
     *    wM2 = (Mzams - Mzams_track1) /  (Mzams_track2 - Mzams_track1)
     *    wZ1 = (Ztrack2 - Z) / (Ztrack2 - Ztrack1)
     *    wZ2 = (Z - Ztrack1) / (Ztrack2 - Ztrack1)
     *    .... etc. etc.
     *
     * If Mzams and Z are smaller or larger than the minimum and maximum values in the table,
     * they are internally set to the minimum/maximum values in the tabe.
     *
     * Some part of the lambda estimate are cached. In particular the coefficient from Klencki+21 depends on the Mzams and Z,
     * therefore these two values are cached and the coefficients (see find_interpolators) are loaded only  if Mzams or Z are different
     * with respect to the last call.
     *
     * @param Mzams Zams mass in Msun
     * @param Z metallicity
     * @param R radius in Rsun
     * @return value of lambda.
     */
    double operator() (double Mzams, double Z, double R) override;

protected:

    /**
     * Read the already loaded table to find all the Mzams and Z in the Klencki tables
     */
    void fill_interpolators_lists();

    /**
     * Find the four interpolating tracks in the Klencki tables for a given Mzams, Z couple.
     * If Mzams and Z are smaller or larger than the minimum and maximum values in the table,
     * they are internally set to the minimum/maximum values in the table.
     * This function set also the weights
     * @param Mzams zams mass in Msun
     * @param Z  metallicity
     */
    void find_interpolators(double Mzams,  double Z);

    static thread_local std::vector<double> Mzams_list;  /*!< list of all zams masses in the Klencki table */
    static thread_local std::vector<std::vector<double>> Z_list;  /*!< list of all metallicity  in the Klencki table */
    static thread_local bool already_loaded_lists;  /*!< flag to check if Mzams_list and Z_list have been already loaded */

    std::vector<std::vector<double>> vector_cache; /*!< 2d vectors with 4 rows (the interpolating tracks) and 17 columns (each one from a given track in the table) */

private:

    double Mzams_interpolators [2] = {0.,0.}; /*!< Mzams of the four interpolator */
    double Z_interpolators [4] = {0.,0.,0.,0.}; /*!< Z of the four interpolator */
    double wM [2] = {0., 0.}; /*!< interpolator weights for the two Mzams */
    double wZ [4] = {0., 0., 0., 0.}; /*!< interpolator weights for the four metallicities */

};

#endif //SEVN_LAMBDA_KLENCKI21_H
