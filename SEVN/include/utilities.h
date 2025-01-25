//
// Created by mario on 04/05/18.
//

#ifndef SEVN_REVISED_UTILITIES_H
#define SEVN_REVISED_UTILITIES_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <sstream>
#include <iomanip>
#include <regex>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <limits>
#include <map>


#include <sevnlog.h>
using sevnstd::SevnLogging;

#define _UNUSED __attribute__ ((unused))
#define openfile(a, b) utilities::_openfile(a, b, __FILE__, __LINE__)

class Star;
class Binstar;

namespace utilities{

    //random number generator (OpenMP thread-safe)
    thread_local extern std::mt19937_64 mtrand;
    //#ifdef _OPENMP
    //    #pragma omp threadprivate(utilities::mtrand)
    //#endif


    //GI 140520 changed this because a firstprivate clause is added in the main files
    //GI 130321 Changed again because a firstprivate clause could not be used with the current implementation
    //this means that this version is not compilable with INTEL compiler  (beacuse mtrand is defined as extern)
    //extern std::uniform_real_distribution<double> rand_unif_0_1(0.0, 1.0);


    //SEVN NAME
    const std::string SEVN_NAME = "SEVN"; //NOTICE THIS HAS TO BE THE NAME OF THE FOLDER CONTAINING THE CODE

    //TODO: Add  cgs -> SEVN units constants

    //constants
    //TODO The timestep and other time property are in Myr, but some costant and the eq. of some process are in yr, should we use only a timescale?

    ///Fundamental quantitis
    //TODO Check all the constants to be consistent with these values
    //GI, 27/07/2022, all the constant updated to match the definition in astropy-4.3.1
    constexpr double G = 3.925125598496094e8; //RSUN^3 YR^-2 MSUN^-1 (astropy: constant.G.to("Rsun^3/yr^2/Msun"))
    constexpr double yr_cgs = 3.1557600e7; //yr in s (astropy: u.yr.to(u.s))
    constexpr double parsec_cgs= 3.085677581491367e+18; //parsec in cm (astropy: u.parsec.to(u.cm))
    constexpr double Rsun_cgs = 6.95700e10; //rsun in cm (astropy constant.R_sun.to("cm"))
    constexpr double Msun_cgs = 1.988409870698051e+33; //msun in g (astropy constant.M_sun.to("cm"))
    constexpr double G_cgs = G*Rsun_cgs*Rsun_cgs*Rsun_cgs/(Msun_cgs*yr_cgs*yr_cgs); //cm^3 s^-2 g^-1


    constexpr double Sigma_StefBoltz = 7.1694165533435e-17; //LSun^3 RSun^-2 K^-4 (astropy: constant.sigma_sb.to('Lsun/(K^4 * Rsun^2)')
    constexpr double Myr_to_yr = 1.0e6;
    constexpr double yr_to_Myr = 1.0e-6;
    constexpr double AU_to_RSun = 215.03215567054764; //(astropy u.au.to(u.Rsun))
    constexpr double kms_to_RSunyr = 45.360931435963785; //(astropy: (u.km/u.s).to(u.Rsun/u.yr))
    constexpr double LSun_to_Solar = 12.500687924182579; //From Lsun to MSun RSun^2 yr^-3 (astropy: u.Lsun.to((u.Msun*u.Rsun**2)/(u.yr**3)))
    constexpr double c = 1.3598865132357053e7; // RSun/yr (astropy: constant.c.to('Rsun/yr'))
    constexpr double km_to_RSun = 1.4374011786689665e-06; //(astropy: u.km.to(u.Rsun))
    constexpr double parsec_to_Rsun= parsec_cgs/Rsun_cgs;
    constexpr double g_to_MSun = 5.029144215870041e-34; //(astropy: u.g.to(u.Msun))
    constexpr double G_over_c2 = G / (c*c);
    constexpr double G3_over_c5 = (G*G*G)/(c*c*c*c*c); //Scaling for GW processes
    constexpr double tH = 13.7*1e3; //Hubble time in Myr
    constexpr double Mchandra = 1.41; //Chandrasekhar mass in Msun




    const std::string PLACEHOLDER="xxx"; //Standard placeholder for input properties
    //SY

    //MAGIC NULL VALUES
    constexpr double NULL_DOUBLE = -9e30;
    constexpr int NULL_INT = -999999999;
    constexpr size_t NULL_SINT = 999999999;
    //constexpr std::string NULL_STR = "FORZAROMA"; //NOT POSSIBLE in C11 (It will be possible in C20)
    const std::string NULL_STR = "FORZAROMA";

    //MAGIC LARGE AN TINY VALUES
    constexpr double DIFF_TOLL = 1e-10; //Tollerance on difference between two values
    constexpr double LARGE = 1e30;
    constexpr double TINY = 1e-15;
    constexpr double DOUBLE_EPS = std::numeric_limits<double>::epsilon();


    //INT CONSTANT TO HANDLE RETURN FROM EVOLUTION
    typedef unsigned int evolution;
    constexpr int SINGLE_STEP_EVOLUTION=0;
    constexpr int REPEATED_EVOLUTION=1;

    //INT CONSTANT TO HANDLE RETURN FROM FUNCTIONS
    typedef unsigned int jump_convergence;
    constexpr int JUMP_CONVERGE=0;
    constexpr int JUMP=1;
    constexpr int NO_JUMP=2;

    //INT CONST FOR SN EXPLOSION
    typedef unsigned int sn_explosion;
    constexpr int SNIA_EXPLODE=1;
    constexpr int SNII_EXPLODE=2;
    constexpr int SN_NOT_EXPLODE=0;

    //INT CONST FOR RLO
    typedef unsigned int rlo;
    constexpr int RLO_FALSE=0; //RLO is happening or happened
    constexpr int RLO_TRUE=1; //RLO is happening or happened

    //bool CONST FOR BINARY EVOLUTION
    typedef bool bse_evolution;
    constexpr int BIN_EV_DONE = 1; //This is the return if the properties of the binary have been evolved with the proper evolve method
    constexpr int BIN_EV_NOT_DONE = 0; //This is the return if the properties of the binary have been evolved without the proper evolve method
    constexpr int BIN_EV_SETBROKEN = 2; //This is the return if the properties of the binary have not been evolved but a set broken has been called

    double maxwellian_cdf(double x, double sigma);
    double maxwellian_pdf(double x, double sigma);
    inline double R_Schwarzschild(double Mass){return 2.0*G_over_c2*Mass; } //rs=2GM/c^2}


    /**
 * @tparam T template type. It should be able to accept anything that can be << to a stream obejct.
 * @param val  value of type T to be converted to a string
 * @param file_input it should be set to  \code{.cpp} __FILE__  \endcode
 * @param line_input it should be set to \code{.cpp} __LINE__  \endcode
 * @param precision  set precision for the scientific format output of numbers
 * @return The \p val converted to a string. If val is a number the string is the scientific format with  precision
 * equal to \p precision.
 */
    template <typename T> const std::string n2s(T val, const char* file_input, const int line_input, const unsigned int precision=6) {
        SevnLogging svlog;

        std::ostringstream stream;
        stream << std::scientific << std::setprecision(precision);
        stream << val;

        if (stream.fail())
            svlog.critical("Cannot convert into a string", file_input, line_input);

        return stream.str();
    }


    ////Phys

    /**
     * Auxiliary struct containing the mass content of a star
     */
    struct MassContainer{
        const double Mass;  //Mass
        const double MHE;  //Mass of the HE core
        const double MCO;  //Mass of the CO core
    };

    /**
     * Estimate the RL follwing the Eggleton formalism (Eq. 53 Hurley02)
     * @param Mass_primary  Mass of the primary in Msun (the star for which we are estimating the RL)
     * @param Mass_secondary  Mass of the secondary in Msun.
     * @param a Semimajor axis in Rsun
     * @return Roche Lobe radius in Rsun
     */
    double roche_lobe_Eg(double Mass_primary, double Mass_secondary, double a);

    /**
     * Estimate the Alfven radius
     * @param s Pointer to the star for which we want to calculate the Alfven radius
     * @param dMdt Accreated mass rate in Msun/yr
     * @param get0 AIf true use the stellar propertie at the beginning of the timestep
     * @return The Alfven radius in Rsun
     */
    double R_Alfven(Star *s, double dMdt, bool get0=false);

    /**
     * Estimate the Hydrogen mass fraction of the star.
     * It is the one used in BSE/MOBSE
     * @param s Pointer to the star
     * @return Hydrogen mass fraction
     */
    double Hfrac(Star *s);

    /**
     * Estimate the maximal accretion rate on an object due to the Eddington limit.
     * It is estimated using Eq. 67 in Hurley+02.
     * In addition it is multiplied by the factor eddfact (>0, default=1) to allow super-eddingont accretion
     * NOTICE the return will be in units of Myr / YR, while SEVN usually assumes Myr as scale, be careful of conversion.
     * @param donor Pointer to the star that is donating mass
     * @param accretor Pointer to the star that is accreting mass
     * @param eddfact Multiplicative factor to allow supereddington accretions
     * @return return the Eddington accretion rate  multiuplied by eddfact in Msun/yr
     */
    double dMdt_Eddington_accretion(Star *donor, Star *accretor, double eddfact=1.0);


    /**
     * Critical angular velocity
     * @param Mass mass in Msun
     * @param Rpolar  polar radius in Rsun
     * @return critical angular velocity in yr^-1
     */
    inline static double omega_crit(double Mass, double Rpolar){
        double Reqc = 1.5 *Rpolar;
        return std::sqrt(utilities::G*Mass/(Reqc*Reqc*Reqc));
    }


    /**
     * Like wait but it works also outside debug mode
     */
    inline void hardwait(){
        std::cout<<"\nWaiting"<<std::endl;
        std::cin.get();
    }
    /**
     * Like wait but it works also outside debug mode
     */
    template<typename T, typename... Tail>
    void hardwait(T head, Tail... tail){
        std::cout << head << " ";
        hardwait(tail...);
    }





    /**
     * Base class for the variadic function wait.
     * @param _message
     */
    inline void wait(){
    #ifdef DEBUG
        std::cout<<"\nWaiting"<<std::endl;
        char _fake;
        std::cin>>_fake;
    #endif
    }

    /**
     * Variadic function. It prints a message in std::cout (from an argument pack) and wait of a cin input from the user.
     * @tparam T template param of the value to be written
     * @tparam Tail  template params of the pack of parameters to iterate over
     * @param head value to be written
     * @param tail rest of the packs
     */
    template<typename T, typename... Tail>
    void wait(_UNUSED T head, _UNUSED Tail... tail){
    #ifdef DEBUG
            std::cout << head << " ";
            wait(tail...);
    #endif
    }

    /********** Template variadic function to be used in the log *******/

    std::string get_name(Star* s);
    long get_ID(Star* s);
    std::string get_name(Binstar* b);
    long get_ID(Binstar* b);
    double get_current_time(Star* s);
    double get_current_time(Binstar* b);




    template <typename T>
    void _log_print_core(std::stringstream &ss, T t){
        ss<<t;
        return;
    }

    template <typename T, typename... ListP>
    void _log_print_core(std::stringstream &ss, T t, ListP... args){

        ss<<t<<":";
        _log_print_core(ss,args...);

        return;
    }

    template <class System, typename... ListP>
    std::string common_log_print(const std::string& label, System* system, ListP... args){

        std::stringstream ss;
        _log_print_core(ss,args...);

        std::string id_str=std::to_string(get_ID(system));
        std::string time = utilities::n2s(get_current_time(system),__FILE__,__LINE__);

        return get_name(system) +";" +id_str + ";" +label  + ";" + time + ";" + ss.str();

    }

    template <class System, typename... ListP>
    std::string common_log_print(const std::string& label, System* system){


        std::string id_str=std::to_string(get_ID(system));
        std::string time = utilities::n2s(get_current_time(system),__FILE__,__LINE__);

        return get_name(system) +";" +id_str + ";" +label  + ";" + time + ";";

    }


    template <typename... ListP>
    std::string log_print(const std::string& label, Star* star, ListP... args){

        return "S;" + common_log_print(label,star,args...);
    }

    template <typename... ListP>
    std::string log_print(const std::string& label, Binstar *binstar, ListP... args){

        return "B;" + common_log_print(label,binstar,args...);
    }

    /**
     * Utility to write Star info in a log format
     * @param s Pointer to star
     * @param oldstep if true get the property of the last step
     * @return a string with ID:Mass:MHE:MCO:Phase:RemnanType
     */
    std::string log_star_info(Star* s, bool oldstep=false);


    /*************************************************************************************/


    unsigned long gen_rseed();
    unsigned long gen_rseed(std::random_device &rd);

    inline const std::string random_keygen(std::mt19937_64 *mtrand){

        //std::string allkeys = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        std::string keys_init = "123456789"; //GI71119: MORE SIMILAT TO A source_id GAIA like (it can easier to load and manage as a pure number table)
        std::string keys = "0123456789";

        std::uniform_int_distribution<int> unifo_keys_init(0, (int)keys_init.size()-1);
        std::uniform_int_distribution<int> unifo_keys(0, (int)keys.size()-1);

        std::string key = "";
        key += keys_init[(unifo_keys_init(*mtrand))];
        for(int i = 0; i < 14; i++){
            int position = unifo_keys(*mtrand);
            key += keys[position];
        }

        return key;

    }

    std::vector<std::string> split(const std::string& s, char delimiter);



    template <class T> inline bool isifstream();
    template <> inline bool isifstream<std::ifstream>() {return true;}

    //general functions
    template <typename T> size_t binary_search(T *array, const size_t left, const size_t right, const T value){


        if (right > left+1) {

            const size_t mid = left + ((right - left)>>1);
            //std::cout<<left<<"   "<<mid<<"   "<<right<<std::endl;
            //std::cout<<"****"<<array[left]<<std::endl;
            //std::cout<<array[mid]<<std::endl;
            //std::cout<<array[right]<<std::endl;

            if(array[mid] > value)
                return binary_search(array, left, mid, value);
            else if (array[mid] == value)
                return mid;
            else
                return binary_search(array, mid, right, value);
        }

        if(array[left] == value) return left;
        else if(array[right] == value) return right;
        else return left;

    }


    template <typename T> bool string_is_number(std::string str) {
        T value;

        std::stringstream stream(str);
        stream >> value;

        return (!stream.fail() and stream.eof()) ;
    }

    /**
     * Solve the Kepler equation to found the Eccentric anomaly starting from eccentricy and Mean anomaly.
     * We use the generalized Newton Raphson’s method  (Eq. 2.6, Nazeer 2016).
     * If in some step we have numerical problem with the second derivative, we revert to a simple NR method for that step
     * @param ecc  Eccentricity
     * @param m Mean anomaly
     * @param tol tollerance
     * @param maxit maximum number of iteration
     * @return Eccentric anomaly
     */
    inline double kepler(const double &ecc, const double &m, const double tol = 1e-6, const int maxit = 50) {  // it solves the Kepler's equation giving the mean anomaly and the eccentricity

        //We use a Generalize Newton-Raphson method that converge cubically (instead of quadratically for the classical  Newton Raphson)
        //Given the Kepler Equation we can have problem if f1 is 0, if f2 is 0 or if f1*f1-2*f*f2 is <=0
        //f1=1-ecc*cos(E) but ecc<1 always so f1<1 as well. When f2 is 0 or f1*f1-2*f*f2 is <=0 we revert to a simple NR for that step

        int k = 1; // to check max iterations
        double em = m; //em = eccentric anomaly (output), m = mean anomaly (input)
        double sinem = sin(em);

        //If eccentricity is 0  we have M=E, i.e. the Eccentricity anomaly  is equal to the Mean anomaly
        if (ecc==0)
            return m;

        double f = em - ecc*sinem - m;
        double f1 = 1.0 - ecc*cos(em);
        double f2 = ecc*sinem;
        double fsquare=f1*f1 - 2*f*f2;
        auto check_problem = [&f2,&fsquare](){ return f2==0 or fsquare<=0;}; //Lmabda expression to check if I have to pass to a NR


        if (check_problem()) //If the second derivatie is 0 make a step using the first order Newthon Methond
            em = em - f/f1;
        else
            em = em + (-f1 + sqrt(fsquare))/f2;

        // iterations
        while(fabs(f) > tol && k < maxit){



            sinem = sin(em);
            f = em - ecc*sinem - m;
            f1 = 1.0 - ecc*cos(em);
            f2 = ecc*sinem;
            fsquare=f1*f1 - 2*f*f2;


            if (check_problem()) //If the second derivatie is 0 make a step using the first order Newthon Methond
                em = em - f/f1;
            else
                em = em + (-f1 + sqrt(fsquare))/f2;



            k++;
        }

        return em;

    }

    //templates
    template <typename T> const T s2n(std::string &str, const char* file_input, const int line_input) {
        SevnLogging svlog;

        T value;


        std::stringstream stream(str);
        stream >> value;


        if (stream.fail())
            svlog.critical("Cannot convert the string " + str + "into a number", file_input, line_input,sevnstd::sevnio_error());

        return value;
    }


    template <typename T> void _openfile(T &in, const std::string f, const char* file_input, const int line_input){
        SevnLogging svlog;


        if(in.is_open()) in.close();
        if(utilities::isifstream<T>())
            in.open(f.c_str(), std::ios::in);
        else
            in.open(f.c_str(), std::ios::out);

        if(!in) svlog.critical("Cannot open file " + f, file_input, line_input,sevnstd::sevnio_error());
    }




    //useful function to sort arrays and indexes
    inline bool wayToSort(int i, int j) { return i > j; }


    template <typename T> T dirname2n(std::string str, const char* file_input, const int line_input){

        std::size_t found = str.find('.');
        if (found == std::string::npos)
            str.insert(1,".");

        return utilities::s2n<T>(str, file_input, line_input);

    }


  //GI
   /**
   *  Functions to print all the elements of a vector in a Python numpy style (GI)
   */
  template <typename T> void print_vector(const std::vector<T>& v){
      std::cout<< "[ ";
      for (auto& element : v)
          std::cout << element << " ";
      std::cout<< "]"<<std::endl;
  }


  //GI
  /**
   *  Functions to generate a complete filename path containing also the directory
   *   @param _folder   A string containing the folder path
   *   @param _fname    Complete filename including also the .extension if any (_fname=_fname_root . _fname_extension)
   *   @param print_threads     If true the number of threads will be added to the filename
   *   @return  a string with  _folder/_fname_root_nthreads._fname_etension if print_threads=true otherwise _folder/_fname_root._fname_extension
  */
  inline std::string gen_filename(const std::string &_folder, const std::string &_fname, bool print_threads=true){

      std::string return_string;
      std::string folder = _folder.back()=='/' ? _folder.substr(0, _folder.length()-1) : _folder; //GI 141219: To avoid to have a double / if the _folder in input already contains a / at the end

      if (print_threads){

          std::size_t found_extension;

          //GI 141219: Simple loop to get only the last occurence of . in the string
          for ( std::size_t pos=0; pos!=std::string::npos; pos=_fname.find('.',pos+1)) found_extension = pos;
          if (found_extension==0) found_extension = _fname.length();

          #ifdef _OPENMP
                  int num_thread=omp_get_thread_num();
          #else
                  int num_thread=0;
          #endif

          return_string = folder + "/" + _fname.substr(0, found_extension) + "_" + std::to_string(num_thread) + _fname.substr(found_extension);
      } else {

          return_string = folder + "/" + _fname;
      }

      return return_string;
  }


  /**
   * Get the absolute path of this SEVN folder
   * @return The absolute SEVN path
   */
  inline std::string get_absolute_SEVN_path(){
        std::string here_token = "include/";
        std::string sevn_path = __FILE__;
        std::size_t found = sevn_path.find("/src/");


        static sevnstd::SevnLogging svlog;
        if (found!=std::string::npos){
            sevn_path = sevn_path.substr(0,found);
        }
        else{
            svlog.critical("Error in the initialisation of parameter myself (path of the SEVN folder)."
                           " The path finding it is based on the fact that the utilities.h file is inside the folder include "
                           "in the main SEVN folder. Did you change it?",__FILE__,__LINE__,sevnstd::sevnio_error());
        }

        return sevn_path;
  }

  /**
   * Find the slope and the intercept of the line passing through (x1,y1), (x2,y2)
   * @param x1 x-value of the first point
   * @param x2 x-value of the second point
   * @param y1 y-value of the first point
   * @param y2 y-value of the second point
   * @param slope this value stores the estimated slope
   * @param intercept this value store the estimated intercept
   * @return
   */
  inline int find_line(const double & x1, const double & x2, const double & y1, const double & y2, double & slope, double & intercept){
      slope = (y2 - y1)/(x2-x1);
      intercept = y2 -slope*x2;

      return EXIT_SUCCESS;
  }

  template<typename T>
  double rel_difference(T val1, T val2){

        return fabs( (val1-val2)/val1 );
    }



  inline void swap_stars(Star* & s1, Star* & s2){
        Star *stmp=s1;
        s1=s2;
        s2=stmp;
    }

    /**
   * Trim a string from all whitespaces
   * @param s
   * @return
   */
  inline std::string trim(const std::string& s) {
        return std::regex_replace(s, std::regex("^[ \\s]+|[ \\s]+$"), std::string(""));
  }

  /**
   * Template function to check if a given element is inside a container.
   * @tparam T typename of the element to be checked
   * @tparam Iter typename of the iterator
   * @param element element to be checked
   * @param it iterator to the begin of the container
   * @param end iteratore to the end of the container
   * @return true if the element is inside the list, otherwise false.
   */
  template<typename T, typename Iter>
  bool isinlist(T element, Iter it, Iter end){

      return std::find(it, end, element)!=end;
  }

  /**
   * Make the string %plife:phase
   * @param plife Life percentage at a given Phase
   * @param Phase Phase to consider
   * @return string %plife:phase
   */
  inline std::string make_pfile_str(const double plife, const size_t Phase, const unsigned int min_precision=6){

      //In the following rows we estimate plife and the transform the number to a string to initialise stars.
      //Hovewer, if plife is very close to 0 or 1 without enough number of digits we can artificially force plife to be 0 or 1
      //Theregore we dynamically set the precision estimating the difference between plife and 0 or 1, taking the exponent of log10
      //and using the next larger integer.
      unsigned int precision; //Set the precision to transform pfile from number to string
      unsigned int digit_0= std::ceil(std::abs(std::log10(std::abs(plife-0))));
      unsigned int digit_1= std::ceil(std::abs(std::log10(std::abs(plife-1))));
      precision=std::max(std::max(digit_0,digit_1),min_precision);

      //Consistency check
      std::stringstream tini_ss; //string containing the plife (see below)
      //Close to plife=1 is necessary to use a large number of digits to avoid to set plife=1 when trasformi to string
      tini_ss << "%" << std::setprecision(precision)  << plife * 100 << ":" << Phase; //Starting time

      return tini_ss.str();
  }


  //Interpolator
  /**
   * Estimate the y value using 1D linear interpolation
   * @tparam T  type of the  of interpolated value
   * @param xp
   * @param x_interp  1D vector containing the x-value of the interpolating tables (in ascending order)
   * @param y_interp  1D vector containing the y-value of the interpolating tables
   * @param equispaced_interval if true the value in x_interp are equi-spaced
   * @param ext_raise if true and xp is out of boundary raise an error otherwise return the extrema of yinterp
   * @return interpolated y value at xp.
   */
  template <typename T> T interpolate_1D(T xp,  std::vector<T>& x_interp,  std::vector<T>& y_interp, bool equispaced_interval=false, bool ext_raise=false){


      if (xp<x_interp[0] and ext_raise)
          throw sevnstd::sevnerr("Error in interpolate_1D in utility.h: xp is out of boundary.");
      else if (xp<=x_interp[0])
          return y_interp[0];

      if (xp>x_interp.back() and ext_raise)
          throw sevnstd::sevnerr("Error in interpolate_1D in utility.h: xp is out of boundary.");
      else if (xp>=x_interp.back())
          return y_interp.back();

      size_t pos;
      if (equispaced_interval){
          double dx   = x_interp[1]-x_interp[0];
          pos = int( (xp-x_interp[0])/dx); //int return floor
      }
      else{
          pos = binary_search(&x_interp[0], 0, x_interp.size()-1, xp);
      }
      return (y_interp[pos+1]-y_interp[pos])/(x_interp[pos+1]-x_interp[pos])*(xp-x_interp[pos])+y_interp[pos];
  }

  /**
   * Get a portion of path name till the given split_string. The split_string is include.
   * @param path  Complete path
   * @param split_string  string patter where to cut the complete path
   * @return the substring from 0 to split_string
   */
  inline std::string get_subpath(std::string path, std::string split_string, bool include_split_string=true){

      size_t tt= path.find(split_string);
      std::string subpath;

      if (include_split_string){
          subpath = path.substr(0,tt+split_string.size());
      } else{
          subpath = path.substr(0,tt);
      }

      return subpath;
  }

  /**
   * Transpose  a matrix
   * @tparam T
   * @param MatrixT Transposed Matrix
   * @param Matrix  Matrix to be transposed
   */
  template <typename T>
  void transpose(std::vector<std::vector<T>>& MatrixT, std::vector<std::vector<T>>& Matrix){

      MatrixT.resize(Matrix.size());

      for (auto& Matrix_row : Matrix){
          for (int j=0; j<(int)Matrix_row.size(); j++)
              MatrixT[j].push_back(Matrix_row[j]);
      }

  }


    /**
     *   Generic function to find an element in vector and also its position.
     * @tparam T
     * @param vecOfElements Vector to look for the element
     * @param element  element to find in the vector
     * @return   It returns an integer
      int : Represents the index of element in vector if its found else -1
     */
    template < typename T>
    int findInVector(const std::vector<T>  & vecOfElements, const T  & element)
    {
        int result ;
        // Find given element in vector
        auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
        if (it != vecOfElements.end())
        {
            result = distance(vecOfElements.begin(), it);

        }
        else
        {
            result=-1;
        }
        return result;
    }


    template <typename Key, typename Value>
    std::map<Value,Key> flip_map(const std::map<Key,Value> &original_map){
        std::map<Value,Key> flipped_map;

        for (auto& pair : original_map){
            flipped_map[pair.second] = pair.first;
        }

        return flipped_map;
    }

    /**
     * Check if two values can be considered equal within the machine precision tolerance.
     * This can be used as a robust alternative to the operator == for doubles
     * Absolute tolerance works well for values<1, relative tolerance for larger values
     * @param x First value to compare
     * @param y Second value to compare
     * @return True if the |x-y|<=epsilon, where
     *          - epsilon=std::numeric_limits<double>::epsilon() if x<=1 and y<=1
     *          - epsilon=std::numeric_limits<double>::epsilon()*max(x,y) otherwise
     */
    bool areEqual(double x, double y);

    /**
     * Estimate the smallest possible step given X so that X+step>X
     * @param x number
     * @return the step is estimated as 1.01*std::numeric_limits<double>::epsilon()*maxXYOne
     * where maxXYOne = std::max( 1.0, fabs(x));
     */
    double smallestSignificativeStep(double x);


    /**
     * Define make_unique for pre-C++14. Notice this is the  the same of the implementation of C++14 standard for single object,
     * in the std there is also another implementation for arrays.
     * @tparam T Type of the unique ptr
     * @tparam Args
     * @param args arg to initialise the object pointed by the unique ptr
     * @return return the unique pointer to the element of type T.
     */
    template<typename T, typename... Args>
    inline std::unique_ptr<T> make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    //!  Class to handle the generation of a list of values
    /*!
      The class usage is based on the methods get(), next() and empty().
      get() return the current value, next() upate the value and empty() return true if the generator reach the ends.
      If empty() is true calling next() or get() throws an out_of_range error.
    */
  class ListGenerator{

  private:

                  double vcurrent;

                  double vstep=std::nan("");
                  double vmin=0.;
                  double vmax=2E30;

                  std::vector<double> vlist;
                  std::vector<double>::iterator begin;
                  std::vector<double>::iterator end;
                  std::vector<double>::iterator  current;
                  bool end_of_list=false;

                  inline void initialise_tlist_iterators(){
                      begin=vlist.begin();
                      end=vlist.end();
                      current=begin;
                  }
  public:
                  ListGenerator(){
                      vlist={};
                      initialise_tlist_iterators();
                  }

                  /**
                   * Constructor based on a step value. The values will be generated from _vstep_min to _vstep_max (can be not included)
                   * with step _vstep. If _vstep_max is not given the generation can last theoretically forever (_vstep_max default is 1E30).
                   * @param _vstep step value to generate the values
                   * @param _vstep_max maximum value.
                   * @param _vstep_min Starting value.
                   */
                  explicit ListGenerator(double _vstep, double _vstep_max=std::nan(""), double _vstep_min=std::nan(""))
                  : vstep{_vstep}, vmin{_vstep_min}, vmax{_vstep_max} {

                      if (vstep<=0)
                          throw std::runtime_error("ListGenerator::vstep cannot be negative or zero");
                      if (std::isnan(vmin))
                          vmin=vstep;
                      if (std::isnan(vmax))
                          vmax=2E30;

                      if(vmin >= vmax)
                          throw std::runtime_error("ListGenerator::vmin cannot be larger than ListGenerator::vmax");

                      vcurrent=vmin;
                  }

                  /**
                   * Constructor based on a vector of doubles. The class will slice all the values in the vector and will be marked as
                   * empty when the end of the vector is reached.
                   * @param _tlist Vector of doubles
                   */
                  explicit ListGenerator(const std::vector<double>& _tlist) : vlist{_tlist} {

                      if (!std::is_sorted(vlist.begin(),vlist.end()))
                          throw std::runtime_error("The input vector in the ListGenerator constructor is not sorted");

                      initialise_tlist_iterators();
                      if (!vlist.empty()){
                          vcurrent=*current;
                          vmin=*begin;
                          vmax=*(end-1);
                      }
                  }

                  static std::unique_ptr<ListGenerator> make_unique(double _vstep, double _vstep_max=std::nan(""), double _vstep_min=std::nan(""));
                  static std::unique_ptr<ListGenerator> make_unique(std::vector<double> _vlist);

                  inline double get() const{
                      if (empty())
                          throw std::out_of_range("The list of times reached the end. The current value is undefined");
                      return vcurrent;
                  }
                  inline double get_max() const {return vmax;}
                  inline double get_min() const {return vmin;}

                  inline bool empty() const {return end_of_list;}
                  inline void next(){
                      //If empty just throw an out_of_range error
                      if (empty())
                          throw std::out_of_range("The list of times reached the end. The next value is undefined");
                      //If we are using the vstep implementation and the next step is  going beyond vmax flag end_of_list
                      else if(!std::isnan(vstep) and vcurrent+vstep>vmax)
                          end_of_list=true;
                      //If we are using the vstep implementation and the next step is  not going beyond vmax just increment vcurrent
                      else if(!std::isnan(vstep))
                          vcurrent+=vstep;
                      //Now start with the vlist implementation check
                      //Update the iterator and check if we   reach the end
                      else if(++current==end){
                          end_of_list=true;
                      }
                      //If we are here We are still inside the list, update the iterator and vcurrent.
                      else{
                          vcurrent=*current;
                      }
                  }
                  inline double operator++(){
                      next();
                      return vcurrent;
                  }
                  inline double operator++(int){
                      double old_t=vcurrent;
                      next();
                      return old_t;
                  }
                  /**
                   * Similar to next, but just return the next value withotu uptading anything
                   * @return The next value if not empty, otherwise std::nan
                   */
                  inline double forecast() const {
                      //If empty just throw an out_of_range error
                      if (empty())
                          return std::nan("");
                      else if(!std::isnan(vstep) and vcurrent+vstep>vmax)
                          return std::nan("");
                      else if(!std::isnan(vstep))
                          return  vcurrent+vstep;
                      else if(current+1==end)
                          return std::nan("");
                      else
                          return *(current+1);
                  }

              };



}


#endif //SEVN_REVISED_UTILITIES_H

