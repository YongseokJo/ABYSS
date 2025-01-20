//
// Created by iorio on 16/03/20.
//

#ifndef SEVN_PARAMS_H
#define SEVN_PARAMS_H

#include <map>
#include <string>
#include <errhand.h>
#include <utilities.h>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>
#include <sevnlog.h>
using sevnstd::SevnLogging;
#include <lookup_and_phases.h>
using namespace Lookup;

/***
 * This class stores and handles the various parameters used in SEVN.
 *
 * The parameters are divided in two big groups:
 *      - params_num: these are numerical parameters stored as doubles (even if they are integer);
 *      - params_str: these are litteral parameters stored as strings.
 *      - params_bool:
 *
 *  Each parameter is composed by three value:
 *      - a name (string),
 *      - a value (double for params_num and string in params_str),
 *      - a short documentation (string).
 *
 * The parameters are also divided in groups depending on their "set" behaviour:
 *      - NOT SETTABLE: these parameters cannot be changed at runtime and have to be hardcoded in the class;
 *      - PRIVATE SETTABLE: these parameters cannot be set directly at runtime, but through auxilary methods (e.g. load);
 *      - PUBLIC SETTABLE: these parameters can be set directly at runtime with the method set(name, value).
 *
 *  In order to retrieve a given parameter there are two get methods:
 *      - get_num(name) to retrieve a numerical parameter,
 *      - get_str(name) to retrieve a str parameter.
 *      - get_str(bool) to retrieve a str parameter.
 *
 *  This class contains also an override version of the operator <<, so that it can print in a fancy way the
 *  current value of the parameters.
 *
 *  @note There is a safety check in the constructors so that only one instance of this class can be present in code
 *  (to avoid to have different istance with different parameters value).
 */

typedef std::pair<double, std::string> num_entry;
typedef std::pair<std::string, std::string> str_entry;
typedef std::pair<bool, std::string> bool_entry;
typedef std::map<std::string, num_entry> map_num;
typedef std::map<std::string, str_entry> map_str;
typedef std::map<std::string, bool_entry> map_bool;

class SEVNpar{

public:

    ///Constructor
    /**
     * Default constructor: just the parameters to the default value.
     */
    SEVNpar(){
        if (++instance_counter<2)
            init();
        else
            svlog.critical("Only one active instance of SEVNpar is allowed",__FILE__,__LINE__,sevnstd::params_error());
    }
    /**
     * Constructor to load parameters from line command
     * @param argc Number of arguments
     * @param argv arguments
     */
    SEVNpar(int argc, char **argv){
        if (++instance_counter<2)
            load(argc,argv,true);
        else
            svlog.critical("Only one active instance of SEVNpar is allowed",__FILE__,__LINE__,sevnstd::params_error());
    }

    ~SEVNpar() {instance_counter--;}

    //load
    /**
     * load parameters from a list of arguments
     * @param n  number of arguments
     * @param val arguments
     * @param initialise if true initialise all the parameters to default before to set them
     * @return EXIT_SUCCESS or throw an error
     *
     * @throws sevnstd::params_error Thrown if one of the parameters is not settable.
     *
     * @note: the function checks only whether a parameter is not settable, it skips all the parameters
     * that are not present in the class.
     *
     */
    int load(int n, char **val, bool initialise=false);


    ///GET and GET like stuff
    inline const map_num& get_num_map(){
        return params_num;
    }

    inline const map_str& get_str_map(){
        return params_str;
    }

    inline const map_bool& get_bool_map(){
        return params_bool;
    }


    /**
     * get a numerical variable
     * @param name  name of the parameter
     * @return the parameter value or throw an error
     *
     * @throws sevnstd::params_error Thrown if the parameter name is not present in the numerical group.
     */
    inline double get_num(std::string name){

        auto it=params_num.find(name); //Check if the key is in the map

        if (it!=params_num.end())
            return params_num[name].first;
        else {
            svlog.critical("key " + name + " not present", __FILE__, __LINE__, sevnstd::params_error());
            return -1.0;
        }
    }
    /**
     * get a literal variable
     * @param name of the parameter
     * @return parameter value or throw an error
     *
     * @throws sevnstd::params_error Thrown if the parameter name is not present in the str group.
     */
    inline std::string get_str(std::string name){
        auto it=params_str.find(name); //Check if the key is in the map
        if (it!=params_str.end())
            return params_str[name].first;
        else {
            svlog.critical("key " + name + " not present", __FILE__, __LINE__, sevnstd::params_error());
            return "error";
        }
    }
    /**
     * get a boolean variable
     * @param name of the parameter
     * @return parameter value or throw an error
     *
     * @throws sevnstd::params_error Thrown if the parameter name is not present in the str group.
     */
    inline bool get_bool(std::string name){

        auto it=params_bool.find(name); //Check if the key is in the map

        if (it!=params_bool.end())
            return params_bool[name].first;
        else
            svlog.critical("key "+name+" not present",__FILE__,__LINE__,sevnstd::params_error());

        return false;
    }
    /**
     * Store directly the value of a numerical parameters in an input variable.
     * @param name  name of the parameter
     * @param val   variable where to store the parameter's value
     * @return EXIT_SUCCESS or thrown an error
     *
     * @throws sevnstd::params_error Thrown if the parameter name is not present.
     */
    inline int get(std::string name, double& val){
        val=get_num(name);
        return EXIT_SUCCESS;
    }
    /**
     * Store directly the value of a literal parameters in an input variable.
     * @overload
     */
    inline int get(std::string name, std::string& val){
        val=get_str(name);
        return EXIT_SUCCESS;
    }
    /**
     * Store directly the value of a boolean parameters in an input variable.
     * @overload
     */
    inline int get(std::string name, bool& val){
        val=get_bool(name);
        return EXIT_SUCCESS;
    }

    /**
     * Store all the names of the parameters in a vector
     * @return vector of string containing of all the parameters names.
     */
    std::vector<std::string> keys(){

        std::vector<std::string> keys_list;
        for (auto& key :  params_num)
            keys_list.push_back(key.first);
        for (auto& key :  params_str)
            keys_list.push_back(key.first);

        return keys_list;
    }

    ///SET and SET like stuff
    //These sets are only available for properties listed in public_settable
    /**
     * Set the value of a given numeric parameter
     * @param name name of the parameter.
     * @param val value to save in the parameter.
     * @return EXIT_SUCCESS of thrown an error
     *
     * @throws sevnstd::params_error Thrown if:
     *  - the parameter is not directly settable,
     *  - the parameter is not present,
     *  - the parameter's value fails the check(s).
     *
     */
    inline int set(std::string name, double val){

        auto it=public_settable.find(name); //Check if the key is in the map
        if (it!=public_settable.end())
            params_num[name].first=val;
        else
            svlog.critical("Trying to set key "+name+". The key is not present or it is not directly settable",__FILE__,__LINE__,sevnstd::params_error());

        if (!check()){
            svlog.critical("Check failed", __FILE__,__LINE__,sevnstd::params_error());
        }

        return EXIT_SUCCESS;
    }
    /**
     * Set the value of a given literal parameter
     * @overload
     */
    inline int set(std::string name, std::string val){
        auto it=public_settable.find(name); //Check if the key is in the map
        if (it!=public_settable.end())
            params_str[name].first=val;
        else
            svlog.critical("Trying to set key "+name+". The key is not present",__FILE__,__LINE__,sevnstd::params_error());

        if (!check()){
            svlog.critical("Check failed", __FILE__,__LINE__,sevnstd::params_error());
        }


        return EXIT_SUCCESS;
    }

    ///SPECIAL SET

    /**
     * Set the min and max Z
     * @param Zlist Vector of doubles containing the list of loaded metallicites (loaded  in IO).
     * @return EXIT_SUCCESS or throw error.
     *
     *
     * @throws sevnstd::params_error Thrown if:
     *  - check_z is false or return an error.
     *
     */
    inline int set_zlimit(const std::vector<double >& Zlist){
        //Direct set
        params_num["min_z"].first=*std::min_element(Zlist.begin(), Zlist.end());
        params_num["max_z"].first=*std::max_element(Zlist.begin(), Zlist.end());
        if (check_z())
            return EXIT_SUCCESS;
        else
            svlog.critical("Paramters check failed while settin min_Z, max_Z",__FILE__,__LINE__,sevnstd::params_error());

        return EXIT_FAILURE;

    }

    /**
     * Set the min and max Z of the pureHE tables
     * @param Zlist_HE Vector of doubles containing the list of loaded metallicites (loaded  in IO).
     * @return EXIT_SUCCESS or throw error.
     *
     *
     * @throws sevnstd::params_error Thrown if:
     *  - check_z is false or return an error.
     *
     */
    inline int set_zlimit_HE(const std::vector<double >& Zlist_HE){
        //Direct set
        params_num["min_z_he"].first=*std::min_element(Zlist_HE.begin(), Zlist_HE.end());
        params_num["max_z_he"].first=*std::max_element(Zlist_HE.begin(), Zlist_HE.end());
        if (check_z_he())
            return EXIT_SUCCESS;
        else
            svlog.critical("Paramters check failed while setting min_Z_he, max_Z_he",__FILE__,__LINE__,sevnstd::params_error());

        return EXIT_FAILURE;

    }
    /**
     * Set the min and max ZAMS Mass.
     * Notice: the max Zams is the minimum between the zams upper limits considering all the metallicites.
     * The min zams is the maximum between the zams lower limits.
     * @param Zams_list  Vector of vector of doubles containing the list of zams Mass at each metallicity (loaded  in IO).
     * @return EXIT_SUCCESS or throw error.
     *
     * @throws sevnstd::params_error Thrown if:
     *  - check_mzams is false or return an error.
     *
     */
    inline int set_zams_limit(const std::vector<std::vector<double>>& Zams_list){


        // Notice: the max Zams is the minimum between the zams upper limits considering all the metallicites.
        // The min zams is the maximum between the zams lower limits.
        double minzams=-1e30, maxzams=1e30;
        for (auto& zams_at_Z : Zams_list){
            size_t last_track = zams_at_Z.size();
            minzams = zams_at_Z[0]>minzams ? zams_at_Z[0] : minzams;
            maxzams = zams_at_Z[last_track-1]<maxzams ? zams_at_Z[last_track-1] : maxzams;
        }

        params_num["min_zams"].first=minzams;
        params_num["max_zams"].first=maxzams;

        if (check_mzams())
            return EXIT_SUCCESS;
        else
            svlog.critical("Paramters check failed while settin min_zams, max_zams",__FILE__,__LINE__,sevnstd::params_error());

        return EXIT_FAILURE;
    }
    /**
     * Set the min and max ZAMS pureHE Mass.
     * Notice: the max Zams is the minimum between the zams upper limits considering all the metallicites.
     * The min zams is the maximum between the zams lower limits.
     * @param Zams_list_HE  Vector of vector of doubles containing the list of zams Mass at each metallicity (loaded  in IO).
     * @return EXIT_SUCCESS or throw error.
     *
     * @throws sevnstd::params_error Thrown if:
     *  - check_mzams is false or return an error.
     *
     */
    inline int set_zams_limit_he(const std::vector<std::vector<double>>& Zams_list_HE){


        // Notice: the max Zams is the minimum between the zams upper limits considering all the metallicites.
        // The min zams is the maximum between the zams lower limits.
        double minzams=-1e30, maxzams=1e30;
        for (auto& zams_at_Z : Zams_list_HE){
            size_t last_track = zams_at_Z.size();
            minzams = zams_at_Z[0]>minzams ? zams_at_Z[0] : minzams;
            maxzams = zams_at_Z[last_track-1]<maxzams ? zams_at_Z[last_track-1] : maxzams;
        }

        params_num["min_zams_he"].first=minzams;
        params_num["max_zams_he"].first=maxzams;

        if (check_mzams_he())
            return EXIT_SUCCESS;
        else
            svlog.critical("Paramters check failed while settin min_zams_he, max_zams_he",__FILE__,__LINE__,sevnstd::params_error());

        return EXIT_FAILURE;

    }
    /**
     * Set the paramter myself containing the PATH of the SEVN folder.
     * @return Exist Success if the __FILE__ contains the /include/ in the string, if /include/ is not found
     * a critical error is raised
     * @Notice this is based on the fact that the folder structure contains the folder include that contains
     * (considering all the possibile subfolders) the file param.h
     */
    inline int set_myself(){
        params_str["myself"].first=utilities::get_absolute_SEVN_path()+"/";
        return EXIT_SUCCESS;
    }
    /**
     * Set the tables path
     * @param tables path to the tables
     * @return EXIT_SUCCESS
     */
    inline int set_tables(std::string tables){
        params_str["tables"].first=tables;

        return EXIT_SUCCESS;
    }
    /**
     * Set the pureHE tables path
     * @param tablesHE  path to the pureHE tables
     * @return EXIT_SUCCESS
     */
    inline int set_tablesHE(std::string tablesHE){
        params_str["tables_HE"].first=tablesHE;

        return EXIT_SUCCESS;
    }


    /**
     * Create a formatted string storing the name, values and comment of each parameter.
     * @param use_default  if true the string will contain the default values independently of the current values.
     * (the current values are not overwritten).
     * @return A string containing the name, values and comment of each parameter.
     *
     * @note In the comment there are two automatic additions inside square brackets:
     *  - The first contains info about the parameter type [N] numeric, [L] literal;
     *  - the second contains info about the set property [S] directle settable, [PS] private settable, [NS] non settable at runtime.
     */
    std::string print(bool use_default=false){

        int col_width = 35;
        std::ostringstream out;
        std::string header;
        map_num  saved_params_num;
        map_str  saved_params_str;
        map_bool saved_params_bool;


        if (use_default){
            header="#DEFAULT PARAMS";
            saved_params_num  = params_num;
            saved_params_str  = params_str;
            saved_params_bool = params_bool;
            init();
        }
        else{
            header="#USED PARAMS";
        }

        //Header
        out<<header<<std::endl;

        std::string info_set;

        for(auto& element :  params_num){

            if (public_settable.find(element.first)!=public_settable.end())
                info_set="[N][S] ";
            else if (private_settable.find(element.first)!=private_settable.end())
                info_set="[N][PS] ";
            else
                info_set="[N][NS] ";

            out<<std::setw(col_width)<<std::left<<element.first+": "<<
               std::setw(col_width)<<element.second.first<<"//"+info_set+element.second.second<<std::endl;
        }

        for(auto& element :  params_str){
            if (public_settable.find(element.first)!=public_settable.end())
                info_set="[L][S] ";
            else if (private_settable.find(element.first)!=private_settable.end())
                info_set="[L][PS] ";
            else
                info_set="[L][NS] ";

            out<<std::setw(col_width)<<std::left<<element.first+": "<<
               std::setw(col_width)<<element.second.first<<"//"+info_set+element.second.second<<std::endl;
        }

        for(auto& element : params_bool){

            if (public_settable.find(element.first)!=public_settable.end())
                info_set="[B][S] ";
            else if (private_settable.find(element.first)!=private_settable.end())
                info_set="[B][PS] ";
            else
                info_set="[B][NS] ";

            std::string val = element.second.first==true ? "true" : "false";

            out<<std::setw(col_width)<<std::left<<element.first+": "<<
               std::setw(col_width)<<val<<"//"+info_set+element.second.second<<std::endl;
        }


        if (use_default){
            params_num  = saved_params_num;
            params_str  = saved_params_str;
            params_bool = saved_params_bool;
        }

        return out.str();
    }

    //NB, since it is a friend function it not a method of the class
    /**
     * Overload of the << operator, it just calls the print method
     *
     * @overloaded
     *
     * @note since it is a friend function it not a method of the class
     */
    friend inline std::ostream& operator<< (std::ostream &os,  SEVNpar &svpar){
        os<<svpar.print();
        return os;
    }

public:
    //TODO Call all the parameters in that way, so that it can help IDE
    ///Parameter function
    inline double ts_maximum_variation() { return get_num("ts_maximum_variation");}


protected:

    SevnLogging svlog; /*!< logging instance */
    static std::size_t instance_counter; /*!< Stores the number of active instance of this class */

private:


    map_num  params_num;   /*!< map containing the numer parameters   {name: {double value, comment}} */
    map_str  params_str;   /*!< map containing the literal parameters {name: {string value, comment}} */
    map_bool params_bool;  /*!< map containing the boolean parameters {name: {boolean value, comment}} */

    //set are faster than vector
    static std::set<std::string> public_settable; /*!< set containing the name of all the properties that can be set directly */
    static std::set<std::string> private_settable; /*!< set containing the name of all the properties that can be set directly */

    ///SET
    /**
     * Set the value of a given numeric parameter
     *
     * @param name name of the parameter.
     * @param val value to save in the parameter.
     * @return EXIT_SUCCESS of throw an error
     *
     * @throws sevnstd::params_error Thrown if:
     *  - the parameter is not present,
     *
     */
    inline int prv_set(std::string name, double val){

        if (!is_num_par(name))
            svlog.critical("Trying to set the num key "+name+". The key is not present",__FILE__,__LINE__,sevnstd::params_error());

        if (is_settable(name))
            params_num[name].first=val;
        else{
            svlog.warning("Trying to set the num key "+name+". The key is non settable ",__FILE__,__LINE__);
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    };
    /**
     * Set the value of a given literal parameter
     *
     * @overload
     */
    inline int prv_set(std::string name, std::string val){

        if (!is_str_par(name))
            svlog.critical("Trying to set the str key "+name+". The key is not present",__FILE__,__LINE__,sevnstd::params_error());

        if (is_settable(name))
            params_str[name].first=val;
        else{
            svlog.warning("Trying to set the str key "+name+". The key is non settable ",__FILE__,__LINE__);
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    }
    /**
     * Set the value of a given boolean parameter
     *
     * @overload
     */
    inline int prv_set(std::string name, bool val){

        if (!is_bool_par(name))
            svlog.critical("Trying to set the boolean key "+name+". The key is not present",__FILE__,__LINE__,sevnstd::params_error());

        if (is_settable(name))
            params_bool[name].first=val;
        else{
            svlog.warning("Trying to set the boolean key "+name+". The key is non settable ",__FILE__,__LINE__);
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    }
    /**
     * Basic method to set value from a string.
     *
     * @param name name of the parameter.
     * @param val value to save in the parameter.
     * @return EXIT_SUCCESS, EXIT_FAILURE or thrown an error.
     *  EXIT_FAILURE means that the paremeter is not in the list
     *
     *
     * @throws sevnstd::params_error Thrown if:
     *  - the parameter is not settable,
     *
     * @note: if the parameter is not in the list, the function return EXIT_FAILURE but it does not throw an error.
     *
     */
    inline int set_from_string(std::string name, std::string val){

        int ret;
        std::string svalue=val;
        //transform to lower letter
        //std::transform(svalue.begin(), svalue.end(), svalue.begin(),
        //               [](unsigned char c){ return std::tolower(c); });
        if (is_num_par(name)) {
            double value;
            value = utilities::s2n<double>(svalue, __FILE__, __LINE__);
            ret = prv_set(name, value); //EXIT_SUCCESS if the value is settable (publicly or private)
            if (ret==EXIT_FAILURE)
                svlog.critical("The num key "+name+" is not settable at runtime");
        } else if (is_str_par(name)) {
            ret = prv_set(name, svalue);  //EXIT_SUCCESS if the value is settable (publicly or private)
            if (ret==EXIT_FAILURE)
                svlog.critical("The str key "+name+" is not settable at runtime");
        } else if (is_bool_par(name)) {
            bool value;
            if (val=="true" || val=="t" || val=="yes" || val=="y") value=true;
            else if (val=="false" || val=="f" || val=="no" || val=="yn") value=false;
            else svlog.critical("A boolean parameter can be set only with a string: (t)rue, (y)es,"
                                "(f)alse, (n)o",__FILE__,__LINE__,sevnstd::params_error());
            ret = prv_set(name, value);
            if (ret==EXIT_FAILURE)
                svlog.critical("The bool key "+name+" is not settable at runtime");
        }
        else{
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    }

    //aux
    /**
     * Check if it is directly settable
     * @param name  name of the parameter
     * @return true or false
     */
    inline bool is_public_settable(std::string name){
        return  public_settable.find(name)!=public_settable.end();
    }
    /**
     * Check if it is  settable
     * @param name name of the parameter
     * @return true or false
     */
    inline bool is_private_settable(std::string name){
        return  private_settable.find(name)!=private_settable.end();
    }
    /**
     * Check if a given parameter is settable
     * @param name name of the parameter
     * @return true or false
     */
    inline bool is_settable(std::string name){
        return is_public_settable(name)||is_private_settable(name);
    }
    /**
     * Check if a given parameter is a numeric parameter.
     * @param name name of the parameter
     * @return true or false
     */
    inline bool is_num_par(std::string name){
        return params_num.find(name)!=params_num.end();
    }
    /**
     * Check if a given parameter is a literal parameter.
     * @param name name of the parameter
     * @return true or false
     */
    inline bool is_str_par(std::string name){
        return params_str.find(name)!=params_str.end();
    }
    /**
     * Check if a giben parameter is a boolean parameter.
     * @param name name of the parameter
     * @return true or false
     */
    inline bool is_bool_par(std::string name){
        return params_bool.find(name)!=params_bool.end();
    }
    /**
     * Check if a given parameter is a  parameter of the class.
     * @param name name of the parameter
     * @return true or false
     */
    inline bool is_par(std::string name){
        return is_num_par(name) || is_str_par(name) || is_bool_par(name);
    }

    ///Initialisation methods
    /**
     * Method to define the numeric parameters and their default value.
     * A parameter can be added as
     * params_num["name"]             = std::make_pair(value,documentation);
     * @return EXIT_SUCCESS
     */
    int default_value_num();
    /**
     * Method to define the literal parameters and their default value.
     * A parameter can be added as
     * params_str["name"]             = std::make_pair(value,documentation);
     * @return EXIT_SUCCESS
     */
    int default_value_str();
    /**
     * Method to define the boolean parameters and their default value.
     * A parameter can be added as
     * params_bool["name"]             = std::make_pair(value,documentation);
     * @return EXIT_SUCCESS
     */
    int default_value_bool();

    /**
     * Initialise all the parameters to default values
     * @return EXIT_SUCCESS or throw an error
     *
     * @throws sevnstd::params_error Thrown if the one of the parameters checks fail.
     */
    int init(){
        default_value_num();
        default_value_str();
        default_value_bool();
        if (!check())
            svlog.critical("Parameters check failed",__FILE__,__LINE__,sevnstd::params_error());
        return EXIT_SUCCESS;
    }


    ///CHECK
    /**
     * Method to check the values of the parameters
     * @return true if all the checks are true, false if at least one check is false or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks raise an error when it is not satisfied
     *
     */
    bool inline check(){
        bool ret=true;
        ret=ret && check_star();
        ret=ret &&check_ns();
        ret=ret &&check_timestep();
        ret=ret &&check_GW();
        ret=ret &&check_mzams();
        ret=ret &&check_z();
        ret=ret &&check_jtrack();
        ret=ret &&check_rlobe();
        ret=ret &&check_accretion();
        ret=ret &&check_sn();
        ret=ret &&check_winds();
        ret=ret &&check_ce();
        ret=ret &&check_hard();
        ret=ret && check_parameters();
        ret=ret &&check_option_mode();
        ret=ret &&check_log();
        ret=ret &&check_ev();
        ret=ret &&check_systems();
        return ret;
    }

    //TODO put together check_star, mzams z
    //Star
    /**
     * check the parameter(s) related to the star properties.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_star(){
        if (get_num("star_lambda")<0. and get_num("star_lambda")!=-1 and get_num("star_lambda")!=-11 and get_num("star_lambda")!=-12 and get_num("star_lambda")!=-13 and get_num("star_lambda")!=-2 and get_num("star_lambda")!=-21
                                                                                                                                                                                                        and get_num("star_lambda")!=-3 and get_num("star_lambda")!=-31
                                                                                                                                                                                                        and get_num("star_lambda")!=-4 and  get_num("star_lambda")!=-41
                                                                                                                                                                                                        and get_num("star_lambda")!=-5 and  get_num("star_lambda")!=-51)
            svlog.critical("star_lambda>0 or star_lambda=(-1,-11,-12,-13,-2,-21,-3,-31,-4,-41,,-5,-51), current value is "+n2s("star_lambda"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("star_lambda_pureHe")<=0)
            svlog.critical("star_lambda_pureHe>0, current value is "+n2s("star_lambda_pureHe"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("star_lambda_fth")<0. or get_num("star_lambda_fth")>1.)
            svlog.critical("0<=star_lambda_fth<=1, current value is "+n2s("star_lambda_fth"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("star_lambda")>0 and get_num("star_lambda_fth")!=0)
            svlog.warning("The parameter star_lambda_fth is larger than 0, but the parameter star_lambda is not -1. The star_lambda_fth will not be used!",
                          __FILE__,__LINE__);
        if (get_num("star_tshold_WR_envelope")<0 or get_num("star_tshold_WR_envelope")>1)
            svlog.critical("0<=star_tshold_WR_envelope<=1, current value is "+n2s("star_tshold_WR_envelope"),__FILE__,__LINE__,sevnstd::params_error());



        return true;
    }

    //NS
    bool inline check_ns(){
        if (get_num("ns_magnetic_tscale")<=0.)
            svlog.critical("ns_magnetic_tscale>0, current value is "+n2s("ns_magnetic_tscale"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("ns_magnetic_mscale")<=0.)
            svlog.critical("ns_magnetic_mscale>0, current value is "+n2s("ns_magnetic_mscale"),__FILE__,__LINE__,sevnstd::params_error());

        return true;
    }

    //SN
    bool inline check_sn(){
        if (get_num("sn_co_lower_sn")<=0)
            svlog.critical("sn_co_lower_sn>0, current value is "+n2s("sn_co_lower_sn"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_co_lower_ecsn")<=0)
            svlog.critical("sn_co_lower_ecsn>0, current value is "+n2s("sn_co_lower_ecsn"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_co_lower_ecsn")>get_num("sn_co_lower_sn"))
            svlog.critical("sn_co_lower_sn>sn_co_lower_ecsn, current value are "+n2s("sn_co_lower_sn")+" and " +n2s("sn_co_lower_ecsn"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_max_ns_mass")<=0)
            svlog.critical("sn_max_ns_mass>0, current value is "+n2s("sn_max_ns_mass"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_Mchandra")<=0.)
            svlog.critical("sn_Mchandra>0, current value is "+n2s("sn_Mchandra"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_compact_csi25_tshold")<0. and get_num("sn_compact_csi25_tshold")!=-1)
            svlog.critical("sn_compact_csi25_tshold>0 or -1, current value is "+n2s("sn_compact_csi25_tshold"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_Mremnant_average_NS")<=0.)
            svlog.critical("sn_Mremnant_average_NS>0, current value is "+n2s("sn_Mremnant_average_NS"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_Mremnant_std_NS")>get_num("sn_max_ns_mass"))
            svlog.critical("M_mean for NS mass distribution is larger than the maximum ns mass, this means that more than 50% of the sample will be rejected. Please use a different the input distribution.",__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_Mremnant_std_NS")<=0.)
            svlog.critical("sn_Mremnant_std_NS>0, current value is "+n2s("sn_Mremnant_std_NS"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_compact_fallback")<0. or get_num("sn_compact_fallback")>1.)
            svlog.critical("0.<sn_compact_fallback<1., current value is "+n2s("sn_compact_fallback"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_Mremnant_average")<=0. and get_num("sn_Mremnant_average")!=-1.)
            svlog.critical("sn_Mremnant_average>0 (or -1), current value is "+n2s("sn_Mremnant_average"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("sn_Mejected_average")<=0. and get_num("sn_Mremnant_average")!=-1.)
            svlog.critical("sn_Mejected_average>0 (or -1), current value is "+n2s("sn_Mejected_average"),__FILE__,__LINE__,sevnstd::params_error());

        return true;
    }
    //Timestep
    /**
     * check the parameter(s) related to the adaptive time step.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_timestep(){

        if (get_num("ts_min_dt")<0. && get_num("ts_min_dt")!=-1.0)
            svlog.critical("ts_min_dt>0 or ts_min_dt=-1, current value is "+n2s("ts_min_dt"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("ts_max_dt")<0. && get_num("ts_max_dt")!=-1.0)
            svlog.critical("ts_max_dt>0 or ts_max_dt=-1, current value is "+n2s("ts_max_dt"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("ts_max_dt")<get_num("ts_min_dt") and get_num("ts_max_dt")!=-1.0)
            svlog.critical("ts_max_dt>ts_min_dt, current values are ts_max_dt= "+n2s("ts_max_dt")+
                    " and ts_min_dt= "+n2s("ts_min_dt"),__FILE__,__LINE__,sevnstd::params_error());

        return true;
    };
    //Mass
    /**
     * check the parameter(s) containing the limit of the ZAMS mass from the tables.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_mzams(){

        const double MIN_ZAMS=get_num("min_zams");
        const double MAX_ZAMS=get_num("max_zams");

        if (MIN_ZAMS<0)
            svlog.critical("MIN_ZAMS>0, current value is "+ n2s(MIN_ZAMS),__FILE__,__LINE__,sevnstd::params_error());
        if (MAX_ZAMS<0)
            svlog.critical("MAX_ZAMS>0, current value is "+ n2s(MAX_ZAMS),__FILE__,__LINE__,sevnstd::params_error());
        if (MAX_ZAMS<MIN_ZAMS)
            svlog.critical("MAX_ZAMS>MIN_ZAMS, current values are: MAX_ZAMS= "+ n2s(MAX_ZAMS)
                    +", MIN_ZAMS= "+n2s(MIN_ZAMS),__FILE__,__LINE__,sevnstd::params_error());

        return true;
    }
    //Mass HE
    /**
     * check the parameter(s) containing the limit of the ZAMS mass from the pureHE tables.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_mzams_he(){

        const double MIN_ZAMS=get_num("min_zams_he");
        const double MAX_ZAMS=get_num("max_zams_he");

        if (MIN_ZAMS<0)
            svlog.critical("MIN_ZAMS_HE>0, current value is "+ n2s(MIN_ZAMS),__FILE__,__LINE__,sevnstd::params_error());
        if (MAX_ZAMS<0)
            svlog.critical("MAX_ZAMS_HE>0, current value is "+ n2s(MAX_ZAMS),__FILE__,__LINE__,sevnstd::params_error());
        if (MAX_ZAMS<MIN_ZAMS)
            svlog.critical("MAX_ZAMS_HE>MIN_ZAMS_HE, current values are: MAX_ZAMS_HE= "+ n2s(MAX_ZAMS)
                           +", MIN_ZAMS_HE= "+n2s(MIN_ZAMS),__FILE__,__LINE__,sevnstd::params_error());

        return true;
    }
    //Z
    /**
     * check the parameter(s) containing the limit of Z   from the tables.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_z(){

        const double MINZ=get_num("min_z");
        const double MAXZ=get_num("max_z");

        if (MINZ<0)
            svlog.critical("MIN_Z>0, current value is "+ n2s(MINZ),__FILE__,__LINE__,sevnstd::params_error());
        if (MAXZ<0)
            svlog.critical("MAX_Z>0, current value is "+ n2s(MAXZ),__FILE__,__LINE__,sevnstd::params_error());
        if (MAXZ<MINZ)
            svlog.critical("MAX_Z>MIN_Z, current values are: MAX_Z= "+ n2s(MAXZ)
                           +", MIN_Z= "+n2s(MINZ),__FILE__,__LINE__,sevnstd::params_error());

        return true;

    }
    //Z
    /**
     * check the parameter(s) containing the limit of Z   from the tables.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_z_he(){

        const double MINZ=get_num("min_z_he");
        const double MAXZ=get_num("max_z_he");

        if (MINZ<0)
            svlog.critical("MIN_Z_HE>0, current value is "+ n2s(MINZ),__FILE__,__LINE__,sevnstd::params_error());
        if (MAXZ<0)
            svlog.critical("MAX_Z_HE>0, current value is "+ n2s(MAXZ),__FILE__,__LINE__,sevnstd::params_error());
        if (MAXZ<MINZ)
            svlog.critical("MAX_Z_HE>MIN_Z_HE, current values are: MAX_Z_HE= "+ n2s(MAXZ)
                           +", MIN_Z_HE= "+n2s(MINZ),__FILE__,__LINE__,sevnstd::params_error());

        return true;

    }
    //jtrack
    /**
     * check the parameter(s) related the tracks jump
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_jtrack(){


        if (get_num("jtrack_tshold_dm_rel")<=0)
            svlog.critical("jtrack_tshold_dm_rel>0, current value is "+ n2s("jtrack_tshold_dm_rel"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("jtrack_h_err_rel_max")<=0)
            svlog.critical("jtrack_h_err_rel_max>0, current value is "+ n2s("jtrack_h_err_rel_max"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("jtrack_max_dm_factor")<0)
            svlog.critical("jtrack_max_dm_factor>=0, current value is "+ n2s("jtrack_max_dm_factor"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("jtrack_min_dm_factor")<0)
            svlog.critical("jtrack_min_dm_factor>=0, current value is "+ n2s("jtrack_min_dm_factor"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("jtrack_max_iteration")<=0)
            svlog.critical("jtrack_max_iteration>=0, current value is "+ n2s("jtrack_max_iteration"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("jtrack_max_dm_factor")<=get_num("jtrack_min_dm_factor"))
            svlog.critical("jtrack_max_dm_factor>jtrack_min_dm_factor, current values are: jtrack_max_dm_factor="+
            n2s("jtrack_max_dm_factor")+" jtrack_min_dm_factor="+n2s("jtrack_min_dm_factor"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("jtrack_dm_step")<=0)
            svlog.critical("jtrack_dm_step>0, currenta value is "+ n2s("jtrack_dm_step"),__FILE__,__LINE__,sevnstd::params_error());

        return true;

    }
    //RLO
    /**
     * Check the parameter(s) related to the Roche Lobe Overflow
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_rlobe(){
        if (get_num("rlo_f_mass_accreted")<0 || get_num("rlo_f_mass_accreted")>1)
            svlog.critical("0<rlo_f_mass_accreted<1, current value is "+n2s("rlo_f_mass_accreted"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("rlo_eps_nova")<0 || get_num("rlo_eps_nova")>1)
            svlog.critical("0<rlo_eps_nova<1, current value is "+n2s("rlo_eps_nova"),__FILE__,__LINE__,sevnstd::params_error());

        bool negative_gamma_ok = (get_num("rlo_gamma_angmom")==-1.0) || (get_num("rlo_gamma_angmom")==-2.0);
        if ( (get_num("rlo_gamma_angmom")<0) && (negative_gamma_ok==false) )
            svlog.critical("rlo_gamma_angmom>0 or rlo_gamma_angmom=-1 or rlo_gamma_angmom=-2, current value is "+n2s("rlo_gamma_angmom"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("rlo_gamma_angmom")>1)
            svlog.critical("rlo_gamma_angmom<=1, current value is "+n2s("rlo_gamma_angmom"),__FILE__,__LINE__,sevnstd::params_error());

        if (get_num("rlo_max_nuclearmt")<0)
            svlog.critical("rlo_max_nuclearmt>0, current value is "+n2s("rlo_max_nuclearmt"),__FILE__,__LINE__,sevnstd::params_error());


        return true;
    }
    //GW
    /**
     * check the parameter(s) related to the GW rad process.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_GW(){

        if (get_num("gw_tshold")<0)
            svlog.critical("gw_tshold>0, current value is "+n2s("gw_tshold"),__FILE__,__LINE__,sevnstd::params_error());
        return true;
    }
    //Wind
    /**
     * check the parameter(s) related to the winds  process.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_winds(){

        if (get_num("w_alpha")<0.)
            svlog.critical("w_alpha>0, current value is "+n2s("w_alpha"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("w_beta")<0.)
            svlog.critical("w_beta>0, current value is "+n2s("w_beta"),__FILE__,__LINE__,sevnstd::params_error());


        return true;
    }
    //General accretion
    bool inline check_accretion(){
        if (get_num("eddington_factor")<0)
            svlog.critical("eddington_factor>0, current value is "+n2s("eddington_factor"),__FILE__,__LINE__,sevnstd::params_error());

        return true;
    }
    //CE
    /**
     * check the parameter(s) related to the common envelope  process.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_ce(){


        if (get_num("ce_alpha")<0.)
            svlog.critical("ce_alpha>0, current value is "+n2s("ce_alpha"),__FILE__,__LINE__,sevnstd::params_error());


        if (get_num("ce_kce")!=-1){
            if (get_num("ce_kce")<0. || get_num("ce_kce")>1.)
                svlog.critical("ce_kce within [0,1] or -1, current value is "+n2s("ce_kce"),__FILE__,__LINE__,sevnstd::params_error());
        }
        if (get_num("ce_knce")!=-1){
            if (get_num("ce_knce")<0. || get_num("ce_kce")>1.)
                svlog.critical("ce_knce within [0,1] or -1, current value is "+n2s("ce_knce"),__FILE__,__LINE__,sevnstd::params_error());
        }


        return true;
    }
    //Hardening
    /**
     * check the paramters related to the hardening process.
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_hard(){

        if (get_num("hard_rhoc")<0.)
            svlog.critical("hard_rhoc>0, current value is "+n2s("hard_rhoc"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("hard_sigma")<0.)
            svlog.critical("hard_sigma>0, current value is "+n2s("hard_sigma"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("hard_xi")<0.)
            svlog.critical("hard_xi>0, current value is "+n2s("hard_xi"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("hard_kappa")<0.)
            svlog.critical("hard_kappa>0, current value is "+n2s("hard_kappa"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("hard_mass_average")<0.)
            svlog.critical("hard_mass_average>0, current value is "+n2s("hard_mass_average"),__FILE__,__LINE__,sevnstd::params_error());

        return true;
    }
    //Parameters
    bool inline check_parameters(){

        if( get_str("spin")!="list"  and !utilities::string_is_number<double>(get_str("spin"))){
            svlog.critical("spin has to be 'list' or a number, current value is "+get_str("spin"),__FILE__,__LINE__,sevnstd::params_error());
        }
        else if(get_str("spin")!="list"){
            std::string s=get_str("spin");
            double spin = utilities::s2n<double>(s,__FILE__,__LINE__);
            if(spin<0 or spin>1){
                svlog.critical("spin (Omega/Omega_crit) has to be in the range [0,1] current value is "+get_str("spin"),__FILE__,__LINE__,sevnstd::params_error());
            }
        }

        return true;
    }
    //Omode
    bool inline check_option_mode(){

        if (rlmap.count(get_str("rlmode"))==0)
            svlog.critical("Roche Lobe mode \"" + get_str("rlmode") + "\" not available", __FILE__, __LINE__,sevnstd::params_error());
        if (tidesmap.count(get_str("tmode"))==0)
            svlog.critical("Tides  mode \"" + get_str("tmode") + "\" not available", __FILE__, __LINE__,sevnstd::params_error());
        if (gwmap.count(get_str("gwmode"))==0)
            svlog.critical("Gravitational Wave mode \"" + get_str("gwmode") + "\" not available", __FILE__, __LINE__,sevnstd::params_error());
        if (mixmap.count(get_str("mixmode"))==0)
            svlog.critical("Mix mode \"" + get_str("mixmode") + "\" not available", __FILE__, __LINE__,sevnstd::params_error());
        if (snkmap.count(get_str("kmode"))==0)
            svlog.critical("SN kick mode \"" + get_str("kmode") + "\" not available", __FILE__, __LINE__,sevnstd::params_error());
        if (cemap.count(get_str("cemode"))==0)
            svlog.critical("CE mode \""+ get_str("cemode") + "\" not available", __FILE__, __LINE__,sevnstd::params_error());
        if (outputmap.count(get_str("omode"))==0)
            svlog.critical("Output mode \"" + get_str("omode") +"\" not available", __FILE__, __LINE__);
        #ifndef H5OUT
        if (outputmap.at(get_str("omode"))==OutputOption::_hdf5)
            svlog.critical("Output mode set to  \"" + get_str("omode") +"\"  but hdf5 is not available (you can make it available with the option -Dh5=ON in cmake)", __FILE__, __LINE__);
        #endif
        if (std::find(Lookup::xspinmodes.begin(), Lookup::xspinmodes.end(), get_str("xspinmode")) == Lookup::xspinmodes.end())
            svlog.critical("Black Hole Xspin mode \"" + get_str("xspinmode") + "\" not available", __FILE__, __LINE__,sevnstd::params_error());

        return true;

    }
    //Evolution
    /**
     * Check the parameter(s) related to the general evolution process
     * @return true or throw an error.
     *
     * @throws sevnstd::params_error Thrown if one of the checks fail.
     */
    bool inline check_ev(){

        if (get_num("ev_max_repetitions")<=0)
            svlog.critical("ev_max_repetitions>0, current value is "+n2s("ev_max_repetitions"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("ev_Nchunk")<=0)
            svlog.critical("ev_Nchunk>0, current value is "+n2s("ev_Nchunk"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("ev_naked_tshold")>1e-4)
            svlog.critical("ev_naked_tshold<=1e-4, current value is "+n2s("ev_naked_tshold"),__FILE__,__LINE__,sevnstd::params_error());

        if(get_str("ev_setwM")!="linear" and get_str("ev_setwM")!="rational" and get_str("ev_setwM")!="log"){
            svlog.critical("ev_setwM option "+get_str("ev_setwM") + " not available",__FILE__,__LINE__,sevnstd::params_error());
        }
        if(get_str("ev_setwM_log")!="linear" and get_str("ev_setwM_log")!="rational" and get_str("ev_setwM_log")!="log"){
            svlog.critical("ev_setwM_log option "+get_str("ev_setwM_log") + " not available",__FILE__,__LINE__,sevnstd::params_error());
        }
        if(get_str("ev_setwM_tphase")!="linear" and get_str("ev_setwM_tphase")!="rational" and get_str("ev_setwM_tphase")!="log"){
            svlog.critical("ev_setwM_tphase option "+get_str("ev_setwM_tphase") + " not available",__FILE__,__LINE__,sevnstd::params_error());
        }


        return true;
    }
    //LOG
    bool inline check_log(){

        std::string allowed_log[5]={"debug","info","warning","error","critical"};

        if( !utilities::isinlist(get_str("log_level"),std::begin(allowed_log),std::end(allowed_log)) )
            svlog.critical("Loglevel "+get_str("log_level")+" not allowed. The possible options are"
                                                           ": debug, info, warning, error",__FILE__,__LINE__,sevnstd::params_error());
        return true;
    }
    //Systems
    bool inline check_systems(){
        if (get_num("nthreads")<=0)
            svlog.critical("nthreads>0, current value is "+n2s("nthreads"),__FILE__,__LINE__,sevnstd::params_error());
        if (get_num("check_stalling_time")<=0)
            svlog.critical("check_stalling_time>0, current value is "+n2s("check_stalling_time"),__FILE__,__LINE__,sevnstd::params_error());

        return true;
    }

    /**
     * Transform the numeric parameter of a given name to a string.
     * @param name name of the parameter
     * @return  a string containing the value of the parameter.
     */
    inline std::string n2s(std::string name){
        return utilities::n2s(get_num(name),__FILE__,__LINE__);
    }
    /**
     * Transform a double to a string.
     * @param value to transform to a string
     * @return  a string containing the \p value.
     */
    inline std::string n2s(double val){
        return utilities::n2s(val,__FILE__,__LINE__);
    }


};






#endif //SEVN_PARAMS_H
