//
// Created by mario on 13/11/18.
//

#ifndef SEVN_REVISED_IO_H
#define SEVN_REVISED_IO_H

#include <iostream>
#include <vector>
#include <set>
#include <dirent.h>
#include <algorithm>
#include <map>
#include <sevnlog.h>
#include <utilities.h>


#include <starparameter.h>
#include <fstream>
#include <numeric> //for iota function
#include <set>

#include <sevnlog.h>
using sevnstd::SevnLogging;

#include <lookup_and_phases.h>
using namespace Lookup;

#include <params.h>

#ifdef H5OUT
#include <H5out.h>

using sevnstd::H5out;
#else
#include <fstream>
#include <iomanip>
#endif

#include <sys/stat.h> //GI271119: moved here from region inside ifdef H5OUT 



class Star;


/*
NB:  Every IO constructor needs to call set_init_variable() at the beginning to properly set the following
- tablesloaded = 0; Otherwise it will never load the tables even if load is called
- nthreads = 1;
-  ntables = _Ntables;
 */
class IO {

public:

    IO(){
        set_init_variables();
    }



    IO(int argc, char **argv){

        set_init_variables();

        load(argc, argv);
        load_stars();

        if (!get_output_folder_name().empty())
            create_folder(get_output_folder_name());

        svlog.debug("Output folder name is "+get_output_folder_name());

    }



    ~IO(){
            if(liststars.is_open())
                liststars.close();
            if(failedstars.is_open())
                failedstars.close();
            if(outputfile.is_open())
                outputfile.close();
            if(logfile.is_open())
                logfile.close();
    }


    /* LOAD METHODS */
    /**
     * Load stars or binary from a list
     */
    void load_stars();
    /**
     * Load the simulation parameters and the tables
     * @param n number of elements in val
     * @param val list of strings
     */
    void load(int n, char** val);



    std::string get_logstring(){ return logstring;};
    inline void  reset_log(){logstring.clear();}

    template <typename T> void fill_matrix_test(std::vector< std::vector<T> > &matrix, std::ifstream &file){

        std::vector<T> vline;
        std::string line, value;
        std::istringstream stream;

        std::cout<<" Filling matrix "<<std::endl;

        for(;;){

            vline.erase(vline.begin(), vline.end());

            getline(file, line);
            if(file.eof() || line.length()==0 ) break; //break if end of file or if a blank line is included.
            stream.clear();
            stream.str(line);

            //If it is a comment do not read
            if(line[0]=='/' or line[0]=='#')
                continue;


            for(;;){ //keep reading the line till its last element

                if(stream >> value){
                    vline.push_back(utilities::s2n<T>(value, __FILE__, __LINE__));
                    stream.clear();
                }
                else
                    break;
            }

            matrix.push_back(vline);

        }

    }


    //tables[ID][i][j][k]:
    //ID = taken from enum lookup_names. It tells what kind of matrix you want to access
    // i = metallicity index
    // j = track number
    // k = time on the j-th track





    std::vector<std::vector<std::vector<std::vector<double>>>> tables,tables_HE;
    std::vector<double> Z, Z_HE;
    std::vector<std::vector<double>> allzams, allzams_HE;
    std::vector< std::vector<std::string> > STARS_MATRIX;

    std::random_device rd;

    //TODO Let the user select the output columns (not just what add to the default one)
    std::vector<std::string> printcolumns_star;   /*!< List of string containing the name of the stellar properties to print in output. It
 *is composed by a default lit of columns (IO::list_cols_star, defined in IO::load) plus some extra columns selected at runtime by the user with the -scol flag
 * see IO::load*/
    std::vector<std::string> printcolumns_binary; /*< As above, but for the properties of the binary system. Use -bcol flag to let the user add more columns
 * to the default defined in IO::list_cols_binary*/

    std::vector<size_t> printIDs_star;
    std::vector<size_t> printIDs_binary;


    size_t ntables;

    int         nthreads;
    std::string output_mode;     /**< Handle the file output type */
    std::string binput_mode;     /**< Handle the  formalism for the input for binaries */
    std::string winds_mode;      /**< Formalism for DA and DE in Winds accretion processes */
    std::string RL_mode;         /**< Formalism for DA and DE in RL  process */
    std::string tides_mode;      /**< Formalism for DA and DE in tidal processes */
    std::string GW_mode;         /**< Formalism for DA and DE in GW rad processes */
    std::string mix_mode;         /**< Formalism for  mixing processes options */
    std::string COLL_mode;         /**< Formalism for  mixing processes options */
    std::string HARD_mode;         /**< Formalism for  hardening processes options */
    std::string CIRC_mode;         /**< Formalism for  hardening processes options */
    std::string SNK_mode;         /**< Formalism for  SN kick processes options */
    std::string CE_mode;         /**< Formalism for  Common envelope  processes options */
    std::string SEVNpath;


    ///Parameters call
    //TODO Transform all the paramters to parameters call
    inline bool rseed_provided(){return svpar.get_bool("rseed");}


    SEVNpar svpar;


    //variable to check if look-up tables have already been loaded
    int tablesloaded;


    void print_output(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint=false);

    inline void print_formatted_output(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint=false){

        //This is not needed since the print_evolved_summary is called inside the call of binstar or star
        //print_evolved_summary(_name, _rseed, _ID);

        //TODO change the way we choose among different outputs. New class structure IO --> PRINT (virtual with instances) --> -- ASCII, HDF5, CSV
        // here we will just call PRINT.do()

        OutputOption out = outputmap.at(output_mode);

        if(out == OutputOption::_ascii)
            print_ascii(printmatrix, _name, _rseed, _ID, binaryprint);
        else if (out == OutputOption::_csv)
            print_csv(printmatrix, _name, _rseed, _ID, binaryprint);
        else if (out == OutputOption::_binary)
            print_bin(printmatrix, _name, _rseed, _ID, binaryprint);
#ifdef H5OUT
            else if (out == OutputOption::_hdf5)
        print_hdf5(printmatrix, _name, binaryprint);
#endif
        else
            svlog.critical("Output option not recognized: [" + output_mode +"] ", __FILE__, __LINE__,sevnstd::sevnio_error());

        return;
    }


    /**
     * Print the evolved summary in ascii format. The info are printed in the file opened as liststars in IO.h. The name of the file
     * is equal to evolved_<NTHREAD>.dat and it will be saved in the ouput folder chosen in input.
     * @param _name     name assigned to the star
     * @param _ID       ID assigen to the star
     */
    void print_evolved_summary(const std::string &_name, const unsigned long &_rseed, const size_t &_ID);

    /**
     * Print the  summary of the failed system in ascii format. The info are printed in the file opened as failedstars in IO.h. The name of the file
     * is equal to failed_<NTHREAD>.dat and it will be saved in the ouput folder chosen in input.
     * @param _name
     * @param _ID
     */
    void print_failed_summary(const std::string &_name, const unsigned long &_rseed, const size_t &_ID);

    void print_failed_initilisation_summary(const size_t &_ID);

    inline void print_params(std::string filename="used_params.svpar"){

        std::string fname = get_output_folder_name()+"/"+filename;

        std::ofstream output_file;
        output_file.open(fname);
        output_file<< svpar.print();
        output_file.close();
    }

    /**
     * Flush the string logstring to the logfile output and clear logstring
     * @param filename  Name of the file where to save the log
     */
     void print_log(std::string filename="logfile.dat");


    /**
     *
     * @param loginfo
     */
    void log_put(std::string& loginfo);


#ifdef DEBUG
    /**
    * Print step by step the info on the times after a single evolution step (see evolvestar in star.h).
    * with respect to the normal time output it prints also the step that in evolve are rejected due to a to fast evolution in Mass (see Timestep in property.cpp)
    * @param s      Pointer to the star object.
    * @return       Nothing, but the info are saved in timelog_x.dat in the output folder
    */
    void print_timelog(Star *s);
#endif


    void create_folder(const std::string &name) {
        #if defined(__linux__)
        int results = mkdir(name.c_str(), 0777);
        #else
        //The above code seems to not work properly for certain path name on  macOS
        int results = system(("mkdir -p "+name).c_str());
        #endif

        if(results){
            svlog.pdebug("Directory exist",__FILE__,__LINE__);
        }

        //GI: this was extremely dangerous because it removes everything inside the folder
        //if it already exists.
        /*
        if(results){
            svlog.pdebug("Directory exist",__FILE__,__LINE__);

            DIR *theFolder = opendir(name.c_str());
            struct dirent *next_file;

            while ( (next_file = readdir(theFolder)) != nullptr ){
                std::string path = name + "/" + next_file->d_name;
                remove(path.c_str());
            }

            closedir(theFolder);
        }
        */
    }

    //GET
    inline std::string get_output_folder_name() { return output_folder_name;}
    inline std::string get_SEVNpath() const { return SEVNpath;}

    //load
    /**
     * Utilities to load auxiliary data tables.
     * The files needs to be all inside the folder SEVN/auxiliary_data
     * @param name  name of the file
     * @return a 2D vector of double storing the given auxiliary table
     * @Note lines starting with # or / are skipped
     */
    std::vector<std::vector<double>> load_auxiliary_table(std::string name) const{

        std::ifstream table_file;
        std::vector<std::vector<double>> Matrix;

        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
            table_file.open(get_SEVNpath() + "auxiliary_data/" + name);
            //std::cout<<get_SEVNpath() + "auxiliary_data/" + name<<std::endl;
            if (table_file.is_open()) {
                fill_matrix(Matrix, table_file);
                table_file.close();
            } else
                svlog.critical("Raised an error reading auxiliary table file " + name +
                               ". Check if the file exist, if you specified the right path in the run script, and if you have reading permission.",
                               __FILE__, __LINE__,
                               sevnstd::sevnio_error());
        }

        return Matrix;
    }

protected:
    inline void set_init_variables(){

        tablesloaded = 0;
        ntables = _Ntables;

    };




private:

#ifdef H5OUT
        static H5out h5; //used by each thread to print the HDF5 dataset
    #pragma omp threadprivate(IO::h5)
#endif
    static std::ofstream liststars;
    static std::ofstream failedstars;
    static std::ofstream failedinits;
    static std::ofstream outputfile; //GI 81119: File where to save the (non h5) output.
    static std::ofstream logfile; //File to output the log
    static std::string logstring; //File to store info to output in log, declared threadprivate

#ifdef _OPENMP
#pragma omp threadprivate(IO::liststars)
#pragma omp threadprivate(IO::failedstars)
#pragma omp threadprivate(IO::failedinits)
#pragma omp threadprivate(IO::outputfile)
#pragma omp threadprivate(IO::logfile)
#pragma omp threadprivate(IO::logstring)
#endif

#ifdef DEBUG
    static std::ofstream timelog;
#pragma omp threadprivate(IO::timelog)
#endif


    bool load_called = false;  /*!< Flag to check if load has been already called for this instance. The idea is that if we want to change parameters we have to instanciate another object*/


    static std::vector<std::string> labels_STARMATRIX; /*!< Vector that stores the label of the STARMATRIX column, it is set with set_STARMATRIX_labels, called by load_txt*/

    std::string output_folder_name="output";


    std::vector<std::vector<double>> printmatrix;

    SevnLogging svlog;


    //main
    std::string tables_dir, tables_dir_HE, list_file; //folder where to find all the look-up tables
    std::string list_cols_star;  /*!< Default list of stellar properties to print in output, defined in IO::load*/
    std::string list_cols_binary; /*!< Default list of binary properties to print in output, defined in IO::load*/


    //auxiliaryload
    std::ifstream in;

    std::vector<std::string> zstring, zstring_HE;

    /**
     * Transfer the name from a string containing the name separated by a : to a vector of string.
     * It checks that each name  in list cols is one of the name included in Property::PrintMap or in BinaryProperty::PrintMAP
     * @param list_cols  A string containing the property names separated with a :
     * @param printcolumns  Vector of sting to fill with names from list_cols.
     */
    void columns_to_print_star(std::string &list_cols, std::vector<std::string> &printcolumns);

    /**
     * Transfer the name from a string containing the name separated by a : to a vector of string.
     * It checks that each name  in list cols is one of the name included in  BinaryProperty::PrintMAP
     * @param list_cols  A string containing the property names separated with a :
     * @param printcolumns  Vector of sting to fill with names from list_cols.
     */
    void columns_to_print_binary(std::string &list_cols, std::vector<std::string> &printcolumns);

    void read(std::vector<std::vector<std::vector<std::vector<double>>>>& tables,
              const std::string& tables_dir, const std::vector<std::string>& zstring,  std::vector<double>& Z);

    void read_tables(){

        read(tables,tables_dir,zstring,Z);
        read(tables_HE,tables_dir_HE,zstring_HE,Z_HE);

    }

    void inspect_tables();
    //It loads all the available tables, at all metallicities
    void load_tables();
    //It loads all the available tables for pure-HE stars (if any), at all metallicities

    /**
     * Fill a matrix (2D vector) taking elements from a file
     * @tparam T matrix type
     * @param matrix 2D vector
     * @param file open file containg the data
     * @param reset if true clear the matrix before to fill it otherwise just append
     */
    template <typename T> void fill_matrix(std::vector< std::vector<T> > &matrix, std::ifstream &file, bool reset=false) const{

        std::vector<T> vline;
        std::string line, value;
        std::istringstream stream;

        if (reset){
            matrix.clear();
        }

        for(;;){

            vline.erase(vline.begin(), vline.end());

            getline(file, line);
            if(file.eof() || line.length()==0 ) break; //break if end of file or if a blank line is included.
            stream.clear();
            stream.str(line);

            //Skip comments
            if (line[0]=='#' or line[0]=='/')
                continue;


            for(;;){ //keep reading the line till its last element
                if(stream >> value and stream.str()[0]!='/' and stream.str()[0]!='#'){
                    vline.push_back(utilities::s2n<T>(value, __FILE__, __LINE__));
                    stream.clear();
                }
                else
                    break;
            }

            matrix.push_back(vline);

        }

    }


    void print_list_summary(std::ofstream& outstream, std::string basename, const std::string &_name, const unsigned long &_rseed, const size_t &_ID);


    //GI 81119: specific functions to be called by print_output depending on the output type
    //void print_evolved_summary(const std::string &_name, const size_t &_ID);
#ifdef H5OUT
    void print_hdf5(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const bool binaryprint);
#endif

    /**
     * Formatted ascii print  of values contained in the a 2D vector of doubles.
     * In addition to the value in printmatrix two columns are added with the name and the ID of the star.
     * @param filename name of the file where to write the data
     * @param printmatrix 2D vector of doubles to print
     * @param _name name of the star (it will be the second column in addition to the data in @p printmatrix)
     * @param _ID ID of the star (it will the first column in addition to the data in @p printmatrix)
     * @param binaryprint if true the system calling print is a binary, so the printmatrix contains the info of the two stars and of their orbit
     * @param separator  field delimiter for the output file
     * @param _w_id  field width of the id column
     * @param _w_name  field width of the name column
     * @param _w_header  field with of all the other columns
     * @param _precision  decimal precision
     * @note  The name of the columns shown in the header are taken from IO::printcolumns_star and IO::printcolumns_binary. The order of the columns name and of the columns in
     * printmatrix has to be the same, but it is not checked here (see IO::columns_to_print for other details on how this is managed in the class).
     */
    void print_formatted_ascii(const std::string &filename,std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint, const std::string &separator,
                                   const size_t &_w_id, const size_t &_w_name, const  size_t &_w_header, const size_t &_precision, const std::string &comment);

    void print_ascii(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint);
    void print_csv(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint);
    void print_bin(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint);


    void inspect_dirs(){
        inspect_dir(tables_dir,zstring,Z);
        inspect_dir(tables_dir_HE,zstring_HE,Z_HE);
    }

    void inspect_dir(const std::string &motherdir, std::vector<std::string>& _zstring, std::vector<double>& _Z,std::string inspectdir="") {

        std::string dir = motherdir + "/" + inspectdir;

        //std::cout<<" inspectig dir = "<<dir<<std::endl;

        if(dir.empty())
            svlog.critical("You are trying to inspect an empty directory: " + dir, __FILE__, __LINE__);

        dir += "/"; //safe instruction.. if the user forgets to add the / at the end of the directory string
        auto directory = opendir(dir.c_str()); //from dirent.h


        if(nullptr == directory)
            svlog.critical("Cannot open directory " + dir, __FILE__, __LINE__);

        dirent *gotobject = readdir(directory); //get the first object inside the directory


        while(gotobject != nullptr) { //inspect all the objects inside the current directory
            inspect_object(gotobject, motherdir,_zstring,_Z);
            gotobject = readdir(directory);
        }



        closedir(directory);

    }

    void inspect_object(dirent *object, const std::string &motherdir, std::vector<std::string>& _zstring, std::vector<double>& _Z) {


        if (object->d_name[0] == '.') { //exclude . and ..
            return;
        }

        if (object->d_type == DT_DIR) {//the object is a directory

            _zstring.emplace_back(object->d_name); //Vector containing the name of all the metallicity folders (GI why not directly casted from Z?)
            //it's a directory, so process it... it should be a metallicity dir
            _Z.push_back(utilities::dirname2n<double>(std::string(object->d_name), __FILE__, __LINE__)); //save the metallicity value... use the string_to_number function

            //_tables_file.resize(_Z.size());

            inspect_dir(motherdir, _zstring, _Z, std::string(object->d_name));//process a directory (a sub_directory)
            return; //done inspection in the current sub-folder
        }


        //for each metallicity we fill the "tables_file vector... it contains all the tables we have for each different metallicity"
        if (object->d_type == DT_REG) {

            std::string filename = motherdir + "/" + _zstring[_zstring.size()-1] + "/" + std::string(object->d_name);

            //GI: Do not read files that are not defined in filemap
            if(filemap.find(std::string(object->d_name)) == filemap.end() and filemap_optional.find(std::string(object->d_name)) == filemap_optional.end()){
                svlog.warning("The unexpected file '"+std::string(object->d_name)+"' has been found  in the directory "+
                                      motherdir + "/" + _zstring[_zstring.size()-1] + ". It will be not loaded. Check if this is a misspelled table."
                                      ,__FILE__,__LINE__);
                return;
            }

            if(filename.empty())
                svlog.critical("This filename should not be empty: " + filename, __FILE__, __LINE__);


            return; //return the name of the file
        }


        svlog.critical("This is nor a file or a directory: " + std::string(object->d_name), __FILE__, __LINE__);

    }


protected:



    /**
     * Check if the tables Z, ZHE, allzams, allzams_HE are sorted
     * @return
     */
    //TODO probably this can be made void (no need to return bool?)
    bool check_sorted(){

        //Z
        if (!std::is_sorted(Z.begin(),Z.end()))
            svlog.critical("Z table is not sorted",__FILE__,__LINE__,sevnstd::sevnio_error());
        //ZHE
        if (!std::is_sorted(Z_HE.begin(),Z_HE.end()))
            svlog.critical("Z_HE table is not sorted",__FILE__,__LINE__,sevnstd::sevnio_error());
        //Mass
        for (auto& zams_table : allzams){
            if (!std::is_sorted(zams_table.begin(),zams_table.end()))
                svlog.critical("zams table is not sorted",__FILE__,__LINE__,sevnstd::sevnio_error());
        }
        //Mass HE
        for (auto& zams_table : allzams_HE){
            if (!std::is_sorted(zams_table.begin(),zams_table.end()))
                svlog.critical("zams_HE table is not sorted",__FILE__,__LINE__,sevnstd::sevnio_error());
        }

        return true;

    }

    /**
     * Set the labels referred to the columns of the STARMATRIX table depending on input type
     * NOTICE: so far we assume that a list of input system contains all binaries or all stars not a mix
     */
    inline void set_STARMATRIX_labels(){
        //TODO We should have a flag to understand if we are dealing with a binary or not.
        //TODO 8 is hardcoded for now;
        unsigned int N_expected = 8;
        bool is_binary = STARS_MATRIX[0].size()>N_expected;
        int inchoice=inputmapbin.at(binput_mode).first;

        if (is_binary and inchoice==InputBinaryOption::_new){
            labels_STARMATRIX={"Mass_0","Z_0","spin_0", "SN_0", "Tstart_0",
                               "Mass_1","Z_1","spin_1", "SN_1", "Tstart_1",
                               "a","e","Tend", "Dtout"};
        }
        else if (is_binary and inchoice==InputBinaryOption::_legacy)
            labels_STARMATRIX={"Mass_0","Mass_1","Z_0","Z_1","spin_0","spin_1",
                               "a","e","Tend", "Tstart", "dt","SN_0","SN_1","Dtout"};
        else if (!is_binary){
            labels_STARMATRIX={"Mass","Z","spin", "SN", "Tstart", "Tend", "Dtout"};
        }
        else
            svlog.critical("Option unkown",__FILE__,__LINE__,sevnstd::sevnio_error());
    }

    //TODO Put a safe get starmatrix table? It check if is has been initialised (i.e. is not empty) or not

};





#endif //SEVN_REVISED_IO_H
