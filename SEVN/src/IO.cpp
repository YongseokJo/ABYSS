//
// Created by mario on 13/11/18.
//

#include <IO.h>
#include <algorithm>
#include <iomanip>
#include <property.h>
#include <BinaryProperty.h>
#include <sevnlog.h>
#include <star.h>
#include <utilities.h>
#include <params.h>
#include <sstream>
#ifdef H5OUT
H5out IO::h5;
#endif

std::ofstream IO::liststars;
std::ofstream IO::failedstars;
std::ofstream IO::failedinits;
std::ofstream IO::outputfile;
std::ofstream IO::logfile;
std::string IO::logstring;


std::vector<std::string> IO::labels_STARMATRIX;

#ifdef DEBUG
std::ofstream IO::timelog;
#endif

#define N2S(val) utilities::n2s(val,__FILE__,__LINE__)


void IO::load(int n, char **val) {

    if (load_called){
        svlog.critical("IO::load cannot be called more then once, if you want to set new parameters, please"
                       "define a new IO object",__FILE__,__LINE__,sevnstd::sevnio_error());
    }


    /***************************************************************
    * Load Parameters
    ****************************************************************/
    svpar.load(n,val);
    svpar.set_myself();

    //Process options
    SEVNpath = svpar.get_str("myself");
    //default tables
    if (svpar.get_str("tables").empty()){
        svpar.set_tables(SEVNpath+"tables/SEVNtracks_parsec_ov05_AGB");
    }
    if (svpar.get_str("tables_HE").empty()){
        svpar.set_tablesHE(SEVNpath+"tables/SEVNtracks_parsec_pureHe36");
    }

    winds_mode   = svpar.get_str("wmode");
    RL_mode      = svpar.get_str("rlmode");
    tides_mode   = svpar.get_str("tmode");
    GW_mode      = svpar.get_str("gwmode");
    mix_mode     = svpar.get_str("mixmode");
    COLL_mode    = svpar.get_str("collmode");
    CIRC_mode    = svpar.get_str("circmode");
    HARD_mode    = svpar.get_str("hardmode");
    SNK_mode     = svpar.get_str("kmode");
    CE_mode      = svpar.get_str("cemode");
    nthreads     = (int)svpar.get_num("nthreads");
    output_mode  = svpar.get_str("omode");
    tables_dir   = svpar.get_str("tables");
    list_cols_star += svpar.get_str("scol");
    list_cols_binary += svpar.get_str("bcol");
    tables_dir_HE = svpar.get_str("tables_HE");
    list_file = svpar.get_str("list");
    binput_mode = svpar.get_str("ibmode");
    output_folder_name = svpar.get_str("o");



    //set loglevel
    svlog.set_level(svpar.get_str("log_level"));
    /****************************************************************/


    /***************************************************************
    * Select output column to print
    ****************************************************************/
    columns_to_print_star(list_cols_star, printcolumns_star);
    columns_to_print_binary(list_cols_binary, printcolumns_binary);
    /****************************************************************/

    /***************************************************************
    * Load tables and stars
    ****************************************************************/
    load_tables();
    svlog.info(" Tables loaded");
    /****************************************************************/

    /***************************************************************
    * Check if tables are sorted
    ****************************************************************/
    check_sorted();
    /****************************************************************/

    //Set load_called to true
    load_called=true;

}

void IO::load_tables() {

    inspect_dirs();
    read_tables();
    inspect_tables();


    tablesloaded = 1;

}

/**
 * Fill  STARS_MATRIX (2D string matrix) reading the line of the file given in input (name of the file is stored in list_file)
 */
void IO::load_stars(){

    list_file = svpar.get_str("list");
    openfile(in, list_file);
    svlog.info("File opened");
    fill_matrix(STARS_MATRIX, in);
    svlog.info("Stars loaded:");
    /*
    for(size_t i = 0; i < STARS_MATRIX.size(); i++){
        for(size_t j = 0; j < STARS_MATRIX[i].size(); j++)
            std::cout<<STARS_MATRIX[i][j]<<"   ";
        std::cout<<std::endl;
    }
     */
    std::cout<<  STARS_MATRIX.size() <<" loaded stars!"<<std::endl;

    if (STARS_MATRIX.empty())
        svlog.critical("The input file " + list_file + " is empty (if this sound strange, check if the first row is accidentally"
                       " empty or check if the last row is not empty)",__FILE__,__LINE__,sevnstd::sevnio_error());


    if (STARS_MATRIX.empty())
        svlog.critical("The input file " + list_file + " is empty (if this sound strange, check if the first row is accidentally"
                                                       " empty or check if the last row is not empty)",__FILE__,__LINE__,sevnstd::sevnio_error());


    set_STARMATRIX_labels();

}

void IO::read(std::vector<std::vector<std::vector<std::vector<double>>>>& tables,
               const std::string& tables_dir, const std::vector<std::string>& zstring,  std::vector<double>& Z){


    std::vector<size_t> oldindex(Z.size());

#ifdef DEBUG
    utilities::print_vector(Z);
#endif

    iota(oldindex.begin(), oldindex.end(), 0); //GI: fill the array with a sequence of number starting from 0

    // sort indexes based on comparing values in Z
    sort(oldindex.begin(), oldindex.end(), [&](size_t i1, size_t i2) {return Z[i1] < Z[i2];}); //Gi: sort the index as a function of the metallicity
    sort(Z.begin(), Z.end());

    DEBUG_LOG("Z sorted");
#ifdef DEBUG
    utilities::print_vector(Z);
#endif



    const size_t nZ = Z.size();

    //now that I know the number of all available metallicities, resize all the matrices
    tables.resize(ntables);

    for(size_t i = 0; i < ntables; i++)
        tables[i].resize(nZ); //now I am ready to collect the matrices for all the metallicities



    for(size_t k = 0; k < nZ; k++){

        size_t i = oldindex[k]; //Gi: transform the unordered metallicity list to the ordered one.

        std::string Z_directory=tables_dir+"/"+zstring[i]+"/";

        //Cycle over required tables
        for (auto& pair : filemap){

            std::string filetoread=Z_directory+pair.first;
            //utilities::hardwait(filetoread,pair.first,pair.second);

            //std::cout<<" Search required file tables = "<<filetoread<<std::endl;
            try{
                openfile(in, filetoread); //it automatically closes and open a file (it could be input or output depending on the definition of "in")
            }
            catch(sevnstd::sevnio_error& e){
                svlog.critical("Cannot open or found required table "+pair.first
                +" in directory "+Z_directory,__FILE__,__LINE__,sevnstd::sevnio_error());
            }
            //std::cout<<" Found and opened "<<std::endl;

            //std::cout<<" read the tables "<<std::endl;
            fill_matrix(tables[pair.second][k], in); //k is the metallicity index here
        }


        //Cycle over optional tables
        for (auto& pair : filemap_optional){

            std::string filetoread=Z_directory+pair.first;
            //utilities::hardwait(filetoread,pair.first,pair.second);
            bool skip=false;

            //std::cout<<" Search optional file tables = "<<filetoread<<std::endl;
            try{
                openfile(in, filetoread); //it automatically closes and open a file (it could be input or output depending on the definition of "in")
            }
            catch(sevnstd::sevnio_error& e){
                skip=true;
            }
            //std::cout<<" Found and opened "<<std::endl;

            if (!skip) {
                fill_matrix(tables[pair.second][k], in); //k is the metallicity index here
            }

        }


    }
}

void IO::inspect_tables(){

    //Resize allzams
    allzams.resize(Z.size());

    //TODO This should become a method whit tables in input so that we can use tableHe and tables without hardcoding them

    for(size_t i = 0; i < Z.size(); i++) {
        for (size_t k = 0; k < tables[_Mass][i].size(); k++)
            allzams[i].push_back(tables[_Mass][i][k][0]);
    }
    //transform luminosity and radii tables in LOG... this is useful for the interpolation
    for(size_t i = 0; i < tables[_Lumi].size(); i++){ //cycle over all the tables at all metallicities
        for(size_t k = 0; k < tables[_Lumi][i].size(); k++){ //cycle over all the tracks in the table
            for(size_t j = 0; j < tables[_Lumi][i][k].size(); j++){ //cycle over all the elements in the single track
                tables[_Lumi][i][k][j] = log10(tables[_Lumi][i][k][j]);
                tables[_Radius][i][k][j] = log10(tables[_Radius][i][k][j]);

                //Inertia is an optional table, check if it has been loaded
                //TODO We shoukd write a function to check if the table has been loaded
                if (!tables[_Inertia][i].empty())
                    tables[_Inertia][i][k][j] = log10(tables[_Inertia][i][k][j]);
            }
        }
    }


    //Resize allzams
    allzams_HE.resize(Z_HE.size());
    for(size_t i = 0; i < Z_HE.size(); i++) {
        for (size_t k = 0; k < tables_HE[_Mass][i].size(); k++)
            allzams_HE[i].push_back(tables_HE[_Mass][i][k][0]);
    }
    //transform luminosity and radii tables in LOG... this is useful for the interpolation
    for(size_t i = 0; i < tables_HE[_Lumi].size(); i++){ //cycle over all the tables at all metallicities
        for(size_t k = 0; k < tables_HE[_Lumi][i].size(); k++){ //cycle over all the tracks in the table
            for(size_t j = 0; j < tables_HE[_Lumi][i][k].size(); j++){ //cycle over all the elements in the single track
                tables_HE[_Lumi][i][k][j] = log10(tables_HE[_Lumi][i][k][j]);
                tables_HE[_Radius][i][k][j] = log10(tables_HE[_Radius][i][k][j]);
                //Inertia is an optional table, check if it has been loaded
                //TODO We shoukd write a function to check if the table has been loaded
                if (!tables_HE[_Inertia][i].empty())
                    tables_HE[_Inertia][i][k][j] = log10(tables_HE[_Inertia][i][k][j]);
            }
        }
    }


    svpar.set_zlimit(Z);
    svpar.set_zlimit_HE(Z_HE);
    svpar.set_zams_limit(allzams);
    svpar.set_zams_limit_he(allzams_HE);

}

void IO::print_output(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint) {



    print_evolved_summary(_name, _rseed, _ID);

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


//TODO: change the _name and _ID to a pointer to the star so that we can call s->name, s->ID
/**
 * Print the evolved summary in ascii format. The info are printed in the file opened as liststars in IO.h. The name of the file
 * is equal to evolved.dat and it will be saved in the ouput folder chosen in input.
 * @param _name     name assigned to the star
 * @param _ID       ID assigen to the star
 */
void IO::print_list_summary(std::ofstream& outstream, std::string basename, const std::string &_name, const unsigned long &_rseed, const size_t &_ID){

    ///PRINT FORMAT PARAMETERS
    std::string separator="    ";
    const unsigned int w_id=7;
    const unsigned int w_short=15;
    const unsigned int w_long=17;



    ///PRINT HEADER
    if(!outstream.is_open()) {
        std::string fname = utilities::gen_filename(get_output_folder_name(),basename);
        //utilities::wait(fname);
        outstream.open(fname.c_str(), std::ios::out);

        //ID HEADER
        outstream<<std::setw(w_id)<<std::left<<"#ID"<<separator<<std::setw(w_long)<<"name";
        //PARAMTERS HEADER
        for (const auto& label : labels_STARMATRIX){
            outstream<<separator<<std::setw(w_short)<<label;
        }
        //SEED HEADER
        outstream<<separator<<std::left<<std::setw(w_long)<<"Seed";
        outstream<<std::endl;

    }

    ///PRINT DATA
    outstream<<std::setw(w_id)<<std::left<< _ID<<separator<<std::setw(w_long)<<_name;
    //Check if STAR_MATRIX contains the seed (in that case does not print it here, see below)
    int  STARS_MATRIX_dimension  = svpar.get_bool("rseed") ? STARS_MATRIX[_ID].size()-1  : STARS_MATRIX[_ID].size();
    for(int j = 0; j < STARS_MATRIX_dimension; j++){


        if ( (labels_STARMATRIX[j]=="Z_0" or labels_STARMATRIX[j]=="Z_1" or labels_STARMATRIX[j]=="Z")
            and svpar.get_str("Z")!="list"){
                outstream<<separator<<std::setw(w_short)<<svpar.get_str("Z");
        }
        else if ( (labels_STARMATRIX[j]=="spin_0" or labels_STARMATRIX[j]=="spin_1" or labels_STARMATRIX[j]=="spin")
             and svpar.get_str("spin")!="list"){
            outstream<<separator<<std::setw(w_short)<<svpar.get_str("spin");
        }
        else if ( (labels_STARMATRIX[j]=="SN_0" or labels_STARMATRIX[j]=="SN_1" or labels_STARMATRIX[j]=="SN")
            and svpar.get_str("snmode")!="list"){
                outstream<<separator<<std::setw(w_short)<<svpar.get_str("snmode");
        }
        else if ( labels_STARMATRIX[j]=="Dtout" and svpar.get_str("dtout")!="list"){
            outstream<<separator<<std::setw(w_short)<<svpar.get_str("dtout");
        }
        else if (labels_STARMATRIX[j]=="Tend" and svpar.get_str("tf")!="list"){
            outstream<<separator<<std::setw(w_short)<<svpar.get_str("tf");
        }
        else if ( (labels_STARMATRIX[j]=="Tstart" or labels_STARMATRIX[j]=="Tstart_0" or labels_STARMATRIX[j]=="Tstart_1") and svpar.get_str("tini")!="list"){
            outstream<<separator<<std::setw(w_short)<<svpar.get_str("tini");
        }
        else
            outstream<<separator<<std::setw(w_short)<<STARS_MATRIX[_ID][j];
    }
    //PRINT SEED
    outstream<<std::setw(w_long)<<_rseed;
    //SPACE FOR THE NEXT OUTPUT
    outstream<<std::endl;
}


//TODO: change the _name and _ID to a pointer to the star so that we can call s->name, s->ID
/**
 * Print the evolved summary in ascii format. The info are printed in the file opened as liststars in IO.h. The name of the file
 * is equal to evolved.dat and it will be saved in the ouput folder chosen in input.
 * @param _name     name assigned to the star
 * @param _ID       ID assigen to the star
 */
void IO::print_evolved_summary(const std::string &_name, const unsigned long &_rseed, const size_t &_ID){
    print_list_summary(liststars, "evolved.dat", _name, _rseed, _ID);
}

void IO::print_failed_summary(const std::string &_name, const unsigned long &_rseed, const size_t &_ID){
    print_list_summary(failedstars, "failed.dat", _name, _rseed, _ID);
}

void IO::print_failed_initilisation_summary(const size_t &_ID) {
    print_list_summary(failedinits, "failed_initialisation.dat", "NOTASSIGNED", 0, _ID);
}

/**
 * Print the output in hdf5 format with header
 */
#ifdef H5OUT
void IO::print_hdf5(std::vector<std::vector<double>> &printmatrix, const std::string &_name,  const bool binaryprint){


    std::string new_f5filename = utilities::gen_filename(get_output_folder_name(),"output.h5");
    h5.set_filename(new_f5filename);


    if(binaryprint) {

        std::vector<std::string> printcolumns = printcolumns_star;
        printcolumns.insert( printcolumns.end(), printcolumns_star.begin(), printcolumns_star.end() );
        printcolumns.insert( printcolumns.end(), printcolumns_binary.begin(), printcolumns_binary.end() );

        h5.print(printcolumns, printmatrix, _name); //header, data and name of the star

    }
    else
        h5.print(printcolumns_star, printmatrix, _name);

}
#endif


void IO::print_formatted_ascii(const std::string &filename,std::vector<std::vector<double>> &printmatrix, const std::string &_name, _UNUSED const unsigned long &_rseed, const size_t &_ID, const bool binaryprint, const std::string &separator,
        const size_t &_w_id, const size_t &_w_name, const  size_t &_w_header, const size_t &_precision, const std::string &comment="#"){


    //Set variable to handle the binaryprint
    std::string star_id[2]={"",""};  //suffix to the column names
    size_t ncomp=1; //number of stars
    if (binaryprint){
        star_id[0]="_0";
        star_id[1]="_1";
        ncomp=2;
    }

    /*******************************************************************************
     * Check if the file is open. If it is not open, open it and write the header.
     ******************************************************************************/
    if (!outputfile.is_open()) {
        std::string fname_out = utilities::gen_filename(get_output_folder_name(), filename);
        outputfile.open(fname_out.c_str(), std::ios::out);
        outputfile << std::setw(_w_id) << std::left << comment << "ID" << separator << std::setw(_w_name) << "name";
        for (size_t k=0; k<ncomp; k++){
            for (size_t i = 0; i < printcolumns_star.size(); i++)
                outputfile << separator << std::setw(_w_header) << printcolumns_star[i]+star_id[k];
        }
        if (binaryprint){
            for (size_t i = 0; i < printcolumns_binary.size(); i++)
                outputfile << separator << std::setw(_w_header) << printcolumns_binary[i];
        }
        outputfile << std::endl;
    }
    /******************************************************************************/

    /*******************************************************************************
    * Print data.
    ******************************************************************************/

    outputfile << std::scientific << std::setprecision(_precision); //set the format for the output

    /********** Print the data in output   *************/
    for (size_t i = 0; i < printmatrix.size(); i++) {
        outputfile << std::setw(_w_id) << _ID << separator << std::setw(_w_name) << _name;
        int col_index;
        for (size_t k = 0; k < ncomp; k++) {
            for (size_t j = 0; j < printcolumns_star.size(); j++) {
                col_index = j + k * printcolumns_star.size();
                if (printcolumns_star[j] == "Phase"){
                    //Below we add 6 to the precision because the Phase is an integer and to mantain the format, we have to account
                    //for the x. at the beginning and the e+xx at the end (2+4 chars). The _precision set the number of decimals
                    if (svpar.get_bool("io_literal_phases"))
                        outputfile << separator << std::setw(_w_header) << Lookup::literal((Lookup::Phases)(int)(printmatrix[i][col_index]));
                    else if (separator.size()==1)
                        outputfile << separator << (int)(printmatrix[i][col_index]);
                    else
                        outputfile << separator << std::setw(_precision + 6) << (int)(printmatrix[i][col_index]);
                }
                else if(printcolumns_star[j] == "RemnantType"){
                    if (svpar.get_bool("io_literal_phases"))
                        outputfile << separator  << std::setw(_w_header) << Lookup::literal((Lookup::Remnants)(int)(printmatrix[i][col_index]));
                    else if (separator.size()==1)
                        outputfile << separator << (int)(printmatrix[i][col_index]);
                    else
                        outputfile << separator << std::setw(_precision + 6) << (int)(printmatrix[i][col_index]);
                }
                else if(printcolumns_star[j] == "Event"){
                    if (svpar.get_bool("io_literal_phases"))
                        outputfile << separator << Lookup::literal((Lookup::EventsList)(int)(printmatrix[i][col_index]));
                    else if (separator.size()==1)
                        outputfile << separator << (int)(printmatrix[i][col_index]);
                    else
                        outputfile << separator << std::setw(_precision + 6) << (int)(printmatrix[i][col_index]);
                }
                else if (std::isnan(printmatrix[i][col_index]))
                    outputfile << separator << std::setw(_w_header) << printmatrix[i][col_index];
                else
                    outputfile << separator << std::setw(_precision) << printmatrix[i][col_index];
            }
        }

        if (binaryprint) {
            for (size_t j = 0; j < printcolumns_binary.size(); j++) {
                col_index = j + 2 * printcolumns_star.size();
                if(printcolumns_binary[j] == "BEvent"){
                    if (svpar.get_bool("io_literal_phases"))
                        outputfile << separator << Lookup::literal((Lookup::EventsList)(int)(printmatrix[i][col_index]));
                    else if (separator.size()==1)
                        outputfile << separator << (int)(printmatrix[i][col_index]);
                    else
                        outputfile << separator << std::setw(_precision + 6) << (int)(printmatrix[i][col_index]);
                }
                else if (std::isnan(printmatrix[i][col_index]))
                    outputfile << separator << std::setw(_w_header) << printmatrix[i][col_index];
                else
                    outputfile << separator << std::setw(_precision) << printmatrix[i][col_index];

            }
        }
        outputfile << std::endl;
    }
   /******************************************************************************/

}

/**
 * Print the output in ascii format with header
 */
void IO::print_ascii(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint){

    std::string separator="    ", filename="output.dat";
    int _w_id=6;
    int _w_name=10;
    int _precision=10;
    int _w_header=_precision + 6; //to take in account the extra characters x. and e+00

    print_formatted_ascii(filename,printmatrix, _name, _rseed, _ID, binaryprint, separator,_w_id,_w_name,_w_header,_precision);
}

/**
 * Print the output in csv format with header
 */
void IO::print_csv(std::vector<std::vector<double>> &printmatrix, const std::string &_name, const unsigned long &_rseed, const size_t &_ID, const bool binaryprint){

    std::string separator=",", filename="output.csv";
    int _w_id=0;
    int _w_name=0;
    int _precision=6;
    int _w_header= 0; //to take in account the extra characters x. and e+00
    std::string comment="";

    print_formatted_ascii(filename,printmatrix, _name, _rseed, _ID, binaryprint, separator,_w_id,_w_name,_w_header,_precision, comment);
}

/**
 * Print the output in binary format with header
 */
void IO::print_bin(_UNUSED std::vector<std::vector<double>> &printmatrix, _UNUSED const std::string &_name, _UNUSED const unsigned long &_rseed, _UNUSED const size_t &_ID, _UNUSED const bool binaryprint){

    svlog.critical("You are trying  to print the output in binary format, but it is not implemented yet: " , __FILE__, __LINE__);

}

void IO::columns_to_print_star(std::string &list_cols, std::vector<std::string> &printcolumns) {

    std::string delimiter = ":";


    size_t pos = 0;
    std::string token;
    //Check the position of delimiters and push back the words within the delimiters to printcolumns
    while ((pos = list_cols.find(delimiter)) != std::string::npos) {
        token = list_cols.substr(0, pos);
        //If token is not already inside printcolumns  add it (in this way we avoid to have a duplicate columns)
        if (std::find(printcolumns.begin(),printcolumns.end(),token)==printcolumns.end()){
            printcolumns.push_back(token);
        }
        list_cols.erase(0, pos + delimiter.length());
    }
    //Check last elements
    if (!list_cols.empty()){
        token = list_cols;
        if (std::find(printcolumns.begin(),printcolumns.end(),token)==printcolumns.end()){
            printcolumns.push_back(token);
        }
    }


    for(size_t i = 0; i < printcolumns.size(); i++)
        svlog.debug("Column to print = "+utilities::n2s(printcolumns[i],__FILE__,__LINE__));


    //Check if the names in printcolumns are present in Property::PrintMap or in Binary::PrintMAP
    for(size_t i = 0; i < printcolumns.size(); i++) {
        auto it = Property::PrintMap.find(printcolumns[i]);


        if (it != Property::PrintMap.end()) { //Does the column you want to print exists in the map?
            printIDs_star.push_back(it->second);
        }
        else {
                std::string sterr;
                sterr = "I cannot recognize the Star parameter you want to print in  output: [" + printcolumns[i] + "]\n";
                sterr += "Possible parameters are the following: \n";
                sterr += "-----> For single stars: \n";
                for (auto &x: Property::PrintMap)
                    sterr += "      " + x.first + ";  ";
                sterr += "]\n";
            svlog.critical(sterr, __FILE__, __LINE__);
            }
        }
    }



void IO::columns_to_print_binary(std::string &list_cols, std::vector<std::string> &printcolumns) {

    std::string delimiter = ":";


    size_t pos = 0;
    std::string token;
    //Check the position of delimiters and push back the words within the delimiters to printcolumns
    while ((pos = list_cols.find(delimiter)) != std::string::npos) {
        token = list_cols.substr(0, pos);
        //If token is not already inside printcolumns  add it (in this way we avoid to have a duplicate columns)
        if (std::find(printcolumns.begin(),printcolumns.end(),token)==printcolumns.end()){
            printcolumns.push_back(token);
        }
        list_cols.erase(0, pos + delimiter.length());
    }
    //Check last elements
    if (!list_cols.empty()){
        token = list_cols;
        if (std::find(printcolumns.begin(),printcolumns.end(),token)==printcolumns.end()){
            printcolumns.push_back(token);
        }
    }


    for(size_t i = 0; i < printcolumns.size(); i++)
        svlog.debug("Column to print = "+utilities::n2s(printcolumns[i],__FILE__,__LINE__));


    //Check if the names in printcolumns are present in Property::PrintMap or in Binary::PrintMAP
    for(size_t i = 0; i < printcolumns.size(); i++) {
        auto it = BinaryProperty::PrintMap.find(printcolumns[i]);

        if (it != BinaryProperty::PrintMap.end()) { //Does the column you want to print exists in the map?
            printIDs_binary.push_back(it->second);
        }
        else {
                std::string sterr;
                sterr = "I cannot recognize the Binary parameter you want to print in  output: [" + printcolumns[i] + "]\n";
                sterr += "Possible parameters are the following: \n";
                sterr += "-----> For binary stars: \n";
                for (auto &x: BinaryProperty::PrintMap)
                    sterr += "      " + x.first + ";  ";
                sterr += "]\n";

                svlog.critical(sterr, __FILE__, __LINE__);
            }
        }
    }



void IO::print_log(std::string filename){
    if (svpar.get_bool("io_logfile")){
        if(!logfile.is_open()) {
            std::string fname = utilities::gen_filename(get_output_folder_name(),filename);
            logfile.open(fname.c_str(), std::ios::out);
        }
        logfile<<logstring;
        //Flush all the data
        logfile.flush();
        //Clear the string
        logstring.clear();
    }
    return;
}

void IO::log_put(std::string& loginfo){

    if (svpar.get_bool("io_logfile")){
        logstring+=loginfo+"\n";
    }
    return;
}



#ifdef DEBUG
void IO::print_timelog(Star *s){

    std::string separator = "    ";

    if (!timelog.is_open()) {
        std::string fname_timelog = utilities::gen_filename(get_output_folder_name(),"timelog.dat");
        timelog.open(fname_timelog.c_str(), std::ios::out);
        timelog <<  "ID" << separator << "name"<<endl;
    }
    else{
        timelog <<  s->get_ID() <<separator<< s->get_name()<<separator<<s->getp(Localtime::ID)<<separator<<s->getp(Worldtime::ID);
        timelog << separator << s->repeatstep;
        FOR4 timelog << separator << s->ttimes[_i];
        timelog << std::endl;
    }
}
#endif




