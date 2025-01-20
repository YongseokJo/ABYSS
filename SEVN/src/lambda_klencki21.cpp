//
// Created by Giuliano Iorio on 07/12/2021.
//

#include <vector>
#include <lambda_klencki21.h>
#include <star.h>
#include <IO.h>

thread_local std::vector<std::vector<double>> Lambda_Klencki::table;
thread_local bool Lambda_Klencki::already_loaded=false;

Lambda_Klencki::Lambda_Klencki(const IO *io){
    if (!already_loaded){
        table = io->load_auxiliary_table("lambda_fit_Klencki21.dat");
        already_loaded = true;
    }
}

Lambda_Klencki::Lambda_Klencki(const Star *s){
    if (!already_loaded){
        table = s->load_auxiliary_table("lambda_fit_Klencki21.dat");
        already_loaded = true;
    }
}

std::vector<double> Lambda_Klencki::find_row(double Mzams,double Z){

    //std::cout<<" Z "<< Z << std::endl;
    Z= Z/Zsun;
    //std::cout<<" Z "<< Z << std::endl;

    Mzams =  table.back()[0] < Mzams ? table.back()[0] : Mzams;

    auto condition = [Mzams,Z](std::vector<double> vv){ return Mzams<=vv[0] and Z>=vv[1]; };
    auto result3 = std::find_if(begin(table), end(table), condition);


    if (result3==table.end()){
        svlog.critical("Mass-Z for Klencki lambda table not found",__FILE__,
                       __LINE__,sevnstd::sanity_error());
    }

    return *result3;
}

double Lambda_Klencki::operator()(const Star *s) {
    int Phase = int(s->getp(PhaseBSE::ID));

    if (Phase==7 or Phase==8 or Phase==9) //WR or naked helium
        return s->get_svpar_num("star_lambda_pureHe");
    return operator()(s->get_zams(),s->get_Z(),s->getp(Radius::ID));
}

double Lambda_Klencki::operator()(double Mzams, double Z, double R) {

    if (Mzams!=Mzams_cache or Z!=Z_cache){
        Mzams_cache=Mzams;
        Z_cache=Z;
        vector_cache= find_row(Mzams,Z);
    }

    return estimate_lambda(R, vector_cache);
}

double Lambda_Klencki::estimate_lambda(double R,  const std::vector<double>& coefficients) {

    const double Rtams = coefficients[3];
    const double R12   = coefficients[4];
    const double R23   = coefficients[5];
    const double Rmax  = coefficients[6];
    // Consider the minimum radius as Rtams and the maximum as Rmax
    double Reff = std::max(Rtams,std::min(R,Rmax));
    double a,b,c,d;

    if (Reff<R12){
        a=coefficients[7];
        b=coefficients[8];
        c=coefficients[9];
        d=coefficients[10];
    } else if (Reff<R23){
        a=coefficients[11];
        b=coefficients[12];
        c=coefficients[13];
        d=coefficients[14];
    } else{
        a=coefficients[15];
        b=coefficients[16];
        c=coefficients[17];
        d=coefficients[18];
    }

    double logReff = std::log10(Reff);
    double logLambda=a*logReff*logReff*logReff + b*logReff*logReff + c*logReff + d;
    double lambda = std::pow(10,logLambda);

    //Now handle one special cases:
    //R<Rtams, in this case we set the Binding energy to the one at the TAMS
    if (R<=Rtams){
        lambda = lambda*Rtams/R; //so that the Binding energy remaints constant for R>Rtams
    }
    // for R>Rmax, instead we do nothing since mantain constant lambda, in this case the
    // Binding energy decreases that is the grend shown by Klencki+21 for expanding stars in the red

    return lambda;
}


thread_local std::vector<double> Lambda_Klencki_interpolator::Mzams_list;
thread_local std::vector<std::vector<double>> Lambda_Klencki_interpolator::Z_list;
thread_local bool Lambda_Klencki_interpolator::already_loaded_lists=false;

Lambda_Klencki_interpolator::Lambda_Klencki_interpolator(const IO *io) :  Lambda_Klencki(io){
    if (!already_loaded_lists){
        fill_interpolators_lists();
        already_loaded_lists = true;
    }
}
Lambda_Klencki_interpolator::Lambda_Klencki_interpolator(const Star *s) :  Lambda_Klencki(s){
    if (!already_loaded_lists){
        fill_interpolators_lists();
        already_loaded_lists = true;
    }
}

void Lambda_Klencki_interpolator::fill_interpolators_lists() {

    double Mzams_old=-9999;
    std::vector<double> set_tmp;
    for (auto& v : table){
        if (v[0]!=Mzams_old){ //Here we have found a new zams
            Mzams_list.push_back(v[0]);
            if (!set_tmp.empty()){
                std::sort(set_tmp.begin(),set_tmp.end());
                Z_list.push_back(set_tmp);
                set_tmp.clear();
            }
            set_tmp.push_back(v[1]);
        } else{
            set_tmp.push_back(v[1]);
        }
        Mzams_old = v[0];
    }
    //Add the last pack
    std::sort(set_tmp.begin(),set_tmp.end());
    Z_list.push_back(set_tmp);
    set_tmp.clear();

    std::sort(Mzams_list.begin(),Mzams_list.end());


}

void Lambda_Klencki_interpolator::find_interpolators(double Mzams,  double Z) {

    const double Znorm = Z/Zsun;

    size_t Mzams_id [2];
    size_t Z_id [4];

    Mzams_id[0] = std::min(utilities::binary_search(&Mzams_list[0], 0, Mzams_list.size()-1, Mzams),Mzams_list.size()-2); //search the metallicity of the star (Z)
    Mzams_id[1] = Mzams_id[0] + 1; //search the metallicity of the star (Z)

    Mzams_interpolators[0] = Mzams_list.at(Mzams_id[0]);
    Mzams_interpolators[1] = Mzams_list.at(Mzams_id[1]);

    for(size_t i = 0, k = 0; i < 2; i++, k+=2){
        size_t dim = Z_list[Mzams_id[i]].size();
        Z_id[k] =  std::min(utilities::binary_search(&Z_list[Mzams_id[i]][0], 0, dim-1, Znorm), dim-2); //search for the ZAMS of the star (mzams)
        Z_id[k+1] = Z_id[k] + 1;
    }

    //In the list, the Z are normalised to Z=0.017, so find, the normalised Z and then multiply  for the normalisation factor
    Z_interpolators[0] = Zsun*Z_list.at(Mzams_id[0]).at(Z_id[0]);
    Z_interpolators[1] = Zsun*Z_list.at(Mzams_id[0]).at(Z_id[1]);
    Z_interpolators[2] = Zsun*Z_list.at(Mzams_id[1]).at(Z_id[2]);
    Z_interpolators[3] = Zsun*Z_list.at(Mzams_id[1]).at(Z_id[3]);

    Mzams = std::max(Mzams_interpolators[0],std::min(Mzams,Mzams_interpolators[1]));
    wM[0] = (Mzams_interpolators[1] - Mzams) / (Mzams_interpolators[1] - Mzams_interpolators[0]);
    wM[1] = (Mzams - Mzams_interpolators[0]) / (Mzams_interpolators[1] - Mzams_interpolators[0]);


    double Zt = std::max(Z_interpolators[0],std::min(Z,Z_interpolators[1]));
    wZ[0] = (Z_interpolators[1] - Zt) / (Z_interpolators[1] - Z_interpolators[0]);
    wZ[1] = (Zt - Z_interpolators[0]) / (Z_interpolators[1] - Z_interpolators[0]);

    Zt = std::max(Z_interpolators[2],std::min(Z,Z_interpolators[3]));
    wZ[2] = (Z_interpolators[3] - Zt) / (Z_interpolators[3] - Z_interpolators[2]);
    wZ[3] = (Zt - Z_interpolators[2]) / (Z_interpolators[3] - Z_interpolators[2]);
}

double Lambda_Klencki_interpolator::operator()(double Mzams, double Z, double R) {

    if (Mzams!=Mzams_cache or Z!=Z_cache){
        Mzams_cache=Mzams;
        Z_cache=Z;
        vector_cache.resize(4); //reset vector_cache

        find_interpolators(Mzams_cache,Z_cache);

        vector_cache[0]=find_row(Mzams_interpolators[0],Z_interpolators[0]);
        vector_cache[1]=find_row(Mzams_interpolators[0],Z_interpolators[1]);
        vector_cache[2]=find_row(Mzams_interpolators[1],Z_interpolators[2]);
        vector_cache[3]=find_row(Mzams_interpolators[1],Z_interpolators[3]);
    }

    //for (auto& v : vector_cache){
    //    std::cout<< "vc" << v[0] << " " << v[1] << std::endl;
    //}


    double M0Z0 = wZ[0]*estimate_lambda(R, vector_cache[0]);
    double M0Z1 = wZ[1]*estimate_lambda(R, vector_cache[1]);
    double M1Z2 = wZ[2]*estimate_lambda(R, vector_cache[2]);
    double M1Z3 = wZ[3]*estimate_lambda(R, vector_cache[3]);


    /*
    for (auto& v : wZ){
        std::cout<< "wZ " << v  << std::endl;
    }

    for (auto& v : wM){
        std::cout<< "wM " << v  << std::endl;
    }
    */


    return  wM[0]*(M0Z0+M0Z1) + wM[1]*(M1Z2+M1Z3);
}