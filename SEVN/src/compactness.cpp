//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <supernova.h>
#include <star.h>
#include <compactness.h>



/****************** Compactness ********************/

compactness::compactness(Star *s) : supernova(s) {
    if ( s == nullptr) {
        Register(this);
    }
    else{
        csi25_explosion_tshold = s->get_svpar_num("sn_compact_csi25_tshold");
        Average_Mremnant_NS    = s->get_svpar_num("sn_Mremnant_average_NS");
        Std_Mremnant_NS        = s->get_svpar_num("sn_Mremnant_std_NS");
        //fallback_frac          = s->get_svpar_num("sn_compact_fallback");  //Fallback fraction in Eq. 3 Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf)
        //TODO it is wrong to assumed that the fH parameter in Mapelli20 is the fallback as in the rapid and delayed models
        //It is a parameter that take into account if the Hydrogen layes are conserved until the very end
       //Therefore we do not set the fallback fraction
        auxiliary_table_name   = "xi25_explprobability_PS20.dat";

        //Use Average_Mremnant_NS as _default_Mremnant_average
        double _default_Mejected_average=10.45; //Same as delayed
        set_Average_for_Unified(s,Average_Mremnant_NS,_default_Mejected_average);

        //Load auxiliary table only if really needed
        if (csi25_explosion_tshold==-1 and csi25_vs_explosion_probability.empty())
            load_table(s);

    }

}

void compactness::explosion(_UNUSED Star *s) {

    double csi25 = csi25_mapelli20(get_preSN_Mass(MCO::ID));

    double mproto = get_preSN_Mass(MHE::ID);
    double final_mass = get_preSN_Mass(Mass::ID);

    ///Following  Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf)
    //EXPLOSION
    if(triggering_explosion(csi25)){
        double Mgau=generate_random_gau(s->get_svpar_num("sn_Mremnant_average_NS"),s->get_svpar_num("sn_Mremnant_std_NS"));
        Mremnant = std::max(1.1,std::min(Mgau,s->get_svpar_num("sn_max_ns_mass")));
        Mremnant = std::min(Mremnant,final_mass); // Set a limit to the final mass
        //fallback_frac = 0.; //Fallback for NS assumed 0.
        //TODO Is this right? what should we assume for fallback_frac?
        //We do not set and we leave the hobbs sn kick model to estimate it
    }
        //IMPLOSION
    else{
        const double& fH = s->get_svpar_num("sn_compact_fallback");
        Mremnant = mproto + fH*(final_mass-mproto);
        //Assume have always fallback_frack=1
        fallback_frac=1; //TODO, maybe it is better to assume that an implosion always has fallback_frac=1, see the TODO in the ctor
        //Yes it is better, because in this case the kick will be zero (if unified or hobbs are used) and the system just receive the Blauw kick do the
        //differene between the preSN mass and the remnant mass
        //if (final_mass==mproto) fallback_frac=1; //In this case we have not an envelope (e.g. naked He or naked CO star)
    }

}

double compactness::csi25_mapelli20(double MCO)  {
    double a=0.55, b=-1.1;
    return std::max(a+b/MCO,0.); //Mapelli+20, Eq. 2 (https://arxiv.org/pdf/1909.01371.pdf)
}

double compactness::csi25_mapelli20(Star *s) const {
    return csi25_mapelli20(s->getp(MCO::ID)); //Mapelli+20, Eq. 2 (https://arxiv.org/pdf/1909.01371.pdf)
}

std::vector<std::vector<double>>  compactness::csi25_vs_explosion_probability={};

void compactness::load_table(Star *s){
    std::vector<std::vector<double>> _temp_Matrix=s->load_auxiliary_table(auxiliary_table_name);
    utilities::transpose(csi25_vs_explosion_probability,_temp_Matrix);
}

bool compactness::triggering_explosion(double csi25){

    if (csi25_explosion_tshold==-1){
        //Explosion probability and given xi25
        double prob_explosion = utilities::interpolate_1D(csi25,csi25_vs_explosion_probability[0],csi25_vs_explosion_probability[1],true);
        //Random number between 0 and 1
        double random_shoot   = rand_unif_0_1(utilities::mtrand);
        //utilities::hardwait("Prob",csi25,prob_explosion,random_shoot,__FILE__,__LINE__);
        //If random number is lower than the Explosion probability trigger explosion
        return random_shoot<=prob_explosion;
    }
    else if(csi25_explosion_tshold>0){
        //If lower than threshold trigger explosion
        return csi25<csi25_explosion_tshold;
    }
    else{
        svlog.critical("csi25 is negative and not -1, this is not allowes",__FILE__,__LINE__,sevnstd::params_error());
    }

    return false;

}

double compactness::generate_random_gau(double mean, double std) {
    return normal_dist(utilities::mtrand)*std + mean;
}

/****************** ************* ********************/