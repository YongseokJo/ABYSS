//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <supernova.h>
#include <star.h>
#include <deathmatrix.h>



/****************** Death matrix ********************/


/// DeathMatrix by Woosley+20 https://arxiv.org/pdf/2001.10492.pdf
std::vector<std::vector<double>> DeathMatrix::death_matrix;

DeathMatrix::DeathMatrix(Star *s) : supernova(s) {
    if(s == nullptr) {
        Register(this);
    } else{
        //Fake values needed for Unified kick
        Average_remnant = 1.;
        Average_ejected = 1.;
        //Check if we need to load the death matrix
        if(death_matrix.empty()){
            load_table(s);
        }
    }
}

void DeathMatrix::load_table(Star *s){
    std::vector<std::vector<double>> _temp_Matrix=s->load_auxiliary_table(auxiliary_table_name);
    utilities::transpose(death_matrix,_temp_Matrix);
}



void DeathMatrix::explosion(_UNUSED Star *s) {

    sn_type = Lookup::SNExplosionType::CoreCollapse;
    double MHE= get_preSN_Mass(MHE::ID);

    //From DeathMatrix Tab.2 in Woosley+20
    if (MHE<=preSN_MHE_min_NS){
        Mremnant=NS_min_mass;
    }
    else if (MHE>preSN_MHE_max_PISN){
        Mremnant=0.;
        pisn_correction=0.;
    } else{
        Mremnant = utilities::interpolate_1D(MHE,death_matrix[0],death_matrix[1],false);
        //If Mremnant is withing 5% of MHE consider it has a direct collapse.
        if ( (MHE-Mremnant)/MHE<0.05 ) Mremnant=get_preSN_Mass(Mass::ID);
    }
    //Rough fallback fraction estimateassuming a mproto 1.1
    //fallback_frac = std::max((Mremnant - 1.1),0.) / std::max((s->getp(Mass::ID) - 1.1),1.1);
    //TODO at the moment if we don't set fallback frac, in the hobbs kick it is se to 0, so no correction
    //is applied
    //The point is that this fallback is just used in the model hobbs and makes sense only for the rapid and delayed model
    //what should we do? Maybe limit this model just to the rapid and delayed sn model? Or use always the simple estimate above if
    //fallback is not set?
}

void DeathMatrix::CCexplosion(Star *s){
    explosion(s);
    if (Mremnant==0.) sn_type=Lookup::SNExplosionType::PISN;
    else sn_type=Lookup::SNExplosionType::CoreCollapse;
    set_remnant_type_after_explosion(s,Mremnant,sn_type);
}

void DeathMatrix::ECSN(Star *s) {
    explosion(s);
    remnant_type = Lookup::Remnants::NS_ECSN;
    //fallback_frac = 0.0; //Unset
    sn_type = Lookup::SNExplosionType::ElectronCapture;
    pisn_correction=1;
}

