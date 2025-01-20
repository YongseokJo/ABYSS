//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <supernova.h>
#include <star.h>
#include <directcollapse.h>


/****************** Direct collapse ********************/

directcollapse::directcollapse(Star *s) : supernova(s) {
    if(s == nullptr) {
        Register(this);
    }
    //Fake values needed for Unified kick, but the kick will be always zero given that we have always a direct collapse
    Average_remnant = 1.;
    Average_ejected = 1.;

}
void directcollapse::explosion(_UNUSED Star *s) {
    fallback_frac = 1;
    Mremnant= get_preSN_Mass(Mass::ID);
    sn_type = Lookup::SNExplosionType::CoreCollapse;
}

/****************** ************* ********************/