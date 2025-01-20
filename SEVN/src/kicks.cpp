//
// Created by spera on 16/06/20.
//


#include <kicks.h>
#include <utilities.h>
#include <star.h>
#include <supernova.h>


Kicks* Kicks::Instance(std::string const &name){
    auto it = GetStaticMap().find(name);
    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

void Kicks::check_and_correct_vkick(Star* s) {

    //Sanity check on the module of the velocity
    if(s->vkick[3]<0){
        svlog.critical("SN kick module is negative: "+utilities::n2s(s->vkick[3],__FILE__,__LINE__),
                       __FILE__,__LINE__, sevnstd::sanity_error(" "));
    }

    //Check if Mremant 0 (this happens after a PPSIN)
    //In this case set the kick to 0
    if(s->get_supernova()->get_Mremnant()==0){
        s->vkick[3] = s->vkick[2] = s->vkick[1] = s->vkick[0] = 0.0;
        return;
    }

    //Check if below Vmin and correct
    if (s->vkick[3]<s->get_svpar_num("sn_min_vkick")){
        double new_kick = s->get_svpar_num("sn_min_vkick");

        //if minimum vkick is 0 just set all the components to 0
        if (new_kick==0){
            s->vkick[3] = s->vkick[2] = s->vkick[1] = s->vkick[0] = 0.0;
        }
        //if minimum vkick>0 just but old vkick=0, set the new module equal to vkick and draw new isotropic direction for the components
        else if(s->vkick[3]==0){
            //Generate random teta angle
            double random = uniformRealDistribution(utilities::mtrand);
            double teta = 2.0 * M_PI * random;

            s->vkick[3] = new_kick; //Total velocity
            s->vkick[2] = 2.0 * new_kick * random - new_kick;  //Vz (random between -Vkick and Vkick
            double vel2d = sqrt(new_kick*new_kick - s->vkick[2]*s->vkick[2]);
            s->vkick[0] = vel2d * cos(teta);  // Vx
            s->vkick[1] = vel2d * sin(teta);  // Vy

        }
        //if minimum vkick>0 just and old vkick>0, set the new module equal to vkick and rescale all the components by a factor new_module/old_module
        else{
            double correction=new_kick/s->vkick[3];
            //Correct so that vkick=vmin
            for (auto& vkick : s->vkick)
                vkick*=correction;
        }
    }
}

void Kicks::apply(_UNUSED Star *s){
    //If star is empty after a SN kick always 0
    _apply(s);
    //Check and correct vkick
    check_and_correct_vkick(s);
    return;
}
