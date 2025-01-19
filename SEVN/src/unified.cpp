//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <unified.h>
#include <star.h>
#include <supernova.h>

void Unified::_apply(Star *s) {

    ///Initial check
    if (s->get_supernova()->get_AverageEjected()<0)
        svlog.critical("Average Ejected Mass in Unified kick is negative. It is likely not initialised  in your chosen SN model ("
                       +s->get_supernova()->name()+").",__FILE__,__LINE__,sevnstd::sn_error());
    else if (s->get_supernova()->get_AverageRemnant()<0)
        svlog.critical("Average remnant Mass in Unified kick is negative. It is likely not initialised  in your chosen SN model ("
                       +s->get_supernova()->name()+").",__FILE__,__LINE__,sevnstd::sn_error());

    double unified_std = s -> get_svpar_num("sn_kick_velocity_stdev");

    // Neutrino mass loss is included in s->get_supernova()->get_Mejected() = M_star_final - Mremnant, i.e.  s->get_supernova()->get_Mejected() is ALWAYS != 0.
    // Note: for ECSNe, Mejected = MCO - Mremnant, i.e. kicks are != 0 because of neutrino mass loss (in principle, kicks != 0 should come from the ejection of tiny envelopes around ONe-WD)

    double ejected_mass = s->get_supernova()->get_fallback_frac() == 1 ? 0.0 : (s->get_supernova()->get_remnant_type() == Lookup::Remnants::NS_ECSN ? s->getp(MCO::ID) - s->get_supernova()->get_Mremnant() : s->get_supernova()->get_Mejected()); //do not consider neutrinos in case of direct collapse
    double correction = (ejected_mass/s->get_supernova()->get_AverageEjected()) * (s->get_supernova()->get_AverageRemnant()/s->get_supernova()->get_Mremnant());

    //Sanity check in the correction factor
    if (correction<0){
        svlog.critical("SN kicks correction factor for the unified model is negative: "+utilities::n2s(correction,__FILE__,__LINE__)+
                       ". This has been estimated with the following properties:"+
                       "\n fallback: "+utilities::n2s(s->get_supernova()->get_fallback_frac(),__FILE__,__LINE__)+
                       "\n Remnant type: "+utilities::n2s(s->get_supernova()->get_remnant_type(),__FILE__,__LINE__)+
                       "\n Mremnant: "+utilities::n2s(s->get_supernova()->get_Mremnant(),__FILE__,__LINE__) +
                       "\n MCO: "+utilities::n2s(s->getp(MCO::ID),__FILE__,__LINE__) +
                       "\n Mejecta: "+utilities::n2s(s->get_supernova()->get_Mejected(),__FILE__,__LINE__)+
                       "\n Finale ejected mass: "+utilities::n2s(ejected_mass,__FILE__,__LINE__)+
                       "\n Mejecta average: "+utilities::n2s(s->get_supernova()->get_AverageEjected(),__FILE__,__LINE__)+
                       "\n Mremnant average: "+utilities::n2s(s->get_supernova()->get_AverageRemnant(),__FILE__,__LINE__),
                       __FILE__,__LINE__, sevnstd::sanity_error(" "));

    }

    s->vkick[0] = draw_from_gaussian(unified_std);
    s->vkick[1] = draw_from_gaussian(unified_std);
    s->vkick[2] = draw_from_gaussian(unified_std);
    s->vkick[3] = sqrt( s->vkick[0]*s->vkick[0] + s->vkick[1]*s->vkick[1] + s->vkick[2]*s->vkick[2]);
    //Velocity kick before correction
    random_velocity_kick = s->vkick[3];
    //Correct
    for (auto& vkick : s->vkick){
        vkick*=correction;
    }
}