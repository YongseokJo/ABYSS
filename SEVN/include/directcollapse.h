//
// Created by Giuliano Iorio on 26/07/2022.
//

#ifndef SEVN_DIRECTCOLLAPSE_H
#define SEVN_DIRECTCOLLAPSE_H

#include <supernova.h>
#include <remnant.h>

class Star;
class PairInstability;

/**
 * Simple direct collapse model where Mrem=Mstar
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: estimated
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: from basic supernova
 * WD: from base supernova class
 * fallback: always set to 1.0
 * ppisn_correction: not applied
 * neutrino_correction: not applied
 */
class directcollapse : virtual public supernova, PisnOFF, NeutrinoMassLossOFF{

    directcollapse(Star *s = nullptr);

    static directcollapse _directcollapse;

    void explosion(Star *s) override;

    directcollapse* instance(Star *s){
        //TODO here we use new but never delete, check
        return (new directcollapse(s));
    }

    inline std::string name() const override { return "directcollapse"; }

};

#endif //SEVN_DIRECTCOLLAPSE_H
