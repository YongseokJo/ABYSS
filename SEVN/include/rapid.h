//
// Created by Giuliano Iorio on 26/07/2022.
//

#ifndef SEVN_RAPID_H
#define SEVN_RAPID_H

#include <supernova.h>
#include <remnant.h>

class Star;
class PairInstability;

/**
 * Rapid SN model by Fryer+12
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: estimated.
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: from base supernova class
 * WD: from base supernova class
 * fallback: estimated
 * ppisn_correction: applied
 * neutrino_correction: applied
 */
class rapid : virtual public supernova, PisnON{
public:
    rapid(Star *s = nullptr);

    static rapid _rapid;

    void explosion(Star *s) override;

    rapid* instance(Star *s){
        return (new rapid(s));
    }

    inline std::string name() const override { return "rapid"; }

};

/**
 * Rapid SN model by Fryer+12 with NS mass drawn from a Maxwellian
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * but always larger than 1.1 Msun
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: NS mass drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * WD: from base supernova class
 * fallback: estimated as in rapid
 * ppisn_correction: applied
 * neutrino_correction: applied
 */
class rapid_gauNS : public rapid, NSfromGau{
public:
    rapid_gauNS(Star *s = nullptr);

    static rapid_gauNS _rapid_gauns;

    void ECSN(Star *s) override {NSfromGau::ECSN(s);}
    void CCexplosion(Star *s) override {NSfromGau::CCexplosion(s);}

    rapid_gauNS* instance(Star *s){
        return (new rapid_gauNS(s));
    }

    inline std::string name() const override { return "rapid_gauNS"; }

};

#endif //SEVN_RAPID_H
