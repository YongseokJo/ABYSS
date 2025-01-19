//
// Created by Giuliano Iorio on 26/07/2022.
//

#ifndef SEVN_DELAYED_H
#define SEVN_DELAYED_H

#include <supernova.h>
#include <remnant.h>

class Star;
class PairInstability;

/**
 * Delayed SN model by Fryer+12
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: estimated.
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: from base supernova class
 * WD: from base supernova class
 * fallback: estimated
 * ppisn_correction: applied
 * neutrino_correction: applied
 */
class delayed : virtual public supernova,PisnON{
public:
    delayed(Star *s = nullptr);

    static delayed _delayed;

    void explosion(Star *s) override;

    delayed* instance(Star *s) {
        return (new delayed(s));
    }

    inline std::string name() const override { return "delayed"; }


};

/**
 * Delayed SN model by Fryer+12 with NS mass drawn from a Maxwellian
 * BH/NS: agnostic, based on the max NS mass
 * NS mass: drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * but always larger than 1.1 Msun
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: NS mass drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * fallback: estimated as in delayed
 * WD: from base supernova class
 * ppisn_correction: applied
 * neutrino_correction: applied
 */
class delayed_gauNS :  public  delayed, NSfromGau{
public:

    delayed_gauNS(Star *s = nullptr);

    static delayed_gauNS _delayed_gauns;

    void ECSN(Star *s) override {NSfromGau::ECSN(s);}
    void CCexplosion(Star *s) override {NSfromGau::CCexplosion(s);}

    delayed_gauNS* instance(Star *s) {
        return (new delayed_gauNS(s));
    }

    inline std::string name() const override { return "delayed_gauNS"; }
};


#endif //SEVN_DELAYED_H
