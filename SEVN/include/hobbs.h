/*
 Created by iorio on 7/26/22.
 it includes:
- Class HobbsPure:  Pure Hobbs kick, kick drawn from a Maxwellian with sigma=sn_kick_velocity_stdev.
- Class Hobbs: Same as HobbePure, but the kick are corrected for the SN fallback (larger fallback, lower velocities)
*/

#ifndef SEVN_HOBBS_H
#define SEVN_HOBBS_H

#include <kicks.h>

/**
 * Pure Hobbs kick, kick drawn from a Maxwellian with sigma=sn_kick_velocity_stdev.
 */
class HobbsPure : public Kicks{
    HobbsPure(bool reg = true){
        if(reg) {
            Register(this, name());
        }
    }

    static HobbsPure _hobbspure;

    void _apply(Star *s) override;

    //when I use the regist function I use the non trivial constructor of the Hermite6th class
    HobbsPure* instance() {
        return (new HobbsPure(false));
    }

    inline std::string name() override { return "hobbs_pure"; }
};

/**
 * Same as HobbePure, but the kick are corrected for the SN fallback (larger fallback, lower velocities)
 */
class Hobbs : virtual public Kicks{

public:
    Hobbs(bool reg = true){
        if(reg) {
            Register(this, name());
        }
    }

    static Hobbs _hobbs;

    void _apply(Star *s) override;

    //when I use the regist function I use the non trivial constructor of the Hermite6th class
    Hobbs* instance() {
        return (new Hobbs(false));
    }

    inline std::string name() override { return "hobbs"; }

protected:

    /**
     * Estimate the fallback f_b to correct the kick, i.e.
     * Vk_new = Vk_old *f_b, where f_b = (1-fallback)
     * The fallback is introduced in the rapid and delayed SN model as
     *
     * fallback = (Mremnant - Mproto) / (MpreSN - Mproto)
     *
     * where Mproto depends on the SN model, in SEVN Mproto = 1.1 for rapid and 1.15 for delayed
     *
     * However, the Sn model can set the fallback, therefore in this method we firstly check if the
     * fallback in the SN model is set (i.e. >=0), in this case we use the value set by the SN model
     * otherwise we use the equation above assuming Mproto=1.1 Msun.
     *
     * @param s Pointer to the exploding star
     * @return the fallback correction term fb
     */
     double estimate_fallback_correction(Star *s);

};


#endif //SEVN_HOBBS_H
