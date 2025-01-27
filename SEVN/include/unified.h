/*
 Created by iorio on 7/26/22.
 it includes:
- Class Unified: Kicks from Giacobbo & Mapelli 2020, https://ui.adsabs.harvard.edu/abs/2020ApJ...891..141G/abstract
*/
#ifndef SEVN_UNIFIED_H
#define SEVN_UNIFIED_H

#include <kicks.h>

/**
 * Kicks from Giacobbo & Mapelli 2020, https://ui.adsabs.harvard.edu/abs/2020ApJ...891..141G/abstract
 * Natal kicks are drawn from a Maxwellina, but then they are rescaled by a factor proportional to Mejecta/Mremnant.
 * Vkick propto Mejecta takes into accoun the SN asymmetries, while Vkick propto 1/MBH takes into account the
 * conservation of linear momentum
 */
class Unified : public Kicks{

    Unified(bool reg = true){
        if(reg) {
            Register(this, name());
        }
    }

    static Unified _unified;

    void _apply(Star *s) override;

    Unified* instance() {
        return (new Unified(false));
    }

    inline std::string name() override { return "unified"; }

};


#endif //SEVN_UNIFIED_H
