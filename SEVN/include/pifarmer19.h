/*
 Created by iorio on 08/02/22.
 it includes:
- Class PIFarmer19:  BH mass in the PPISN following Farmer19
*/

#ifndef SEVN_FARMER19_H
#define SEVN_FARMER19_H

#include <pairinstability.h>

/**
 * farmer19 option
 * Fitting equations from  Farmer+19: https://ui.adsabs.harvard.edu/abs/2019ApJ...887...53F/abstract
 *
 *  Notice: The Mbh derived in farmer have been obtained starting from pureHe model, therefore here
 * we are assuming the Hydrogen envelope is always easily removed after the first pulse.
 *
 * Modification with respect to the standard
 * - apply_beforeSN:  no modifications
 * - apply_afterSN:
 *          The mass of the BH in the PPISN is derived using the fitting equations
 *          A1 from Farmer19, it the minimum between the Mremnant in input and the one
 *          derived from Eq. A1 using the star MCO
 *
 *
 */
class PIFarmer19 : public  PairInstability{

public:
    //Ctor
    explicit PIFarmer19(bool reg=true);

    //Static instance
    static PIFarmer19 _pifarmer19;

    //name
    inline std::string name() const override { return "farmer19";}

    //return an heap allocated instance
    PIFarmer19* instance() override {
        return (new PIFarmer19(false));
    }

    /**
     * PISN/PPISN implementation from Farmer+19   https://ui.adsabs.harvard.edu/abs/2019ApJ...887...53F/abstract
     * PISN for MCO>60 Msun and MHE<135 Msun
     * PPISN for MCO>=38 Msun and MCO=60 Msun
     *
     * The final mass of the PPISN is the minimum between the mass coming from the SN prescriptions
     * and the BH mass in equation A1 of Farmer+19
     *
     * Notice1: We limit the MBH from Eq. A1 to 11 Msun, the minimum remnant mass shown in the Farmer+19 plot
     * Notice2: The Mbh derived in farmer have been obtained starting from pureHe model, therefore here
     * we are assuming the Hydrogen envelope is always easily removed after the first pulse.
     *
     * @param s
     * @param mremnant
     * @return
     */
    PISNreturn apply_afterSN(_UNUSED Star* s, _UNUSED double mremnant) const override;

};


#endif //SEVN_FARMER19_H
