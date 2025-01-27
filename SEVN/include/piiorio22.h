/*
 Created by iorio on 7/26/22.
 it includes:
- Class PIIorio22:  Simple prescription presented in Iorio+22 (in prep.) based on the Farmer+19 results
- Class PIIorio22Limited: Same as PIIorio22, but the BH mass is limited using Eq. 2 in Farmer+19
*/

#ifndef SEVN_PIIORIO22_H
#define SEVN_PIIORIO22_H

#include <pairinstability.h>

/**
 * iorio22 option
 *
 * Simple prescription presented in Iorio+22 (in prep.) based on the Farmer+19 results
 * see:
 * - Farmer+19: https://ui.adsabs.harvard.edu/abs/2019ApJ...887...53F/abstract
 *
 *
 * Modification with respect to the standard
 * - apply_beforeSN: modify the preSN mass
 *         - MCO<38 Msun or MHE>135 Msun -> not modifications
 *         - MCO>=38 Msun and (MHE<=135 Msun)-> Mass=MHE=MCO
 * - apply_afterSN:
 *         - MCO<38 Msun or MHE>135 Msun, no mass correction ->Core Collapse
 *         - 38 Msun<=MCO<=60 Msun, no mass correction (but apply_beforeSN already set Mass=MHE=MCO) -> PPISN
 *         - MCO>60 Msun and (MHE<135 Msun), set Mremnant to 0 -> PISN
 *
 */
class PIIorio22 : public  PairInstability{

public:
    //Ctor
    explicit PIIorio22(bool reg=true);

    //Static instance
    static PIIorio22 _piiorio22;

    //name
    inline std::string name() const override { return "iorio22";}

    //return an heap allocated instance
    PIIorio22* instance() override {
        return (new PIIorio22(false));
    }

    utilities::MassContainer apply_beforeSN(_UNUSED Star* s) const override;

    PISNreturn apply_afterSN(_UNUSED Star* s, _UNUSED double mremnant) const override;

protected:

    /**
     * Apply a mass limit to the BH mass after ppisn
     * By default not mass limit is applied
     * @param Mremnant Remnant mass before the application of the mass limit
     * @param s Pointer to the exploding star
     * @return New value of the remnant mass after the application of the limits
     */
    virtual double ppisn_mass_limitation(_UNUSED double Mremnant, _UNUSED Star* s) const {return Mremnant;}


};


/**
 * iorio22_limited option
 *
 * Same as iorio22 but with an upper limit on the PPISN BH mass following Eq. 2 in Farmer+19
 * see:
 * - Farmer+19: https://ui.adsabs.harvard.edu/abs/2019ApJ...887...53F/abstract
 *
 * Modification with respect to the standard
 * - apply_beforeSN: modify the preSN mass
 *         - MCO<38 Msun or MHE>135 Msun -> not modifications
 *         - MCO>=38 Msun and (MHE<=135 Msun)-> Mass=MHE=MCO
 * - apply_afterSN:
 *         - MCO<38 Msun or MHE>135 Msun, no mass correction ->Core Collapse
 *         - 38 Msun<=MCO<=60 Msun, limit BH mass using the *MBH_max from Equation 2 in Farmer+19  -> PPISN
 *                                  *MBH_max=35.1 -3.9 log10(Z) -0.31 log10(Z)*log10(Z), where Z is the metallicity
 *         - MCO>60 Msun and (MHE<135 Msun), set Mremnant to 0 -> PISN
 *
 */
class PIIorio22Limited : public PIIorio22{

public:
    //Ctor
    explicit PIIorio22Limited(bool reg=true);

    //Static instance
    static PIIorio22Limited _piiorio22limited;

    //name
    inline std::string name() const override { return "iorio22_limited";}

    //return an heap allocated instance
    PIIorio22Limited* instance() override {
        return (new PIIorio22Limited(false));
    }

protected:

    /**
     * Apply a mass limit to the BH mass after ppisn
     * Limit the BH mass to the MBH_max obtained with Eq. 2 in Farmer+19:
     * MBH_max=35.1 -3.9 log10(Z) -0.31 log10(Z)*log10(Z), where Z is the metallicity
     * @param Mremnant Remnant mass before the application of the mass limit
     * @param s Pointer to the exploding star
     * @return New value of the remnant mass after the application of the limits
     */
     double ppisn_mass_limitation(_UNUSED double Mremnant, _UNUSED Star* s) const override;

};

#endif //SEVN_PIIORIO22_H
