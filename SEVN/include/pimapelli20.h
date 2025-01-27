/*
 Created by iorio on 7/25/22.
 it includes:
- Class PIMapelli20:  Pair instability model from Appendix A of Mapelli+20, based on the fit (Spera+17) of the Woosley+17 tables
*/

#ifndef SEVN_PIMAPELLI20_H
#define SEVN_PIMAPELLI20_H

#include <pairinstability.h>

/**
 * mapelli20 option
 *
 * Fit (based on the preSN MHE core mass) to the Woosley17 tables as reported in the Appendix of Mapelli+20 (see also Spera+17)
 * see:
 * - Mapelli+20: https://ui.adsabs.harvard.edu/abs/2020ApJ...888...76M/abstract
 * - Spera+17: https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.4739S/abstract
 * - Woosley+17: https://ui.adsabs.harvard.edu/abs/2017ApJ...836..244W/abstract
 *
 * Modification with respect to the standard
 * - apply_beforeSN: no modifications
 * - apply_afterSN: apply corrections based on finale core He mass, correct Mremnant and set sn_type
 *   Approximately: 32<MHE<64 PPISN SN, 64<MHE<135 PISN (no remnant) (the fit depends also on MHE/MTOT)
 *
 */
class PIMapelli20 : public PairInstability{

public:
    //Ctor
    explicit PIMapelli20(bool reg=true);

    //Static instance
    static PIMapelli20 _pimapelli20;

    //name
    inline std::string name() const override { return "mapelli20";}

    //return an heap allocated instance
    PIMapelli20* instance() override {
        return (new PIMapelli20(false));
    }

    //Main function to call
    /**
     * PISN correction based on the
     * Fit (based on the preSN MHE core mass) to the Woosley17 tables as reported in the Appendix of Mapelli+20 (see also Spera+17).
     * Approximately:
     * 32<MHE<64 PPISN SN
     * 64<MHE<135 PISN (no remnant)
     * @param s Pointer to the exploding star to which we want to apply the correction
     * @param mremnant Mass of the remnant before the pair isntability correction
     * @return an Instance of the struct PISNeturn with only two members:
     *          - Mremnant_after_pi, Mass of the remnant after the application of pi
     *          - SNType: SN explosion type an object of the enum Lookup::SNExplosionType (ElectronCapture, CoreCollapse, PPISN,PISN)
     */
    PISNreturn apply_afterSN(_UNUSED Star* s, _UNUSED double mremnant) const override;

protected:

    /**
     * PISN correction based on the
     * Fit (based on the preSN MHE core mass) to the Woosley17 tables as reported in the Appendix of Mapelli+20 (see also Spera+17).
     * Approximately:
     * 32<MHE<64 PPISN SN
     * 64<MHE<135 PISN (no remnant)
     * @param mass  Mass of the remnant before the pair isntability correction
     * @param s Pointer to the exploding star to which we want to apply the correction
     * @param sntype  reference to a variable of type Lookup::SNExplosionType, after the correction, this variable
     * stores the type of triggered supernova.
     * @return
     */
    static double pisn_correction(const double mass, Star *s, Lookup::SNExplosionType& sntype) ;

};

#endif //SEVN_PIMAPELLI20_H
