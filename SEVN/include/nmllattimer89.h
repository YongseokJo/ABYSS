/*
 Created by iorio on 7/25/22.
it includes:
- Class NMLLattimer89:  Mass loss correction following Lattimer&Yahil89. https://ui.adsabs.harvard.edu/abs/1989ApJ...340..426L/abstract
*/

#ifndef SEVN_NMLLATTIMER89_H
#define SEVN_NMLLATTIMER89_H

#include <neutrinomassloss.h>
#include <star.h>

/**
 * option lattimer89
 *
 * Apply the neutrino mass loss correction from Lattimer&Yahil89. https://ui.adsabs.harvard.edu/abs/1989ApJ...340..426L/abstract
 * as reported in Eq. 4 in Giacobbo & Mapelli 2018 (https://arxiv.org/pdf/1805.11100.pdf).
 * In addition, the neutrino mass loss is limited to 0.5 Msun.
 */
class NMLLattimer89 : public NeutrinoMassLoss{
public:
    //Ctor
    NMLLattimer89(bool reg=true);

    //name
    inline std::string name() const override { return "lattimer89";}

    //Static instance
    static NMLLattimer89 _nmllattimer89;

    //return an heap allocated instance
    NMLLattimer89* instance() override {
        return (new NMLLattimer89(false));
    }

    /**
     * Estimate the neutrino mass loss correction following Lattimer&Yahil89. https://ui.adsabs.harvard.edu/abs/1989ApJ...340..426L/abstract
     * as reported in Eq. 4 in Giacobbo & Mapelli 2018 (https://arxiv.org/pdf/1805.11100.pdf).
     * In addition, the neutrino mass loss is limited to 0.5 Msun.
     * @param mremnant Mass of the remnant before the pair instability correction
     * @param s  Pointer to the exploding star to which we want to apply the correction
     * @return mremnant
     */
    double apply(const double mremnant, _UNUSED Star *s) const override;
};


#endif //SEVN_NMLLATTIMER89_H
