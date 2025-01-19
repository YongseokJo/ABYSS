//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <nmllattimer89.h>

NMLLattimer89::NMLLattimer89(bool reg) : NeutrinoMassLoss(false) {
    if (reg){
        Register(this);
    }
}


/**
 * Estimate the neutrino mass loss correction following Lattimer&Yahil89. https://ui.adsabs.harvard.edu/abs/1989ApJ...340..426L/abstract
 * as reported in Eq. 4 in Giacobbo & Mapelli 2018 (https://arxiv.org/pdf/1805.11100.pdf).
 * In addition, the neutrino mass loss is limited to 0.5 Msun.
 * @param mremnant Mass of the remnant before the pair instability correction
 * @param s  Pointer to the exploding star to which we want to apply the correction
 * @return mremnant
 */
double NMLLattimer89::apply(const double mremnant, _UNUSED Star *s) const {
    //Lattimer & Yahil89 mass loss, valid for NSs...
    //Eq. 4 in Giacobbo & Mapelli, 2018 (https://arxiv.org/pdf/1805.11100.pdf)
    double mremnant_after = 6.6667*(sqrt(1.0+0.3*mremnant)-1.0);

    double Mass_neutrinos = mremnant - mremnant_after;
    return (Mass_neutrinos < 0.5 ? mremnant_after : mremnant-0.5);
}

