//
// Created by Giuliano Iorio on 26/07/2022.
//

#ifndef SEVN_COMPACTNESS_H
#define SEVN_COMPACTNESS_H

#include <supernova.h>
#include <remnant.h>

class Star;
class PairInstability;

/**
 * Compactness Model, bases on Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf).
 * BH/NS: estimated
 * NS mass: drawn from a Gaussian with average and std set by sn_Mremnant_average_NS and sn_Mremnant_std_NS.
 * but always larger than 1.1 Msun
 * max NS mass: from SEVN option sn_max_ns_mass
 * ECSN: from base supernova class
 * WD: from base supernova class
 * fallback: unset for explosion/NS (estimate in Hobbs Snkick if needed), set to 1 for implosion/BH
 * ppisn_correction: applied
 * neutrino_correction: applied
 */
class compactness : virtual public supernova, PisnON{

    compactness(Star *s = nullptr);

    static compactness _compactness;

    void explosion(Star *s) override ;

    compactness* instance(Star *s){
        return (new compactness(s));
    }

    inline std::string name() const  override { return "compact";}

protected:

    /**
     * Estimate the compactenss from the properties at the onset of the collpase
     * from Eq. 2 in Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf).
     * @param MCO mass of the CO core
     * @return  Compactness parameter (see Mapelli+20)
     */
    static double csi25_mapelli20(double MCO);

    /**
     * Estimate the compactenss from the properties at the onset of the collpase
     * from Eq. 2 in Mapelli+20 (https://arxiv.org/pdf/1909.01371.pdf).
     * @param s  Pointer to star
     * @return  Compactness parameter (see Mapelli+20)
     */
    double csi25_mapelli20(Star *s) const;


    SevnLogging svlog;
    double csi25_explosion_tshold;
    double Average_Mremnant_NS;
    double Std_Mremnant_NS;
    std::string auxiliary_table_name;

    std::uniform_real_distribution<double> rand_unif_0_1{0.,1.};
    std::normal_distribution<double> normal_dist{0.,1.};
    inline double generate_random_gau(double mean, double std);

    static std::vector<std::vector<double>> csi25_vs_explosion_probability;

    /**
     * Estimate if the SN trigger an explosion or an implosion
     * @param csi25 Compacteness parameters csi25
     * @return true if an explosion is trigerred, false if an implosion is triggered
     * @Note the estimate method depends on the parameter csi25_explosion_tshold
     * if it is positive we just compare csi25 with csi25_explosion_tshold. If
     * csi25_explosion_tshold=-1, we estimate the likelihood of explosion at given csi25
     * (from Patton&Sukhbold20) and then we draw a random number between 0 aand 1, if it is lower
     * than per probability lkelihood a explosion is triggered.
     */
    bool triggering_explosion(double csi25);

    /**
     * Load the auxiliary table to esimate the explosion likelihood (from Patton&Sukhbold20),
     * @param s Pointer to the star (from which we get the io object)
     */
    void load_table(Star *s);


};


#endif //SEVN_COMPACTNESS_H
