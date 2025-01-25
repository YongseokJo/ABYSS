//
// Created by Giuliano Iorio on 26/07/2022.
//

#ifndef SEVN_DEATHMATRIX_H
#define SEVN_DEATHMATRIX_H

#include <supernova.h>
#include <remnant.h>

class Star;
class PairInstability;

/**
 * SN model based on the Death Matrix by Woosley+20 (Tab.2, https://arxiv.org/pdf/2001.10492.pdf)
 * BH/NS: based on the max NS mass
 * NS mass: estimated
 * max NS mass: set to 2.3
 * ECSN: same as other NS
 * WD: from base supernova class
 * fallback: unset
 * ppisn_correction: not applied (already included in the model)
 * neutrino_correction: not applied (already included in the model)
*/
class DeathMatrix : virtual  public supernova, PisnOFF, NeutrinoMassLossOFF{
public:
    DeathMatrix(Star *s = nullptr);

    static DeathMatrix _deathmatrix;

    void explosion(Star *s) override;
    void ECSN(Star *s) override;
    void CCexplosion(Star *s) override;


    DeathMatrix* instance(Star *s){
        return (new DeathMatrix(s));
    }

    inline std::string name() const  override { return "deathmatrix";}

protected:

    //Some specific varaibles
    const std::string auxiliary_table_name="deathmatrix_Woosley+20.dat";  /**< name of the input file containing the death matrix values */
    const double NS_max_mass=2.30;   /**< Max NS mass (see Woosley+20) */
    const double NS_min_mass=1.24;   /**< Min NS mass (see Woosley+20) */
    const double preSN_MHE_min_NS=2.07;  /**< below this value of MHE set the mass of the remnant to NS_min_mass  */
    const double preSN_MHE_max_PISN=60.12;  /**< below this value of MHE the supernova leaves no remnant  */


    static std::vector<std::vector<double>> death_matrix;

    /**
     * Load the auxiliary table to estimate the Mrem given the preSN MHE (deathmatrix_Woosley+20.dat)
     * @param s Pointer to the star (from which we get the io object)
     */
    void load_table(Star *s);
};






#endif //SEVN_DEATHMATRIX_H
