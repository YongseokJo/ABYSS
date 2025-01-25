/*
 Created by iorio on 7/26/22.
 it includes special "ad-hoc" kick models
- Class CC15:  Kicks as in model CC15 by Giacobbo&Mapelli2018  https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2011G/abstract
- Class EC15CC265: Kicks as in models alpha in Giacobbo&Mapelli2018  https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2011G/abstract
- Class ECUS30: Kicks similar to Broekgaarden+21 (https://arxiv.org/pdf/2103.02608.pdf)
*/

#ifndef SEVN_SPECIALKICKS_H
#define SEVN_SPECIALKICKS_H

#include <kicks.h>
#include <hobbs.h>


/**
 * Kick as in model CC15 by Giacobbo&Mapelli2018  https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2011G/abstract
 * Both ECSN and CCSN receive a kick drawn from a Maxwellian with sigma= 15 km/s,
 * the BH kick is then corrected for the fallback fraction vkick_corrected =(1−ffb)vkick_maxwellian,
 */
class CC15 : public Hobbs {

public:

    CC15(bool reg = true) : gaussian15(0,sigma){
        if(reg) {
            Register(this, name());
        }
    }

    CC15* instance() {
        return (new CC15(false));
    }

    static CC15 _cc15;
    void _apply(Star *s) override;
    inline std::string name() override { return "cc15";}

protected:

    const double sigma=15.0; /*!< dispersion of the Maxwellian used to draw the natal kick*/
    std::normal_distribution<> gaussian15;
};


/**
 * Kick as in models alpha in Giacobbo&Mapelli2018  https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2011G/abstract
 * The ECSN receives a kick from a Maxwellian with sigma= 15 km/s the CCSN from a Maxwellian with sigma= 265 km/s,
 * the BH kick is then corrected for the fallback fraction vkick_corrected =(1−ffb)vkick_maxwellian,
 */
class EC15CC265 : public Hobbs{

    EC15CC265(bool reg = true) : gaussian_ecsn(0,sigma_ecsn), gaussian_ccsn(0,sigma_ccsn){
        if(reg) {
            Register(this, name());
        }
    }

    static EC15CC265 _ec15cc265;

    void _apply(Star *s) override;

    //when I use the regist function I use the non trivial constructor of the Hermite6th class
    EC15CC265* instance() {
        return (new EC15CC265(false));
    }

    inline std::string name() override { return "ec15cc265"; }

protected:

    const double sigma_ecsn=15.0; /*!< dispersion of the Maxwellian used to draw the natal kick of ECSN remnant*/
    const double sigma_ccsn=256.0; /*!< dispersion of the Maxwellian used to draw the natal kick of CCSN remnant*/
    //Make them static? Maybe
    std::normal_distribution<> gaussian_ecsn;
    std::normal_distribution<> gaussian_ccsn;
};


/**
 * Kicks similar to Broekgaarden+21 (https://arxiv.org/pdf/2103.02608.pdf)
 * The ECSN and Ultra stripped SN receives a kick from a Maxwellian with sigma= 15 km/s
 * Ultra stripped SN are nakedCO or star stripped by the HE layers while they are purHe stars (case BB),
 * we flag a SN as Ultra stripped if it is a nakedCO or a nakedHe with MHE-MCO<0.1 (see Broekgaarden+21)
 * the CCSN from a Maxwellian with sigma set by the parameter sn_kick_velocity_stdev.
 * The BH kick is NOT corrected for the fallback fraction.
 */
class ECUS30 : public Kicks{

    ECUS30(bool reg = true) {
        if(reg) {
            Register(this, name());
        }
    }

    static ECUS30 _ecus30;

    void _apply(Star *s) override;

    //when I use the regist function I use the non trivial constructor of the Hermite6th class
    ECUS30* instance() {
        return (new ECUS30(false));
    }

    inline std::string name() override { return "ecus30"; }

protected:

    const double sigma_ecsn=30.0; /*!< dispersion of the Maxwellian used to draw the natal kick of ECSN and Ultra stipped SN remnants*/
};


#endif //SEVN_SPECIALKICKS_H