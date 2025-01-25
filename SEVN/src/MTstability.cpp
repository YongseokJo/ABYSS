//
// Created by Giuliano Iorio on 17/11/2021.
//

#include <star.h>
#include <MTstability.h>
#include <utilities.h>
#include <binstar.h>

MTstability* MTstability::Instance(std::string const &name){
    auto it = GetStaticMap().find(name);
    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    return it == GetStaticMap().end() ? nullptr : (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

bool MT_Qcrit::mt_unstable(Star *donor, Star *accretor, _UNUSED Binstar *binary) const {
    return q(donor,accretor) > qcrit(donor,accretor);
}

double MT_Qcrit::q(Star *donor, Star *accretor) const {
    double Md = donor->getp_0(Mass::ID);
    double Ma = accretor->getp_0(Mass::ID);
    return Md/Ma;
}



double Qcrit_Hurley::qcrit(Star *donor, Star *accretor) const {


    int bse_type = donor->get_bse_phase_0();
    double qc;


    if(bse_type==0)
        //fully convective low mass star (Sec. 2.6.1 in Hurley+02)
        qc = 0.695;
    else if(bse_type==2)
        qc = 4;
    else if(bse_type==3 or bse_type==5 or bse_type==6){
        qc = qcrit_giant(donor,accretor);
    }
    else if (bse_type==8 or bse_type==9)
        //This value is used in Hurley+2002 for Naked Helium  stars
        qc = 0.784;
    else if (bse_type>=10 and bse_type<=12)
        qc = 0.628;
    else
        //This means that stars with phase 1,4,7,>13 have qc=3.
        //This is what used in BSE, but this is not reported in Hurley+2002
        qc = 3.0;

    // Note in SEVN1 the condition     else if (bse_type>=10 && bse_type<=12) is not present (the qc is 3), but then
    // inside the RLO it is hardcoded as 0.628 that is what is reported in Sec. 2.6.1 of Hurley+2002 for White Dwarf.

    return qc;
}

double Qcrit_Hurley::qcrit_giant(Star *donor,  _UNUSED Star *accretor) const {

        double MHE = donor->getp_0(MHE::ID);
        double Mass = donor->getp_0(Mass::ID);

        //qc = (1.67-zpars(7)+2.0*pow(mcore1/m1),5.))/2.13 Eq.57 in Hurely+02

        //Eq. 47 in Hurley+00
        double _Z = std::log10(donor->get_Z()/0.02);
        double zpar = 0.30406 + 0.0805*_Z + 0.0897*_Z*_Z + 0.0878*_Z*_Z*_Z + 0.0222*_Z*_Z*_Z*_Z;


        double Mcore_fraction = MHE/Mass;
        return  (1.67-zpar + 2.0*pow(Mcore_fraction,5))/2.13;

}

double Qcrit_Hurley_Webbink::qcrit_giant(Star *donor, _UNUSED  Star *accretor) const {
    double MHE = donor->getp_0(MHE::ID);
    double Mass = donor->getp_0(Mass::ID);

    return 0.362 + 1.0/(3.0*(1.0 - MHE/Mass));
}

double Qcrit_Hurley_Webbink_Shao::qcrit(Star *donor, Star *accretor) const {

    //In the Shao+21 paper (https://arxiv.org/pdf/2107.03565.pdf) there is a special treatment
    //for the accretion on a BH. Therefore if the accretor is not a BH we use the classical qc from Hurley.
    //In Shao+21 there is not really a constant value for qc, but there are a number of criteria that can we
    //have to take into account to estimate the stability of the mass transfer, then if the mass trasnfer is stable we
    //use a very large value for qcrit otherwise a very small one.
    //The criteria are:
    //q<qmin always stable (qmin ~ 1.5 - 2, we use 2)
    //q>qmax always unstable (qmax ~ 2.1 +0.8*MBH)
    //Between qmin and qmax the conditions for unstability are on the donor radius:
    //  - Rd< Rs~6.6-26.1*q +11.4q*q or
    //  - Rd>Ru ~ -173.8+45.5*Md -0.18*Md*Md
    // If at least one of the condition is true, the mass transfer is unstable.
    double qc;
    bool wasIBH = accretor->getp_0(RemnantType::ID)==Lookup::BH;

    if (wasIBH) {
        double current_q = donor->getp_0(Mass::ID) / accretor->getp_0(Mass::ID);
        double qmin = 2.0;
        double qmax = 2.1 + 0.8 * accretor->getp_0(Mass::ID);

        if (current_q < qmin) {
            qc = utilities::LARGE; //ALWAYS STABLE
        } else if (current_q > qmax) {
            qc = 0; //ALWAYS UNSTABLE
        } else {
            double Rs = 6.6 - 26.1 * current_q + 11.4 * current_q * current_q; //Eq.9 Shao+21
            double Ru = -173.8 + 45.5 * donor->getp_0(Mass::ID) - 0.18 * donor->getp_0(Mass::ID) * donor->getp_0(Mass::ID); //Eq.9 Shao+21

            if (donor->getp_0(Radius::ID) < Rs or donor->getp_0(Radius::ID) > Ru) {
                qc = 0; //ALWAYS UNSTABLE
            } else {
                qc = utilities::LARGE; //ALWAYS STABLE
            }
        }
    }
    else{
        qc = Qcrit_Hurley_Webbink::qcrit(donor,accretor);
    }

    return qc;
}

double Qcrit_COSMIC_Neijssel::qcrit(Star *donor, Star *accretor) const {

    //Use the COMPAS prescriptions taken from Neijssel+2020, section 2.3
    //This implementation is taken directly from COSMIC (https://github.com/COSMIC-PopSynth/COSMIC/blob/develop/cosmic/src/evolv2.f
    //In the documentation they report:
    //"We convert from radial response to qcrit for MS and HG,
    //which assumes conservative mass transfer,
    //Stable MT is always assumed for stripped stars,
    //Assume standard qcrit from BSE for kstar>=10"


    int bse_type = donor->get_bse_phase_0();
    double qc;


    if(bse_type<=1)
        qc = 1.717;
    else if(bse_type==2)
        qc = 3.825;
    else if(bse_type==3 or bse_type==5 or bse_type==6){
        //As in Hurley+02, second option (See _Webbink_giant implementation)
        qc = Qcrit_Hurley_Webbink::qcrit(donor, accretor);
    }
    else if (bse_type>=7 and bse_type<=9)
        //Mass transfer is always stable for pureHE.
        qc = 1000.0;
    else if (bse_type>=10)
        qc = 0.628;
    else
        //This means that stars with phase 4 have qc=3.
        qc = 3.0;

    return qc;

}

double Qcrit_COSMIC_Claeys::qcrit(Star *donor, Star *accretor) const {

    int bse_type = donor->get_bse_phase_0();
    double qc;
    //TODO Are BH and NS degenarete in the defintion of Clays+14?
    bool is_accretor_degenerate = accretor->amiremnant() or accretor->aminakedco();

    if(bse_type==0)
        //fully convective low mass star (Sec. 2.6.1 in Hurley+02)
        qc = is_accretor_degenerate? 1.0 : 0.695;
    else if(bse_type==1)
        qc = is_accretor_degenerate? 1.0 : 1.6;
    else if(bse_type==2)
        qc = is_accretor_degenerate? 4.7619 : 4.0;
    else if(bse_type==3 or bse_type==5 or bse_type==6){
        qc = is_accretor_degenerate? 1.15 : qcrit_giant(donor, accretor);
    }
    else if(bse_type==4)
        qc = 3.0;
    else if(bse_type==7)
        qc = 3.0;
    else if(bse_type==8)
        qc = is_accretor_degenerate? 4.762 : 4.0;
    else if(bse_type==9)
        qc = is_accretor_degenerate? 1.15 : 0.784;
    else if (bse_type>=10)
        qc = is_accretor_degenerate? 0.625 : 3.0;
    else
        qc = 0;

    return qc;
}

double Qcrit_StarTrack::qcrit(Star *donor, _UNUSED Star *accretor) const {
    int bse_type = donor->get_bse_phase_0();
    double qc;

    if(bse_type<7) //All Hydrogen stars
        qc = 3.0;
    else if(bse_type==7) //MS WR stars
        qc=1.7;
    else if(bse_type<9) //evolved WR stars
        qc=3.5;
    else             //all the degenerate objects
        qc=0.628;

    return qc;
}

double Qcrit_Radiative_Stable::qcrit(Star *donor, Star *accretor) const {

    //     * Just a wrapper of the Qcrit_Hurley_Webbink qcrit with
    //     * the difference that qc for MS (1), HG(2) and pureHE stars(7-9) is set to 1000

    int bse_type = donor->get_bse_phase_0();
    double qc;

    if(bse_type==0)
        //fully convective low mass star (Sec. 2.6.1 in Hurley+02)
        qc = 0.695;
    else if(bse_type<=2 or (bse_type>=7 and bse_type<=9))
        //Mass transfer always stable for MS and HG stars (radiative envelope)
        //and pureHe stars (case BB)
        qc = 1000.0;
    else if(bse_type==3 or bse_type==5 or bse_type==6){
        //As in Hurley+02, second option (See _Webbink_giant implementation)
        qc = Qcrit_Hurley_Webbink::qcrit(donor, accretor);
    }
    else if (bse_type>=10)
        qc = 0.628;
    else
        //This means that stars with phase 4 have qc=3.
        qc = 3.0;

    return qc;

}

double Qcrit_HRadiative_Stable::qcrit(Star *donor, Star *accretor) const {

    //     * Just a wrapper of the Qcrit_Hurley_Webbink qcrit with
    //     * the difference that qc for MS (1), HG(2) and pureHE stars(7-9) is set to 1000

    int bse_type = donor->get_bse_phase_0();
    double qc;

    if(bse_type==0)
        //fully convective low mass star (Sec. 2.6.1 in Hurley+02)
        qc = 0.695;
    else if(bse_type==1 or bse_type==2)
        qc = 1000.;
    else if(bse_type==3 or bse_type==5 or bse_type==6){
        qc = qcrit_giant(donor,accretor);
    }
    else if (bse_type==8 or bse_type==9)
        //This value is used in Hurley+2002 for Naked Helium  stars
        qc = 0.784;
    else if (bse_type>=10 and bse_type<=12)
        qc = 0.628;
    else
        //This means that stars with phase 1,4,7,>13 have qc=3.
        //This is what used in BSE, but this is not reported in Hurley+2002
        qc = 3.0;

    return qc;

}