//
// Created by Giuliano Iorio on 14/12/2021.
//

#include <lambda_nanjing.h>
#include <lookup_and_phases.h>
#include <star.h>

Lambda_Nanjing::coeff_status Lambda_Nanjing::set_coeff_Giant(double Mzams, double Z, double Radius) {

    /*
    bool check_phases= Phase_cache==Lookup::Phases::TerminalMainSequence or Phase_cache==Lookup::Phases::ShellHBurning;
    //Now check if Mzams and Z are the same of the cached value,
    //check also the phase, if the last phase_0 was not the only for which this function is for, we have to
    //estimate the coefficient even if Mzams and Z have not changed
    if(check_cached(Mzams,Z) and check_phases){
        return COEF_NOT_SET;
    }
    reset_cache(Mzams,Z,Phase);
    */
    reset();


    if (Z>ZpopI) {                     // Z>0.5 Zsun: popI
        if (Mzams < 1.5) {
            maxBG = { 2.5, 1.5 };
            if (Radius > 200) lambdaBG = { 0.05, 0.05 };
            else if (Radius > 2.7) lambdaBG = { 2.33 - (Radius * 9.18E-03), 1.12 - (Radius * 4.59E-03) };
            else {
                coeff_lambdab = {  8.35897, -18.89048, 10.47651, 0.99352, 0.0, 0.0 };
                coeff_lambdag = { 17.58328, -34.84355, 10.70536, 8.49042, 0.0, 0.0 };
            }
        }
        else if (Mzams < 2.5) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 340) lambdaBG = { 3.589970, 0.514132 };
            else {
                coeff_lambdab = { 2.05363, -0.00685, -3.42739E-04, 3.93987E-06, -1.18237E-08, 0.0 };
                coeff_lambdag = { 1.07658, -0.01041, -4.90553E-05, 1.13528E-06, -3.91609E-09, 0.0 };
            }
        }
        else if (Mzams < 3.5) {
            maxBG = { 500.0, 10.0 };
            if (Radius > 400) lambdaBG = { 116.935557, 0.848808 };
            else {
                maxBG = { 2.5, 1.5 };
                coeff_lambdab     = { 2.40831, -0.42459, 0.03431, -9.26879E-04, 8.24522E-06, 0.0 };
                coeff_lambdag     = { 1.30705, -0.22924, 0.01847, -5.06216E-04, 4.57098E-06, 0.0 };
            }
        }
        else if (Mzams < 4.5) {
            maxBG = { 1000.0, 8.0 };
            if (Radius > 410) lambdaBG = { 52.980056, 1.109736 };
            else {
                maxBG = { 2.5, 1.5 };
                coeff_lambdab     = { 1.8186 , -0.17464, 0.00828, -1.31727E-04, 7.08329E-07, 0.0 };
                coeff_lambdag     = { 1.02183, -0.1024 , 0.00493, -8.16343E-05, 4.55426E-07, 0.0 };
            }
        }
        else if (Mzams < 5.5) {
            maxBG = { 1000.0, 8.0 };
            if (Radius > 430.0) lambdaBG = { 109.593522, 1.324248 };
            else {

                coeff_lambdab = { 1.52581, -0.08125, 0.00219, -2.0527E-05 , 6.79169E-08, 0.0 };
                coeff_lambdag = { 0.85723, -0.04922, 0.00137, -1.36163E-05, 4.68683E-08, 0.0 };
            }
        }
        else if (Mzams < 6.5) {
            maxBG = { 25.5, 5.0 };
            if (Radius > 440.0) lambdaBG = { 16.279603, 1.352166 };
            else {
                coeff_lambdab = { 1.41601, -0.04965, 8.51527E-04, -5.54384E-06, 1.32336E-08, 0.0 };
                coeff_lambdag = { 0.78428, -0.02959, 5.2013E-04 , -3.45172E-06, 8.17248E-09, 0.0 };
            }
        }
        else if (Mzams < 7.5) {
            maxBG = { 9.0, 3.0 };
            if (Radius > 420.0) lambdaBG = { 5.133959, 1.004036 };
            else {
                coeff_lambdab = { 1.38344, -0.04093, 5.78952E-04, -3.19227E-06, 6.40902E-09, 0.0 };
                coeff_lambdag = { 0.76009, -0.02412, 3.47104E-04, -1.92347E-06, 3.79609E-09, 0.0 };
            }
        }
        else if (Mzams < 8.5) {
            maxBG = { 7.0, 3.0 };
            if (Radius > 490.0) lambdaBG = { 4.342985, 0.934659 };
            else {
                coeff_lambdab = { 1.35516, -0.03414, 4.02065E-04, -1.85931E-06, 3.08832E-09, 0.0 };
                coeff_lambdag = { 0.73826, -0.01995, 2.37842E-04, -1.09803E-06, 1.79044E-09, 0.0 };
            }
        }
        else if (Mzams < 9.5) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 530.0) lambdaBG = { 2.441672, 0.702310 };
            else {
                coeff_lambdab  = { 1.32549, -0.02845, 2.79097E-04, -1.07254E-06, 1.46801E-09, 0.0 };
                coeff_lambdag = { 0.71571, -0.01657, 1.64607E-04, -6.31935E-07, 8.52082E-10, 0.0 };
            }
        }
        else if (Mzams < 11.0) {
            maxBG = { 3.0, 1.5 };
            if (Radius > 600.0) lambdaBG = { 1.842314, 0.593854 };
            else {
                maxBG = { 1.0, 0.6 };
                coeff_lambdab     = { 1.29312, -0.02371, 1.93764E-04, -6.19576E-07, 7.04227E-10, 0.0 };
                coeff_lambdag     = { 0.69245, -0.01398, 1.17256E-04, -3.81487E-07, 4.35818E-10, 0.0 };
            }
        }
        else if (Mzams < 13.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 850.0) lambdaBG = { 0.392470, 0.176660 };
            else if (Radius  > 0 && Radius <= 350.0) {
                coeff_lambdab = { 1.28593, -0.02209, 1.79764E-04, -6.21556E-07, 7.59444E-10, 0.0 };
                coeff_lambdag  = { 0.68544, -0.01394, 1.20845E-04, -4.29071E-07, 5.29169E-10, 0.0 };
            }
            else if (Radius > 350.0 && Radius <= 600.0) {
                coeff_lambdab = { -11.99537,  0.0992, -2.8981E-04,  3.62751E-07, -1.65585E-10, 0.0 };
                coeff_lambdag = {   0.46156, -0.0066,  3.9625E-05, -9.98667E-08, -8.84134E-11, 0.0 };
            }
            else {
                coeff_lambdab = { -58.03732, 0.23633, -3.20535E-04, 1.45129E-07, 0.0, 0.0 };
                coeff_lambdag = { -15.11672, 0.06331, -8.81542E-05, 4.0982E-08 , 0.0, 0.0 };
            }
        }
        else if (Mzams < 15) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1000.0)                                       lambdaBG = { 0.414200, 0.189008 };
            else if (Radius > 190.0 && Radius < 600.0) lambdaBG = { 0.15, 0.15 };
            else {
                coeff_lambdab = { 1.39332, -0.0318 , 3.95917E-04, -2.23132E-06, 4.50831E-09, 0.0 };
                coeff_lambdag = { 0.78215, -0.02326, 3.25984E-04, -1.94991E-06, 4.08044E-09, 0.0 };
            }
        }
        else if (Mzams < 18.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1050.0)    lambdaBG = { 0.2, 0.1 };
            else if (Radius > 120.0 && Radius < 170.0) lambdaBG = { 0.2, 0.2 };
            else {
                coeff_lambdab = { 1.43177, -0.03533, 5.11128E-04, -3.57633E-06, 9.36778E-09, 0.0 };
                coeff_lambdag = { 0.85384, -0.03086, 5.50878E-04, -4.37671E-06, 1.25075E-08, 0.0 };
            }
        }
        else if (Mzams < 35.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1200.0) lambdaBG = { 0.05, 0.05 };
            else  lambdaBG = { 1.2 * exp(-Radius / 90.0), 0.55 * exp(-Radius / 160.0) };
        }
        else if (Mzams < 75.0) {
            maxBG = { 1.0, 0.5 };
            coeff_lambdab     = { 0.31321, -7.50384E-04, 5.38545E-07, -1.16946E-10, 0.0, 0.0 };
            coeff_lambdag     = { 0.159  , -3.94451E-04, 2.88452E-07, -6.35132E-11, 0.0, 0.0 };
        }
        else {
            maxBG = { 1.0, 0.5 };
            coeff_lambdab     = { 0.376 , -0.0018 , 2.81083E-06, -1.67386E-09, 3.35056E-13, 0.0 };
            coeff_lambdag     = { 0.2466, -0.00121, 1.89029E-06, -1.12066E-09, 2.2258E-13 , 0.0 };
        }
    }
    else {                                                                  //  popII
        if (Mzams < 1.5) {
            maxBG = { 2.0, 1.5 };
            if (Radius > 160.0) lambdaBG = { 0.05, 0.05 };
            else if (Radius  > 12.0) lambdaBG = { 1.8 * exp(-Radius / 80.0), exp(-Radius / 45.0) };
            else {
                coeff_lambdab = { 0.24012, -0.01907, 6.09529E-04, -8.17819E-06, 4.83789E-08, -1.04568E-10 };
                coeff_lambdag = { 0.15504, -0.01238, 3.96633E-04, -5.3329E-06 , 3.16052E-08, -6.84288E-11 };
            }
        }
        else if (Mzams < 2.50) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 350.0)         lambdaBG = { 2.868539, 0.389991 };
            else if (Radius > 22.0 && Radius < 87.0) lambdaBG = { 1.95, 0.85 };
            else {
                coeff_lambdab = { 2.56108, -0.75562, 0.1027 , -0.00495, 8.05436E-05, 0.0 };
                coeff_lambdag = { 1.41896, -0.4266 , 0.05792, -0.00281, 4.61E-05   , 0.0 };
            }
        }
        else if (Mzams < 3.5) {
            maxBG = { 600.0, 2.0 };
            if (Radius > 400.0) lambdaBG = { 398.126442, 0.648560 };
            else {
                coeff_lambdab = { 1.7814 , -0.17138, 0.00754, -9.02652E-05, 0.0, 0.0 };
                coeff_lambdag = { 0.99218, -0.10082, 0.00451, -5.53632E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 4.5) {
            maxBG = { 600.0, 2.0 };
            if (Radius > 410.0) lambdaBG = { 91.579093, 1.032432 };
            else {
                coeff_lambdab = { 1.65914, -0.10398, 0.0029 , -2.24862E-05, 0.0, 0.0 };
                coeff_lambdag = { 0.92172, -0.06187, 0.00177, -1.42677E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 5.5) {
            maxBG = { 10.0, 3.0 };
            if (Radius > 320.0) lambdaBG = { 7.618019, 1.257919 };
            else {
                coeff_lambdab = { 1.58701, -0.06897, 0.00129    , -6.99399E-06, 0.0, 0.0 };
                coeff_lambdag = { 0.87647, -0.04103, 7.91444E-04, -4.41644E-06, 0.0, 0.0 };
            }
        }
        else if (Mzams < 6.5) {
            maxBG = { 4.0, 1.5 };
            if (Radius > 330.0) lambdaBG = { 2.390575, 0.772091 };
            else {
                coeff_lambdab = { 1.527  , -0.04738, 6.1373E-04 , -2.36835E-06, 0.0, 0.0 };
                coeff_lambdag = { 0.83636, -0.02806, 3.73346E-04, -1.47016E-06, 0.0, 0.0 };
            }
        }
        else if (Mzams < 7.5) {
            maxBG = { 2.5, 1.0 };
            if (Radius > 360.0) lambdaBG = { 1.878174, 0.646353 };
            else {
                coeff_lambdab = { 1.49995, -0.03921, 4.2327E-04, -1.37646E-06, 0.0, 0.0 };
                coeff_lambdag = { 0.81688, -0.02324, 2.5804E-04, -8.54696E-07, 0.0, 0.0 };
            }
        }
        else if (Mzams < 8.5) {
            maxBG = { 2.0, 1.0 };
            if (Radius > 400.0) lambdaBG = { 1.517662, 0.553169 };
            else {
                coeff_lambdab = { 1.46826, -0.03184, 2.85622E-04, -7.91228E-07, 0.0, 0.0 };
                coeff_lambdag = { 0.79396, -0.01903, 1.77574E-04, -5.04262E-07, 0.0, 0.0 };
            }
        }
        else if (Mzams < 9.5) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 440.0) lambdaBG = { 1.136394, 0.478963 };
            else {
                coeff_lambdab = { 1.49196, -0.03247, 3.08066E-04, -9.53247E-07, 0.0, 0.0 };
                coeff_lambdag = { 0.805  , -0.02   , 2.01872E-04, -6.4295E-07 , 0.0, 0.0 };
            }
        }
        else if (Mzams < 11.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 500.0) lambdaBG = { 1.068300, 0.424706 };
            else   lambdaBG = { 1.75 * exp(-Radius / 35.0), 0.9 * exp(-Radius /35.0) };
        }
        else if (Mzams < 13.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 600.0) lambdaBG = { 0.537155, 0.211105 };
            else {
                coeff_lambdab = { 1.63634, -0.04646, 7.49351E-04, -5.23622E-06, 0.0, 0.0 };
                coeff_lambdag = { 1.17934, -0.08481, 0.00329    , -4.69096E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 15.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 650.0) lambdaBG = { 0.3, 0.160696 };
            else {
                coeff_lambdab = { 1.45573, -0.00937, -0.00131,  3.07004E-05, 0.0, 0.0 };
                coeff_lambdag = { 1.19526, -0.08503,  0.00324, -4.58919E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 18.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 750.0) lambdaBG = { 0.5, 0.204092 };
            else {
                coeff_lambdab = { 1.33378,  0.01274, -0.00234,  4.6036E-05 , 0.0, 0.0 };
                coeff_lambdag = { 1.17731, -0.07834,  0.00275, -3.58108E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 35.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 900.0) lambdaBG = { 0.2, 0.107914 };
            else {
                coeff_lambdab = { 1.27138,  0.00538, -0.0012 ,  1.80776E-05, 0.0, 0.0 };
                coeff_lambdag = { 1.07496, -0.05737,  0.00153, -1.49005E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 75.0) {
            maxBG = { 20.0, 3.0 };
            coeff_lambdab     = { 0.821  , -0.00669, 1.57665E-05, -1.3427E-08 , 3.74204E-12, 0.0 };
            coeff_lambdag     = { 0.49287, -0.00439, 1.06766E-05, -9.22015E-09, 2.58926E-12, 0.0 };
        }
        else {
            maxBG = { 4.0, 2.0 };
            coeff_lambdab     = { 1.25332, -0.02065, 1.3107E-04 , -3.67006E-07, 4.58792E-10, -2.09069E-13 };
            coeff_lambdag     = { 0.81716, -0.01436, 9.31143E-05, -2.6539E-07 , 3.30773E-10, -1.51207E-13 };
        }
    }

    return COEF_SET;
}
Lambda_Nanjing::coeff_status Lambda_Nanjing::set_coeff_Cheb(double Mzams, double Z, double Radius) {

    /*
    bool check_phases = Phase_cache==Lookup::Phases::CoreHeBurning;
    //Now check if Mzams and Z are the same of the cached value,
    //check also the phase, if the last phase_0 was not the only for which this function is for, we have to
    //estimate the coefficient even if Mzams and Z have not changed
    if(check_cached(Mzams,Z) and check_phases){
        return COEF_NOT_SET;
    }
    reset_cache(Mzams,Z,Phase);
    */
    reset();

    if (Z > ZpopI) {                 // Z>0.5 Zsun: popI
        if (Mzams < 1.5) {
            maxBG = { 2.5, 1.5 };
            if (Radius >200.00) lambdaBG = { 0.05, 0.05 };
            else {
                coeff_lambdab = { 46.00978, -298.64993, 727.40936, -607.66797, 0.0, 0.0 };
                coeff_lambdag = { 63.61259, -399.89494, 959.62055, -795.20699, 0.0, 0.0 };
            }
        }
        else if (Mzams < 2.5) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 340.0)  lambdaBG = { 3.589970, 0.514132 };
            else if (Radius > 8.5 && Radius < 60.0) lambdaBG = { 3.0, 1.2 };
            else {
                coeff_lambdab = { 34.41826, -6.65259, 0.43823, -0.00953, 0.0, 0.0 };
                coeff_lambdag = { 13.66058, -2.48031, 0.15275, -0.00303, 0.0, 0.0 };
            }
        }
        else if (Mzams < 3.5) {
            maxBG = { 500.0, 10.0 };
            if (Radius > 400.0) lambdaBG = { 116.935557, 0.848808 };
            else {
                maxBG = { 2.5, 1.5 };
                coeff_lambdab     = { -42.98513, 7.90134, -0.54646, 0.01863,  3.13101E-04, 2.07468E-06 };
                coeff_lambdag     = { -6.73842 , 1.06656, -0.05344, 0.00116, -9.34446E-06, 0.0 };
            }
        }
        else if (Mzams < 4.5) {
            maxBG = { 1000.0, 8.0 };
            if (Radius > 410.0) lambdaBG = { 52.980056, 1.109736 };
            else {
                maxBG = { 2.5, 1.5 };
                coeff_lambdab     = { -7.3098 , 0.56647, -0.01176, 7.90112E-05, 0.0, 0.0 };
                coeff_lambdag     = { -3.80455, 0.29308, -0.00603, 4.00471E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 5.5) {
            maxBG = { 1000.0, 8.0 };
            if (Radius > 430.0) lambdaBG = { 109.593522, 1.324248 };
            else {
                coeff_lambdab = { -9.93647, 0.42831, -0.00544, 2.25848E-05, 0.0, 0.0 };
                coeff_lambdag = { -5.33279, 0.22728, -0.00285, 1.16408E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 6.5) {
            maxBG = { 25.5, 5.0 };
            if (Radius > 440.0) lambdaBG = { 16.279603, 1.352166 };
            else {
                coeff_lambdab = { 13.91465, -0.55579, 0.00809, -4.94872E-05, 1.08899E-07, 0.0 };
                coeff_lambdag = {  7.68768, -0.30723, 0.00445, -2.70449E-05, 5.89712E-08, 0.0 };
            }
        }
        else if (Mzams < 7.5) {
            maxBG = { 9.0, 3.0 };
            if (Radius > 420.0) lambdaBG = { 5.133959, 1.004036 };
            else {
                coeff_lambdab = { 4.12387, -0.12979, 0.00153    , -7.43227E-06, 1.29418E-08, 0.0 };
                coeff_lambdag = { 2.18952, -0.06892, 8.00936E-04, -3.78092E-06, 6.3482E-09 , 0.0 };
            }
        }
        else if (Mzams < 8.5) {
            maxBG = { 7.0, 3.0 };
            if (Radius > 490.0) lambdaBG = { 4.342985, 0.934659 };
            else {
                maxBG = { 1.0, 0.5 };
                coeff_lambdab     = { -3.89189, 0.19378, -0.0032 , 2.39504E-05, -8.28959E-08, 1.07843E-10 };
                coeff_lambdag     = { -2.24354, 0.10918, -0.00179, 1.33244E-05, -4.57829E-08, 5.90313E-11 };
            }
        }
        else if (Mzams < 9.5) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 530.0) lambdaBG = { 2.441672, 0.702310 };
            else {
                coeff_lambdab = { 0.86369, -0.00995,  4.80837E-05, -6.10454E-08, -2.79504E-12, 0.0 };
                coeff_lambdag = { -0.7299,  0.0391 , -5.78132E-04,  3.7072E-06 , -1.07036E-08, 1.14833E-11 };
            }
        }
        else if (Mzams < 11.0) {
            maxBG = { 3.0, 1.5 };
            if (Radius > 600.0) lambdaBG = { 1.842314, 0.593854 };
            else {
                coeff_lambdab = { 0.74233, -0.00623, 2.04197E-05, -1.30388E-08, 0.0, 0.0 };
                coeff_lambdag = { 0.36742, -0.00344, 1.27838E-05, -1.0722E-08 , 0.0, 0.0 };
            }
        }
        else if (Mzams < 13.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 850.0) lambdaBG = { 0.392470, 0.176660 };
            else if (Radius > 0. && Radius <= 350.0) {
                coeff_lambdab = { 1.28593, -0.02209, 1.79764E-04, -6.21556E-07, 7.59444E-10, 0.0 };
                coeff_lambdag = { 0.68544, -0.01394, 1.20845E-04, -4.29071E-07, 5.29169E-10, 0.0 };
            }
            else if (Radius > 350.0 && Radius <= 600.0) {
                coeff_lambdab = { -11.99537,  0.0992, -2.8981E-04,  3.62751E-07, -1.65585E-10, 0.0 };
                coeff_lambdag = {   0.46156, -0.0066,  3.9625E-05, -9.98667E-08, -8.84134E-11, 0.0 };
            }
            else {
                coeff_lambdab = { -58.03732, 0.23633, -3.20535E-04, 1.45129E-07, 0.0, 0.0 };
                coeff_lambdag = { -15.11672, 0.06331, -8.81542E-05, 4.0982E-08 , 0.0, 0.0 };
            }
        }
        else if (Mzams < 15.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1000.0)  lambdaBG = { 0.414200, 0.189008 };
            else if (Radius > 69.0 && Radius < 126.0) lambdaBG = { 0.5 - (Radius * 8.77E-04), 0.18 };
            else {
                coeff_lambdab = { 1.12889, -0.00901, 3.04077E-05, -4.31964E-08, 2.14545E-11, 0.0 };
                coeff_lambdag = { 0.568  , -0.0047 , 1.57818E-05, -2.21207E-08, 1.08472E-11, 0.0 };
            }
        }
        else if (Mzams < 18.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1050.0) lambdaBG = { 0.2, 0.1 };
            else {
                coeff_lambdab = { 0.84143, -0.00576, 1.68854E-05, -2.0827E-08 , 8.97813E-12, 0.0 };
                coeff_lambdag = { 0.36014, -0.00254, 7.49639E-06, -9.20103E-09, 3.93828E-12, 0.0 };
            }
        }
        else if (Mzams < 35.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1200.0) lambdaBG = { 0.05, 0.05 };
            else {
                coeff_lambdab = { 0.48724, -0.00177   , 2.60254E-06, -1.25824E-09, 0.0, 0.0 };
                coeff_lambdag = { 0.22693, -8.7678E-04, 1.28852E-06, -6.12912E-10, 0.0, 0.0 };
            }
        }
        else if (Mzams < 75.0) {
            maxBG = { 1.0, 0.5 };
            coeff_lambdab     = { 0.31321, -7.50384E-04, 5.38545E-07, -1.16946E-10, 0.0, 0.0 };
            coeff_lambdag     = { 0.159  , -3.94451E-04, 2.88452E-07, -6.35132E-11, 0.0, 0.0 };
        }
        else {
            maxBG = { 1.0, 0.5 };
            coeff_lambdab     = { 0.376 , -0.0018 , 2.81083E-06, -1.67386E-09, 3.35056E-13, 0.0 };
            coeff_lambdag     = { 0.2466, -0.00121, 1.89029E-06, -1.12066E-09, 2.2258E-13 , 0.0 };
        }
    }
    else {                                                                  // Z<=0.5 Zsun: popI and popII
        if (Mzams < 1.5) {
            maxBG = { 2.0, 1.5 };
            if (Radius > 160.0) lambdaBG = { 0.05, 0.05 };
            else {
                coeff_lambdab = { 0.37294, -0.05825, 0.00375, -7.59191E-05, 0.0, 0.0 };
                coeff_lambdag = { 0.24816, -0.04102, 0.0028 , -6.20419E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 2.5) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 350.0)  lambdaBG = { 2.868539, 0.389991 };
            else if (Radius > 6.0 && Radius < 50.0) lambdaBG = { 0.8, 0.35 };
            else {
                coeff_lambdab = { -103.92538, 25.37325, -2.03273, 0.0543 , 0.0, 0.0 };
                coeff_lambdag = {  -56.03478, 13.6749 , -1.09533, 0.02925, 0.0, 0.0 };
            }
        }
        else if (Mzams < 3.5) {
            maxBG = { 600.0, 2.0 };
            if (Radius > 400.0)   lambdaBG = { 398.126442, 0.648560 };
            else if (Radius > 36.0 && Radius < 53.0) lambdaBG = { 1.0, 1.0 };
            else {
                coeff_lambdab = { -12.40832, 1.59021, -0.06494, 8.69587E-04, 0.0, 0.0 };
                coeff_lambdag = { -6.47476 , 0.8328 , -0.03412, 4.58399E-04, 0.0, 0.0 };
            }
        }
        else if (Mzams < 4.5) {
            maxBG = { 600.0, 2.0 };
            if (Radius > 410.0)  lambdaBG = { 91.579093, 1.032432 };
            else if (Radius > 19.0 && Radius < 85.0) lambdaBG = { 0.255, 0.115 };
            else {
                coeff_lambdab = { -5.89253, 0.54296, -0.01527, 1.38354E-04, 0.0, 0.0 };
                coeff_lambdag = { -3.21299, 0.29583, -0.00833, 7.55646E-05, 0.0, 0.0 };
            }
        }
        else if (Mzams < 5.5) {
            maxBG = { 10.0, 3.0 };
            if (Radius > 320.0)  lambdaBG = { 7.618019, 1.257919 };
            else if (Radius > 85.0 && Radius < 120.0) lambdaBG = { 0.4, 0.1 };
            else {
                coeff_lambdab = { -0.67176, 0.07708, -0.00175   , 1.1991E-05 , 0.0, 0.0 };
                coeff_lambdag = { -0.38561, 0.0427 , -9.6948E-04, 6.64455E-06, 0.0, 0.0 };
            }
        }
        else if (Mzams < 6.5) {
            maxBG = { 4.0, 1.5 };
            if (Radius > 330.0)  lambdaBG = { 2.390575, 0.772091 };
            else if (Radius > 115.0 && Radius < 165.0) lambdaBG = { 0.2, 0.1 };
            else {
                coeff_lambdab = { 0.30941, 0.00965, -2.31975E-04, 1.26273E-06, 0.0, 0.0 };
                coeff_lambdag = { 0.14576, 0.00562, -1.30273E-04, 7.06459E-07, 0.0, 0.0 };
            }
        }
        else if (Mzams < 7.5) {
            maxBG = { 2.5, 1.0 };
            if (Radius > 360.0)   lambdaBG = { 1.878174, 0.646353 };
            else if (Radius > 150.0 && Radius < 210.0) lambdaBG = { 0.2, 0.1 };
            else {
                coeff_lambdab = { 0.44862, 0.00234, -9.23152E-05, 4.67797E-07, 0.0, 0.0 };
                coeff_lambdag = { 0.21873, 0.00154, -5.18806E-05, 2.60283E-07, 0.0, 0.0 };
            }
        }
        else if (Radius < 8.5) {
            maxBG = { 2.0, 1.0 };
            if (Radius > 400.0) lambdaBG = { 1.517662, 0.553169 };
            else if (Radius > 190.0 && Radius < 260.0) lambdaBG = { 0.2, 0.1 };
            else {
                coeff_lambdab = { 0.50221, -3.19021E-04, -3.81717E-05, 1.80726E-07, 0.0, 0.0 };
                coeff_lambdag = { 0.24748, -9.9338E-05 , -1.99272E-05, 9.47504E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 9.5) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 440.0)  lambdaBG = { 1.136394, 0.478963 };
            else if (Radius > 180.0 && Radius < 300.0) lambdaBG = { 0.15, 0.08 };
            else {
                coeff_lambdab = { 0.39342, 0.00259    , -4.97778E-05, 1.69533E-07, 0.0, 0.0 };
                coeff_lambdag = { 0.20796, 6.62921E-04, -1.84663E-05, 6.58983E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 11.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 500.0) lambdaBG = { 1.068300, 0.424706 };
            else {
                coeff_lambdab = { 0.75746, -0.00852, 3.51646E-05, -4.57725E-08, 0.0, 0.0 };
                coeff_lambdag = { 0.35355, -0.00388, 1.56573E-05, -1.98173E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 13.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 600.0)   lambdaBG = { 0.537155, 0.211105 };
            else if (Radius > 200.0 && Radius < 410.0) lambdaBG = { 0.08, 0.05 };
            else {
                coeff_lambdab = { 0.85249, -0.00861, 2.99246E-05, -3.21416E-08, 0.0, 0.0 };
                coeff_lambdag = { 0.37188, -0.00365, 1.24944E-05, -1.32388E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 15.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 650.0) lambdaBG = { 0.3, 0.160696 };
            else if (Radius > 250.0 && Radius< 490.0) lambdaBG = { 0.06, 0.05 };
            else {
                coeff_lambdab = { 0.85271, -0.00793, 2.5174E-05 , -2.4456E-08 , 0.0, 0.0 };
                coeff_lambdag = { 0.36163, -0.00328, 1.03119E-05, -9.92712E-09, 0.0, 0.0 };
            }
        }
        else if (Mzams < 18.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 750.0) lambdaBG = { 0.5, 0.204092 };
            else if (Radius > 200.0 && Radius < 570.0) lambdaBG = { 0.1, 0.05 };
            else {
                coeff_lambdab = { 0.83254, -0.00696, 1.9597E-05 , -1.67985E-08, 0.0, 0.0 };
                coeff_lambdag = { 0.34196, -0.0028 , 7.82865E-06, -6.66684E-09, 0.0, 0.0 };
            }
        }
        else if (Mzams < 35.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 900.0)   lambdaBG = { 0.2, 0.107914 };
            else if (Radius > 230.0 && Radius < 755.0) lambdaBG = { 0.1, 0.05 };
            else {
                coeff_lambdab = { 0.69746, -0.0043 , 8.97312E-06, -5.83402E-09, 0.0, 0.0 };
                coeff_lambdag = { 0.26691, -0.00161, 3.3378E-06 , -2.1555E-09 , 0.0, 0.0 };
            }
        }
        else if (Mzams < 75.0) {
            maxBG = { 20.0, 3.0 };
            coeff_lambdab     = { 0.821  , -0.00669, 1.57665E-05, -1.3427E-08 , 3.74204E-12, 0.0 };
            coeff_lambdag     = { 0.49287, -0.00439, 1.06766E-05, -9.22015E-09, 2.58926E-12, 0.0 };
        }
        else {
            maxBG = { 4.0, 2.0 };
            coeff_lambdab     = { 1.25332, -0.02065, 1.3107E-04 , -3.67006E-07, 4.58792E-10, -2.09069E-13 };
            coeff_lambdag     = { 0.81716, -0.01436, 9.31143E-05, -2.6539E-07 , 3.30773E-10, -1.51207E-13 };
        }
    }


    return COEF_SET;
}
Lambda_Nanjing::coeff_status Lambda_Nanjing::set_coeff_AGB(double Mzams, double Z, double Radius) {

    /*
    bool check_phases= Phase_cache==Lookup::Phases::TerminalCoreHeBurning or Phase_cache==Lookup::Phases::ShellHeBurning;
    //Now check if Mzams and Z are the same of the cached value,
    //check also the phase, if the last phase_0 was not the only for which this function is for, we have to
    //estimate the coefficient even if Mzams and Z have not changed
    if(check_cached(Mzams,Z) and check_phases){
        return COEF_NOT_SET;
    }
    reset_cache(Mzams,Z,Phase);
    */
    reset();

    if (Z > ZpopI) {                 // Z>0.5 Zsun: popI
        if (Mzams < 1.5) {
            maxBG = { 2.5, 1.5 };
            if (Radius > 200.0) lambdaBG = { 0.05, 0.05 };
            else  {
                double tmp = 0.1 - (Radius * 3.57E-04);
                lambdaBG   = { tmp, tmp };
            }
        }
        else if (Mzams < 2.5) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 340.0) lambdaBG = { 3.589970, 0.514132 };
            else {
                coeff_lambdab = { 0.88954, 0.0098 , -3.1411E-05 , 7.66979E-08,  0.0       , 0.0 };
                coeff_lambdag = { 0.48271, 0.00584, -6.22051E-05, 2.41531E-07, -3.1872E-10, 0.0 };
            }
        }
        else if (Mzams < 3.5) {
            maxBG = { 500.0, 10.0 };
            if (Radius > 400.0) lambdaBG = { 116.935557, 0.848808 };
            else {
                coeff_lambdab = { -0.04669, 0.00764, -4.32726E-05, 9.31942E-08, 0.0        ,  0.0 };
                coeff_lambdag = {  0.44889, 0.01102, -6.46629E-05, 5.66857E-09, 7.21818E-10, -1.2201E-12 };
            }
        }
        else if (Mzams < 4.5) {
            maxBG = { 1000.0, 8.0 };
            if (Radius > 410.0) lambdaBG = { 52.980056, 1.109736 };
            else {
                coeff_lambdab = { -0.37322, 0.00943, -3.26033E-05, 5.37823E-08, 0.0, 0.0 };
                coeff_lambdag = {  0.13153, 0.00984, -2.89832E-05, 2.63519E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 5.5) {
            maxBG = { 1000.0, 8.0 };
            if (Radius > 430.0) lambdaBG = { 109.593522, 1.324248 };
            else {
                coeff_lambdab = { -0.80011, 0.00992, -3.03247E-05,  5.26235E-08, 0.0, 0.0 };
                coeff_lambdag = { -0.00456, 0.00426,  4.71117E-06, -1.72858E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 6.5) {
            maxBG = { 25.5, 5.0 };
            if (Radius > 440.0) lambdaBG = { 16.279603, 1.352166 };
            else {
                coeff_lambdab = { -2.7714 ,  0.06467, -4.01537E-04,  7.98466E-07, 0.0, 0.0 };
                coeff_lambdag = {  0.23083, -0.00266,  2.21788E-05, -2.35696E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 7.5) {
            maxBG = { 9.0, 3.0 };
            if (Radius > 420.0) lambdaBG = { 5.133959, 1.004036 };
            else {
                coeff_lambdab = { -0.63266,  0.02054, -1.3646E-04 ,  2.8661E-07 , 0.0, 0.0 };
                coeff_lambdag = {  0.26294, -0.00253,  1.32272E-05, -7.12205E-09, 0.0, 0.0 };
            }
        }
        else if (Mzams < 8.5) {
            maxBG = { 7.0, 3.0 };
            if (Radius > 490.0) lambdaBG = { 4.342985, 0.934659 };
            else {
                coeff_lambdab = { -0.1288 ,  0.0099 , -6.71455E-05,  1.33568E-07, 0.0, 0.0 };
                coeff_lambdag = {  0.26956, -0.00219,  7.97743E-06, -1.53296E-09, 0.0, 0.0 };
            }
        }
        else if (Mzams < 9.5) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 530.0) lambdaBG = { 2.441672, 0.702310 };
            else {
                coeff_lambdab = { 1.19804, -0.01961, 1.28222E-04, -3.41278E-07, 3.35614E-10, 0.0 };
                coeff_lambdag = { 0.40587, -0.0051 , 2.73866E-05, -5.74476E-08, 4.90218E-11, 0.0 };
            }
        }
        else if (Mzams < 11.0) {
            maxBG = { 3.0, 1.5 };
            if (Radius >600.0) lambdaBG = { 1.842314, 0.593854 };
            else {
                coeff_lambdab = { 0.3707 ,  2.67221E-04, -9.86464E-06, 2.26185E-08, 0.0, 0.0 };
                coeff_lambdag = { 0.25549, -0.00152    ,  3.35239E-06, 2.24224E-10, 0.0, 0.0 };
            }
        }
        else if (Mzams < 13.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 850.0) lambdaBG = { 0.392470, 0.176660 };
            else if (Radius > 0. && Radius <= 350.0) {
                coeff_lambdab = { 1.28593, -0.02209, 1.79764E-04, -6.21556E-07, 7.59444E-10, 0.0 };
                coeff_lambdag = { 0.68544, -0.01394, 1.20845E-04, -4.29071E-07, 5.29169E-10, 0.0 };
            }
            else if (Radius > 350.0 && Radius <= 600.0) {
                coeff_lambdab = { -11.99537,  0.0992, -2.8981E-04,  3.62751E-07, -1.65585E-10, 0.0 };
                coeff_lambdag = {   0.46156, -0.0066,  3.9625E-05, -9.98667E-08, -8.84134E-11, 0.0 };
            }
            else {
                coeff_lambdab = { -58.03732, 0.23633, -3.20535E-04, 1.45129E-07, 0.0, 0.0 };
                coeff_lambdag = { -15.11672, 0.06331, -8.81542E-05, 4.0982E-08 , 0.0, 0.0 };
            }
        }
        else if (Mzams < 15.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1000.0) lambdaBG = { 0.414200, 0.189008 };
            else {
                coeff_lambdab = { -106.90553, 0.36469, -4.1472E-04 , 1.57349E-07, 0.0, 0.0 };
                coeff_lambdag = {  -39.93089, 0.13667, -1.55958E-04, 5.94076E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 18.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1050.0) lambdaBG = { 0.2, 0.1 };
            else {
                coeff_lambdab = { -154.70559, 0.46718, -4.70169E-04, 1.57773E-07, 0.0, 0.0 };
                coeff_lambdag = {  -65.39602, 0.19763, -1.99078E-04, 6.68766E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 35.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius > 1200.0) lambdaBG = { 0.05, 0.05 };
            else {
                coeff_lambdab = { -260484.85724, 4.26759E+06, -2.33016E+07, 4.24102E+07, 0.0, 0.0 };
                coeff_lambdag = { -480055.67991, 7.87484E+06, -4.30546E+07, 7.84699E+07, 0.0, 0.0 };
            }
        }
        else if (Mzams < 75.0) {
            maxBG = { 1.0, 0.5 };
            coeff_lambdab     = { 0.31321, -7.50384E-04, 5.38545E-07, -1.16946E-10, 0.0, 0.0 };
            coeff_lambdag     = { 0.159  , -3.94451E-04, 2.88452E-07, -6.35132E-11, 0.0, 0.0 };
        }
        else {
            maxBG = { 1.0, 0.5 };
            coeff_lambdab     = { 0.376 , -0.0018 , 2.81083E-06, -1.67386E-09, 3.35056E-13, 0.0 };
            coeff_lambdag     = { 0.2466, -0.00121, 1.89029E-06, -1.12066E-09, 2.2258E-13 , 0.0 };
        }
    }
    else {                                                                          // popII
        if (Mzams < 1.5) {
            maxBG = { 2.0, 1.5 };
            if (Radius > 160.0) lambdaBG = { 0.05, 0.05 };
            else {
                coeff_lambdab = { 0.24012, -0.01907, 6.09529E-04, -8.17819E-06, 4.83789E-08, -1.04568e-10 };
                coeff_lambdag = { 0.15504, -0.01238, 3.96633E-04, -5.3329E-06 , 3.16052E-08, -6.84288e-11 };
            }
        }
        else if (Mzams < 2.5) {
            maxBG = { 4.0, 2.0 };
            if (Radius > 350.0) lambdaBG = { 2.868539, 0.389991 };
            else {
                coeff_lambdab = { 0.5452 ,  0.00212    , 6.42941E-05, -1.46783E-07, 0.0       ,  0.0 };
                coeff_lambdag = { 0.30594, -9.58858E-04, 1.12174E-04, -1.04079E-06, 3.4564E-09, -3.91536e-12 };
            }
        }
        else if (Mzams < 3.5) {
            maxBG = { 600.0, 2.0 };
            if (Radius > 400.0)  lambdaBG = { 398.126442, 0.648560 };
            else if (Radius > 36.0 && Radius < 53.0) lambdaBG = { 1.0, 1.0 };
            else {
                coeff_lambdab = { -0.475  , -0.00328, 1.31101E-04, -6.03669E-07, 8.49549E-10, 0.0 };
                coeff_lambdag = {  0.05434,  0.0039 , 9.44609E-06, -3.87278E-08, 0.0        , 0.0 };
            }
        }
        else if (Mzams < 4.5) {
            maxBG = { 600.0, 2.0 };
            if (Radius > 410.0) lambdaBG = { 91.579093, 1.032432 };
            else {
                coeff_lambdab = { -0.2106 , -0.01574, 2.01107E-04, -6.90334E-07, 7.92713E-10, 0.0 };
                coeff_lambdag = {  0.36779, -0.00991, 1.19411E-04, -3.59574E-07, 3.33957E-10, 0.0 };
            }
        }
        else if (Mzams < 5.5) {
            maxBG = { 10.0, 3.0 };
            if (Radius > 320.0) lambdaBG = { 7.618019, 1.257919 };
            else {
                coeff_lambdab  = { -0.12027,  0.01981, -2.27908E-04,  7.55556E-07, 0.0, 0.0 };
                coeff_lambdag = {  0.31252, -0.00527,  3.60348E-05, -3.22445E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 6.5) {
            maxBG = { 4.0, 1.5 };
            if (Radius > 330.0) lambdaBG = { 2.390575, 0.772091 };
            else {
                coeff_lambdab = { 0.26578,  0.00494, -7.02203E-05, 2.25289E-07, 0.0, 0.0 };
                coeff_lambdag = { 0.26802, -0.00248,  6.45229E-06, 1.69609E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 7.5) {
            maxBG = { 2.5, 1.0 };
            if (Radius > 360.0) lambdaBG = { 1.878174, 0.646353 };
            else {
                coeff_lambdab = { 0.8158 , -0.01633, 1.46552E-04, -5.75308E-07, 8.77711E-10, 0.0 };
                coeff_lambdag = { 0.26883, -0.00219, 4.12941E-06,  1.33138E-08, 0.0        , 0.0 };
            }
        }
        else if (Mzams < 8.5) {
            maxBG = { 2.0, 1.0 };
            if (Radius > 400.0) lambdaBG = { 1.517662, 0.553169 };
            else {
                coeff_lambdab = { 0.74924, -0.01233, 9.55715E-05, -3.37117E-07, 4.67367E-10, 0.0 };
                coeff_lambdag = { 0.25249, -0.00161, 8.35478E-07,  1.25999E-08, 0.0        , 0.0 };
            }
        }
        else if (Mzams < 9.5) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 440.0) lambdaBG = { 1.136394, 0.478963 };
            else {
                coeff_lambdab = { 0.73147, -0.01076, 7.54308E-05, -2.4114E-07 , 2.95543E-10, 0.0 };
                coeff_lambdag = { 0.31951, -0.00392, 2.31815E-05, -6.59418E-08, 7.99575E-11, 0.0 };
            }
        }
        else if (Mzams < 11.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 500.0) lambdaBG = { 1.068300, 0.424706 };
            else {
                coeff_lambdab = { -9.26519,  0.08064, -2.30952E-04, 2.21986E-07, 0.0, 0.0 };
                coeff_lambdag = {  0.81491, -0.00161, -8.13352E-06, 1.95775E-08, 0.0, 0.0 };
            }
        }
        else if (Mzams < 13.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius > 600.0)  lambdaBG = { 0.537155, 0.211105 };
            else if (Radius > 390.0 && Radius < 460.0) lambdaBG = { 0.08, 0.05 };
            else {
                coeff_lambdab = { -51.15252, 0.30238, -5.95397E-04, 3.91798E-07, 0.0, 0.0 };
                coeff_lambdag = { -13.44   , 0.08141, -1.641E-04  , 1.106E-07  , 0.0, 0.0 };
            }
        }
        else if (Mzams < 15.0) {
            maxBG = { 1.6, 1.0 };
            if (Radius >  650.0)  lambdaBG = { 0.3, 0.160696 };
            else if (Radius >  480.0 && Radius <  540.0) lambdaBG = { 0.06, 0.05 };
            else if (Radius >= 540.0 && Radius <= 650.0) lambdaBG = { (Radius * 1.8E-03) - 0.88, (Radius * 9.1E-04) - 0.43 };
            else {
                coeff_lambdab = { -140.0   , 0.7126 , -0.00121    , 6.846E-07  , 0.0, 0.0 };
                coeff_lambdag = {  -44.1964, 0.22592, -3.85124E-04, 2.19324E-07, 0.0, 0.0 };
            }
        }
        else if (Mzams < 18.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius >  750.0)  lambdaBG = { 0.5, 0.204092 };
            else if (Radius >  560.0 && Radius <  650.0) lambdaBG = { 0.1, 0.05 };
            else if (Radius >= 650.0 && Radius <= 750.0) lambdaBG = { (Radius * 4.0E-03) - 2.5, (Radius * 1.5E-03) - 0.93 };
            else {
                coeff_lambdab = { -358.4    , 1.599  , -0.00238   , 1.178E-06  , 0.0, 0.0 };
                coeff_lambdag = { -118.13757, 0.52737, -7.8479E-04, 3.89585E-07, 0.0, 0.0 };
            }
        }
        else if (Mzams < 35.0) {
            maxBG = { 1.5, 1.0 };
            if (Radius >  900.0)  lambdaBG = { 0.2, 0.107914 };
            else if (Radius >  725.0 && Radius <  850.0) lambdaBG = { 0.1, 0.05 };
            else if (Radius >= 850.0 && Radius <= 900.0) lambdaBG = { (Radius * 2.0E-03) - 1.6, (Radius * 1.0E-03) - 0.8 };
            else {
                coeff_lambdab = { -436.00777, 1.41375, -0.00153    , 5.47573E-07, 0.0, 0.0 };
                coeff_lambdag = { -144.53456, 0.46579, -4.99197E-04, 1.78027E-07, 0.0, 0.0 };
            }
        }
        else if (Mzams < 75.0) {
            maxBG = { 20.0, 3.0 };
            coeff_lambdab     = { 0.821  , -0.00669, 1.57665E-05, -1.3427E-08 , 3.74204E-12, 0.0 };
            coeff_lambdag     = { 0.49287, -0.00439, 1.06766E-05, -9.22015E-09, 2.58926E-12, 0.0 };
        }
        else {
            maxBG = { 4.0, 2.0 };
            coeff_lambdab     = { 1.25332, -0.02065, 1.3107E-04 , -3.67006E-07, 4.58792E-10, -2.09069E-13 };
            coeff_lambdag     = { 0.81716, -0.01436, 9.31143E-05, -2.6539E-07 , 3.30773E-10, -1.51207E-13 };
        }
    }

    return COEF_SET;
}

double Lambda_Nanjing::lambda_Giant(double Mzams,   double Z, double Radius, double Menv, double lambda_th) {


    set_coeff_Giant(Mzams,Z,Radius);


    if (lambdaBG.empty()){
        if (Z>ZpopI and Mzams<1.5){
            auto y = fitting_equation(Menv);
            lambdaBG = {1.0/y[0], 1.0/y[1]};
        } else{
            lambdaBG = fitting_equation(Radius);
        }
    }

    lambdaBG[0] = std::min(std::max(0.05, lambdaBG[0]), maxBG[0]);
    lambdaBG[1] = std::min(std::max(0.05, lambdaBG[1]), maxBG[1]);

    return lambda_estimate(lambda_th);
}
double Lambda_Nanjing::lambda_Cheb(double Mzams,   double Z, double Radius, double Menv, double lambda_th) {


    set_coeff_Cheb(Mzams,Z,Radius);

    if (lambdaBG.empty()){
        if (Z>ZpopI and Mzams<1.5){
            auto y = fitting_equation(Menv);
            lambdaBG = {1.0/y[0], 1.0/y[1]};
        } else{
            lambdaBG = fitting_equation(Radius);
        }
    } else{
        lambdaBG = lambdaBG;
    }

    lambdaBG[0] = std::min(std::max(0.05, lambdaBG[0]), maxBG[0]);
    lambdaBG[1] = std::min(std::max(0.05, lambdaBG[1]), maxBG[1]);

    return lambda_estimate(lambda_th);
}
double Lambda_Nanjing::lambda_AGB(double Mzams,   double Z, double Radius, double Menv, double lambda_th) {

    set_coeff_AGB(Mzams,Z,Radius);

    if (lambdaBG.empty()){
        if (Z>ZpopI and (Mzams<1.5 or (Mzams>=18 and Mzams<25.0))){
            auto y = fitting_equation(Menv);
            //std::cout<< " MW "<< Menv << " "<< 1/y[0] << " " << 1/y[1] << std::endl;
                     lambdaBG = {1.0/y[0], 1.0/y[1]};
        } else if ( (Z>ZpopI and Mzams>=2.5 and Mzams<5.5) or
                    (Z<=ZpopI and Mzams>=2.5 and Mzams<4.5)){
            auto y = fitting_equation(Radius);
            lambdaBG = {std::pow(10.0,y[0]), y[1]};
        }else{
            lambdaBG = fitting_equation(Radius);
        }
    } else{
        lambdaBG = lambdaBG;
    }

    lambdaBG[0] = std::min(std::max(0.05, lambdaBG[0]), maxBG[0]);          // clamp lambda B to [0.05, maxB]
    lambdaBG[1] = std::min(std::max(0.05, lambdaBG[1]), maxBG[1]);          // clamp lambda G to [0.05, maxG]

    return lambda_estimate(lambda_th);
}
double Lambda_Nanjing::lambda_pureHe(_UNUSED double Mzams, _UNUSED   double Z, double Radius, _UNUSED double Menv, _UNUSED double lambda_th) {

    double rMin = 0.25;                              // minimum considered radius
    double rMax = 120.0;                             // maximum considered radius
    Radius = std::min(std::max(rMin,Radius),rMax);

    return 0.3 * std::pow(Radius,-0.8);
}


double Lambda_Nanjing::lambda_Giant(const Star *s) {

    const double& Z         = s->get_Z();
    const double& Mzams     = s->get_zams();
    const double& Radius    = s->getp(Radius::ID);
    const double& Menv      = (s->getp(Mass::ID)-s->getp(MHE::ID))/s->getp(Mass::ID);
    const double& lambda_th = s->get_svpar_num("star_lambda_fth");

    return lambda_Giant(Mzams,Z,Radius,Menv,lambda_th);
}
double Lambda_Nanjing::lambda_Cheb(const Star *s) {

    const double& Z         = s->get_Z();
    const double& Mzams     = s->get_zams();
    const double& Radius    = s->getp(Radius::ID);
    const double& Menv      = (s->getp(Mass::ID)-s->getp(MHE::ID))/s->getp(Mass::ID);
    const double& lambda_th = s->get_svpar_num("star_lambda_fth");

    return lambda_Cheb(Mzams,Z,Radius,Menv,lambda_th);
}
double Lambda_Nanjing::lambda_AGB(const Star *s) {

    const double& Z         = s->get_Z();
    const double& Mzams     = s->get_zams();
    const double& Radius    = s->getp(Radius::ID);
    const double& Menv      = (s->getp(Mass::ID)-s->getp(MHE::ID))/s->getp(Mass::ID);
    const double& lambda_th = s->get_svpar_num("star_lambda_fth");

    return lambda_AGB(Mzams,Z,Radius,Menv,lambda_th);
}
double Lambda_Nanjing::lambda_pureHe(const Star *s) {

    const double& Z         = s->get_Z();
    const double& Mzams     = s->get_zams();
    const double& Radius    = s->getp(Radius::ID);
    const double& Menv      = (s->getp(Mass::ID)-s->getp(MHE::ID))/s->getp(Mass::ID);
    const double& lambda_th = s->get_svpar_num("star_lambda_fth");

    return lambda_pureHe(Mzams,Z,Radius,Menv,lambda_th);
}

double Lambda_Nanjing::operator()(const Star *s) {



    int Phase = int(s->getp(PhaseBSE::ID));
    if (Phase==7 or Phase==8 or Phase==9) //WR or naked helium
        return lambda_pureHe(s);
    else if (Phase==1 or Phase==2 or Phase==3){
        return lambda_Giant(s);
    }
    else if (Phase==4)
        return lambda_Cheb(s);
    else if (Phase==5 or Phase==6)
        return lambda_AGB(s);

    return std::nan("");

}

//Mzams list from the update version of Xu&Li+10 - errata corrige https://iopscience.iop.org/article/10.1088/0004-637X/722/2/1985/pdf
std::vector<double> Lambda_Nanjing_interpolator::Mzams_list = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,20.,50.,100.};

void Lambda_Nanjing_interpolator::find_interpolators(double Mzams) {

    size_t Mzams_id [2];

    Mzams_id[0] = std::min(utilities::binary_search(&Mzams_list[0], 0, Mzams_list.size()-1, Mzams),Mzams_list.size()-2); //search the metallicity of the star (Z)
    Mzams_id[1] = Mzams_id[0] + 1; //search the metallicity of the star (Z)

    Mzams_interpolators[0] = Mzams_list.at(Mzams_id[0]);
    Mzams_interpolators[1] = Mzams_list.at(Mzams_id[1]);

    Mzams = std::max(Mzams_interpolators[0],std::min(Mzams,Mzams_interpolators[1]));
    wM[0] = (Mzams_interpolators[1] - Mzams) / (Mzams_interpolators[1] - Mzams_interpolators[0]);
    wM[1] = (Mzams - Mzams_interpolators[0]) / (Mzams_interpolators[1] - Mzams_interpolators[0]);

}

double Lambda_Nanjing_interpolator::operator()(const Star *s) {

    const int& Phase = int(s->getp(PhaseBSE::ID));
    //If naked helium/WR we do not need to use the interpolator, lambda depends only on radius
    if (Phase==7 or Phase==8 or Phase==9){
        return lambda_pureHe(s);
    }

    const double& Z         = s->get_Z();
    const double& Mzams     = s->get_zams();
    const double& Radius    = s->getp(Radius::ID);
    const double& Menv      = (s->getp(Mass::ID)-s->getp(MHE::ID))/s->getp(Mass::ID);
    const double& lambda_th = s->get_svpar_num("star_lambda_fth");

    if (Mzams!=Mzams_cache or Z!=Z_cache){
        Mzams_cache=Mzams;
        Z_cache=Z;
        find_interpolators(Mzams_cache);
    }

    double lambda0 = operator()(Phase,Mzams_interpolators[0],Z,Radius,Menv,lambda_th);
    double lambda1 = operator()(Phase,Mzams_interpolators[1],Z,Radius,Menv,lambda_th);


    return wM[0]*lambda0 + wM[1]*lambda1;
}

double Lambda_Nanjing_interpolator::operator()(unsigned int Phase, double Mzams, double Z, double Radius, double Menv,
                                               double lambda_th) {


    if (Phase==1 or Phase==2 or Phase==3){
        return lambda_Giant(Mzams, Z, Radius, Menv, lambda_th);
    }
    else if (Phase==4)
        return lambda_Cheb(Mzams, Z, Radius, Menv, lambda_th);
    else if (Phase==5 or Phase==6)
        return lambda_AGB(Mzams, Z, Radius, Menv, lambda_th);

    return std::nan("");
}
