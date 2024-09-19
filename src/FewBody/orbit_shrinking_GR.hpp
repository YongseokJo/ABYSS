#pragma once

#include <cmath>

#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "Particle.h"


// made 2024.09.19 by Eunwoo Chung

// reference: tides3.f from Nbody6++GPU (Rizzuto et al. 2020, Arca Sedda et al. 2023)
// Gravitational wave energy loss of hard binaries according to the orbit averaged approximation of Peters & Mathews 1963.
// Calculate average change of energy and angular momentum per orbit.
// Modify semi-major axis and eccentricity per time step in SDAR.
inline void orbit_shrinking_GR(AR::BinaryTree<Particle> _bin, REAL dtime) { // AR::BinaryTree<Particle>& _bin

    REAL c = 299752.458/(velocity_unit/yr*pc/1e5); // speed of light

    REAL e = _bin.ecc;
    REAL semi = _bin.semi;
    REAL m1 = _bin.m1;
    REAL m2 = _bin.m2;

    REAL e2 = e*e;
    REAL e4 = e2*e2;
    REAL mtot = m1 + m2;
    REAL cost = pow(c, -5)*m1*m2*mtot;

    REAL FE = 0;
    REAL FE1 = 0;
    REAL FE2 = 0;
    REAL de = 0;
    REAL dsemi = 0;

// Eliptic case
    if (e < 1) {
        FE = (1.0 + 73./24.*e2 + 37./96.*e4);
        FE1 = pow((1 - e2), -3.5)*FE;
        FE2 = e*pow((1-e2), -2.5)*(1. + 121./304.*e2);

        de = 304./15.*cost*pow(semi, -4.)*FE2*dtime;
        dsemi = 64./5.*cost*pow(semi, -3.)*FE1*dtime;
    }
// Hyperbolic case
    else if (e > 1) {
        FE1 = 3.*(96. + 292.*e2 + 37.*e4)*acos(-1/e);
        FE1 = FE1 + (673.*e2 + 602.)*sqrt(e2-1);
        FE1 = 1./M_PI*pow((e2-1), -3.5)*FE1;

        FE2 = 3.*e2*(304. + 121.*e2)*acos(-1/e);
        FE2 = FE2 + (134. + 1069.*e2 + 72.*e4)*sqrt(e2-1);
        FE2 = 1./M_PI/e*pow((e2-1), -2.5)*FE2;

        de = 1./45.*cost*pow(abs(semi), -4.)*FE2*dtime;
        dsemi = 2./45.*cost*pow(abs(semi), -3.)*FE1*dtime;
    }
    _bin.semi -= de;
    _bin.semi -= dsemi;
}

// Kick BH when BH merger happens
// void GWmergerKick(COMM::BinaryTree<Particle,COMM::Binary> &_bin) {

// }