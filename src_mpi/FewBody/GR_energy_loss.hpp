#pragma once

#include <cmath>

#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "Particle.h"
#include "ar_interaction.hpp"

extern REAL min_timestep;


// made 2024.09.19 by Eunwoo Chung

// reference: tides3.f from Nbody6++GPU (Rizzuto et al. 2020, Arca Sedda et al. 2023)
// Gravitational wave energy loss of hard binaries according to the orbit averaged approximation of Peters & Mathews 1963.
// Calculate average change of energy and angular momentum per orbit.
// Modify semi-major axis, eccentricity, omega(argument of periapsis) per time step in SDAR.
// Orbit shrinking by PN2.5
// Precession by PN1.0 & PN2.0
/*
inline void GR_energy_loss(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, REAL current_time, REAL next_time) {


    REAL c = 299752.458/(velocity_unit/yr*pc/1e5); // speed of light

    REAL m1 = _bin.m1;
    REAL m2 = _bin.m2;
    REAL mtot = m1 + m2;

    REAL e = 0;
    REAL semi = 0;
    REAL e2 = 0;
    REAL e4 = 0;
    REAL cost = 0;

    REAL FE = 0;
    REAL FE1 = 0;
    REAL FE2 = 0;
    REAL de = 0;
    REAL dsemi = 0;
    REAL domega = 0;

    REAL dt = min_timestep*EnzoTimeStep;
    REAL dtime = next_time - current_time;
    int iteration_num = static_cast<int>(round(dtime/min_timestep));

    // if (iteration_num != 1 && iteration_num % 2 != 0)
    //     fprintf(stderr, "iteration_num: %d, _dtime: %e", iteration_num, _dtime);
    // assert(iteration_num == 1 || iteration_num % 2 == 0);

    int num = 0;
    while (num < iteration_num) {

        e = _bin.ecc;
        semi = _bin.semi;
        e2 = e*e;
        e4 = e2*e2;
        cost = pow(c, -5)*m1*m2*mtot;

        FE = (1.0 + 73./24.*e2 + 37./96.*e4);
        FE1 = pow((1 - e2), -3.5)*FE;
        FE2 = e*pow((1-e2), -2.5)*(1. + 121./304.*e2);

        de = 304./15.*cost*pow(semi, -4.)*FE2*dt;
        dsemi = 64./5.*cost*pow(semi, -3.)*FE1*dt;
        domega = (6.*M_PI/(c*c*_bin.period)*mtot/(semi*(1-e2)) + 3.*(18.+e2)*M_PI/(2*pow(c, 4.)*_bin.period)*pow((mtot/(semi*(1-e2))), 2.))*dt;

        _bin.ecc -= de;
        _bin.semi -= dsemi;
        _bin.rot_self += domega;

        num += 1;

        // if (_bin.getMember(0)->PID == 716 || _bin.getMember(0)->PID == 44)
        //     fprintf(mergerout, "dt: %e Myr, de: %e, dsemi %e pc\n", dt*1e4, de, dsemi*position_unit);

        if (_bin.ecc < 0 || _bin.semi < 0) {

            Interaction interaction;

            _bin_interrupt.time_now = current_time*EnzoTimeStep + dt*num;
            interaction.merge(_bin_interrupt, _bin);

            return;
        }
    }
    _bin.calcParticles(REAL(1.0));
    for (int dim=0; dim<Dim; dim++) {
        _bin.getLeftMember()->Position[dim] += _bin.Position[dim];
        _bin.getLeftMember()->Velocity[dim] += _bin.Velocity[dim];
        _bin.getRightMember()->Position[dim] += _bin.Position[dim];
        _bin.getRightMember()->Velocity[dim] += _bin.Velocity[dim];
    }

    // // Eliptic case
    //     if (e < 1) {
    //         FE = (1.0 + 73./24.*e2 + 37./96.*e4);
    //         FE1 = pow((1 - e2), -3.5)*FE;
    //         FE2 = e*pow((1-e2), -2.5)*(1. + 121./304.*e2);

    //         de = 304./15.*cost*pow(semi, -4.)*FE2*_dtime;
    //         dsemi = 64./5.*cost*pow(semi, -3.)*FE1*_dtime;
    //     // fprintf(mergerout, "dt: %e, de: %e, dsemi %e (pc)\n", _dtime*1e4, de, dsemi*position_unit);
    //     }
    // // Hyperbolic case

    //     else if (e > 1) {
    //         FE1 = 3.*(96. + 292.*e2 + 37.*e4)*acos(-1/e);
    //         FE1 = FE1 + (673.*e2 + 602.)*sqrt(e2-1);
    //         FE1 = 1./M_PI*pow((e2-1), -3.5)*FE1;

    //         FE2 = 3.*e2*(304. + 121.*e2)*acos(-1/e);
    //         FE2 = FE2 + (134. + 1069.*e2 + 72.*e4)*sqrt(e2-1);
    //         FE2 = 1./M_PI/e*pow((e2-1), -2.5)*FE2;

    //         de = 1./45.*cost*pow(abs(semi), -4.)*FE2*_dtime;
    //         dsemi = 2./45.*cost*pow(abs(semi), -3.)*FE1*_dtime;
    //     }
}
*/

inline void GR_energy_loss(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, REAL current_time, REAL next_time) {

    const REAL c = 299752.458 / (velocity_unit / yr * pc / 1e5); // speed of light in your units
    const REAL m1 = _bin.m1;
    const REAL m2 = _bin.m2;
    const REAL mtot = m1 + m2;
    const REAL cost = pow(c, -5) * m1 * m2 * mtot;

    // REAL dt = min_timestep * EnzoTimeStep;
    // REAL dtime = next_time - current_time;
    // int iteration_num = static_cast<int>(round(dtime / min_timestep));

    REAL dt = (next_time - current_time) * EnzoTimeStep;

    // for (int num = 0; num < iteration_num; num++) {
    // Current state variables
    REAL e = _bin.ecc;
    REAL semi = _bin.semi;
    
    // Define the derivative functions for `e`, `semi`, and `rot_self`
    auto de_dt = [&](REAL e, REAL semi) {
        REAL e2 = e * e;
        REAL FE2 = e * pow((1 - e2), -2.5) * (1 + 121.0 / 304.0 * e2);
        return 304.0 / 15.0 * cost * pow(semi, -4.0) * FE2;
    };
    auto dsemi_dt = [&](REAL e, REAL semi) {
        REAL e2 = e * e;
        REAL e4 = e2 * e2;
        REAL FE1 = pow((1 - e2), -3.5) * (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e4);
        return 64.0 / 5.0 * cost * pow(semi, -3.0) * FE1;
    };
    auto domega_dt = [&](REAL e, REAL semi) {
        REAL e2 = e * e;
        return (6.0 * M_PI / (c * c * _bin.period) * mtot / (semi * (1 - e2)) +
                3.0 * (18.0 + e2) * M_PI / (2 * pow(c, 4) * _bin.period) * pow((mtot / (semi * (1 - e2))), 2.0));
    };

    // Runge-Kutta 4th Order Method
    // k1
    REAL k1_de = de_dt(e, semi) * dt;
    REAL k1_dsemi = dsemi_dt(e, semi) * dt;
    REAL k1_domega = domega_dt(e, semi) * dt;

    // if (_bin.getLeftMember()->PID == 35) {
    //     fprintf(mergerout, "k1_de: %e\n", k1_de);
    //     fprintf(mergerout, "k1_dsemi: %e\n", k1_dsemi);
    //     fprintf(mergerout, "k1_domega: %e\n", k1_domega);
    //     fflush(mergerout);
    // } 

    // k2
    REAL k2_de = de_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;
    REAL k2_dsemi = dsemi_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;
    REAL k2_domega = domega_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;

    // if (_bin.getLeftMember()->PID == 35) {
    //     fprintf(mergerout, "k2_de: %e\n", k2_de);
    //     fprintf(mergerout, "k2_dsemi: %e\n", k2_dsemi);
    //     fprintf(mergerout, "k2_domega: %e\n", k2_domega);
    //     fflush(mergerout);
    // } 

    // k3
    REAL k3_de = de_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;
    REAL k3_dsemi = dsemi_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;
    REAL k3_domega = domega_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;

    // if (_bin.getLeftMember()->PID == 35) {
    //     fprintf(mergerout, "k3_de: %e\n", k3_de);
    //     fprintf(mergerout, "k3_dsemi: %e\n", k3_dsemi);
    //     fprintf(mergerout, "k3_domega: %e\n", k3_domega);
    //     fflush(mergerout);
    // } 

    // k4
    REAL k4_de = de_dt(e - k3_de, semi - k3_dsemi) * dt;
    REAL k4_dsemi = dsemi_dt(e - k3_de, semi - k3_dsemi) * dt;
    REAL k4_domega = domega_dt(e - k3_de, semi - k3_dsemi) * dt;

    // if (_bin.getLeftMember()->PID == 35) {
    //     fprintf(mergerout, "k4_de: %e\n", k4_de);
    //     fprintf(mergerout, "k4_dsemi: %e\n", k4_dsemi);
    //     fprintf(mergerout, "k4_domega: %e\n", k4_domega);
    //     fflush(mergerout);
    // } 

    // Update variables using weighted sum of Runge-Kutta increments
    _bin.ecc -= (k1_de + 2 * k2_de + 2 * k3_de + k4_de) / 6.0;
    _bin.semi -= (k1_dsemi + 2 * k2_dsemi + 2 * k3_dsemi + k4_dsemi) / 6.0;
    _bin.rot_self += (k1_domega + 2 * k2_domega + 2 * k3_domega + k4_domega) / 6.0;  

    // Check for invalid state
    if (!(_bin.ecc > 0) || !(_bin.semi > 0)) {
        fprintf(mergerout, "GW driven Merger happened! (a < da)\n");
        fprintf(mergerout, "ecc: %e, semi: %e pc, dsemi: %e pc, timestep: %e Myr\n", _bin.ecc, _bin.semi*position_unit, (k1_dsemi + 2 * k2_dsemi + 2 * k3_dsemi + k4_dsemi) / 6.0*position_unit, dt*1e4);
        fflush(mergerout);
        Interaction interaction;
        // _bin_interrupt.time_now = current_time * EnzoTimeStep + dt * num;
        _bin_interrupt.time_now = next_time * EnzoTimeStep;
        interaction.merge(_bin_interrupt, _bin);
        return;
    }

    _bin.calcParticles(REAL(1.0));
    for (int dim = 0; dim < Dim; dim++) {
        _bin.getLeftMember()->Position[dim] += _bin.Position[dim];
        _bin.getLeftMember()->Velocity[dim] += _bin.Velocity[dim];
        _bin.getRightMember()->Position[dim] += _bin.Position[dim];
        _bin.getRightMember()->Velocity[dim] += _bin.Velocity[dim];
    }
}


inline void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, REAL current_time, REAL next_time) {
    if (_bin.getMemberN() == 2) {
        if (_bin.ecc < 1)
            GR_energy_loss(_bin_interrupt, _bin, current_time, next_time);
    }
    else {
        for (int k=0; k<2; k++) {
            if (_bin.isMemberTree(k)) {
                auto memberTree = _bin.getMemberAsTree(k);

                memberTree->calcOrbit(REAL(1.0));
                if (memberTree->ecc < 1)
                    GR_energy_loss_iter(_bin_interrupt, *memberTree, current_time, next_time);
            }
        }
    }
}