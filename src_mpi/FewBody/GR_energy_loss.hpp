#ifndef GR_ENERGY_LOSS
#define GR_ENERGY_LOSS
#ifdef FEWBODY
#include <cmath>

#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "../particle.h"
#include "ar_interaction.hpp"

// extern REAL min_timestep;


// made 2024.09.19 by Eunwoo Chung

// reference: tides3.f from Nbody6++GPU (Rizzuto et al. 2020, Arca Sedda et al. 2023)
// Gravitational wave energy loss of hard binaries according to the orbit averaged approximation of Peters & Mathews 1963.
// Calculate average change of energy and angular momentum per orbit.
// Modify semi-major axis, eccentricity, omega(argument of periapsis) per time step in SDAR.
// Orbit shrinking by PN2.5
// Precession by PN1.0 & PN2.0

inline void GR_energy_loss(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time) {

    const double c = 299752.458 / (velocity_unit / yr * pc / 1e5); // speed of light in code unit
    const double m1 = _bin.m1;
    const double m2 = _bin.m2;
    const double mtot = m1 + m2;
    const double cost = pow(c, -5) * m1 * m2 * mtot;

    double dt = (next_time - current_time) * EnzoTimeStep;

    double e = _bin.ecc;
    double semi = _bin.semi;
    
    // Define the derivative functions for `e`, `semi`, and `rot_self`
    auto de_dt = [&](double e, double semi) {
        double e2 = e * e;
        double FE2 = e * pow((1 - e2), -2.5) * (1 + 121.0 / 304.0 * e2);
        return 304.0 / 15.0 * cost * pow(semi, -4.0) * FE2;
    };
    auto dsemi_dt = [&](double e, double semi) {
        double e2 = e * e;
        double e4 = e2 * e2;
        double FE1 = pow((1 - e2), -3.5) * (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e4);
        return 64.0 / 5.0 * cost * pow(semi, -3.0) * FE1;
    };
    auto domega_dt = [&](double e, double semi) {
        double e2 = e * e;
        return (6.0 * M_PI / (c * c * _bin.period) * mtot / (semi * (1 - e2)) +
                3.0 * (18.0 + e2) * M_PI / (2 * pow(c, 4) * _bin.period) * pow((mtot / (semi * (1 - e2))), 2.0));
    };

    // Runge-Kutta 4th Order Method
    // k1
    double k1_de = de_dt(e, semi) * dt;
    double k1_dsemi = dsemi_dt(e, semi) * dt;
    double k1_domega = domega_dt(e, semi) * dt;

    // k2
    double k2_de = de_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;
    double k2_dsemi = dsemi_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;
    double k2_domega = domega_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;

    // k3
    double k3_de = de_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;
    double k3_dsemi = dsemi_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;
    double k3_domega = domega_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;

    // k4
    double k4_de = de_dt(e - k3_de, semi - k3_dsemi) * dt;
    double k4_dsemi = dsemi_dt(e - k3_de, semi - k3_dsemi) * dt;
    double k4_domega = domega_dt(e - k3_de, semi - k3_dsemi) * dt;

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
        // interaction.merge(_bin_interrupt, _bin);

        auto* p1 = _bin.getLeftMember();
        auto* p2 = _bin.getRightMember();

        p1->setBinaryInterruptState(BinaryInterruptState::collision);
        p2->setBinaryInterruptState(BinaryInterruptState::collision);
        p1->setBinaryPairID(p2->ParticleIndex);
        p2->setBinaryPairID(p1->ParticleIndex);
        _bin_interrupt.status = AR::InterruptStatus::merge;
        _bin_interrupt.adr = &_bin;
        return;
    }

    _bin.calcParticles(double(1.0));
    for (int dim = 0; dim < Dim; dim++) {
        _bin.getLeftMember()->Position[dim] += _bin.Position[dim];
        _bin.getLeftMember()->Velocity[dim] += _bin.Velocity[dim];
        _bin.getRightMember()->Position[dim] += _bin.Position[dim];
        _bin.getRightMember()->Velocity[dim] += _bin.Velocity[dim];
    }
}


inline void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time) {
    if (_bin.getMemberN() == 2) {
        if (_bin.ecc < 1)
            GR_energy_loss(_bin_interrupt, _bin, current_time, next_time);
    }
    else {
        for (int k=0; k<2; k++) {
            if (_bin.isMemberTree(k)) {
                auto memberTree = _bin.getMemberAsTree(k);

                memberTree->calcOrbit(double(1.0));
                if (memberTree->ecc < 1)
                    GR_energy_loss_iter(_bin_interrupt, *memberTree, current_time, next_time);
            }
        }
    }
}
#endif
#endif