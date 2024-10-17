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
// Modify semi-major axis and eccentricity per time step in SDAR.
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
        _bin.ecc -= de;
        _bin.semi -= dsemi;

        num += 1;

        // if (_bin.getMember(0)->PID == 716 || _bin.getMember(0)->PID == 44)
        //     fprintf(mergerout, "dt: %e Myr, de: %e, dsemi %e pc\n", dt*1e4, de, dsemi*position_unit);

        if (_bin.ecc < 0 || _bin.semi < 0) {

            Interaction interaction;

            if (_bin.getMember(0)->ParticleType == (Blackhole+SingleParticle) && _bin.getMember(1)->ParticleType == (Blackhole+SingleParticle)) {
                _bin_interrupt.time_now = current_time*EnzoTimeStep + dt*num; // correct
                interaction.GWmerge(_bin_interrupt, _bin);
                // _bin_interrupt.time_now = _bin_interrupt.time_now - dtime*EnzoTimeStep + dt*num; // test
            }
            else {
                _bin_interrupt.time_now = current_time*EnzoTimeStep + dt*num; // correct
                interaction.TDE(_bin_interrupt, _bin);
                // _bin_interrupt.time_now = _bin_interrupt.time_now - dtime*EnzoTimeStep + dt*num; // test
            }
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

    /*
    // Eliptic case
        if (e < 1) {
            FE = (1.0 + 73./24.*e2 + 37./96.*e4);
            FE1 = pow((1 - e2), -3.5)*FE;
            FE2 = e*pow((1-e2), -2.5)*(1. + 121./304.*e2);

            de = 304./15.*cost*pow(semi, -4.)*FE2*_dtime;
            dsemi = 64./5.*cost*pow(semi, -3.)*FE1*_dtime;
        // fprintf(mergerout, "dt: %e, de: %e, dsemi %e (pc)\n", _dtime*1e4, de, dsemi*position_unit);
        }
    // Hyperbolic case

        else if (e > 1) {
            FE1 = 3.*(96. + 292.*e2 + 37.*e4)*acos(-1/e);
            FE1 = FE1 + (673.*e2 + 602.)*sqrt(e2-1);
            FE1 = 1./M_PI*pow((e2-1), -3.5)*FE1;

            FE2 = 3.*e2*(304. + 121.*e2)*acos(-1/e);
            FE2 = FE2 + (134. + 1069.*e2 + 72.*e4)*sqrt(e2-1);
            FE2 = 1./M_PI/e*pow((e2-1), -2.5)*FE2;

            de = 1./45.*cost*pow(abs(semi), -4.)*FE2*_dtime;
            dsemi = 2./45.*cost*pow(abs(semi), -3.)*FE1*_dtime;
        }
    */
}

inline void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, REAL current_time, REAL next_time) {
    for (int k=0; k<2; k++) {
        if (_bin.isMemberTree(k)) {
            auto memberTree = _bin.getMemberAsTree(k);

            memberTree->calcOrbit(REAL(1.0));
            if (memberTree->ecc < 1)
                GR_energy_loss(_bin_interrupt, *memberTree, current_time, next_time);
        }
    }
}


// Kick BH when BH merger happens
// void GWmergerKick(COMM::BinaryTree<Particle,COMM::Binary> &_bin) {

// }