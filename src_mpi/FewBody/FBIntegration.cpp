#define ASSERT(x) assert(x)

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "stdio.h"
#include <vector>
#include <algorithm>
#include "../global.h"
#include "../def.h"


#include "Common/binary_tree.h"
#include "Common/Float.h"
#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "ar_interaction.hpp"
#include "ar_perturber.hpp"
#include "GR_energy_loss.hpp"

// #define SEVN
#ifdef SEVN
void UpdateEvolution(Particle* ptcl);
#endif

void GR_energy_loss(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time);
void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time);

void NewFBInitialization3(Group* group);


// made 2024.08.12 by Eunwoo Chung

// Reference: SDAR/sample/AR/ar.cxx & PeTar/src/hard.hpp

// true: integrated normally, false: terminated by stellar merger, TDE, GW merger, etc.
// If Intererrupt_mode != none, then bin_termination = true;
bool Group::ARIntegration(double next_time){

    // fprintf(stderr, "PID: %d\n", groupCM->PID);
    // if (groupCM->PID == -1035)
    //     fprintf(stderr, "NumOfAC: %d\n", groupCM->NumberOfAC);

    // fprintf(mergerout, "PID: %d\n", groupCM->PID);
    // fprintf(mergerout, "before posx: %e pc\n", sym_int.particles[0].Position[0]*position_unit);

    auto bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);

    // AR::InterruptBinary<Particle> bin_interrupt;
    // if (groupCM->PID != -1025)
    //     bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);
    // else {
    //     auto& bin_root = sym_int.info.getBinaryTreeRoot();
    //     bin_root.calcOrbit(REAL(1.0));

    //     bin_root.evolve((next_time - CurrentTime)*EnzoTimeStep);
    //     bin_root.calcParticles(REAL(1.0));
    // }


    // if (groupCM->PID == -1035) {
    //     // fprintf(mergerout, "current_time: %e Myr\n", CurrentTime*EnzoTimeStep*1e4);
    //     // fprintf(mergerout, "next_time: %e Myr\n", next_time*EnzoTimeStep*1e4);
    //     fprintf(mergerout, "Time: %e Myr\n", next_time*EnzoTimeStep*1e4);
    //     fprintf(mergerout, "Time step: %e Myr\n", (next_time - CurrentTime)*1e4);
    //     // fprintf(mergerout, "bin_interrupt.time_now: %e Myr\n", bin_interrupt.time_now*1e4);
    //     // fprintf(mergerout, "bin_interrupt.time_end: %e Myr\n", bin_interrupt.time_end*1e4);
    //     fprintf(mergerout, "PID: %d. x: %e, y: %e, z: %e\n", sym_int.particles[0].PID, sym_int.particles[0].Position[0]*position_unit, sym_int.particles[0].Position[1]*position_unit, sym_int.particles[0].Position[2]*position_unit);
    //     fprintf(mergerout, "PID: %d. x: %e, y: %e, z: %e\n\n", sym_int.particles[1].PID, sym_int.particles[1].Position[0]*position_unit, sym_int.particles[1].Position[1]*position_unit, sym_int.particles[1].Position[2]*position_unit);
    //     fflush(mergerout);
    // }

    // if (next_time*EnzoTimeStep == bin_interrupt.time_now) {
    //     fprintf(mergerout, "Same. PID: %d\n", groupCM->PID);
    /*
    fprintf(mergerout, "after posx: %e pc\n", sym_int.particles[0].Position[0]);
    fprintf(mergerout, "current time: %e\n", CurrentTime*EnzoTimeStep);
    fprintf(mergerout, "next_time: %e\n", next_time*EnzoTimeStep);
    fprintf(mergerout, "bin_interrupt.time_now: %e\n", bin_interrupt.time_now);
    fprintf(mergerout, "bin_interrupt.time_end: %e\n", bin_interrupt.time_end);
    fprintf(mergerout, "\n");
    fflush(mergerout);
    */
    // }
    // else {
    //     fprintf(mergerout, "Diff. PID: %d\n", groupCM->PID);
    //     fflush(mergerout);
    // }

/* // Precession & Kepler solver test with constant dt
    AR::InterruptBinary<Particle> bin_interrupt;
    if (groupCM->PID == -1025) {
        REAL time = 0;
        fprintf(mergerout, "Iteration num: %d\n", static_cast<int>(round((next_time - CurrentTime)/min_timestep)));
        while (time <= (next_time - CurrentTime)) {
            bin_interrupt.status = AR::InterruptStatus::none;

            auto& bin_root = sym_int.info.getBinaryTreeRoot();
            bin_root.calcOrbit(REAL(1.0));

            bin_root.evolve(min_timestep*EnzoTimeStep);
            bin_root.calcParticles(REAL(1.0));
            GR_energy_loss(bin_interrupt, bin_root, 0, min_timestep);

            auto* p1 = bin_root.getLeftMember();
            auto* p2 = bin_root.getRightMember();
            fprintf(mergerout, "Time: %e Myr\n", (CurrentTime + time)*EnzoTimeStep*1e4);
            fprintf(mergerout, "PID: %d. x: %e, y: %e, z: %e\n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            fprintf(mergerout, "PID: %d. x: %e, y: %e, z: %e\n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
            fprintf(mergerout, "Ecc: %e\n", bin_root.ecc);
            fprintf(mergerout, "Semi: %e pc\n", bin_root.semi*position_unit);
            fprintf(mergerout, "omega: %e rad\n", bin_root.rot_self);
            fflush(mergerout);

            time += min_timestep;
        }
        fprintf(mergerout, "Done!\n");
        fflush(mergerout);
    }
    else
        bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);
*/ // Precession test

/* // Widely used version, but let's just use only SDAR, not Kepler integration
    AR::InterruptBinary<Particle> bin_interrupt;
    if (Members.size() == 2 && groupCM->NumberOfAC == 0) { // for unperturbed binary
    // if (Members.size() == 2) { // test

        bin_interrupt.status = AR::InterruptStatus::none;

        auto& bin_root = sym_int.info.getBinaryTreeRoot();
        bin_root.calcOrbit(REAL(1.0));

        bin_root.evolve((next_time - CurrentTime)*EnzoTimeStep);
        bin_root.calcParticles(REAL(1.0));
        bin_interrupt.time_now = next_time*EnzoTimeStep;

        if (manager.interrupt_detection_option > 0) {
            Interaction interaction;
            interaction.modifyAndInterruptKepler(bin_interrupt, bin_root, (next_time - CurrentTime)*EnzoTimeStep);
        }
        if (bin_interrupt.status == AR::InterruptStatus::none) {
        // if (!(sym_int.info.getBinaryTreeRoot().Velocity[0]*sym_int.info.getBinaryTreeRoot().Velocity[0]<1e-10)) { 
        //     fprintf(stderr, "Kepler\n");
        //     fprintf(stderr, "pos: (%e, %e, %e)\n", sym_int.info.getBinaryTreeRoot().Position[0], sym_int.info.getBinaryTreeRoot().Position[1], sym_int.info.getBinaryTreeRoot().Position[2]);
        //     fprintf(stderr, "vel: (%e, %e, %e)\n", sym_int.info.getBinaryTreeRoot().Velocity[0], sym_int.info.getBinaryTreeRoot().Velocity[1], sym_int.info.getBinaryTreeRoot().Velocity[2]);
        //     for (Particle* members: Members) {
        //         fprintf(stderr, "PID: %d\n", members->PID);
        //     }
        //     fflush(stderr);
        // }
            // sym_int.particles.calcCenterOfMass();        
            sym_int.initialIntegration(next_time*EnzoTimeStep);
        }
    }
    else {
        bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);
    }
*/

// /* PN corrections
    // if (PNon && bin_interrupt.status == AR::InterruptStatus::none) { // Only particles with BH
    if (bin_interrupt.status == AR::InterruptStatus::none) { // Every bound orbit
        
        auto& bin_root = sym_int.info.getBinaryTreeRoot();
        // if (groupCM->PID == -1025) {
        //     fprintf(mergerout, "next_time: %e, real_time: %e\n", next_time*EnzoTimeStep*1e4, bin_interrupt.time_now*1e4);
        //     fflush(mergerout);
        // }

        if (bin_root.getMemberN() == 2) {
            bin_root.calcOrbit(double(1.0));
            if (bin_root.ecc < 1) {
                GR_energy_loss(bin_interrupt, bin_root, CurrentTime, next_time);

                // if (groupCM->PID == -1035) {
                //     auto* p1 = bin_root.getLeftMember();
                //     auto* p2 = bin_root.getRightMember();
                //     fprintf(mergerout, "Time: %e Myr\n", next_time*EnzoTimeStep*1e4);
                //     fprintf(mergerout, "PID: %d. x: %e, y: %e, z: %e\n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
                //     fprintf(mergerout, "PID: %d. x: %e, y: %e, z: %e\n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
                //     fprintf(mergerout, "Ecc: %e\n", bin_root.ecc);
                //     fprintf(mergerout, "Semi: %e pc\n", bin_root.semi*position_unit);
                //     fprintf(mergerout, "omega: %e rad\n", bin_root.rot_self);
                //     fflush(mergerout);
                // }     
            }
        }
        else {
            GR_energy_loss_iter(bin_interrupt, bin_root, CurrentTime, next_time);
        }
        if (bin_interrupt.status == AR::InterruptStatus::none) {
            // if (!(sym_int.info.getBinaryTreeRoot().Velocity[0]*sym_int.info.getBinaryTreeRoot().Velocity[0]<1e-10)) { 
            //     fprintf(stderr, "PN correction\n");
            //     // fprintf(stderr, "pos: (%e, %e, %e)\n", sym_int.info.getBinaryTreeRoot().Position[0], sym_int.info.getBinaryTreeRoot().Position[1], sym_int.info.getBinaryTreeRoot().Position[2]);
            //     // fprintf(stderr, "vel: (%e, %e, %e)\n", sym_int.info.getBinaryTreeRoot().Velocity[0], sym_int.info.getBinaryTreeRoot().Velocity[1], sym_int.info.getBinaryTreeRoot().Velocity[2]);
            //     for (Particle* members: Members) {
            //         fprintf(stderr, "PID: %d\n", members->PID);
            //         // fprintf(stderr, "pos: %e, %e, %e\n", members->Position[0], members->Position[1], members->Position[2]);
            //         // fprintf(stderr, "vel: %e, %e, %e\n", members->Velocity[0], members->Velocity[1], members->Velocity[2]);
            //     }
            //     fflush(stderr);
            // }
            // if (groupCM->PID == -10917) {
            //     fprintf(stderr, "CurrentTime: %e Myr\n", next_time*EnzoTimeStep*1e4);
            //     fprintf(stderr, "treepos: %e, %e, %e\n", bin_root.Position[0], bin_root.Position[1], bin_root.Position[2]);
            //     fprintf(stderr, "treevel: %e, %e, %e\n", bin_root.Velocity[0], bin_root.Velocity[1], bin_root.Velocity[2]);
            //     fflush(stderr);
            // }
            sym_int.initialIntegration(next_time*EnzoTimeStep); // Eunwoo: this should be fixed later // Eunwoo: I don't think so!
        }
    }    
// */

    if (bin_interrupt.status != AR::InterruptStatus::none) {

        isMerger = true;

        if (sym_int.particles.getSize() == 2) {

            double pos[Dim], vel[Dim];

            groupCM->predictParticleSecondOrder(bin_interrupt.time_now/EnzoTimeStep - CurrentTime, pos, vel);

            for (int dim=0; dim<Dim; dim++) {
                groupCM->Position[dim] = pos[dim];
                groupCM->Velocity[dim] = vel[dim];
                sym_int.particles.cm.Position[dim] = pos[dim];
                sym_int.particles.cm.Velocity[dim] = vel[dim];
            }
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<Particle>();

            isTerminate = true;

            CurrentTime = bin_interrupt.time_now/EnzoTimeStep;

            return false;   
        }
        else {

            CurrentTime = bin_interrupt.time_now/EnzoTimeStep;

            for (int dim=0; dim<Dim; dim++) {
                sym_int.particles.cm.Position[dim] = groupCM->Position[dim];
                sym_int.particles.cm.Velocity[dim] = groupCM->Velocity[dim];
            }
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<Particle>();

            for (int i = 0; i < sym_int.particles.getSize(); i++) {
                if (sym_int.particles[i].Mass != 0)
                    sym_int.particles[i].CurrentTimeIrr  = next_time;
            }

            NewFBInitialization3(this);
            return false; // Eunwoo: return true makes an SIGMENTATION FAULT error when doing checkBreak!
        }
    }
#ifdef SEVN

    if (useSEVN) {

        bool breakEvolution = false;
        bool evolved = false;
        while (EvolutionTime + EvolutionTimeStep < next_time*EnzoTimeStep*1e4) {

            evolved = true;

            // if (groupCM->PID == -10917) {
            //     fprintf(stderr, "Stellar evolution!\n");
            //     fflush(stderr);
            // }

            EvolutionTime += EvolutionTimeStep;

            REAL dt_evolve_next = NUMERIC_FLOAT_MAX; // Myr

            for (int i=0; i < sym_int.particles.getSize(); i++) {
                Particle* members = &sym_int.particles[i];
                if (members->star == nullptr || members->star->amiremnant())
                    continue;
                members->star->sync_with(EvolutionTimeStep);
                members->EvolutionTime += members->star->getp(Timestep::ID);
                // fprintf(SEVNout, "INT. PID: %d. GT: %e, ET: %e\n", members->PID, CurrentTime*EnzoTimeStep*1e4, members->EvolutionTime);
                // fflush(SEVNout);
                members->star->evolve();
                // fprintf(SEVNout, "INT. After. PID: %d, ET: %e, WT: %e\n", members->PID, members->EvolutionTime, members->star->getp(Worldtime::ID));
                // fflush(SEVNout);
                UpdateEvolution(members);
                // fprintf(SEVNout, "Int. PID: %d, Mass: %e, Radius: %e, EvolutionTime: %e\n", members->PID, members->Mass*mass_unit, members->radius*position_unit, members->EvolutionTime);
                // fflush(SEVNout);
                if (members->star->amiempty() || members->star->vkick[3] > 0.0) {
                    breakEvolution = true;
                    break;
                }
                if (members->star->amiBH())
                    PNon = true;
                if (!members->star->amiremnant() && members->star->getp(Timestep::ID) < dt_evolve_next)
                    dt_evolve_next = members->star->getp(Timestep::ID);
            }
            EvolutionTimeStep = dt_evolve_next;
            // fprintf(SEVNout, "INT. EvolutionTimeStep: %e\n", EvolutionTimeStep);
            // fflush(SEVNout);
        }

        useSEVN = false;
        for (int i=0; i < sym_int.particles.getSize(); i++) {
            Particle* members = &sym_int.particles[i];
            if (members->star != nullptr && !members->star->amiremnant()) {
                useSEVN = true;
                break;
            }
        }
        if (breakEvolution) { // break SDAR
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<Particle>();

            FBTermination2(groupCM, next_time, particle);
            return false;
        }
        if (evolved) { // Eunwoo: orbital parameters should be re-calculated due to mass changes during stellar evolution!
            sym_int.particles.shiftToOriginFrame();
            sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);
            sym_int.initialIntegration(next_time*EnzoTimeStep);
        }
    }
    
#endif

    CurrentTime = next_time;
    return true; 
}