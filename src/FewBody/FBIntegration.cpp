#define ASSERT(x) assert(x)

#include "FB_defs.h"

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
#include "../defs.h"


#include "Common/binary_tree.h"
#include "Common/Float.h"
#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "ar_interaction.hpp"
#include "ar_perturber.hpp"
#include "GR_energy_loss.hpp"


void GR_energy_loss(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, REAL current_time, REAL next_time);
void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, REAL current_time, REAL next_time);

void FBTermination2(Particle* ptclCM, REAL current_time, std::vector<Particle*> &particle);


// made 2024.08.12 by Eunwoo Chung

// Reference: SDAR/sample/AR/ar.cxx & PeTar/src/hard.hpp

// true: integrated normally, false: terminated by stellar merger, TDE, GW merger, etc.
// If Intererrupt_mode != none, then bin_termination = true;
bool Group::ARIntegration(REAL next_time, std::vector<Particle*> &particle){

    // fprintf(stderr, "PID: %d\n", groupCM->PID);
    // if (groupCM->PID == -1035)
    //     fprintf(stderr, "NumOfAC: %d\n", groupCM->NumberOfAC);

    // fprintf(mergerout, "PID: %d\n", groupCM->PID);
    // fprintf(mergerout, "before posx: %e pc\n", sym_int.particles[0].Position[0]*position_unit);

    // auto bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);

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
        sym_int.initialIntegration(next_time*EnzoTimeStep);
    }
    else {
        bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);
    }

// /* PN corrections
    // if (PNon && bin_interrupt.status == AR::InterruptStatus::none) { // Only particles with BH
    if (bin_interrupt.status == AR::InterruptStatus::none) { // Every bound orbit
        
        auto& bin_root = sym_int.info.getBinaryTreeRoot();
        // if (groupCM->PID == -1025) {
        //     fprintf(mergerout, "next_time: %e, real_time: %e\n", next_time*EnzoTimeStep*1e4, bin_interrupt.time_now*1e4);
        //     fflush(mergerout);
        // }

        if (bin_root.getMemberN() == 2) {
            bin_root.calcOrbit(REAL(1.0));
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
        if (bin_interrupt.status == AR::InterruptStatus::none)
            sym_int.initialIntegration(next_time*EnzoTimeStep); // Eunwoo: this should be fixed later // Eunwoo: I don't think so!
    }    
// */

    if (bin_interrupt.status != AR::InterruptStatus::none) {

        if (Members.size() == 2) {

            groupCM->predictParticleSecondOrderIrr(bin_interrupt.time_now/EnzoTimeStep);
            for (int dim=0; dim<Dim; dim++) {
                groupCM->Position[dim] = groupCM->PredPosition[dim];
                groupCM->Velocity[dim] = groupCM->PredVelocity[dim];
            }

            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<Particle>();

            FBTermination2(groupCM, bin_interrupt.time_now/EnzoTimeStep, particle);
            fflush(mergerout);
            return false;   
        }
        else {
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<Particle>();

            Group *newGroup = new Group();

            for (int i = 0; i < Members.size(); i++) {
                if (Members[i]->Mass == 0)
                    Members[i]->isErase = true;
                else
                    newGroup->Members.push_back(Members[i]);
            }

            Members.erase(
                std::remove_if(
                    Members.begin(), Members.end(),
                    [](Particle* p) {
                        bool to_remove = p->isErase;
                        if (to_remove) delete p; // Delete the memory of zero mass particles.
                        return to_remove;
                    }
                ),
                Members.end()
            );

            newGroup->initialManager();
            newGroup->initialIntegrator();
            newGroup->sym_int.particles.cm = *groupCM;
            newGroup->sym_int.initialIntegration(bin_interrupt.time_now);
            newGroup->sym_int.info.calcDsAndStepOption(newGroup->manager.step.getOrder(), newGroup->manager.interaction.gravitational_constant, newGroup->manager.ds_scale);

            if (bin_interrupt.time_now != next_time*EnzoTimeStep)
                newGroup->ARIntegration(next_time*EnzoTimeStep, particle);
            // /* // This is correct
            auto& bin_root = sym_int.info.getBinaryTreeRoot();
            if (bin_root.semi>0.0)
                newGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, groupCM->RadiusOfAC);
            else
                newGroup->sym_int.info.r_break_crit = 2*bin_root.semi*(1-bin_root.ecc); // r_break_crit = 2*peri
            // */ // This is correct
            CurrentTime = next_time;
            sym_int = newGroup->sym_int;
            manager = newGroup->manager;
            return true;
        }
    }
    CurrentTime = next_time;
    return true; 
}