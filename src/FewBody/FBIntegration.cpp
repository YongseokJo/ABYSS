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


    auto bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);


    if (PNon && bin_interrupt.status == AR::InterruptStatus::none) {

        auto& bin_root = sym_int.info.getBinaryTreeRoot();

        if (bin_root.getMemberN() == 2) {
            bin_root.calcOrbit(REAL(1.0));
            if (bin_root.ecc < 1) {
                GR_energy_loss(bin_interrupt, bin_root, CurrentTime, next_time);
                if (bin_root.getMember(0)->PID == 716 && bin_root.getMember(1)->PID == 44) {
                // if (bin_root.getMember(0)->PID == 45 && bin_root.getMember(1)->PID == 44) {
                    auto* p1 = bin_root.getLeftMember();
                    auto* p2 = bin_root.getRightMember();
                    fprintf(mergerout, "Time: %e Myr\n", next_time*EnzoTimeStep*1e4);
                    fprintf(mergerout, "PID: %d. x: %e, y: %e, z: %e\n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
                    fprintf(mergerout, "PID: %d. x: %e, y: %e, z: %e\n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
                }
            }
        }
        else {
            GR_energy_loss_iter(bin_interrupt, bin_root, CurrentTime, next_time);
        }
// /* // correct
        if (bin_interrupt.status == AR::InterruptStatus::none)
            sym_int.initialIntegration(next_time*EnzoTimeStep);
// */ // correct
        // sym_int.initialIntegration(next_time*EnzoTimeStep); // test
    }    


    if (bin_interrupt.status != AR::InterruptStatus::none) {

        if (Members.size() == 2) {

            // groupCM->CurrentTimeIrr = CurrentTime;

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
            int j = 0;
            for (int i=0; i<sym_int.particles.getSize(); i++) {

                Particle* member = &sym_int.particles[i];

                if (member->Mass == 0) {
                    delete member;
                    sym_int.particles.removeMember(i, false);
                    assert(Members[j]->Mass == 0);
                    Members.erase(Members.begin() + j);
                    continue;
                }
                j++;
            }

            this->ARIntegration(next_time, particle);

            groupCM->updateParticle();
            CurrentTime = next_time;
            fflush(mergerout);
            return true;
        }
    }

    CurrentTime = next_time;
    return true; 
}