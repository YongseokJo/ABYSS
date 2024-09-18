#define ASSERT(x) assert(x)

#include "FB_defs.h"

#include <iostream>
#include <fstream>
//#include <unistd.h>
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

void InitializeFBParticle(Particle* FBParticle, std::vector<Particle*> &particle);
void UpdateNextRegTime(std::vector<Particle*> &particle);




// made 2024.08.12 by Eunwoo Chung

// Reference: SDAR/sample/AR/ar.cxx & PeTar/src/hard.hpp
// No debugging yet

// I think it is better to change void function to bool function when we consider the group termination!
// If Intererrupt_mode != none, then bin_termination = true;
void Group::ARIntegration(REAL next_time, std::vector<Particle*> &particle){

    auto bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);

    sym_int.particles.shiftToOriginFrame();
    sym_int.particles.template writeBackMemberAll<Particle>(); // Eunwoo: I'm not sure
    // for (Particle* member : Members) {
    //     fprintf(binout, "PID: %d. posx: %e, posy: %e, posz: %e\n", member->PID, member->Position[0], member->Position[1], member->Position[2]);
    // }
    // sym_int.particles.writeBackMemberAll<Particle>();
    // for (Particle* member : Members) {
    //     fprintf(binout, "PID: %d. posx: %e, posy: %e, posz: %e\n", member->PID, member->Position[0], member->Position[1], member->Position[2]);
    // }

    if (bin_interrupt.status != AR::InterruptStatus::none) {

        // Almost the same process as FBTermination
        // However it is interrupted and we have to predict the position to the next_time!

        groupCM->isErase = true;
        particle.erase(
                std::remove_if(particle.begin(), particle.end(),
                    [](Particle* p) {
                    bool to_remove = p->isErase;
                    //if (to_remove) delete p;
                    return to_remove;
                    }),
                particle.end());

        for (int i=0; i<particle.size(); i++) {
            particle[i]->ParticleOrder = i;
        }
        

        for (int i=0; i<sym_int.particles.getSize(); i++) {
            if (sym_int.particles[i].Mass == 0) {
                fprintf(binout, "PID: %d, Zero mass particle is removed from the simulation!!!\n", sym_int.particles[i].PID);
                delete &sym_int.particles[i];
                sym_int.particles.removeMember(i, true);    // true: shift last member to current position (defaulted);
                                                            // This might call some memory error so should be checked later!!!
                i--;                                            
                continue;
            }
            sym_int.particles[i].ParticleOrder = particle.size();
            particle.push_back(&sym_int.particles[i]);
        }


        for (Particle* member : Members) {
            assert(member->Mass > 0); // Zero mass particles should be removed already!
            InitializeFBParticle(member, particle);

            member->CurrentTimeIrr = CurrentTime;
            member->predictParticleSecondOrderIrr(bin_interrupt.time_now/EnzoTimeStep);
            member->correctParticleFourthOrder(bin_interrupt.time_now/EnzoTimeStep, next_time, member->a_irr); // Eunwoo: Is this really necessary?
            member->calculateTimeStepReg();
            
            // /* Eunwoo: just for a while
            if (member->TimeLevelReg <= groupCM->TimeLevelReg-1 
                    && member->TimeBlockReg/2+member->CurrentBlockReg > groupCM->CurrentBlockIrr+groupCM->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
                member->TimeLevelReg = groupCM->TimeLevelReg-1;
            }
            else if  (member->TimeLevelReg >= groupCM->TimeLevelReg+1) {
                member->TimeLevelReg = groupCM->TimeLevelReg+1;
            }
            else 
                member->TimeLevelReg = groupCM->TimeLevelReg;
            member->TimeStepReg  = static_cast<REAL>(pow(2, member->TimeLevelReg)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()
            member->TimeBlockReg = static_cast<ULL>(pow(2, member->TimeLevelReg-time_block)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()

            member->calculateTimeStepIrr(member->a_tot, member->a_irr);
        }

        int index = 0;
        for (Particle* ptcl: particle) {

            index = 0;
            for (Particle* neighbor: ptcl->ACList) {
                if (neighbor->PID == groupCM->PID) {
                    ptcl->ACList.erase(ptcl->ACList.begin() + index);
                    ptcl->ACList.insert(ptcl->ACList.end(), Members.begin(), Members.end());
                    ptcl->NumberOfAC = ptcl->ACList.size();
                    break; // Eunwoo check
                }
                index++;
            }
        }


        UpdateNextRegTime(particle);

        for (Particle* member : Members) {
            member->isErase		= false;
            member->isGroup		= false;
            // member->GroupInfo	= nullptr;
        }

        // // we also need to delete ptclGroup from the group list
        // // fprintf(binout,"deleting binary information from the GroupList \n");
        // this->isErase = true;
        // GroupList.erase(
        //         std::remove_if(GroupList.begin(), GroupList.end(),
        //             [](Group* p) {
        //             bool to_remove = p->isErase;
        //             //if (to_remove) delete p;
        //             return to_remove;
        //             }),
        //         GroupList.end());

        delete this;
        // delete ptclCM; // It is deleted automatically when delete ptclGroup!!!
        // this = nullptr;
        groupCM  = nullptr;

        bin_termination = true;

        return;
    }

    sym_int.particles.shiftToCenterOfMassFrame();
    CurrentTime = next_time;
    // fprintf(binout, "AR done! time: %e\n", next_time*EnzoTimeStep*1e4);
}