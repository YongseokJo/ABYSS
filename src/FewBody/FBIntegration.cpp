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
#include "interaction.h"
#include "perturber.h"




// made 2024.08.12 by Eunwoo Chung

// Reference: SDAR/sample/AR/ar.cxx & PeTar/src/hard.hpp
// No debugging yet


void Group::ARIntegration(REAL next_time){

    sym_int.initialIntegration(CurrentTime*EnzoTimeStep);
    sym_int.integrateToTime(next_time*EnzoTimeStep);

    sym_int.particles.shiftToOriginFrame();
    sym_int.particles.template writeBackMemberAll<Particle>(); // Eunwoo: I'm not sure
    sym_int.particles.shiftToCenterOfMassFrame();

    CurrentTime = next_time;



    
//     // integration loop
//     const int n_particle = sym_int.particles.getSize();
//     if (!synch_flag) {
//         float time_out = time_zero.value + dt_out.value;
//         Float time_table[manager.step.getCDPairSize()];
//         sym_int.profile.step_count = 1;
//         auto IntegrateOneStep = [&] (){
// #if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)
//             sym_int.updateSlowDownAndCorrectEnergy(true, false);
// #endif
//             if(n_particle==2) sym_int.integrateTwoOneStep(sym_int.info.ds, time_table);
//             else sym_int.integrateOneStep(sym_int.info.ds, time_table);
//             if (sym_int.getTime()>=time_out) {
//                 // sym_int.printColumn(std::cout, print_width.value, n_sd);
//                 std::cout<<std::endl;
//                 time_out += dt_out.value;
//             }
//             sym_int.profile.step_count_sum++;
//         };

//         if (nstep.value>0) for (int i=0; i<nstep.value; i++) IntegrateOneStep();
//         else while (sym_int.getTime()<time_end.value) IntegrateOneStep();
//     }
//     else {
//         if (dt_out.value>0.0) nstep.value = int(time_end.value/dt_out.value+0.5);
//         else if (nstep.value>0) dt_out.value = time_end.value/nstep.value;
//         for (int i=1; i<=nstep.value; i++) {
//             auto bin_interrupt = sym_int.integrateToTime(dt_out.value*i);
//             if (bin_interrupt.status!=AR::InterruptStatus::none) {
//                 std::cerr<<"Interrupt condition triggered! ";
//                 Particle* p1 = bin_interrupt.adr->getLeftMember();
//                 Particle* p2 = bin_interrupt.adr->getRightMember();
//                 switch (bin_interrupt.status) {
//                 case AR::InterruptStatus::change:
//                     std::cerr<<" Change";
//                     break;
//                 case AR::InterruptStatus::merge:
//                     std::cerr<<" merge";
//                     break;
//                 case AR::InterruptStatus::destroy:
//                     std::cerr<<" Destroy";
//                     break;
//                 case AR::InterruptStatus::none:
//                     break;
//                 }
//                 std::cerr<<" Time: "<<bin_interrupt.time_now<<std::endl;
//                 bin_interrupt.adr->printColumnTitle(std::cerr);
//                 std::cerr<<std::endl;
//                 bin_interrupt.adr->printColumn(std::cerr);
//                 std::cerr<<std::endl;
//                 // Particle::printColumnTitle(std::cerr);
//                 std::cerr<<std::endl;
//                 for (int j=0; j<2; j++) {
//                     // bin_interrupt.adr->getMember(j)->printColumn(std::cerr);
//                     std::cerr<<std::endl;
//                 }

//                 // merger case, quit integration
//                 if (n_particle==2&&(p1->Mass==0||p2->Mass==0)) {
//                     // sym_int.printColumn(std::cout, print_width.value, n_sd);
//                     std::cout<<std::endl;
//                     break;
//                 }
//             }
//             sym_int.info.generateBinaryTree(sym_int.particles, manager.interaction.gravitational_constant);
//             // sym_int.printColumn(std::cout, print_width.value, n_sd);
//             std::cout<<std::endl;
//         }
//     }


}