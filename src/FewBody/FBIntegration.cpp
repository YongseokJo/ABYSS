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

	// for (Particle* ptcl : this->Members) {
    //     fprintf(stdout, "Time: %e\n", this->CurrentTime*EnzoTimeStep);
	// 	fprintf(stdout, "PID: %d, Position - x:%e, y:%e, z:%e\n", ptcl->PID, ptcl->Position[0], ptcl->Position[1], ptcl->Position[2]);
    //     fprintf(stdout, "PID: %d, Velocity - x:%e, y:%e, z:%e\n", ptcl->PID, ptcl->Velocity[0], ptcl->Velocity[1], ptcl->Velocity[2]);
    //     fprintf(stdout, "PID: %d, Mass: %e \n", ptcl->PID, ptcl->Mass);
	// }

    // fprintf(binout, "AR Integration started!\n"); // Eunwoo debug
    // fprintf(binout, "GroupCM PID: %d, current time: %e\n", this->groupCM->PID, this->CurrentTime); // Eunwoo debug

    AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, AR::Information<Particle,Particle>> sym_int;
    AR::TimeTransformedSymplecticManager<Interaction> manager;

    // COMM::IOParamsContainer input_par_store;

    // COMM::IOParams<double> gravitational_constant   (input_par_store, 1.0, "gravitational constant"); // gravitational constant

    manager.interaction.gravitational_constant = 1.0;
    manager.time_step_min = 1e-13; // minimum physical time step // reference: ar.cxx
    manager.ds_scale = 1.0; // step size scaling factor // reference: ar.cxx
    manager.time_error_max = 0.25*1e-13; // time synchronization absolute error limit for AR, default is 0.25*dt-min
    // reference: ar.cxx
    manager.energy_error_relative_max = 1e-10; // relative energy error limit for AR, phase error requirement
    // reference: ar.cxx
    // 1e-8 in PeTar
    manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX; // maximum timescale for maximum slowdown factor, time-end
    // if (slowdown_timescale_max.value>0.0) manager.slowdown_timescale_max = slowdown_timescale_max.value;
    // else if (time_end.value>0.0) manager.slowdown_timescale_max = time_end.value;
    // else manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX;
    // should be positive
    manager.slowdown_pert_ratio_ref = 1e-6; // slowdown perturbation ratio reference
    // reference: ar.cxx
    // 1e-4 in PeTar
    manager.step_count_max = 1000000; // number of maximum (integrate/output) step for AR integration // set symplectic order
    // 1000000 in PeTar & ar.cxx
    manager.step.initialSymplecticCofficients(-6); // Symplectic integrator order, should be even number
    // -6 in PeTar & ar.cxx
    manager.interrupt_detection_option = 0; // modify orbit or check interruption using modifyAndInterruptIter function
                                            // 0: turn off
                                            // 1: modify the binary orbits based on detetion criterion
                                            // 2. modify and also interrupt integrations

    // Eunwoo: it is turned off now but I will turn it on later.
    // Eunwoo: It can be used for merging star (dr < sum of radius) or destroy.

    sym_int.info.r_break_crit = 1e-3/position_unit; // distance criterion for checking stability
    // more information in symplectic_integrator.h
    // Eunwoo set this to 1e-3 pc. reference: ar.cxx
    // check whether the system is stable for 10000 out period and the apo-center is below break criterion
    // PeTar (hard.hpp): sym_int.info.r_break_crit = std::max(sym_int.info.r_break_crit,ptcl_origin[i].getRGroup());
    int fix_step_option = 1; // Eunwoo set this it 'later'.

    //! Fix step options for integration with adjusted step (not for time sychronizatio phase)
    /*! always: use the given step without change
        later: fix step after a few adjustment of initial steps due to energy error
        none: don't fix step
     */

    double s = 0.0; // step size, not physical time step, auto
    // Eunwoo: automatically set during the calculation.
    // bool synch_flag=false; // if true, switch on time synchronization
    // Eunwoo: it seems unnecessary now.


    
    // integrator

    
    sym_int.manager = &manager;

    sym_int.particles.setMode(COMM::ListMode::copy);
    sym_int.particles.reserveMem(Members.size());

    for (size_t i = 0; i < Members.size(); ++i) {
        sym_int.particles.addMemberAndAddress(*this->Members[i]);
    }

    sym_int.particles.calcCenterOfMass();
    sym_int.reserveIntegratorMem();

#ifdef AR_SLOWDOWN_MASSRATIO
    Float m_ave = sym_int.particles.cm.mass/sym_int.particles.getSize();
    if (slowdown_mass_ref.value<=0.0) manager.slowdown_mass_ref = m_ave;
    else manager.slowdown_mass_ref = slowdown_mass_ref.value;
#endif
    // manager.print(std::cerr); // Eunwoo deleted

    sym_int.info.reserveMem(sym_int.particles.getSize());
    sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);


    // initialization 
    sym_int.initialIntegration(this->CurrentTime*EnzoTimeStep); // Eunwoo: Is this the real current time?
    sym_int.info.calcDsAndStepOption(manager.step.getOrder(), manager.interaction.gravitational_constant, manager.ds_scale);


    // use input fix step option
    if (fix_step_option>=0) {
        switch (fix_step_option) {
        case 2:
            sym_int.info.fix_step_option = AR::FixStepOption::none;
            break;
        case 0:
            sym_int.info.fix_step_option = AR::FixStepOption::always;
            break;
        case 1:
            sym_int.info.fix_step_option = AR::FixStepOption::later;
            break;
        }
    }

    // use input ds
    if (s>0.0) sym_int.info.ds = s;

    // precision
    // std::cout<<std::setprecision(print_precision.value);

#ifdef AR_SLOWDOWN_ARRAY
    int n_sd = sym_int.binary_slowdown.getSize();
#elif defined(AR_SLOWDOWN_TREE)
    int n_sd = sym_int.info.binarytree.getSize();
#else
    int n_sd = 0;
#endif
    //print column title
    // sym_int.printColumnTitle(std::cout, print_width.value, n_sd);
    // std::cout<<std::endl;

    //print initial data
    // sym_int.printColumn(std::cout, print_width.value, n_sd);
    // std::cout<<std::endl;

    sym_int.integrateToTime(next_time*EnzoTimeStep); // Eunwoo: sym_int.updateSlowDownAndCorrectEnergy is already built-in.
    sym_int.particles.shiftToOriginFrame(); // Eunwoo: is this the reason?
    sym_int.particles.template writeBackMemberAll<Particle>(); // Eunwoo: I'm not sure

    this->CurrentTime = next_time; // Eunwoo added but I'm not sure...

    // for (Particle* ptcl : this->Members) {
    //     fprintf(stdout, "Time: %e\n", this->CurrentTime*EnzoTimeStep);
	// 	fprintf(stdout, "PID: %d, Position - x:%e, y:%e, z:%e\n", ptcl->PID, ptcl->Position[0], ptcl->Position[1], ptcl->Position[2]);
	// 	fprintf(stdout, "PID: %d, Velocity - x:%e, y:%e, z:%e\n", ptcl->PID, ptcl->Velocity[0], ptcl->Velocity[1], ptcl->Velocity[2]);
    //     fprintf(stdout, "PID: %d, Mass: %e \n", ptcl->PID, ptcl->Mass);
    // }

    sym_int.clear(); // Delocate memory

    // fprintf(binout, "GroupCM PID: %d, current time: %e\n", this->groupCM->PID, this->CurrentTime); // Eunwoo debug
    // fprintf(binout, "AR Integration is done!\n"); // Eunwoo debug

    
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