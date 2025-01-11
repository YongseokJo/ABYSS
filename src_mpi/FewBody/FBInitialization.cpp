#ifdef FEWBODY
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../def.h"

#ifdef SEVN
void UpdateEvolution(Particle* ptcl);
#endif

// REAL getNewTimeStepIrr(REAL f[3][4], REAL df[3][4]);
// void getBlockTimeStep(REAL dt, int& TimeLevel, ULL &TimeBlock, REAL &TimeStep);
// bool UpdateComputationChain(Particle* ptcl);
// void UpdateNextRegTime(std::vector<Particle*> &particle);
// bool CreateComputationList(Particle* ptcl);
// bool CreateComputationChain(std::vector<Particle*> &particle);
void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);

void deleteGroup(Particle* ptclCM) {

	Group* ptclGroup = ptclCM->GroupInfo;

	ptclCM->NewNumberOfNeighbor = ptclGroup->sym_int.particles.getSize();

	assert(!ptclGroup->sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.4

	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		ptclCM->NewNeighbors[i] = members->ParticleIndex;

		for (int dim=0; dim<Dim; dim++) {
			particles[members->ParticleIndex].Position[dim] = ptclCM->Position[dim] + members->Position[dim];
			particles[members->ParticleIndex].Velocity[dim] = ptclCM->Velocity[dim] + members->Velocity[dim];
		}
		particles[members->ParticleIndex].Mass = members->Mass;
	}

/* // original version; this might take so much time by EW 2025.1.4
	for (int dim=0; dim<Dim; dim++) {
		ptclGroup->sym_int.particles.cm.Position[dim] = ptclCM->Position[dim];
		ptclGroup->sym_int.particles.cm.Velocity[dim] = ptclCM->Velocity[dim];
	}
	ptclGroup->sym_int.particles.shiftToOriginFrame();
	group2->sym_int.particles.template writeBackMemberAll<Particle>();
*/

	delete ptclGroup;
}

void Group::initialManager() {

	manager.interaction.gravitational_constant = 1.0;
	manager.time_step_min = 1e-13; // minimum physical time step // 1e-13 in ar.cxx
	manager.ds_scale = 1.0; // step size scaling factor // reference: ar.cxx
	manager.time_error_max = 0.25*manager.time_step_min; // time synchronization absolute error limit for AR, default is 0.25*dt-min
	// reference: ar.cxx
	manager.energy_error_relative_max = 1e-10; // relative energy error limit for AR, phase error requirement
	// 1e-10 in ar.cxx
	// 1e-8 in PeTar
	manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX; // maximum timescale for maximum slowdown factor, time-end
	// if (slowdown_timescale_max.value>0.0) manager.slowdown_timescale_max = slowdown_timescale_max.value;
	// else if (time_end.value>0.0) manager.slowdown_timescale_max = time_end.value;
	// else manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX;
	// should be positive
	manager.slowdown_pert_ratio_ref = 1e-6; // slowdown perturbation ratio reference
	// 1e-6 in ar.cxx
	// 1e-4 in PeTar
	manager.step_count_max = 1000000; // number of maximum (integrate/output) step for AR integration // set symplectic order
	// 1000000 in PeTar & ar.cxx
	manager.step.initialSymplecticCofficients(-6); // Symplectic integrator order, should be even number
	// -6 in PeTar & ar.cxx
	manager.interrupt_detection_option = 2; // modify orbit or check interruption using modifyAndInterruptIter function
											// 0: turn off
											// 1: modify the binary orbits based on detetion criterion
											// 2. modify and also interrupt integrations

	// Eunwoo: it is turned off now but I will turn it on later.
	// Eunwoo: It can be used for merging star (dr < sum of radius) or destroy.
}

void Group::initialIntegrator(int NumMembers) {

	sym_int.manager = &manager;

	sym_int.particles.setMode(COMM::ListMode::copy);
    sym_int.particles.reserveMem(NumMembers);
	sym_int.info.reserveMem(NumMembers);

	fprintf(binout, "Mem PID:");
    for (int i = 0; i < groupCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[groupCM->NewNeighbors[i]];
		if (!members->isCMptcl) {
			members->CMPtclIndex = groupCM->ParticleIndex; // added for write_out_group function by EW 2025.1.6
			sym_int.particles.addMemberAndAddress(*members);
			fprintf(binout, " %d", sym_int.particles[i].PID);
		}
		else {
			for (int j=0; j < members->NewNumberOfNeighbor; j++) {
				Particle* members_members = &particles[members->NewNeighbors[j]];
				members_members->CMPtclIndex = groupCM->ParticleIndex; // added for write_out_group function by EW 2025.1.6
				sym_int.particles.addMemberAndAddress(*members_members);
				fprintf(binout, " %d", sym_int.particles[i].PID);
			}
		}
    }
	fprintf(binout, "\n");
	fflush(binout);

	sym_int.info.r_break_crit = rbin/position_unit; // distance criterion for checking stability
	// more information in symplectic_integrator.h
	// ar.cxx: 1e-3 pc
	// check whether the system is stable for 10000 out period and the apo-center is below break criterion
	// PeTar (hard.hpp): sym_int.info.r_break_crit = std::max(sym_int.info.r_break_crit,ptcl_origin[i].getRGroup());

	// // manager.print(std::cerr); // Eunwoo deleted

    sym_int.reserveIntegratorMem();
	sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);


	sym_int.particles.calcCenterOfMass();

	// sym_int.initialIntegration(CurrentTime*EnzoTimeStep);
    // sym_int.info.calcDsAndStepOption(manager.step.getOrder(), manager.interaction.gravitational_constant, manager.ds_scale);

	//! Fix step options for integration with adjusted step (not for time sychronizatio phase)
	// PeTar doesn't set this value explicitly!
	// sym_int.info.fix_step_option = AR::FixStepOption::none; // none: don't fix step
	// sym_int.info.fix_step_option = AR::FixStepOption::always; // always: use the given step without change
	// sym_int.info.fix_step_option = AR::FixStepOption::later; // later: fix step after a few adjustment of initial steps due to energy error

}

	//        ///////////////////        //
	//        ///////////////////        //
	//        NEWFBINITIALIZATION        //
	//        ///////////////////        //
	//        ///////////////////        //



// Initialize new Few body group
void NewFBInitialization(Particle* ptclCM) {

	Group* ptclGroup = new Group();

	ptclCM->GroupInfo = ptclGroup;
	ptclGroup->groupCM = ptclCM;

	// Find member particle with the biggest CurrentTimeIrr
	Particle* ptcl = &particles[ptclCM->NewNeighbors[0]];
	int NumberOfMembers=0;

	for (int i = 0; i < ptclCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];
		if (members->CurrentTimeIrr > ptcl->CurrentTimeIrr) {
        	ptcl = members;
    	}
		if (!members->isCMptcl)
			NumberOfMembers++;
		else
			NumberOfMembers += members->NewNumberOfNeighbor;
    }

	fprintf(binout, "Starting time CurrentTimeIrr (Myr): %e\n", ptcl->CurrentTimeIrr*EnzoTimeStep*1e4);

#ifdef SEVN
	ptclGroup->useSEVN = false;
	REAL dt_evolve_next = NUMERIC_FLOAT_MAX; // Myr
	for (Particle* members : ptclGroup->Members) {
		if (members->star == nullptr || members->star->amiremnant())
			continue;
		REAL dt_evolve = ptcl->CurrentTimeIrr*EnzoTimeStep*1e4 - members->EvolutionTime; // Myr
		
    	if (dt_evolve > 0) {
			members->star->sync_with(dt_evolve);
			members->EvolutionTime += members->star->getp(Timestep::ID);
			// assert(members->EvolutionTime == ptcl->CurrentTimeIrr*EnzoTimeStep*1e4);
			members->star->evolve();
			// fprintf(SEVNout, "INI. After. PID: %d, ET: %e, WT: %e\n", members->PID, members->EvolutionTime, members->star->getp(Worldtime::ID));
			// fflush(SEVNout);
			// fprintf(SEVNout, "CurrentTimeIrr: %e, PID: %d, EvolutionTime: %e\n", ptcl->CurrentTimeIrr*EnzoTimeStep*1e4, members->PID, members->EvolutionTime);
			UpdateEvolution(members);
			if (members->star->vkick[3] > 0.0) {
				for (Particle* members : group->Members)
					members->isGroup = false;
				delete ptclGroup;
				ptclGroup = nullptr;
				delete group;
				group = nullptr;
				return;
			}
			else if (members->star->amiempty()) {
				for (Particle* members : group->Members)
					members->isGroup = false;
				delete ptclGroup;
				ptclGroup = nullptr;
				delete group;
				group = nullptr;

				members->isErase = true;

				particle.erase(
					std::remove_if(particle.begin(), particle.end(),
						[](Particle* p) {
						bool to_remove = p->isErase;
						//if (to_remove) delete p;
						return to_remove;
						}),
					particle.end());

				RegularList.erase(
					std::remove_if(RegularList.begin(), RegularList.end(),
						[](Particle* p) {
						bool to_remove = p->isErase;
						//if (to_remove) delete p;
						return to_remove;
						}),
					RegularList.end());

				for (int i=0; i<particle.size(); i++) {
					particle[i]->ParticleIndex = i;

					// This might be not necessary because we're moving to the regular acceleration routine, and re-set neighbors.
					// Should be checked again later.
					particle[i]->ACList.erase(
							std::remove_if(particle[i]->ACList.begin(), particle[i]->ACList.end(),
								[](Particle* p) {
								return p->isErase; }),
							particle[i]->ACList.end());

					particle[i]->NumberOfAC = particle[i]->ACList.size();
				}
				return;
			}
		}
		if (!members->star->amiremnant()) {
			ptclGroup->useSEVN = true;
			if (members->star->getp(Timestep::ID) < dt_evolve_next)
				dt_evolve_next = members->star->getp(Timestep::ID);
		}
		
	}
	if (ptclGroup->useSEVN) {
		ptclGroup->EvolutionTime = ptcl->CurrentTimeIrr*EnzoTimeStep*1e4;
		ptclGroup->EvolutionTimeStep = dt_evolve_next;
		// fprintf(SEVNout, "INI. ETS: %e\n", ptclGroup->EvolutionTimeStep);
		// fflush(SEVNout);
		// fprintf(binout, "EvolutionTimeStep: %e Myr\n", ptclGroup->EvolutionTimeStep);
		// fflush(binout);
	}
#endif

	for (int i = 0; i < ptclCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];

		double dt = ptcl->CurrentTimeIrr - members->CurrentTimeIrr;
		double pos[Dim], vel[Dim];
		members->predictParticleSecondOrder(dt, pos, vel);
		
		if (members->isCMptcl) {
			Particle* members_members;

			for (int j=0; j<members->NewNumberOfNeighbor; j++) {
				members_members = &particles[members->NewNeighbors[j]];

				for (int dim=0; dim<Dim; dim++) {
					members_members->Position[dim] += pos[dim] - members->Position[dim];
					members_members->Velocity[dim] += vel[dim] - members->Velocity[dim];
				}
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				members->Position[dim] = pos[dim];
				members->Velocity[dim] = vel[dim];
			}
		}
	}

	// Let's link CM particle with the cm particles made in the binary tree (SDAR).

	std::cout << "2" << std::endl;
	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(NumberOfMembers); // Binary tree is made and CM particle is made automatically.

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->Position[dim] = ptclGroup->sym_int.particles.cm.Position[dim];
		ptclCM->Velocity[dim] = ptclGroup->sym_int.particles.cm.Velocity[dim];
		ptclCM->Mass = ptclGroup->sym_int.particles.cm.Mass;
	}

	// Set ptcl information like time, PID, etc.

	ptclCM->CurrentTimeIrr  = ptcl->CurrentTimeIrr;
	ptclCM->CurrentTimeReg  = ptcl->CurrentTimeReg;
	ptclCM->CurrentBlockIrr = ptcl->CurrentBlockIrr; 
	ptclCM->CurrentBlockReg = ptcl->CurrentBlockReg;
	ptclCM->NewCurrentBlockIrr = ptcl->NewCurrentBlockIrr;

	ptclCM->TimeStepIrr     = ptcl->TimeStepIrr;
	ptclCM->TimeBlockIrr    = ptcl->TimeBlockIrr;
	ptclCM->TimeLevelIrr    = ptcl->TimeLevelIrr;

	ptclCM->TimeStepReg     = ptcl->TimeStepReg;
	ptclCM->TimeBlockReg    = ptcl->TimeBlockReg;
	ptclCM->TimeLevelReg    = ptcl->TimeLevelReg;


	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);

	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		fprintf(binout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(binout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(binout, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
		fprintf(binout, "PID: %d. Particle type: %d\n", members->PID, members->ParticleType);
        fprintf(binout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(binout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(binout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(binout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(binout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "PID: %d. Time Steps (Myr) - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(binout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(binout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);

		// fprintf(binout, "PID: %d. Position - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0], members->Position[1], members->Position[2]);
		// fprintf(binout, "PID: %d. Velocity - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0], members->Velocity[1], members->Velocity[2]);
		// fprintf(binout, "PID: %d. Mass - %e, \n", members->PID, members->Mass);
		// fprintf(binout, "PID: %d. Particle type: %d\n", members->PID, members->ParticleType);
        // fprintf(binout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(binout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(binout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		// fprintf(binout, "PID: %d. Time Steps - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*EnzoTimeStep, members->TimeStepReg*EnzoTimeStep);
		// fprintf(binout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(binout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);
    }

	ptclGroup->CurrentTime	= ptclCM->CurrentTimeIrr;

	// if (!(ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[0]*ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[0]<1e-10)) { 
	// 	fprintf(stderr, "NewFBInitialization\n");
	// 	fprintf(stderr, "pos: (%e, %e, %e)\n", ptclGroup->sym_int.info.getBinaryTreeRoot().Position[0], ptclGroup->sym_int.info.getBinaryTreeRoot().Position[1], ptclGroup->sym_int.info.getBinaryTreeRoot().Position[2]);
	// 	fprintf(stderr, "vel: (%e, %e, %e)\n", ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[0], ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[1], ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[2]);
	// 	for (Particle* members: ptclGroup->Members) {
	// 		fprintf(stderr, "PID: %d\n", members->PID);
	// 	}
	// 	fflush(stderr);
	// }

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	for (int i = 0; i < ptclCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];

		if (members->isCMptcl)
			members->clear();
	}

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	CalculateAcceleration01(ptclCM);
	CalculateAcceleration23(ptclCM);

	ptclCM->calculateTimeStepReg();
	if (ptclCM->TimeLevelReg <= ptcl->TimeLevelReg-1 
			&& ptcl->TimeBlockReg/2+ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr+ptcl->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg-1;
	}
	else if  (ptclCM->TimeLevelReg >= ptcl->TimeLevelReg+1) {
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg+1;
	}
	else 
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg;

	ptclCM->TimeStepReg  = static_cast<double>(pow(2, ptclCM->TimeLevelReg));
	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));

	ptclCM->calculateTimeStepIrr();
	// ptclCM->calculateTimeStepIrr2(); // by EW 2025.1.4
/*
	while (ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr <= global_time_irr 
			&& ptclCM->TimeLevelIrr <= ptcl->TimeLevelIrr) { //first condition guarantees that ptclcm is small than ptcl
		ptclCM->TimeLevelIrr++;
		ptclCM->TimeStepIrr  = static_cast<double>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}
*/

/* // Eunwoo added
	while (ptclCM->TimeStepIrr*EnzoTimeStep*1e4 < 1e-7) {
		ptclCM->TimeLevelIrr += 1;
		ptclCM->TimeStepIrr = static_cast<REAL>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}
*/

/* // Eunwoo test
	if (ptclGroup->sym_int.info.getBinaryTreeRoot().semi < 0) { // Only for hyperbolic case
		while (ptclCM->TimeStepIrr*EnzoTimeStep > abs(ptclGroup->sym_int.info.getBinaryTreeRoot().t_peri)) {
			ptclCM->TimeLevelIrr -= 1;
			ptclCM->TimeStepIrr = static_cast<REAL>(pow(2, ptclCM->TimeLevelIrr));
			ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
		}
	}
*/ // Eunwoo test

// /* // Eunwoo test
	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->RadiusOfNeighbor));
		fprintf(binout, "Bound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(binout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
		fprintf(binout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(binout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(binout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
	else {
		ptclGroup->sym_int.info.r_break_crit = 2*bin_root.semi*(1-bin_root.ecc); // r_break_crit = 2*peri
		fprintf(binout, "Unbound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(binout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(binout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(binout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
// */ // Eunwoo test


	fprintf(binout, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization\n");

	fprintf(binout, "Position (pc) - x:%e, y:%e, z:%e, \n", ptclCM->Position[0]*position_unit, ptclCM->Position[1]*position_unit, ptclCM->Position[2]*position_unit);
	fprintf(binout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[1]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[2]*velocity_unit/yr*pc/1e5);
	fprintf(binout, "Mass (Msol) - %e, \n", ptclCM->Mass*mass_unit);
	fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	fprintf(binout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_reg[0][0], ptclCM->a_reg[1][0], ptclCM->a_reg[2][0]);
	// fprintf(binout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_reg[0][1], ptclCM->a_reg[1][1], ptclCM->a_reg[2][1]);
	// fprintf(binout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_reg[0][2], ptclCM->a_reg[1][2], ptclCM->a_reg[2][2]);
	// fprintf(binout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_reg[0][3], ptclCM->a_reg[1][3], ptclCM->a_reg[2][3]);
	fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_irr[0][0], ptclCM->a_irr[1][0], ptclCM->a_irr[2][0]);
	// fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_irr[0][1], ptclCM->a_irr[1][1], ptclCM->a_irr[2][1]);
	// fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_irr[0][2], ptclCM->a_irr[1][2], ptclCM->a_irr[2][2]);
	// fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_irr[0][3], ptclCM->a_irr[1][3], ptclCM->a_irr[2][3]);
	fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr*EnzoTimeStep*1e4, ptclCM->TimeStepReg*EnzoTimeStep*1e4);
	// fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", ptclCM->TimeBlockIrr, ptclCM->TimeBlockReg);
	// fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", ptclCM->CurrentBlockIrr, ptclCM->CurrentBlockReg);

	fprintf(binout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(binout);
}

// Use this function when many-body (>3) group breaks during SDAR integration.
void NewFBInitialization3(Group* group) {

	Group* ptclGroup = new Group();

	std::cout <<"\n\n\nStarting Routine NewFBInitialization3" << std::endl;

	// Set ptclGroup members first; this will be very useful

	Particle* ptclCM = group->groupCM;

	ptclGroup->groupCM = ptclCM;
	ptclGroup->isTerminate = group->isTerminate;
	ptclGroup->isMerger = group->isMerger;
	ptclGroup->CurrentTime = group->CurrentTime;
#ifdef SEVN
	ptclGroup->useSEVN = group->useSEVN;
	ptclGroup->EvolutionTime = group->EvolutionTime;
	ptclGroup->EvolutionTimeStep = group->EvolutionTimeStep;
#endif

	ptclCM->NewNumberOfNeighbor = 0;
	for (int i = 0; i < group->sym_int.particles.getSize(); i++) {
		if (group->sym_int.particles[i].Mass != 0) {
			ptclCM->NewNeighbors[ptclCM->NewNumberOfNeighbor] = group->sym_int.particles[i].ParticleIndex;
			ptclCM->NewNumberOfNeighbor++;
		}
	}
	// Let's link CM particle with the cm particles made in the binary tree (SDAR).

	std::cout << "3" << std::endl;
	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(ptclCM->NewNumberOfNeighbor); // Binary tree is made and CM particle is made automatically.

	delete group;
	ptclCM->GroupInfo = ptclGroup;

	fprintf(binout, "The ID of CM is %d.\n", ptclCM->PID);


	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];

		fprintf(binout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(binout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(binout, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
        fprintf(binout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(binout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(binout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(binout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(binout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(binout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "PID: %d. Time Steps (Myr) - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(binout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(binout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);
    }

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

// /* // Eunwoo test
	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->RadiusOfNeighbor));
		fprintf(binout, "Bound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(binout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
		fprintf(binout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(binout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(binout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
	else {
		ptclGroup->sym_int.info.r_break_crit = 2*bin_root.semi*(1-bin_root.ecc); // r_break_crit = 2*peri
		fprintf(binout, "Unbound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(binout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(binout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(binout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
// */ // Eunwoo test

	fprintf(binout, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization3\n");

	fprintf(binout, "Position (pc) - x:%e, y:%e, z:%e, \n", ptclCM->Position[0]*position_unit, ptclCM->Position[1]*position_unit, ptclCM->Position[2]*position_unit);
	fprintf(binout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[1]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[2]*velocity_unit/yr*pc/1e5);
	fprintf(binout, "Mass (Msol) - %e, \n", ptclCM->Mass*mass_unit);
	fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	fprintf(binout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_reg[0][0], ptclCM->a_reg[1][0], ptclCM->a_reg[2][0]);
	// fprintf(binout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_reg[0][1], ptclCM->a_reg[1][1], ptclCM->a_reg[2][1]);
	// fprintf(binout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_reg[0][2], ptclCM->a_reg[1][2], ptclCM->a_reg[2][2]);
	// fprintf(binout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_reg[0][3], ptclCM->a_reg[1][3], ptclCM->a_reg[2][3]);
	fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_irr[0][0], ptclCM->a_irr[1][0], ptclCM->a_irr[2][0]);
	// fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_irr[0][1], ptclCM->a_irr[1][1], ptclCM->a_irr[2][1]);
	// fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_irr[0][2], ptclCM->a_irr[1][2], ptclCM->a_irr[2][2]);
	// fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_irr[0][3], ptclCM->a_irr[1][3], ptclCM->a_irr[2][3]);
	fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr*EnzoTimeStep*1e4, ptclCM->TimeStepReg*EnzoTimeStep*1e4);
	// fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", ptclCM->TimeBlockIrr, ptclCM->TimeBlockReg);
	// fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", ptclCM->CurrentBlockIrr, ptclCM->CurrentBlockReg);

	fprintf(binout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(binout);
}
#endif