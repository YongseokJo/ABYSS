#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"
// #define SEVN
#ifdef SEVN
void UpdateEvolution(Particle* ptcl);
#endif

REAL getNewTimeStepIrr(REAL f[3][4], REAL df[3][4]);
void getBlockTimeStep(REAL dt, int& TimeLevel, ULL &TimeBlock, REAL &TimeStep);
bool UpdateComputationChain(Particle* ptcl);
void UpdateNextRegTime(std::vector<Particle*> &particle);
bool CreateComputationList(Particle* ptcl);
bool CreateComputationChain(std::vector<Particle*> &particle);
void InitializeFBParticle(Particle* FBParticle, std::vector<Particle*> &particle);
void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle);


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

void Group::initialIntegrator() {

	sym_int.manager = &manager;

	sym_int.particles.setMode(COMM::ListMode::copy);
    sym_int.particles.reserveMem(Members.size());
	sym_int.info.reserveMem(Members.size());

	// bool BHinside = false;

    for (size_t i = 0; i < Members.size(); ++i) {
        sym_int.particles.addMemberAndAddress(*Members[i]);
		fprintf(binout, "Mem PID: %d\n", Members[i]->PID);
		if (Members[i]->ParticleType == (Blackhole+SingleParticle)) 
			PNon = true;
    }

	// if (BHinside) sym_int.manager->setBHinside();

	// for (Particle* member: Members) {
	// 	sym_int.particles.addMemberAndAddress(*member);
	// 	fprintf(binout, "Mem PID: %d\n", member->PID);
	// }

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
void NewFBInitialization(Group* group, std::vector<Particle*> &particle) {


	Particle *ptclCM;
	Group *ptclGroup;

	std::cout <<"\n\n\nStarting Routine NewFBInitialization" << std::endl;

	// Set ptclGroup members first; this will be very useful

	ptclGroup = new Group();
	for (Particle* members : group->Members) {
		if (members->isCMptcl) 
			ptclGroup->Members.insert(ptclGroup->Members.end(), members->GroupInfo->Members.begin(), members->GroupInfo->Members.end());
		
		else 
			ptclGroup->Members.push_back(members);
	}

	// Find member particle with the biggest CurrentTimeIrr
	Particle* ptcl = group->Members[0];
	for (Particle* members : group->Members) {
    	if (members->CurrentTimeIrr > ptcl->CurrentTimeIrr) {
        	ptcl = members;
    	}
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
					particle[i]->ParticleOrder = i;

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

	for (Particle* members: group->Members) {
		members->predictParticleSecondOrderIrr(ptcl->CurrentTimeIrr);
		for (int dim=0; dim<Dim; dim++) {
			members->Position[dim] = members->PredPosition[dim];
			members->Velocity[dim] = members->PredVelocity[dim];
		}
		if (members->isCMptcl) {
			members->GroupInfo->sym_int.particles.shiftToOriginFrame();
			members->GroupInfo->sym_int.particles.template writeBackMemberAll<Particle>();
		}
	}

	// Let's link CM particle with the cm particles made in the binary tree (SDAR).

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(); // Binary tree is made and CM particle is made automatically.

	ptclCM = &ptclGroup->sym_int.particles.cm;

	// Set ptcl information like time, PID, etc.

	ptclCM->CurrentTimeIrr  = ptcl->CurrentTimeIrr;
	ptclCM->CurrentTimeReg  = ptcl->CurrentTimeReg;
	ptclCM->CurrentBlockIrr = ptcl->CurrentBlockIrr; 
	ptclCM->CurrentBlockReg = ptcl->CurrentBlockReg;

	ptclCM->TimeStepIrr     = ptcl->TimeStepIrr;
	ptclCM->TimeBlockIrr    = ptcl->TimeBlockIrr;
	ptclCM->TimeLevelIrr    = ptcl->TimeLevelIrr;

	ptclCM->TimeStepReg     = ptcl->TimeStepReg;
	ptclCM->TimeBlockReg    = ptcl->TimeBlockReg;
	ptclCM->TimeLevelReg    = ptcl->TimeLevelReg;

	ptclCM->PID             = -(std::abs(ptcl->PID) + NNB);
	ptclCM->GroupInfo		= ptclGroup;
	ptclCM->isCMptcl        = true;

	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->PredPosition[dim] = ptclCM->Position[dim];
		ptclCM->PredVelocity[dim] = ptclCM->Velocity[dim];
	}

	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (Particle* members: ptclGroup->Members) {
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

	ptclGroup->groupCM		= ptclCM;
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

	// Erase group particles from particle vector because now they should be replaced with CM particle

	for (Particle* members : group->Members) {
		members->isErase = true;
	}

	particle.erase(
			std::remove_if(particle.begin(), particle.end(),
				[](Particle* p) {
				bool to_remove = p->isErase;
				return to_remove;
				}),
			particle.end());

	// Update particle order and put CM particle into the last of particle vector
	for (int i=0; i<particle.size(); i++) {
		particle[i]->ParticleOrder = i;
	}
	ptclCM->ParticleOrder = particle.size();
	particle.push_back(ptclCM);

/* // New test starts

	REAL radiusOfAC = 0.0; // Initialize variable to store the maximum distance

	ptclCM->ACList.clear(); // Clear ptclCM's ACList

	for (Particle* member : group->Members) {
		for (Particle* neighbor : member->ACList) {
			// Check if neighbor is not inside ptclCM->ACList and not in ptclGroup->Members
			if (std::find(ptclCM->ACList.begin(), ptclCM->ACList.end(), neighbor) == ptclCM->ACList.end() &&
				std::find(group->Members.begin(), group->Members.end(), neighbor) == group->Members.end()) {
				
				ptclCM->ACList.push_back(neighbor); // Insert neighbor into ptclCM->ACList

				REAL distance = dist(ptclCM->Position, neighbor->Position);
				
				// Update maxDistance if the current distance is greater
				if (distance > radiusOfAC) {
					radiusOfAC = distance;
				}
			}
		}
	}

	ptclCM->NumberOfAC = ptclCM->ACList.size();
	ptclCM->RadiusOfAC = radiusOfAC;

	std::cout << ptclCM->PID << "("<< ptclCM->NumberOfAC <<")" << "=" <<std::flush;
	for (Particle * nn:ptclCM->ACList) {
		std::cout << nn->PID << ", ";
	}
	std::cout << std::endl;

	CalculateAcceleration01(ptclCM, particle);
	CalculateAcceleration23(ptclCM, particle);

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->PredPosition[dim] =  ptclCM->Position[dim];
		ptclCM->PredVelocity[dim] =  ptclCM->Velocity[dim];
		ptclCM->NewPosition[dim]  =  ptclCM->Position[dim];
		ptclCM->NewVelocity[dim]  =  ptclCM->Velocity[dim];
	}

	ptclCM->calculateTimeStepReg();
	// while (ptclCM->TimeStepReg < ptclCM->CurrentTimeIrr - ptclCM->CurrentTimeReg) {
	// 	ptclCM->TimeLevelReg++;
	// 	ptclCM->TimeStepReg  = static_cast<REAL>(pow(2, ptclCM->TimeLevelReg));
	// 	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));
	// }
	if (ptclCM->TimeLevelReg <= ptcl->TimeLevelReg-1 
			&& ptcl->TimeBlockReg/2+ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr+ptcl->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg-1;
	}
	else if  (ptclCM->TimeLevelReg >= ptcl->TimeLevelReg+1) {
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg+1;
	}
	else 
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg;

	ptclCM->TimeStepReg  = static_cast<REAL>(pow(2, ptclCM->TimeLevelReg));
	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));

	ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr);
	// /* // Eunwoo: just for a while
	while (ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr <= global_time_irr 
			&& ptclCM->TimeLevelIrr <= ptcl->TimeLevelIrr) { //first condition guarantees that ptclcm is small than ptcl
		ptclCM->TimeLevelIrr++;
		ptclCM->TimeStepIrr  = static_cast<REAL>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}

*/ // New test ends

// /* // Test starts
	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	InitializeFBParticle(ptclCM, particle);

	std::cout << "Calculating Time steps for the CM particle" << std::endl;
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

	ptclCM->TimeStepReg  = static_cast<REAL>(pow(2, ptclCM->TimeLevelReg));
	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));

	// ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr); // fiducial
	ptclCM->calculateTimeStepIrr2(ptclCM->a_tot, ptclCM->a_irr); // test_1e5_3 & 5
	while (ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr <= global_time_irr 
			&& ptclCM->TimeLevelIrr <= ptcl->TimeLevelIrr) { //first condition guarantees that ptclcm is small than ptcl
		ptclCM->TimeLevelIrr++;
		ptclCM->TimeStepIrr  = static_cast<REAL>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}
// */ // Test end

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
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, ptclCM->RadiusOfAC);
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


	int size = 0;
	for (Particle* ptcl: particle) {
		size = ptcl->ACList.size();

		ptcl->ACList.erase(
				std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					return to_remove; }),
				ptcl->ACList.end());
		ptcl->NumberOfAC = ptcl->ACList.size();

		if (size != ptcl->NumberOfAC) {
			ptcl->ACList.push_back(ptclCM);
			ptcl->NumberOfAC++;
			// InitializeFBParticle(ptcl, particle); // Eunwoo added
			// ptcl->calculateTimeStepReg(); // Eunwoo added
			// ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr); // Eunwoo added
		}
	}

	for (Particle* members : group->Members) {
		members->isErase = false;
		if (members->isCMptcl) {
			// members->GroupInfo->isErase = true;
			// GroupList.erase(
			// 		std::remove_if(GroupList.begin(), GroupList.end(),
			// 			[](Group* p) {
			// 			bool to_remove = p->isErase;
			// 			//if (to_remove) delete p;
			// 			return to_remove;
			// 			}),
			// 		GroupList.end());
			delete members->GroupInfo;
			members->GroupInfo = nullptr;
			members  = nullptr;
		}
	}


	UpdateNextRegTime(particle);
	
	CreateComputationChain(particle);
	CreateComputationList(FirstComputation);

	// // Add the binary to binary integration list
	// GroupList.push_back(ptclGroup);

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

	delete group;
	group = nullptr;

	fprintf(binout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(binout);
}

// Use this function when many-body (>3) group splits up into small group.
void NewFBInitialization2(Group* group, std::vector<Particle*> &particle) {

	Particle *ptclCM;
	Group *ptclGroup;

	std::cout <<"\n\n\nStarting Routine NewFBInitialization2" << std::endl;

	// Set ptclGroup members first; this will be very useful

	ptclGroup = new Group();
	for (Particle* members : group->Members) {
		if (members->isCMptcl) 
			ptclGroup->Members.insert(ptclGroup->Members.end(), members->GroupInfo->Members.begin(), members->GroupInfo->Members.end());
		
		else 
			ptclGroup->Members.push_back(members);
	}

	// Eunwoo added: to prevent too small irregular timesteps when many-body (> 3) repeatedly forms and breaks
	// for (Particle* members : ptclGroup->Members) {
	// 	members->TimeLevelIrr++;
	// 	members->TimeStepIrr = static_cast<REAL>(pow(2, members->TimeLevelIrr));
	// 	members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
	// }

	// Find member particle with the biggest CurrentTimeIrr
	Particle* ptcl = group->Members[0];
	for (Particle* members : group->Members) {
    	if (members->CurrentTimeIrr > ptcl->CurrentTimeIrr) {
        	ptcl = members;
    	}
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
			// fprintf(SEVNout, "INI. After. PID: %d, ET: %e, WT: %e, GT: %e\n", members->PID, members->EvolutionTime, members->star->getp(Worldtime::ID), ptcl->CurrentTimeIrr*EnzoTimeStep*1e4);
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
					particle[i]->ParticleOrder = i;

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

	for (Particle* members: group->Members) {
		members->predictParticleSecondOrderIrr(ptcl->CurrentTimeIrr);
		for (int dim=0; dim<Dim; dim++) {
			members->Position[dim] = members->PredPosition[dim];
			members->Velocity[dim] = members->PredVelocity[dim];
		}
		if (members->isCMptcl) {
			members->GroupInfo->sym_int.particles.shiftToOriginFrame();
			members->GroupInfo->sym_int.particles.template writeBackMemberAll<Particle>();
		}
	}

	// Let's link CM particle with the cm particles made in the binary tree (SDAR).

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(); // Binary tree is made and CM particle is made automatically.

	ptclCM = &ptclGroup->sym_int.particles.cm;

	// Set ptcl information like time, PID, etc.

	ptclCM->CurrentTimeIrr  = ptcl->CurrentTimeIrr;
	ptclCM->CurrentTimeReg  = ptcl->CurrentTimeReg;
	ptclCM->CurrentBlockIrr = ptcl->CurrentBlockIrr; 
	ptclCM->CurrentBlockReg = ptcl->CurrentBlockReg;

	ptclCM->TimeStepIrr     = ptcl->TimeStepIrr;
	ptclCM->TimeBlockIrr    = ptcl->TimeBlockIrr;
	ptclCM->TimeLevelIrr    = ptcl->TimeLevelIrr;

	ptclCM->TimeStepReg     = ptcl->TimeStepReg;
	ptclCM->TimeBlockReg    = ptcl->TimeBlockReg;
	ptclCM->TimeLevelReg    = ptcl->TimeLevelReg;

	ptclCM->PID             = -(std::abs(ptcl->PID) + NNB);
	ptclCM->GroupInfo		= ptclGroup;
	ptclCM->isCMptcl        = true;

	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->PredPosition[dim] = ptclCM->Position[dim];
		ptclCM->PredVelocity[dim] = ptclCM->Velocity[dim];
	}

	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (Particle* members: ptclGroup->Members) {
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

	ptclGroup->groupCM		= ptclCM;
	// ptclGroup->StartTime	= ptclCM->CurrentTimeIrr;
	ptclGroup->CurrentTime	= ptclCM->CurrentTimeIrr;

	// if (!(ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[0]*ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[0]<1e-10)) { 
	// 	fprintf(stderr, "NewFBInitialization2\n");
	// 	fprintf(stderr, "pos: (%e, %e, %e)\n", ptclGroup->sym_int.info.getBinaryTreeRoot().Position[0], ptclGroup->sym_int.info.getBinaryTreeRoot().Position[1], ptclGroup->sym_int.info.getBinaryTreeRoot().Position[2]);
	// 	fprintf(stderr, "vel: (%e, %e, %e)\n", ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[0], ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[1], ptclGroup->sym_int.info.getBinaryTreeRoot().Velocity[2]);
	// 	for (Particle* members: ptclGroup->Members) {
	// 		fprintf(stderr, "PID: %d\n", members->PID);
	// 	}
	// 	fflush(stderr);
	// }

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
	// ptclGroup->sym_int.initialIntegration(0);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	// Erase group particles from particle vector because now they should be replaced with CM particle

	for (Particle* members : group->Members) {
		members->isErase = true;
	}

	particle.erase(
			std::remove_if(particle.begin(), particle.end(),
				[](Particle* p) {
				bool to_remove = p->isErase;
				return to_remove;
				}),
			particle.end());

	// Update particle order and put CM particle into the last of particle vector
	for (int i=0; i<particle.size(); i++) {
		particle[i]->ParticleOrder = i;
	}
	ptclCM->ParticleOrder = particle.size();
	particle.push_back(ptclCM);

/* // New test starts

	REAL radiusOfAC = 0.0; // Initialize variable to store the maximum distance

	ptclCM->ACList.clear(); // Clear ptclCM's ACList

	for (Particle* member : group->Members) {
		for (Particle* neighbor : member->ACList) {
			// Check if neighbor is not inside ptclCM->ACList and not in ptclGroup->Members
			if (std::find(ptclCM->ACList.begin(), ptclCM->ACList.end(), neighbor) == ptclCM->ACList.end() &&
				std::find(group->Members.begin(), group->Members.end(), neighbor) == group->Members.end()) {
				
				ptclCM->ACList.push_back(neighbor); // Insert neighbor into ptclCM->ACList

				REAL distance = dist(ptclCM->Position, neighbor->Position);
				
				// Update maxDistance if the current distance is greater
				if (distance > radiusOfAC) {
					radiusOfAC = distance;
				}
			}
		}
	}

	ptclCM->NumberOfAC = ptclCM->ACList.size();
	ptclCM->RadiusOfAC = radiusOfAC;

	std::cout << ptclCM->PID << "("<< ptclCM->NumberOfAC <<")" << "=" <<std::flush;
	for (Particle * nn:ptclCM->ACList) {
		std::cout << nn->PID << ", ";
	}
	std::cout << std::endl;

	CalculateAcceleration01(ptclCM, particle);
	CalculateAcceleration23(ptclCM, particle);

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->PredPosition[dim] =  ptclCM->Position[dim];
		ptclCM->PredVelocity[dim] =  ptclCM->Velocity[dim];
		ptclCM->NewPosition[dim]  =  ptclCM->Position[dim];
		ptclCM->NewVelocity[dim]  =  ptclCM->Velocity[dim];
	}

	ptclCM->calculateTimeStepReg();
	// while (ptclCM->TimeStepReg < ptclCM->CurrentTimeIrr - ptclCM->CurrentTimeReg) {
	// 	ptclCM->TimeLevelReg++;
	// 	ptclCM->TimeStepReg  = static_cast<REAL>(pow(2, ptclCM->TimeLevelReg));
	// 	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));
	// }
	if (ptclCM->TimeLevelReg <= ptcl->TimeLevelReg-1 
			&& ptcl->TimeBlockReg/2+ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr+ptcl->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg-1;
	}
	else if  (ptclCM->TimeLevelReg >= ptcl->TimeLevelReg+1) {
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg+1;
	}
	else 
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg;

	ptclCM->TimeStepReg  = static_cast<REAL>(pow(2, ptclCM->TimeLevelReg));
	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));

	ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr);
	// /* // Eunwoo: just for a while
	while (ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr <= global_time_irr 
			&& ptclCM->TimeLevelIrr <= ptcl->TimeLevelIrr) { //first condition guarantees that ptclcm is small than ptcl
		ptclCM->TimeLevelIrr++;
		ptclCM->TimeStepIrr  = static_cast<REAL>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}

*/ // New test ends

// /* // Test starts
	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	InitializeFBParticle(ptclCM, particle);

	std::cout << "Calculating Time steps for the CM particle" << std::endl;
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

	ptclCM->TimeStepReg  = static_cast<REAL>(pow(2, ptclCM->TimeLevelReg));
	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));

	ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr);
	while (ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr <= global_time_irr 
			&& ptclCM->TimeLevelIrr <= ptcl->TimeLevelIrr) { //first condition guarantees that ptclcm is small than ptcl
		ptclCM->TimeLevelIrr++;
		ptclCM->TimeStepIrr  = static_cast<REAL>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}
// */ // Test end

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
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, ptclCM->RadiusOfAC);
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


	int size = 0;
	for (Particle* ptcl: particle) {
		size = ptcl->ACList.size();

		ptcl->ACList.erase(
				std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					return to_remove; }),
				ptcl->ACList.end());
		ptcl->NumberOfAC = ptcl->ACList.size();

		if (size != ptcl->NumberOfAC) {
			ptcl->ACList.push_back(ptclCM);
			ptcl->NumberOfAC++;
			// InitializeFBParticle(ptcl, particle); // Eunwoo added
			// ptcl->calculateTimeStepReg(); // Eunwoo added
			// ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr); // Eunwoo added
		}
	}

	for (Particle* members : group->Members) {
		members->isErase = false;
		if (members->isCMptcl) {
			// members->GroupInfo->isErase = true;
			// GroupList.erase(
			// 		std::remove_if(GroupList.begin(), GroupList.end(),
			// 			[](Group* p) {
			// 			bool to_remove = p->isErase;
			// 			//if (to_remove) delete p;
			// 			return to_remove;
			// 			}),
			// 		GroupList.end());
			delete members->GroupInfo;
			members->GroupInfo = nullptr;
			members  = nullptr;
		}
	}


	UpdateNextRegTime(particle);
//	CreateComputationChain(particle);
// 	CreateComputationList(FirstComputation);

	// // Add the binary to binary integration list
	// GroupList.push_back(ptclGroup);

	fprintf(binout, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization2\n");

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

	delete group;
	group = nullptr;

	fprintf(binout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(binout);
}

// Use this function when many-body (>3) group breaks during SDAR integration.
void NewFBInitialization3(Group* group, std::vector<Particle*> &particle) {

	Particle *ptclCM;
	Group *ptclGroup;

	std::cout <<"\n\n\nStarting Routine NewFBInitialization3" << std::endl;

	// Set ptclGroup members first; this will be very useful

	ptclGroup = new Group();

	ptclGroup->Members = group->Members;

#ifdef SEVN
	ptclGroup->useSEVN = group->useSEVN;
	ptclGroup->EvolutionTime = group->EvolutionTime;
	ptclGroup->EvolutionTimeStep = group->EvolutionTimeStep;
#endif

	// Let's link CM particle with the cm particles made in the binary tree (SDAR).

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(); // Binary tree is made and CM particle is made automatically.

	ptclCM = &ptclGroup->sym_int.particles.cm;

	// Set ptcl information like time, PID, etc.

	ptclCM->CurrentTimeIrr  = group->groupCM->CurrentTimeIrr;
	ptclCM->CurrentTimeReg  = group->groupCM->CurrentTimeReg;
	ptclCM->CurrentBlockIrr = group->groupCM->CurrentBlockIrr; 
	ptclCM->CurrentBlockReg = group->groupCM->CurrentBlockReg;

	ptclCM->TimeStepIrr     = group->groupCM->TimeStepIrr;
	ptclCM->TimeBlockIrr    = group->groupCM->TimeBlockIrr;
	ptclCM->TimeLevelIrr    = group->groupCM->TimeLevelIrr;

	ptclCM->TimeStepReg     = group->groupCM->TimeStepReg;
	ptclCM->TimeBlockReg    = group->groupCM->TimeBlockReg;
	ptclCM->TimeLevelReg    = group->groupCM->TimeLevelReg;

	ptclCM->PID             = group->groupCM->PID;
	ptclCM->GroupInfo		= ptclGroup;
	ptclCM->isCMptcl        = true;

	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->PredPosition[dim] = ptclCM->Position[dim];
		ptclCM->PredVelocity[dim] = ptclCM->Velocity[dim];
		ptclCM->NewPosition[dim] = group->groupCM->NewPosition[dim];
		ptclCM->NewVelocity[dim] = group->groupCM->NewVelocity[dim];
	}

	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (Particle* members: ptclGroup->Members) {
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

	ptclGroup->groupCM		= ptclCM;
	ptclGroup->CurrentTime	= group->CurrentTime;


	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	// Erase group particles from particle vector because now they should be replaced with CM particle

	group->groupCM->isErase = true;

	particle.erase(
			std::remove_if(particle.begin(), particle.end(),
				[](Particle* p) {
				bool to_remove = p->isErase;
				return to_remove;
				}),
			particle.end());

	// Update particle order and put CM particle into the last of particle vector
	for (int i=0; i<particle.size(); i++) {
		particle[i]->ParticleOrder = i;
	}
	ptclCM->ParticleOrder = particle.size();
	particle.push_back(ptclCM);

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	ptclCM->ACList = group->groupCM->ACList;
	ptclCM->NumberOfAC = group->groupCM->NumberOfAC;
	ptclCM->RadiusOfAC = group->groupCM->RadiusOfAC;

	for (int i=0; i<Dim; i++) {
		ptclCM->PredPosition[i] = group->groupCM->PredPosition[i];
		ptclCM->PredVelocity[i] = group->groupCM->PredPosition[i];
		ptclCM->BackgroundAcceleration[i] = group->groupCM->BackgroundAcceleration[i];
		for (int j=0; j<HERMITE_ORDER; j++) {
			ptclCM->a_tot[i][j] = group->groupCM->a_tot[i][j];
			ptclCM->a_reg[i][j] = group->groupCM->a_reg[i][j];
			ptclCM->a_irr[i][j] = group->groupCM->a_irr[i][j];
		}
	}
	ptclCM->LocalDensity = group->groupCM->LocalDensity;
	ptclCM->NextParticleInEnzo = group->groupCM->NextParticleInEnzo;
	ptclCM->NextParticleForComputation = group->groupCM->NextParticleForComputation;

// /* // Eunwoo test
	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, ptclCM->RadiusOfAC);
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


	int size = 0;
	for (Particle* ptcl: particle) {
		size = ptcl->ACList.size();

		ptcl->ACList.erase(
				std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					return to_remove; }),
				ptcl->ACList.end());
		ptcl->NumberOfAC = ptcl->ACList.size();

		if (size != ptcl->NumberOfAC) {
			ptcl->ACList.push_back(ptclCM);
			ptcl->NumberOfAC++;
			// InitializeFBParticle(ptcl, particle); // Eunwoo added
			// ptcl->calculateTimeStepReg(); // Eunwoo added
			// ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr); // Eunwoo added
		}
	}

	// UpdateNextRegTime(particle);
//	CreateComputationChain(particle);
// 	CreateComputationList(FirstComputation);

	// // Add the binary to binary integration list
	// GroupList.push_back(ptclGroup);

	if (ptclCM->NumberOfAC != 0) // IAR modified
		ptclCM->updateParticle();
	ptclCM->CurrentBlockIrr += ptclCM->TimeBlockIrr;
	ptclCM->CurrentTimeIrr   = ptclCM->CurrentBlockIrr*time_step;

	ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr);

	ptclCM->NextBlockIrr = ptclCM->CurrentBlockIrr + ptclCM->TimeBlockIrr; // of this particle

	UpdateComputationChain(ptclCM);
	// UpdateNextRegTime(particle);
	for (int i=0; i<RegularList.size(); i++) {
		if (RegularList[i]->PID == ptclCM->PID) {
			RegularList[i] = ptclCM;
			break;
		}
	}

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

	delete group;
	group = nullptr;

	fprintf(binout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(binout);
}