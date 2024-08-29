#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"
#include "Group.h"

REAL getNewTimeStepIrr(REAL f[3][4], REAL df[3][4]);
void getBlockTimeStep(REAL dt, int& TimeLevel, ULL &TimeBlock, REAL &TimeStep);
bool UpdateComputationChain(Particle* ptcl);
void UpdateNextRegTime(std::vector<Particle*> &particle);
bool CreateComputationList(Particle* ptcl);
bool CreateComputationChain(std::vector<Particle*> &particle);
void InitializeFBParticle(Particle* FBParticle, std::vector<Particle*> &particle);


//        ///////////////////        //
//        ///////////////////        //
//           ISFBCANDIDATE           //
//        ///////////////////        //
//        ///////////////////        //



// made 2024.08.08 by Eunwoo Chung

// from particles with shortest time steps...
// Group detection - by distance for the simplest case
// reference - not yet but PeTar/search_group_candidate.hpp might be refered laler

void Particle::isFBCandidate() {

	REAL x[Dim];
	REAL v[Dim];
	REAL xv[Dim];
	REAL current_time;

	REAL r2;
	REAL v2;
	REAL energy; // specific orbital energy
	REAL semi; // semi-major axis

	REAL r_group = 1e-3/position_unit; // group detecting radius: 1e-3 pc (tentative) for FewBody physics

	std::vector<Particle*> groupParticles; // Vector to store pointers of particles in the detected group candidate

	// fprintf(binout, "Function isFBCandidate started\n"); // Eunwoo debug

	for (Particle* ptcl: ACList) {

		r2 = 0.0;
		v2 = 0.0;
		energy = 0.0;

		// find out what the paired particle is
		current_time = this->CurrentTimeIrr > ptcl->CurrentTimeIrr ? \
									 this->CurrentTimeIrr : ptcl->CurrentTimeIrr;
		this->predictParticleSecondOrderIrr(current_time);
		ptcl->predictParticleSecondOrderIrr(current_time);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences
			x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
			v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];
			xv[dim] = x[dim]*v[dim];

			// calculate the square of relative position and velocity
			r2 += x[dim]*x[dim];
			v2 += v[dim]*v[dim];
		}

		// determine they are bound or not
		energy = v2/2 - (this->Mass + ptcl->Mass)/std::sqrt(r2);
		semi = - (this->Mass + ptcl->Mass)/2/energy; // semi-major axis
		REAL p = 1.0 - sqrt(r2)/semi;
		REAL ecc = sqrt(p*p + (xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2])/semi/(this->Mass + ptcl->Mass));

		if (energy < 0) { // bound object

			assert(semi > 0);

			if ((r2 < r_group*r_group) || (semi < r_group)) { // distant binary near periapsis OR close binary
				groupParticles.push_back(ptcl);
				fprintf(binout, "Bound group member added!\n"); // Eunwoo debug
				fprintf(binout, "PID: %d, periapsis dist: %e pc\n", ptcl->PID, semi*(1-ecc)*position_unit);
			}

		}
		else if (r2 < r_group*r_group) { //close fly-by
			groupParticles.push_back(ptcl);
			fprintf(binout, "Unbound group member added!\n"); // Eunwoo debug
			fprintf(binout, "PID: %d, periapsis dist: %e pc\n", ptcl->PID, semi*(1-ecc)*position_unit);
		} // Eunwoo didn't consider distant fly-by because it can be captured later irregular time step!
	}

	if (!groupParticles.empty()) {

		Group *groupCandidate;
		groupCandidate = new Group();

		groupParticles.push_back(this);
		this->isGroup = true; // Assuming isGroup is a member of Particle to mark group status
		for (Particle* ptcl : groupParticles) {
			ptcl->isGroup = true;
		}
		groupCandidate->Members = groupParticles;
		GroupCandidateList.push_back(groupCandidate);
	}
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
	manager.interrupt_detection_option = 0; // modify orbit or check interruption using modifyAndInterruptIter function
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

    for (size_t i = 0; i < Members.size(); ++i) {
        sym_int.particles.addMemberAndAddress(*Members[i]);
		fprintf(binout, "Mem PID: %d\n", Members[i]->PID);
    }

	sym_int.info.r_break_crit = 1e-3/position_unit; // distance criterion for checking stability
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



	// initialize conditions for new Few body group

void NewFBInitialization(Group* group, std::vector<Particle*> &particle) {

	Particle *ptclCM;
	Group *ptclGroup;

	std::cout <<"\n\n\nStarting Routine NewFBInitialization" << std::endl;

	// Set ptclGroup members first; this will be very useful

	ptclGroup = new Group();
	ptclGroup->Members = group->Members;

	// Find member particle with the biggest CurrentTimeIrr
	Particle* ptcl = ptclGroup->Members[0];
	for (Particle* members : ptclGroup->Members) {
    	if (members->CurrentTimeIrr > ptcl->CurrentTimeIrr) {
        	ptcl = members;
    	}
	}

	for (Particle* members : ptclGroup->Members){
		members->predictParticleSecondOrderIrr(ptcl->CurrentTimeIrr);
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

	ptclCM->PID             = -(ptcl->PID + NNB);
	ptclCM->GroupInfo		= ptclGroup;
	ptclCM->isCMptcl        = true;

	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->PredPosition[dim] = ptclCM->Position[dim];
		ptclCM->PredVelocity[dim] = ptclCM->Velocity[dim];
	}

	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (Particle* members: ptclGroup->Members) {
		// fprintf(binout, "PID: %d. Position - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0], members->Position[1], members->Position[2]);
		// fprintf(binout, "PID: %d. Velocity - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0], members->Velocity[1], members->Velocity[2]);
		// fprintf(binout, "PID: %d. Mass - %e, \n", members->PID, members->Mass);
        // fprintf(binout, "PID: %d. Distance from mother particle (pc): %e\n", members->PID, dist(ptclI->Position, members->Position)*position_unit);
        // fprintf(binout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(binout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "PID: %d. Time Steps (Myr) - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(binout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(binout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);
    }

	ptclGroup->groupCM		= ptclCM;
	ptclGroup->CurrentTime	= ptclCM->CurrentTimeIrr;

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	// Erase group particles from particle vector because now they should be replaced with CM particle

	for (Particle* members : ptclGroup->Members) {
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

// /* Eunwoo: How about chainging neighbor list of particles here?

	std::cout << "Changing the Neighbor List of Particles" << std::endl;

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
// */


	UpdateNextRegTime(particle);

  CreateComputationChain(particle);
	CreateComputationList(FirstComputation);

	// Add the binary to binary integration list
	GroupList.push_back(ptclGroup);

	fprintf(binout, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization\n");

	// fprintf(binout, "Position - x:%e, y:%e, z:%e, \n", ptclCM->Position[0], ptclCM->Position[1], ptclCM->Position[2]);
	// fprintf(binout, "Velocity - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0], ptclCM->Velocity[1], ptclCM->Velocity[2]);
	// fprintf(binout, "Mass - %e, \n", ptclCM->Mass);

	// fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	// fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_irr[0][0], ptclCM->a_irr[1][0], ptclCM->a_irr[2][0]);
	// fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_irr[0][1], ptclCM->a_irr[1][1], ptclCM->a_irr[2][1]);
	// fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_irr[0][2], ptclCM->a_irr[1][2], ptclCM->a_irr[2][2]);
	// fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_irr[0][3], ptclCM->a_irr[1][3], ptclCM->a_irr[2][3]);
	fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr*EnzoTimeStep*1e4, ptclCM->TimeStepReg*EnzoTimeStep*1e4);
	// fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", ptclCM->TimeBlockIrr, ptclCM->TimeBlockReg);
	// fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", ptclCM->CurrentBlockIrr, ptclCM->CurrentBlockReg);

	// we also need to change the neighbor list of Particles if they are containing group particles
	// assuming that all the neighbors are bidirectional
	// may need to update later if the radius for neighbor differs depending on the particle

	// std::cout << "Changing the Neighbor List of Particles" << std::endl;

	// // ptclI->isErase = true;
	// // for (Particle* ptclJ : ptclI->GroupParticles) {
	// // 	ptclJ->isErase = true;
	// // }

	// int size = 0;
	// for (Particle* ptcl: particle) {
	// 	size = ptcl->ACList.size();

	// 	ptcl->ACList.erase(
	// 			std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
	// 				[](Particle* p) {
	// 				bool to_remove = p->isErase;
	// 				return to_remove; }),
	// 			ptcl->ACList.end());
	// 	ptcl->NumberOfAC = ptcl->ACList.size();

	// 	if (size != ptcl->NumberOfAC) {
	// 		ptcl->ACList.push_back(ptclCM);
	// 		ptcl->NumberOfAC++;
	// 	}
	// }

	for (Particle* members : ptclGroup->Members) {
		members->isErase = false;
	}

	fprintf(binout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	// fprintf(stdout, "\n---------------------END-OF-NEW-GROUP---------------------\n\n"); // Eunwoo check
	fflush(binout);
	// fflush(stdout); // Eunwoo check
}



// We will use already existing group and CM particle.
void FBModification(std::vector<Particle*> addMembers, std::vector<Particle*> &particle) {

	Particle *ptclCM;
	Group *ptclGroup;

	std::cout <<"\n\n\nStarting Routine FBModification" << std::endl;

	// Set ptclGroup members first; this will be very useful

	ptclGroup = new Group();
	ptclGroup->Members = addMembers;

	// Find member particle with the biggest CurrentTimeIrr
	Particle* ptcl = ptclGroup->Members[0]; // CurrentTimeIrr are not set for existing group particles
	for (Particle* members : ptclGroup->Members) {
    	if (members->CurrentTimeIrr > ptcl->CurrentTimeIrr) {
        	ptcl = members;
    	}
	}

	for (Particle* members : ptclGroup->Members){
		members->predictParticleSecondOrderIrr(ptcl->CurrentTimeIrr);
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

	ptclCM->PID             = -(ptcl->PID + NNB);
	ptclCM->GroupInfo		= ptclGroup;
	ptclCM->isCMptcl        = true;

	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->PredPosition[dim] = ptclCM->Position[dim];
		ptclCM->PredVelocity[dim] = ptclCM->Velocity[dim];
	}

	fprintf(binout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (Particle* members: ptclGroup->Members) {
		// fprintf(binout, "PID: %d. Position - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0], members->Position[1], members->Position[2]);
		// fprintf(binout, "PID: %d. Velocity - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0], members->Velocity[1], members->Velocity[2]);
		// fprintf(binout, "PID: %d. Mass - %e, \n", members->PID, members->Mass);
        // fprintf(binout, "PID: %d. Distance from mother particle (pc): %e\n", members->PID, dist(ptclI->Position, members->Position)*position_unit);
        // fprintf(binout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(binout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "PID: %d. Time Steps (Myr) - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(binout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(binout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);
    }

	ptclGroup->groupCM		= ptclCM;
	ptclGroup->CurrentTime	= ptclCM->CurrentTimeIrr;

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	// Erase group particles from particle vector because now they should be replaced with CM particle

	for (Particle* members : ptclGroup->Members) {
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

// /* Eunwoo: How about chainging neighbor list of particles here?

	std::cout << "Changing the Neighbor List of Particles" << std::endl;

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
// */


	UpdateNextRegTime(particle);

  CreateComputationChain(particle);
	CreateComputationList(FirstComputation);

	// Add the binary to binary integration list
	GroupList.push_back(ptclGroup);

	fprintf(binout, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization\n");

	// fprintf(binout, "Position - x:%e, y:%e, z:%e, \n", ptclCM->Position[0], ptclCM->Position[1], ptclCM->Position[2]);
	// fprintf(binout, "Velocity - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0], ptclCM->Velocity[1], ptclCM->Velocity[2]);
	// fprintf(binout, "Mass - %e, \n", ptclCM->Mass);

	// fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	// fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_irr[0][0], ptclCM->a_irr[1][0], ptclCM->a_irr[2][0]);
	// fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_irr[0][1], ptclCM->a_irr[1][1], ptclCM->a_irr[2][1]);
	// fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_irr[0][2], ptclCM->a_irr[1][2], ptclCM->a_irr[2][2]);
	// fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_irr[0][3], ptclCM->a_irr[1][3], ptclCM->a_irr[2][3]);
	fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr*EnzoTimeStep*1e4, ptclCM->TimeStepReg*EnzoTimeStep*1e4);
	// fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", ptclCM->TimeBlockIrr, ptclCM->TimeBlockReg);
	// fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", ptclCM->CurrentBlockIrr, ptclCM->CurrentBlockReg);

	for (Particle* members : ptclGroup->Members) {
		members->isErase = false;
	}

	fprintf(binout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	// fprintf(stdout, "\n---------------------END-OF-NEW-GROUP---------------------\n\n"); // Eunwoo check
	fflush(binout);
	// fflush(stdout); // Eunwoo check

}