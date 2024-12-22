#include "stdio.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"

// #define SEVN

bool CreateComputationList(Particle* ptcl);
bool CreateComputationChain(std::vector<Particle*> &particle);
void generate_Matrix(REAL a[3], REAL (&A)[3][4]);
void InitializeFBParticle(Particle* FBParticle, std::vector<Particle*> &particle);
void UpdateNextRegTime(std::vector<Particle*> &particle);
bool UpdateComputationChain(Particle* ptcl);

void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle){

	Group* ptclGroup;

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination.cpp... (CM PID: %d)\n\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->CurrentTimeIrr*EnzoTimeStep*1e4);

	ptclGroup	= ptclCM->GroupInfo;


	// Update ptclCM pos & velocity first. Then, update pos & vel of group members to the original frame.
	if (ptclCM->NumberOfAC != 0)
		ptclCM->updateParticle();
	ptclGroup->sym_int.particles.shiftToOriginFrame();
	ptclGroup->sym_int.particles.template writeBackMemberAll<Particle>();
	ptclCM->CurrentBlockIrr += ptclCM->TimeBlockIrr;
	ptclCM->CurrentTimeIrr   = ptclCM->CurrentBlockIrr*time_step;
	ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr);
	
	// Eunwoo: Set CurrentBlock and CurrentTime for group particles.
	for (Particle* members : ptclGroup->Members) {
		members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr;
		members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
		members->CurrentTimeIrr		= ptclCM->CurrentTimeIrr;
		members->CurrentTimeReg		= ptclCM->CurrentTimeReg;

		members->TimeLevelIrr		= ptclCM->TimeLevelIrr; // Eunwoo: I think this is redundant.
		members->TimeLevelReg		= ptclCM->TimeLevelReg; // Eunwoo: I think this is redundant.
	}
	// fprintf(binout, "CM Time Steps (Myr) - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr*EnzoTimeStep*1e4, ptclCM->TimeStepReg*EnzoTimeStep*1e4);

	// Update particle vector: erase ptclCM and add group particles
	ptclCM->isErase = true;
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

	for (Particle* members : ptclGroup->Members) {
		members->ParticleOrder = particle.size();
		particle.push_back(members);
	}

	// Find neighbors and calculate the 0th, 1st, 2nd, 3rd derivative of accleration for group particles 
	fprintf(binout,"initialize group particles \n");
	for (Particle* members : ptclGroup->Members) {
		InitializeFBParticle(members, particle);
		// members->calculateTimeStepReg2();
		members->calculateTimeStepReg();
		// while (members->TimeStepReg < members->CurrentTimeIrr - members->CurrentTimeReg) {
		// 	members->TimeLevelReg++;
		// 	members->TimeStepReg  = static_cast<REAL>(pow(2, members->TimeLevelReg));
		// 	members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block));
		// }

		// /* // Eunwoo: just for a while
		if (members->TimeLevelReg <= ptclCM->TimeLevelReg-1 
				&& members->TimeBlockReg/2+members->CurrentBlockReg > ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
			members->TimeLevelReg = ptclCM->TimeLevelReg-1;
		}
		else if  (members->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
			members->TimeLevelReg = ptclCM->TimeLevelReg+1;
		}
		else 
			members->TimeLevelReg = ptclCM->TimeLevelReg;
		members->TimeStepReg  = static_cast<REAL>(pow(2, members->TimeLevelReg)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()
		members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()
		// /* // Eunwoo: just for a while

		/* // Eunwoo added
		while (members->TimeStepReg*EnzoTimeStep*std::sqrt(members->Velocity[0]*members->Velocity[0]+members->Velocity[1]*members->Velocity[1]+members->Velocity[2]*members->Velocity[2]) > members->RadiusOfAC) {
			members->TimeLevelReg--;
			members->TimeStepReg = static_cast<REAL>(pow(2, members->TimeLevelReg));
			members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block));
		}
		// */ // Eunwoo added

		// members->calculateTimeStepIrr2(members->a_tot, members->a_irr);
		members->calculateTimeStepIrr(members->a_tot, members->a_irr);
		// while (members->CurrentBlockIrr+members->TimeBlockIrr <= global_time_irr) { //first condition guarantees that members is small than ptcl
		// 	members->TimeLevelIrr++;
		// 	members->TimeStepIrr  = static_cast<REAL>(pow(2, members->TimeLevelIrr));
		// 	members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
		// }

		// members->TimeLevelIrr--; // Eunwoo test
		// members->TimeStepIrr = static_cast<REAL>(pow(2, members->TimeLevelIrr));
		// members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
		if (ptclGroup->sym_int.particles.getSize() > 2) {
			members->TimeLevelIrr--;
			members->TimeStepIrr = static_cast<REAL>(pow(2, members->TimeLevelIrr));
			members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
		}
	}


// /* Eunwoo: Is this the problem?
	fprintf(binout,"replacing CM particle in neighbor list to component particles \n");

	int index = 0;
	//fprintf(stderr, "neighbor of %d, ", ptclCM->PID);
	for (Particle* ptcl: particle) {

		index = 0;
		for (Particle* neighbor: ptcl->ACList) {
			if (neighbor->PID == ptclCM->PID) {
				ptcl->ACList.erase(ptcl->ACList.begin() + index);
				ptcl->ACList.insert(ptcl->ACList.end(), ptclGroup->Members.begin(), ptclGroup->Members.end());
				ptcl->NumberOfAC = ptcl->ACList.size();
				break; // Eunwoo check
			}
			index++;
		}
		// Eunwoo added
		// if (dist(ptcl->Position, ptclCM->Position) > ptcl->RadiusOfAC 
		// 		&& dist(ptcl->Position, ptclCM->Position) < 2*ptcl->RadiusOfAC) {
		// 	ptcl->ACList.insert(ptcl->ACList.end(), ptclGroup->Members.begin(), ptclGroup->Members.end());
		// 	ptcl->NumberOfAC = ptcl->ACList.size();
		// }
		// Eunwoo added
	}
// */



	//re-do UpdateNextRegTime
	UpdateNextRegTime(particle);

	// fprintf(binout,"total number of particles = %lu, total number of groups = %lu \n", particle.size(), GroupList.size()-1); // GroupList.size() - 1 because we are terminating it.
	fprintf(binout,"total number of ComputationList = %lu\n", ComputationList.size());

	for (Particle* members : ptclGroup->Members) {
		fprintf(binout,"PID: %d\n", members->PID);
		fprintf(binout, "Position (pc) - x:%e, y:%e, z:%e, \n", members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(binout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(binout, "Mass (Msol) - %e, \n", members->Mass*mass_unit);

		fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
		// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(binout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(binout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(binout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(binout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", members->CurrentBlockIrr, members->CurrentBlockReg);
	}


	// Change the booleans, pointers, and vectors about group information for the group members

	for (Particle* members : ptclGroup->Members) {
		members->isErase		= false;
		members->isGroup		= false;
		// members->GroupInfo		= nullptr; // GroupInfo is only assigned to the CMptcl, so it is needless
	}

	// // we also need to delete ptclGroup from the group list
	// // fprintf(binout,"deleting binary information from the GroupList \n");
	// ptclGroup->isErase = true;
	// GroupList.erase(
	// 		std::remove_if(GroupList.begin(), GroupList.end(),
	// 			[](Group* p) {
	// 			bool to_remove = p->isErase;
	// 			//if (to_remove) delete p;
	// 			return to_remove;
	// 			}),
	// 		GroupList.end());

	delete ptclGroup;
	ptclGroup = nullptr;
	// delete ptclCM; // It is deleted automatically when delete ptclGroup!!!
	// ptclCM  = nullptr;

	fprintf(binout,"end of Few Body Termination\n");
	fprintf(binout,"--------------------------------------\n");
	fflush(binout);

}

// Use this function when SDAR (2-body) is interrupted in the middle of its integration by stellar merger, TDE, GW merger, etc.
// current_time: interrupted time in SDAR
// next_time: intended time to integrate
void FBTermination2(Particle* ptclCM, REAL current_time, std::vector<Particle*> &particle){

	Group* ptclGroup;

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination2.cpp... (CM PID: %d)\n\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of interrupted group (Myr): %e\n", current_time*EnzoTimeStep*1e4);

	ptclGroup	= ptclCM->GroupInfo;

	ptclCM->CurrentBlockIrr += ptclCM->TimeBlockIrr;
	ptclCM->CurrentTimeIrr   = ptclCM->CurrentBlockIrr*time_step;
	ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr);

	// Set CurrentBlock and CurrentTime for group particles.

	for (Particle* members : ptclGroup->Members) {
		members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr; // Block to be integrated
		members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
		members->CurrentTimeIrr		= current_time; // This will be updated later.
		members->CurrentTimeReg		= ptclCM->CurrentTimeReg;

		members->TimeLevelIrr		= ptclCM->TimeLevelIrr; // Eunwoo: I think this is redundant.

		members->TimeLevelReg		= ptclCM->TimeLevelReg; // Eunwoo: I think this is redundant.
	}
	

	// Erase ptclCM from particle vector

	ptclCM->isErase = true;
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

	// Add group members to the particle vector except zero mass particles!

	for (Particle* members : ptclGroup->Members) {
		if (members->Mass == 0)
			members->isErase = true;
		else {
			members->ParticleOrder = particle.size();
			particle.push_back(members);
		}
	}

	// Erase zero mass particles from the ptclGroup->Members!
#ifdef SEVN
	ptclGroup->Members.erase(
		std::remove_if(
			ptclGroup->Members.begin(), ptclGroup->Members.end(),
			[](Particle* p) {
				bool to_remove = p->isErase;
				return to_remove;
			}
		),
		ptclGroup->Members.end()
	);
#else
	ptclGroup->Members.erase(
		std::remove_if(
			ptclGroup->Members.begin(), ptclGroup->Members.end(),
			[](Particle* p) {
				bool to_remove = p->isErase;
				if (to_remove) {
					delete p; // Delete the memory of zero mass particles.
					p = nullptr;
				}
				return to_remove;
			}
		),
		ptclGroup->Members.end()
	);
#endif

	for (Particle* members : ptclGroup->Members) {

		assert(members->Mass > 0); // Zero mass particles should be removed in the above!

		// fprintf(mergerout, "Here I am\n");

		InitializeFBParticle(members, particle); // Find neighbors and set accelerations
		// fprintf(mergerout, "Here I am\n");

		// fprintf(mergerout,"A. PID: %d\n", members->PID);
		// fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", members->Position[0], members->Position[1], members->Position[2]);
		// fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->Velocity[0], members->Velocity[1], members->Velocity[2]);
		// fprintf(mergerout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
		// fprintf(mergerout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(mergerout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(mergerout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		// fprintf(mergerout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(mergerout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(mergerout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(mergerout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		// fprintf(mergerout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(mergerout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(mergerout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(mergerout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);

		// fprintf(mergerout,"N. PID: %d\n", members->ACList[0]->PID);
		// fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", members->ACList[0]->Position[0], members->ACList[0]->Position[1], members->ACList[0]->Position[2]);
		// fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->ACList[0]->Velocity[0], members->ACList[0]->Velocity[1], members->ACList[0]->Velocity[2]);
		// fprintf(mergerout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", members->ACList[0]->a_tot[0][0], members->ACList[0]->a_tot[1][0], members->ACList[0]->a_tot[2][0]);
		// fprintf(mergerout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->ACList[0]->a_tot[0][1], members->ACList[0]->a_tot[1][1], members->ACList[0]->a_tot[2][1]);
		// fprintf(mergerout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->ACList[0]->a_tot[0][2], members->ACList[0]->a_tot[1][2], members->ACList[0]->a_tot[2][2]);
		// fprintf(mergerout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->ACList[0]->a_tot[0][3], members->ACList[0]->a_tot[1][3], members->ACList[0]->a_tot[2][3]);
		// fprintf(mergerout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->ACList[0]->a_reg[0][0], members->ACList[0]->a_reg[1][0], members->ACList[0]->a_reg[2][0]);
		// fprintf(mergerout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->ACList[0]->a_reg[0][1], members->ACList[0]->a_reg[1][1], members->ACList[0]->a_reg[2][1]);
		// fprintf(mergerout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->ACList[0]->a_reg[0][2], members->ACList[0]->a_reg[1][2], members->ACList[0]->a_reg[2][2]);
		// fprintf(mergerout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->ACList[0]->a_reg[0][3], members->ACList[0]->a_reg[1][3], members->ACList[0]->a_reg[2][3]);
		// fprintf(mergerout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->ACList[0]->a_irr[0][0], members->ACList[0]->a_irr[1][0], members->ACList[0]->a_irr[2][0]);
		// fprintf(mergerout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->ACList[0]->a_irr[0][1], members->ACList[0]->a_irr[1][1], members->ACList[0]->a_irr[2][1]);
		// fprintf(mergerout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->ACList[0]->a_irr[0][2], members->ACList[0]->a_irr[1][2], members->ACList[0]->a_irr[2][2]);
		// fprintf(mergerout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->ACList[0]->a_irr[0][3], members->ACList[0]->a_irr[1][3], members->ACList[0]->a_irr[2][3]);

		// fprintf(mergerout, "members current time: %e Myr", members->CurrentTimeIrr);
		// fprintf(mergerout, "CM current time: %e Myr", ptclCM->CurrentTimeIrr);

		// fprintf(mergerout, "a_reg_x: %e\n", members->ACList[0]->Mass/pow(dist(members->Position, members->ACList[0]->Position), 3)*(members->Position[0], members->ACList[0]->Position[0]));

		members->predictParticleSecondOrderIrr(ptclCM->CurrentTimeIrr); // Integrate particle to the intended irregular time
		members->correctParticleFourthOrder(members->CurrentTimeIrr, ptclCM->CurrentTimeIrr, members->a_tot);
		members->CurrentTimeIrr = members->CurrentTimeIrr;

		for (int dim=0; dim<Dim; dim++) {
			members->Position[dim] =  members->PredPosition[dim];
			members->Velocity[dim] =  members->PredVelocity[dim];
			members->NewPosition[dim]  =  members->Position[dim];
			members->NewVelocity[dim]  =  members->Velocity[dim];
		}

		// fprintf(mergerout,"B. PID: %d\n", members->PID);
		// fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", members->Position[0], members->Position[1], members->Position[2]);
		// fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->Velocity[0], members->Velocity[1], members->Velocity[2]);
		// fprintf(mergerout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
		// fprintf(mergerout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(mergerout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(mergerout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		// fprintf(mergerout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(mergerout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(mergerout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(mergerout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		// fprintf(mergerout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(mergerout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(mergerout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(mergerout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);

		// members->calculateTimeStepReg2();
		members->calculateTimeStepReg();
		
		// /* Eunwoo: just for a while
		if (members->TimeLevelReg <= ptclCM->TimeLevelReg-1 
				&& members->TimeBlockReg/2+members->CurrentBlockReg > ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
			members->TimeLevelReg = ptclCM->TimeLevelReg-1;
		}
		else if  (members->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
			members->TimeLevelReg = ptclCM->TimeLevelReg+1;
		}
		else 
			members->TimeLevelReg = ptclCM->TimeLevelReg;
		members->TimeStepReg  = static_cast<REAL>(pow(2, members->TimeLevelReg)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()
		members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()

		// members->calculateTimeStepIrr2(members->a_tot, members->a_irr);
		members->calculateTimeStepIrr(members->a_tot, members->a_irr);
		// members->TimeLevelIrr--;
		// members->TimeStepIrr = static_cast<REAL>(pow(2, members->TimeLevelIrr));
		// members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
	}


// /* Eunwoo: Is this the problem?
	fprintf(binout,"replacing CM particle in neighbor list to component particles \n");

	int index = 0;
	for (Particle* ptcl: particle) {

		index = 0;
		for (Particle* neighbor: ptcl->ACList) {
			if (neighbor->PID == ptclCM->PID) {
				ptcl->ACList.erase(ptcl->ACList.begin() + index);
				ptcl->ACList.insert(ptcl->ACList.end(), ptclGroup->Members.begin(), ptclGroup->Members.end());
				ptcl->NumberOfAC = ptcl->ACList.size();
				// InitializeFBParticle(ptcl, particle); // Eunwoo added
				// ptcl->calculateTimeStepReg(); // Eunwoo added
				// ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr); // Eunwoo added
				break; // Eunwoo check
			}
			index++;
		}
	}
// */


	//re-do UpdateNextRegTime
	UpdateNextRegTime(particle);

	// fprintf(binout,"total number of particles = %lu, total number of groups = %lu \n", particle.size(), GroupList.size()-1); // GroupList.size() - 1 because we are terminating it.
	fprintf(binout,"total number of ComputationList = %lu\n", ComputationList.size());

	for (Particle* members : ptclGroup->Members) {
		fprintf(binout,"PID: %d\n", members->PID);
		fprintf(binout, "Position (pc) - x:%e, y:%e, z:%e, \n", members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(binout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(binout, "Mass (Msol) - %e, \n", members->Mass*mass_unit);

		fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
		// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(binout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(binout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(binout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(binout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", members->CurrentBlockIrr, members->CurrentBlockReg);
	}


	// Change the booleans, pointers, and vectors about group information for the group members

	for (Particle* members : ptclGroup->Members) {
		members->isErase		= false;
		members->isGroup		= false;
		// members->GroupInfo		= nullptr; // GroupInfo is only assigned to the CMptcl, so it is needless
	}

	delete ptclGroup;
	ptclGroup = nullptr;
	// delete ptclCM; // It is deleted automatically when delete ptclGroup!!!
	// ptclCM  = nullptr;

	fprintf(binout,"end of Few Body Termination2\n");
	fprintf(binout,"--------------------------------------\n");
	fflush(binout);
}