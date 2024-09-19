#include "stdio.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"

bool CreateComputationList(Particle* ptcl);
bool CreateComputationChain(std::vector<Particle*> &particle);
void generate_Matrix(REAL a[3], REAL (&A)[3][4]);
void InitializeFBParticle(Particle* FBParticle, std::vector<Particle*> &particle);
void UpdateNextRegTime(std::vector<Particle*> &particle);
bool UpdateComputationChain(Particle* ptcl);

void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle){


	Group* ptclGroup;

	// fprintf(stdout,"--------------------------------------\n");
	// fprintf(stdout,"In FBTermination.cpp...\n\n");
	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination.cpp... (CM PID: %d)\n\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->CurrentTimeIrr*EnzoTimeStep*1e4);

	ptclGroup	= ptclCM->GroupInfo;


	// Eunwoo: Set CurrentBlock and CurrentTime for group particles. This was originally from convertBinaryCoordinatesCartesian().

	for (Particle* members : ptclGroup->Members) {
		members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr;
		members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
		members->CurrentTimeIrr		= ptclCM->CurrentTimeIrr;
		members->CurrentTimeReg		= ptclCM->CurrentTimeReg;

		// members->TimeStepIrr		= ptclCM->TimeStepIrr; // Eunwoo: I think this is redundant.
		// members->TimeBlockIrr		= ptclCM->TimeBlockIrr; // Eunwoo: I think this is redundant.
		// members->TimeLevelIrr		= ptclCM->TimeLevelIrr; // Eunwoo: I think this is redundant.

		// members->TimeStepReg		= ptclCM->TimeStepReg; // Eunwoo: I think this is redundant.
		// members->TimeBlockReg		= ptclCM->TimeBlockReg; // Eunwoo: I think this is redundant.
		// members->TimeLevelReg		= ptclCM->TimeLevelReg; // Eunwoo: I think this is redundant.
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

	// fprintf(stdout,"add the group components to particle list (to be included neighbor search)\n");

	for (Particle* members : ptclGroup->Members) {
		members->ParticleOrder = particle.size();
		particle.push_back(members);
	}

	// Find neighbors and calculate the 0th, 1st, 2nd, 3rd derivative of accleration for group particles 

	fprintf(binout,"initialize group particles \n");
	for (Particle* members : ptclGroup->Members) {
		InitializeFBParticle(members, particle);
		members->calculateTimeStepReg2();
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

		members->calculateTimeStepIrr2(members->a_tot, members->a_irr);
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
		fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(binout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		fprintf(binout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		fprintf(binout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		fprintf(binout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", members->TimeBlockIrr, members->TimeBlockReg);
		fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", members->CurrentBlockIrr, members->CurrentBlockReg);
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
	// delete ptclCM; // It is deleted automatically when delete ptclGroup!!!
	ptclGroup = nullptr;
	ptclCM  = nullptr;

	// fprintf(stdout,"end of Few Body Termination\n");
	fprintf(binout,"end of Few Body Termination\n");

	// fprintf(stdout,"--------------------------------------\n");
	// fflush(stdout);
	fprintf(binout,"--------------------------------------\n");
	fflush(binout);

}

// Use this function when SDAR is interrupted in the middle of its integration by stellar merger, TDE, GW merger, etc.
// current_time: interrupted time in SDAR
// next_time: intended time to integrate
void FBTermination2(Particle* ptclCM, REAL current_time, std::vector<Particle*> &particle){


	Group* ptclGroup;

	// fprintf(stdout,"--------------------------------------\n");
	// fprintf(stdout,"In FBTermination.cpp...\n\n");
	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination2.cpp... (CM PID: %d)\n\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->CurrentTimeIrr*EnzoTimeStep*1e4);

	ptclGroup	= ptclCM->GroupInfo;


	// Set CurrentBlock and CurrentTime for group particles.

	for (Particle* members : ptclGroup->Members) {
		members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr;
		members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
		members->CurrentTimeIrr		= current_time; // This will be updated later.
		members->CurrentTimeReg		= ptclCM->CurrentTimeReg;
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

	ptclGroup->Members.erase(
		std::remove_if(
			ptclGroup->Members.begin(), ptclGroup->Members.end(),
			[](Particle* p) {
				bool to_remove = p->isErase;
				if (to_remove) delete p; // Delete the memory of zero mass particles.
				return to_remove;
			}
		),
		ptclGroup->Members.end()
	);


	for (Particle* members : ptclGroup->Members) {

		assert(members->Mass > 0); // Zero mass particles should be removed in the above!

		InitializeFBParticle(members, particle); // Find neighbors and set accelerations

		members->predictParticleSecondOrderIrr(ptclCM->CurrentTimeIrr); // Integrate particle to the intended irregular time
		members->CurrentTimeIrr = ptclCM->CurrentTimeIrr;

		for (int dim=0; dim<Dim; dim++) {
			members->Position[dim] =  members->PredPosition[dim];
			members->Velocity[dim] =  members->PredVelocity[dim];
			members->NewPosition[dim]  =  members->Position[dim];
			members->NewVelocity[dim]  =  members->Velocity[dim];
		}
		members->calculateTimeStepReg2();
		
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

		members->calculateTimeStepIrr2(members->a_tot, members->a_irr);
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
		fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(binout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		fprintf(binout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		fprintf(binout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		fprintf(binout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", members->TimeBlockIrr, members->TimeBlockReg);
		fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", members->CurrentBlockIrr, members->CurrentBlockReg);
	}


	// Change the booleans, pointers, and vectors about group information for the group members

	for (Particle* members : ptclGroup->Members) {
		members->isErase		= false;
		members->isGroup		= false;
		// members->GroupInfo		= nullptr; // GroupInfo is only assigned to the CMptcl, so it is needless
	}

	delete ptclGroup;
	// delete ptclCM; // It is deleted automatically when delete ptclGroup!!!
	ptclGroup = nullptr;
	ptclCM  = nullptr;

	// fprintf(stdout,"end of Few Body Termination\n");
	fprintf(binout,"end of Few Body Termination\n");

	// fprintf(stdout,"--------------------------------------\n");
	// fflush(stdout);
	fprintf(binout,"--------------------------------------\n");
	fflush(binout);
}