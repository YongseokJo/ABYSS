#include "stdio.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../def.h"

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);

void FBTermination(int Order){

	assert(Order <= global_variable.NumberOfParticle-1);

	Particle* ptclCM = &particles[Order];

	Group* ptclGroup;
	ptclGroup	= ptclCM->GroupInfo;

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination.cpp... (CM PID: %d)\n\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->CurrentTimeIrr*EnzoTimeStep*1e4);


	// Update ptclCM pos & velocity first. Then, update pos & vel of group members to the original frame.
	if (ptclCM->NumberOfNeighbor != 0)
		ptclCM->updateParticle();

	for (int dim=0; dim<Dim; dim++) {
		ptclGroup->sym_int.particles.cm.Position[dim] = ptclCM->Position[dim];
		ptclGroup->sym_int.particles.cm.Velocity[dim] = ptclCM->Velocity[dim];
	}
	ptclGroup->sym_int.particles.shiftToOriginFrame();
	ptclGroup->sym_int.particles.template writeBackMemberAll<Particle>();

	ptclCM->CurrentBlockIrr = ptclCM->NewCurrentBlockIrr;
	ptclCM->CurrentTimeIrr  = ptclCM->CurrentBlockIrr*time_step;
	
	// Eunwoo: Set CurrentBlock and CurrentTime for group particles.
	for (int i=0; i<ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &particles[ptclGroup->sym_int.particles[i].ParticleOrder];

		members->isActive = true;

		members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr;
		members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
		members->CurrentTimeIrr		= ptclCM->CurrentTimeIrr;
		members->CurrentTimeReg		= ptclCM->CurrentTimeReg;

		members->TimeLevelIrr		= ptclCM->TimeLevelIrr; // Eunwoo: I think this is redundant.
		members->TimeLevelReg		= ptclCM->TimeLevelReg; // Eunwoo: I think this is redundant.
	}
	// fprintf(binout, "CM Time Steps (Myr) - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr*EnzoTimeStep*1e4, ptclCM->TimeStepReg*EnzoTimeStep*1e4);

	// Update particle vector: erase ptclCM and add group particles
	ptclCM->isActive = false;

	for (int i=0; i<ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &particles[ptclGroup->sym_int.particles[i].ParticleOrder];
		
		CalculateAcceleration01(members);
		CalculateAcceleration23(members);

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
		members->TimeStepReg  = static_cast<double>(pow(2, members->TimeLevelReg)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()
		members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()
		// /* // Eunwoo: just for a while


		// members->calculateTimeStepIrr2(members->a_tot, members->a_irr);
		members->calculateTimeStepIrr();
		// while (members->CurrentBlockIrr+members->TimeBlockIrr <= global_time_irr) { //first condition guarantees that members is small than ptcl
		// 	members->TimeLevelIrr++;
		// 	members->TimeStepIrr  = static_cast<REAL>(pow(2, members->TimeLevelIrr));
		// 	members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
		// }

		/* // This might be used later
		if (ptclGroup->sym_int.particles.getSize() > 2) {
			members->TimeLevelIrr--;
			members->TimeStepIrr = static_cast<REAL>(pow(2, members->TimeLevelIrr));
			members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
		}
		*/
	}

	for (int i=0; i<global_variable.NumberOfParticle; i++) {
		Particle* ptcl = &particles[i];
		auto newEnd = std::remove_if(
			ptcl->Neighbors, 
			ptcl->Neighbors + ptcl->NumberOfNeighbor, 
			[ptclCM](int j) {
				return j == ptclCM->ParticleOrder;
			}
		);

		if (newEnd != ptcl->Neighbors + ptcl->NumberOfNeighbor) {
			ptcl->NumberOfNeighbor = newEnd - ptcl->Neighbors;
			for (int j=0; i<ptclGroup->sym_int.particles.getSize(); i++) {
				Particle* members = &particles[ptclGroup->sym_int.particles[j].ParticleOrder];
				ptcl->Neighbors[ptcl->NumberOfNeighbor] = members->ParticleOrder;
				ptcl->NumberOfNeighbor++;
			}
		}
	}



	//re-do UpdateNextRegTime
	// UpdateNextRegTime(particle);

	// fprintf(binout,"total number of particles = %lu, total number of groups = %lu \n", particle.size(), GroupList.size()-1); // GroupList.size() - 1 because we are terminating it.
	// fprintf(binout,"total number of ComputationList = %lu\n", ComputationList.size());

	for (int i=0; i<ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &particles[ptclGroup->sym_int.particles[i].ParticleOrder];
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

	delete ptclCM->GroupInfo;
	ptclCM->clear();

	if (Order != global_variable.NumberOfParticle-1) {
		int beforeOrder = particles[global_variable.NumberOfParticle-1].ParticleOrder;
		int afterOrder = particles[Order].ParticleOrder;

		std::swap(particles[Order], particles[global_variable.NumberOfParticle-1]);
		particles[Order].ParticleOrder = afterOrder;
		global_variable.NumberOfParticle--;
		for (int i=0; i<global_variable.NumberOfParticle; i++) {
			Particle* ptcl = &particles[i];
			
			std::replace(
				ptcl->Neighbors,
				ptcl->Neighbors + ptcl->NumberOfNeighbor,
				beforeOrder,
				afterOrder
			);
		}
	}

	fprintf(binout,"end of Few Body Termination\n");
	fprintf(binout,"--------------------------------------\n");
	fflush(binout);

}

// Use this function when SDAR (2-body) is interrupted in the middle of its integration by stellar merger, TDE, GW merger, etc.
// current_time: interrupted time in SDAR
// next_time: intended time to integrate
void FBTermination2(int Order){

	assert(Order <= global_variable.NumberOfParticle-1);

	Particle* ptclCM = &particles[Order];

	Group* ptclGroup;
	ptclGroup	= ptclCM->GroupInfo;

	double current_time = ptclGroup->CurrentTime;

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination2.cpp... (CM PID: %d)\n\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of interrupted group (Myr): %e\n", current_time*EnzoTimeStep*1e4);

	ptclCM->CurrentBlockIrr	= ptclCM->NewCurrentBlockIrr;
	ptclCM->CurrentTimeIrr	= ptclCM->CurrentBlockIrr*time_step;
	

	// Erase ptclCM from particle vector

	ptclCM->isActive = false;

	// Activate group members again except zero mass particles!

	int ptclOrder;

	for (int i=0; i<ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* ptcl = &particles[ptclGroup->sym_int.particles[i].ParticleOrder];
		if (ptcl->Mass != 0.0) {
			ptcl->isActive = true;
			ptclOrder = ptcl->ParticleOrder;

			// Set CurrentBlock and CurrentTime for group particles.
			ptcl->CurrentBlockIrr	= ptclCM->CurrentBlockIrr; // Block to be integrated
			ptcl->CurrentBlockReg	= ptclCM->CurrentBlockReg;
			ptcl->CurrentTimeIrr	= current_time; // This will be updated later.
			ptcl->CurrentTimeReg	= ptclCM->CurrentTimeReg;

			ptcl->TimeLevelIrr		= ptclCM->TimeLevelIrr;
			ptcl->TimeLevelReg		= ptclCM->TimeLevelReg;

			CalculateAcceleration01(ptcl);
			CalculateAcceleration23(ptcl);

			double pos[Dim], vel[Dim];
			ptcl->predictParticleSecondOrder(ptclCM->CurrentTimeIrr - current_time, pos, vel);

			for (int dim=0; dim<Dim; dim++) {
				ptcl->Position[dim] =  pos[dim];
				ptcl->Velocity[dim] =  vel[dim];
			}
			ptcl->CurrentTimeIrr = ptclCM->CurrentTimeIrr;

			ptcl->calculateTimeStepReg();
		
			if (ptcl->TimeLevelReg <= ptclCM->TimeLevelReg-1 
					&& ptcl->TimeBlockReg/2+ptcl->CurrentBlockReg > ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
				ptcl->TimeLevelReg = ptclCM->TimeLevelReg-1;
			}
			else if  (ptcl->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
				ptcl->TimeLevelReg = ptclCM->TimeLevelReg+1;
			}
			else 
				ptcl->TimeLevelReg = ptclCM->TimeLevelReg;
			ptcl->TimeStepReg  = static_cast<double>(pow(2, ptcl->TimeLevelReg)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()
			ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg-time_block)); // Eunwoo: I think this was already calculated in calculateTimeStepReg()

			// ptcl->calculateTimeStepIrr2(ptcl->a_tot, ptcl->a_irr);
			ptcl->calculateTimeStepIrr();
		}
	}

	for (int i=0; i<global_variable.NumberOfParticle; i++) {
		Particle* ptcl = &particles[i];
		
		std::replace(
			ptcl->Neighbors,
			ptcl->Neighbors + ptcl->NumberOfNeighbor,
			ptclCM->ParticleOrder,
			ptclOrder
		);
	}

	//re-do UpdateNextRegTime
	// UpdateNextRegTime(particle);

	for (int i=0; i<ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &particles[ptclGroup->sym_int.particles[i].ParticleOrder];
		if (members->Mass == 0.0) continue;

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

	delete ptclCM->GroupInfo;
	ptclCM->clear();

	// Swap particle orders; last ptclCM <-> terminated ptclCM

	if (Order != global_variable.NumberOfParticle-1) {
		int beforeOrder = particles[global_variable.NumberOfParticle-1].ParticleOrder;
		int afterOrder = particles[Order].ParticleOrder;

		// std::swap(&particles[Order], &particles[global_variable.NumberOfParticle-1]);
		std::swap(particles[Order], particles[global_variable.NumberOfParticle-1]); // Eunwoo: is this right?
		particles[Order].ParticleOrder = afterOrder;
		global_variable.NumberOfParticle--;
		for (int i=0; i<global_variable.NumberOfParticle; i++) {
			Particle* ptcl = &particles[i];
			
			std::replace(
				ptcl->Neighbors,
				ptcl->Neighbors + ptcl->NumberOfNeighbor,
				beforeOrder,
				afterOrder
			);
		}
	}
	
	fprintf(binout,"end of Few Body Termination2\n");
	fprintf(binout,"--------------------------------------\n");
	fflush(binout);
}
