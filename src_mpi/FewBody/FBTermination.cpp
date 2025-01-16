#ifdef FEWBODY
#include "stdio.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../def.h"

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);

void FBTermination(Group* group){

	Particle* ptclCM = group->groupCM;

	ptclCM->CurrentBlockIrr	= ptclCM->NewCurrentBlockIrr;
	ptclCM->CurrentTimeIrr	= ptclCM->CurrentBlockIrr*time_step;

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination.cpp... (CM PID: %d)\n\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->CurrentTimeIrr*EnzoTimeStep*1e4);

	// I don't use the following code anymore because it might take a long time by EW 2025.1.6
	/*
	if (ptclCM->NumberOfNeighbor != 0)
		ptclCM->updateParticle();

	for (int dim=0; dim<Dim; dim++) {
		ptclGroup->sym_int.particles.cm.Position[dim] = ptclCM->Position[dim];
		ptclGroup->sym_int.particles.cm.Velocity[dim] = ptclCM->Velocity[dim];
	}
	ptclGroup->sym_int.particles.shiftToOriginFrame();
	ptclGroup->sym_int.particles.template writeBackMemberAll<Particle>();
	*/

	assert(!group->sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
	if (ptclCM->NumberOfNeighbor != 0) {
		for (int i = 0; i < group->sym_int.particles.getSize(); i++) {
			Particle* members = &group->sym_int.particles[i];
			
			particles[members->ParticleIndex].isActive = true;
			for (int dim=0; dim<Dim; dim++) {
				particles[members->ParticleIndex].Position[dim] = ptclCM->NewPosition[dim] + members->Position[dim];
				particles[members->ParticleIndex].Velocity[dim] = ptclCM->NewVelocity[dim] + members->Velocity[dim];
			}
			particles[members->ParticleIndex].Mass = members->Mass;
		}
	}
	else {
		for (int i = 0; i < group->sym_int.particles.getSize(); i++) {
			Particle* members = &group->sym_int.particles[i];

			particles[members->ParticleIndex].isActive = true;
			for (int dim=0; dim<Dim; dim++) {
				particles[members->ParticleIndex].Position[dim] = ptclCM->Position[dim] + members->Position[dim];
				particles[members->ParticleIndex].Velocity[dim] = ptclCM->Velocity[dim] + members->Velocity[dim];
			}
			particles[members->ParticleIndex].Mass = members->Mass;
		}
	}

	ptclCM->isActive = false;
	ptclCM->NewNumberOfNeighbor = 0;
	
	for (int i=0; i<group->sym_int.particles.getSize(); i++) {
		Particle* members = &particles[group->sym_int.particles[i].ParticleIndex];

		ptclCM->NewNeighbors[ptclCM->NewNumberOfNeighbor] = members->ParticleIndex;
		ptclCM->NewNumberOfNeighbor++;

		members->CMPtclIndex = -1; // added for write_out_group function by EW 2025.1.6

		// For 3-body & 4-body termination case by EW 2025.1.6
		if (group->sym_int.particles.getSize() > 2)
			members->binary_state = -1;

		members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr;
		members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
		members->CurrentTimeIrr		= ptclCM->CurrentTimeIrr;
		members->CurrentTimeReg		= ptclCM->CurrentTimeReg;
		members->NewCurrentBlockIrr	= ptclCM->NewCurrentBlockIrr;

		members->TimeLevelIrr		= ptclCM->TimeLevelIrr;
		members->TimeLevelReg		= ptclCM->TimeLevelReg;

		members->RadiusOfNeighbor = ACRadius*ACRadius; // added by EW 2025.1.16

		CalculateAcceleration01(members);
		CalculateAcceleration23(members);

		members->calculateTimeStepReg();

		if (members->TimeLevelReg <= ptclCM->TimeLevelReg-1 
				&& members->TimeBlockReg/2+members->CurrentBlockReg > ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
			members->TimeLevelReg = ptclCM->TimeLevelReg-1;
		}
		else if  (members->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
			members->TimeLevelReg = ptclCM->TimeLevelReg+1;
		}
		else 
			members->TimeLevelReg = ptclCM->TimeLevelReg;
		members->TimeStepReg  = static_cast<double>(pow(2, members->TimeLevelReg));
		members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block));

		// members->calculateTimeStepIrr2(members->a_tot, members->a_irr);
		members->calculateTimeStepIrr();
		members->NewCurrentBlockIrr = members->CurrentBlockIrr + members->TimeBlockIrr;
		members->NextBlockIrr = members->CurrentBlockIrr + members->TimeBlockIrr;
	}


	for (int i=0; i<ptclCM->NewNumberOfNeighbor; i++) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];

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

	delete group;

	fprintf(binout,"end of Few Body Termination\n");
	fprintf(binout,"--------------------------------------\n");
	fflush(binout);
}

// Use this function when SDAR (2-body) is interrupted in the middle of its integration by stellar merger, TDE, GW merger, etc.
// current_time: interrupted time in SDAR
// next_time: intended time to integrate
void FBTermination2(Group* group){

	Particle* ptclCM = group->groupCM;

	double current_time = group->CurrentTime;

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination2.cpp... (CM PID: %d)\n\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of interrupted group (Myr): %e\n", current_time*EnzoTimeStep*1e4);

	ptclCM->CurrentBlockIrr	= ptclCM->NewCurrentBlockIrr;
	ptclCM->CurrentTimeIrr	= ptclCM->CurrentBlockIrr*time_step;
	

	ptclCM->isActive = false;
	ptclCM->NewNumberOfNeighbor = 0;

	for (int i=0; i<group->sym_int.particles.getSize(); i++) {
		Particle* members = &particles[group->sym_int.particles[i].ParticleIndex];
		if (members->Mass != 0.0) {

			ptclCM->NewNeighbors[ptclCM->NewNumberOfNeighbor] = members->ParticleIndex;
			ptclCM->NewNumberOfNeighbor++;

			members->isActive = true;

			members->CMPtclIndex = -1; // added for write_out_group function by EW 2025.1.6

			// Set CurrentBlock and CurrentTime for group particles.
			members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr; // Block to be integrated
			members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
			members->CurrentTimeIrr		= current_time; // This will be updated later.
			members->CurrentTimeReg		= ptclCM->CurrentTimeReg;
			members->NewCurrentBlockIrr = ptclCM->NewCurrentBlockIrr;

			members->TimeLevelIrr		= ptclCM->TimeLevelIrr;
			members->TimeLevelReg		= ptclCM->TimeLevelReg;

			members->RadiusOfNeighbor = ACRadius*ACRadius; // added by EW 2025.1.16

			CalculateAcceleration01(members);
			CalculateAcceleration23(members);

			double pos[Dim], vel[Dim];
			members->predictParticleSecondOrder(ptclCM->CurrentTimeIrr - current_time, pos, vel);

			for (int dim=0; dim<Dim; dim++) {
				members->Position[dim] =  pos[dim];
				members->Velocity[dim] =  vel[dim];
			}
			members->CurrentTimeIrr = ptclCM->CurrentTimeIrr;

			members->calculateTimeStepReg();
		
			if (members->TimeLevelReg <= ptclCM->TimeLevelReg-1 
					&& members->TimeBlockReg/2+members->CurrentBlockReg > ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
				members->TimeLevelReg = ptclCM->TimeLevelReg-1;
			}
			else if  (members->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
				members->TimeLevelReg = ptclCM->TimeLevelReg+1;
			}
			else 
				members->TimeLevelReg = ptclCM->TimeLevelReg;
			members->TimeStepReg  = static_cast<double>(pow(2, members->TimeLevelReg));
			members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block));

			// members->calculateTimeStepIrr2(ptcl->a_tot, ptcl->a_irr);
			members->calculateTimeStepIrr();
			members->NewCurrentBlockIrr = members->CurrentBlockIrr + members->TimeBlockIrr;
			members->NextBlockIrr = members->CurrentBlockIrr + members->TimeBlockIrr;
		}
	}

	for (int i=0; i<ptclCM->NewNumberOfNeighbor; i++) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];
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

	delete group;
	
	fprintf(binout,"end of Few Body Termination2\n");
	fprintf(binout,"--------------------------------------\n");
	fflush(binout);
}

void insertNeighbors(int Order) {

	Particle* ptclCM = &particles[Order];

	for (int i=0; i<NumberOfParticle; i++) {
		Particle* ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;

		int BeforeNumberOfNeighbor = ptcl->NumberOfNeighbor;
		
		for (int j = 0; j < ptcl->NumberOfNeighbor; ++j) {
			if (ptcl->Neighbors[j] == Order) {
				ptcl->Neighbors[j] = ptcl->Neighbors[ptcl->NumberOfNeighbor - 1];
				ptcl->NumberOfNeighbor--;

				for (int k=0; i<ptclCM->NewNumberOfNeighbor; k++) {
					ptcl->Neighbors[ptcl->NumberOfNeighbor] = ptclCM->NewNeighbors[k];
					ptcl->NumberOfNeighbor++;
				}
				break;
			}
		}
	}
}
#endif