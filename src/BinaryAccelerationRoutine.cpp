#include <iostream>
#include "global.h"
#include "defs.h"

// Eunwoo edited

void NewFBInitialization(Particle* ptclI, std::vector<Particle*> &particle, std::vector<Particle*> &ComputationList);
void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle, REAL current_time, ULL current_block);

bool AddNewBinariesToList(std::vector<Particle*> &particle) {

	//fprintf(binout, "Finding new binaries ... Candidates = %lu\n", BinaryCandidateList.size());
	// add new binaries
	for (Particle *ptcl : GroupCandidateList) {
		//fprintf(binout, "next candidate\n");
		//fprintf(binout, "PID = %d, timestep=%e\n", ptcl->PID, ptcl->TimeStepIrr*EnzoTimeStep*1e4);
		// if the irregular time step is too short, check if it is binary
		//if ((ptcl->TimeStepIrr*EnzoTimeStep*1e4<KSTime) && ( (ptcl->isBinary == false) && (ptcl->isCMptcl == false) )) {
		//fprintf(binout, "BinaryAccelerationRoutine.cpp: new binary particle found! timestep=%e\n",
				//ptcl->TimeStepIrr*EnzoTimeStep*1e4);
		if (!ptcl->isGroup) {

			// fprintf(binout, "Let's check ptcl inside GroupCandidateList is really group!\n"); // Eunwoo debug

			ptcl->isFBCandidate();
			if (ptcl->isGroup) {
				std::cout << "AddNewGroups ... new group found" << std::endl;
				fprintf(binout, "BinaryAccelerationRoutine.cpp: new binary particle found!\n");
				// the initialization of the binary counterpart will be done together in the following function.
				NewFBInitialization(ptcl,particle,ComputationList);
				//fprintf(stdout, "New binary of (%d, %d) initialization finished ...\n",ptcl->PID, ptcl->BinaryPairParticle->PID);
				fprintf(binout,"After group addition, the number of particles are... %d \n",int(particle.size()));
				fflush(binout);
			}
		}
	}
	GroupCandidateList.clear();
	//fprintf(binout,"\n After binary addition, the number of particles are... %d \n",int(particle.size()));
	return true;
}

void BinaryAccelerationRoutine(REAL next_time, std::vector<Particle*> &particle) {

	int count;
	int bincount = 0;

	count = 0;

	if (next_time == 0) {
		return;
	}

	for (Group* ptclGroup: GroupList) {

		// fprintf(stderr, "time error max: %e\n", ptclGroup->manager.time_error_max); // Eunwoo debug

		ptclGroup->ARIntegration(next_time);

		count += 1;

		//fprintf(binout, "\nBinaryAccelerationRoutine.cpp: After KS Integration of %dth binary....\n", count);
		//fprintf(binout, "The ID of ith particle is %d \n",ptclBin->ptclCM->BinaryParticleI->PID);
		//fprintf(binout, "The ID of jth particle is %d \n",ptclBin->ptclCM->BinaryParticleJ->PID);
		//fflush(binout);	

		//if (bincount>0) {
		//	std::cout << "Integrating Binary ..." << std::endl;

		//	fprintf(binout, "KS coordinates - u1:%e, u2:%e, u3:%e, u4:%e\n", ptclBin->u[0], ptclBin->u[1], ptclBin->u[2], ptclBin->u[3]);
		//	fprintf(binout, "KS coordinates - udot1:%e, udot2:%e, udot3:%e, udot4:%e\n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
		//	fprintf(binout, "KS coordinates - udot1:%e, udot2:%e, udot3:%e, udot4:%e\n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
		//	fprintf(binout, "Other important KS variables - r:%e, h:%e, gamma: %e, tau:%e, step:%e, currentTime: %e \n", ptclBin->r, ptclBin->h, ptclBin->gamma, ptclBin->dTau, ptclBin->TimeStep, ptclBin->CurrentTime);
		//	fprintf(binout, "loop number = %d \n", bincount);
		//	fflush(binout);
		//}

	}
	// fprintf(stdout, "BinaryAccelerationRoutine ended, current time: %e\n", next_time); // Eunwoo added for debug
}
