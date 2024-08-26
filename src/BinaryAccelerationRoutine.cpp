#include <iostream>
#include "global.h"
#include "defs.h"

// Eunwoo edited

void NewFBInitialization(Particle* ptclI, std::vector<Particle*> &particle, std::vector<Particle*> &ComputationList);
void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle);

bool AddNewBinariesToList(std::vector<Particle*> &particle) {

	//fprintf(binout, "Finding new binaries ... Candidates = %lu\n", BinaryCandidateList.size());
	// add new binaries
	for (Particle *ptcl : GroupCandidateList) {
		//fprintf(binout, "next candidate\n");
		//fprintf(binout, "PID = %d, timestep=%e\n", ptcl->PID, ptcl->TimeStepIrr*EnzoTimeStep*1e4);

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

void BinaryAccelerationRoutine(REAL next_time) {


	if (next_time == 0) {
		return;
	}

	for (Group* ptclGroup: GroupList) {

		ptclGroup->ARIntegration(next_time);

	}
	// fprintf(stdout, "BinaryAccelerationRoutine ended, current time: %e\n", next_time); // Eunwoo added for debug
}
