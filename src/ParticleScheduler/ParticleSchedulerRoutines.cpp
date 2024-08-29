#include "ParticleScheduler.h"
#include <algorithm>
#include <omp.h>
//#include <tbb/parallel_sort.h>

bool ParticleScheduler::create_level() {

	if (this->debug) {
		fprintf(stdout, "create level starts!\n");
		fflush(stdout);
	}

#ifdef time_trace
	_time.irr_chain.markStart();
#endif
	//int num_threads = omp_get_max_threads();

	Particle* ptcl;
	int size=0;

	for (int i=0; i<__particle->size(); i++) {
		ptcl =  (*__particle)[i];

		if ((ptcl->NumberOfAC != 0) && (ptcl->NextBlockIrr <= NextRegTimeBlock)) {
			//fprintf(stdout, "PID=%d, NBI=%llu\n", ptcl->PID, ptcl->NextBlockIrr);
			if (!skiplist->search(ptcl->NextBlockIrr, ptcl))
				skiplist->insert(ptcl->NextBlockIrr, ptcl);
			size++;
		}
	}


#ifdef time_trace
	_time.irr_chain.markEnd();
	_time.irr_chain.getDuration();
#endif

	if (this->debug) {
		skiplist->display();
	}

	if (skiplist->getFirstNode() == nullptr)
		return true;
	else
		return false;
}






void ParticleScheduler::run(void) {

#ifdef time_trace
	_time.irr_force.markStart();
#endif

	//fprintf(stdout, "ParticleScheduler run starts!\n");
	//fflush(stdout);

	//Particle* ptcl;
	Node* ThisLevelNode = skiplist->getFirstNode();

#ifdef OMP
#pragma omp parallel
	{
		Particle* ptcl;
#pragma omp for
		for (int i=0; i<ThisLevelNode->particle_list.size(); i++) {
			ptcl = ThisLevelNode->particle_list[i];
			ptcl->calculateIrrForce();
		}

#pragma barrier

#pragma omp for
		for (int i=0; i<ThisLevelNode->particle_list.size(); i++) {
			ptcl = ThisLevelNode->particle_list[i];
			ptcl->updateParticle();
			ptcl->CurrentBlockIrr += ptcl->TimeBlockIrr;
			ptcl->CurrentTimeIrr   = ptcl->CurrentBlockIrr*time_step;
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
			ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
		}
	}
#else
	Particle* ptcl;
	for (int i=0; i<ThisLevelNode->particle_list.size(); i++) {
		ptcl = ThisLevelNode->particle_list[i];
		ptcl->calculateIrrForce();
	}

	for (int i=0; i<ThisLevelNode->particle_list.size(); i++) {
		ptcl = ThisLevelNode->particle_list[i];
		ptcl->updateParticle();
		ptcl->CurrentBlockIrr += ptcl->TimeBlockIrr;
		ptcl->CurrentTimeIrr   = ptcl->CurrentBlockIrr*time_step;
		ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
		ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
	}
#endif

	//fprintf(stdout, "all done!\n");
	//fflush(stdout);

#ifdef time_trace
	_time.irr_force.markEnd();
	_time.irr_force.getDuration();
#endif
}





bool ParticleScheduler::update_level(void) {

	if (this->debug) {
		fprintf(stdout, "update level starts!\n");
		fflush(stdout);
	}

#ifdef time_trace
	_time.irr_sort.markStart();
#endif


	/* Update New Time Steps */

	Particle* ptcl;
	Node* ThisLevelNode = skiplist->getFirstNode();

	if (this->debug) {
		fprintf(stdout, "update size = %lu\n", ThisLevelNode->particle_list.size());
		fflush(stdout);
	}

	for (int i=0; i<ThisLevelNode->particle_list.size(); i++) {
		ptcl = ThisLevelNode->particle_list[i];

		/*
		if (this->debug) {
		fprintf(stdout, "PID=%d, NBI=%llu, size=%lu\n", ptcl->PID, ptcl->NextBlockIrr, ThisLevelNode->particle_list.size());
			fprintf(stdout, "NextBlockIrr=%llu\n",ptcl->NextBlockIrr);
			fflush(stdout);
		}
		*/

		if (ptcl->NextBlockIrr > NextRegTimeBlock)
			continue;

		//if (!skiplist->parallel_search(ptcl->NextBlockIrr, ptcl))
		if (!skiplist->search(ptcl->NextBlockIrr, ptcl))
			skiplist->insert(ptcl->NextBlockIrr, ptcl);
	}

	if (this->debug) {
		fprintf(stdout, "For loop done!\n");
		fflush(stdout);
	}


	skiplist->deleteFirstNode();


	if (this->debug) {
		fprintf(stdout, "update level in the middle!\n");
		fflush(stdout);
	}

#ifdef time_trace
	_time.irr_sort.markEnd();
	_time.irr_sort.getDuration();
#endif

	if (skiplist->getFirstNode() == nullptr)
		return false;



	if (this->debug) {
		skiplist->display();
		fprintf(stderr, "This is it.\n\n\n");
	}

	return true;
}




