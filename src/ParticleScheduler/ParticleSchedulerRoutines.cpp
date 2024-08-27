#include "ParticleScheduler.h"
#include <algorithm>
#include <omp.h>
//#include <tbb/parallel_sort.h>

bool ParticleScheduler::create_level() {

	if (this->debug) {
		fprintf(stdout, "create level starts!\n");
	}

#ifdef time_trace
	_time.irr_chain.markStart();
#endif
	int num_threads = omp_get_max_threads();
	std::vector<std::vector<Level*>> local_LevelList(num_threads);
	std::vector<std::vector<ULL>> local_LevelTime(num_threads);

#pragma omp parallel 
	{
		int tid = omp_get_thread_num();
		Particle* ptcl;
		int success = 0;

#pragma omp for 
		for (int i=0; i<__particle->size(); i++) {
			ptcl =  (*__particle)[i];
			/*
			ptcl->NextBlockIrr =
			ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
			*/

			if ((ptcl->NumberOfAC != 0)
					&& (ptcl->NextBlockIrr <= NextRegTimeBlock)) {
				success = 1;
				for (int j=0; j<local_LevelList[tid].size(); j++) {
					if (local_LevelList[tid][j]->LevelTime == ptcl->NextBlockIrr) {
						local_LevelList[tid][j]->ParticleThisLevel.push_back(ptcl);
						success = 2;
					}
				}

				if (success == 1) {
					Level *level = new Level();
					level->LevelTime = ptcl->NextBlockIrr;
					level->ParticleThisLevel.push_back(ptcl);
					local_LevelList[tid].push_back(level);
					local_LevelTime[tid].push_back(ptcl->NextBlockIrr);
				}
				success = 0;
			}
		}
	}





	size_t total_size = 0;
	for (const auto &vec :local_LevelTime) {
		total_size += vec.size();
	}

	if (total_size == 0) {
#ifdef time_trace
		_time.irr_chain.markEnd();
		_time.irr_chain.getDuration();
#endif

		local_LevelList.clear();
		local_LevelList.shrink_to_fit();
		local_LevelTime.clear();
		local_LevelTime.shrink_to_fit();
		return true;
	}


	LevelTimes.resize(total_size);
	/*******************************************************
		Concatenate local_LevelTime vectors into LevelTimes
	 *******************************************************/
#pragma omp parallel
	{
#pragma omp for schedule(static) nowait
		for (size_t i = 0; i < local_LevelTime.size(); ++i) {
			size_t local_offset = 0;

			for (size_t j = 0; j < i; ++j) {
				local_offset += local_LevelTime[j].size();
			}
			std::copy(local_LevelTime[i].begin(), local_LevelTime[i].end(), LevelTimes.begin() + local_offset);
		}
	}
	/*******************************************************/



	//Step 1: Sort the vector
	std::sort(LevelTimes.begin(), LevelTimes.end());

	// Step 2: Remove duplicates
	auto last = std::unique(LevelTimes.begin(), LevelTimes.end());

	// Step 3: Erase the removed elements
	LevelTimes.erase(last, LevelTimes.end());

	

	if (this->debug) {
		fprintf(stderr, "LevelTimes: ");
		for (const auto &vec:LevelTimes) {
			fprintf(stderr, "%e, ", vec*time_step*EnzoTimeStep*1e4);
		}
		fprintf(stderr, "\n");
	}

	for (const auto &vec:LevelTimes) {
		Level *level = new Level();
		level->LevelTime = vec;
		LevelList.push_back(level);
	}

		/*******************************************************
		  Concatenate local_LevelList vectors into LevelList 
		*******************************************************/
#pragma omp parallel
	{
		int tid = omp_get_thread_num();

#pragma omp for schedule(static) nowait
		for (size_t i = 0; i < LevelTimes.size(); i++) {

			for (int j=0; j<local_LevelList.size(); j++) {
				for (int k=0; k<local_LevelList[j].size(); k++) {
					if (local_LevelList[j][k]->LevelTime == LevelList[i]->LevelTime) {

						// Reserve space in LevelList[i] to avoid reallocation
						LevelList[i]->ParticleThisLevel.reserve(
								LevelList[i]->ParticleThisLevel.size() + local_LevelList[j][k]->ParticleThisLevel.size());

						// Insert 
						LevelList[i]->ParticleThisLevel.insert(
								LevelList[i]->ParticleThisLevel.end(),
								local_LevelList[j][k]->ParticleThisLevel.begin(),
								local_LevelList[j][k]->ParticleThisLevel.end());
					}
				}
			}
		}
	}
	/*******************************************************/


	//ttb::parallel_sort(LevelList.begin(), LevelList.end(), [](Level *l1, Level *l2){
			//return l1->LevelTime < l2-> LevelTime;});
	std::sort(LevelList.begin(), LevelList.end(), [](Level *l1, Level *l2){
			return l1->LevelTime < l2-> LevelTime;});
	std::sort(LevelTimes.begin(), LevelTimes.end());
#ifdef time_trace
	_time.irr_chain.markEnd();
	_time.irr_chain.getDuration();
#endif

	if (this->debug) {
		fprintf(stderr, "NextRegTimeBlock=%.4e Myr, size of LevelList=%lu\n",
				NextRegTimeBlock*time_step*EnzoTimeStep*1e4, LevelList.size());
		for (Level* lv: LevelList)
			lv->print_level();
		fprintf(stderr, "This is it.\n\n\n");
	}

		local_LevelList.clear();
		local_LevelList.shrink_to_fit();
		local_LevelTime.clear();
		local_LevelTime.shrink_to_fit();
	return false;
}






void ParticleScheduler::run(void) {

#ifdef time_trace
	_time.irr_force.markStart();
#endif

#pragma omp parallel
	{
		Particle* ptcl;  
#pragma omp for 
		for (int i=0; i<LevelList[0]->ParticleThisLevel.size(); i++) {

			ptcl = LevelList[0]->ParticleThisLevel[i];
			ptcl->calculateIrrForce(); 
		}

#pragma omp barrier

#pragma omp for
		for (int i=0; i<LevelList[0]->ParticleThisLevel.size(); i++) {
			ptcl = LevelList[0]->ParticleThisLevel[i];
			ptcl->updateParticle();
			ptcl->CurrentBlockIrr += ptcl->TimeBlockIrr;
			ptcl->CurrentTimeIrr   = ptcl->CurrentBlockIrr*time_step;
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
			ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
		}
	}

#ifdef time_trace
	_time.irr_force.markEnd();
	_time.irr_force.getDuration();
#endif
}





bool ParticleScheduler::update_level(void) {

#ifdef time_trace
		_time.irr_sort.markStart();
#endif

	int num_threads = omp_get_max_threads();
	std::vector<std::vector<Level*>> local_LevelList(num_threads);
	std::vector<std::vector<ULL>> local_LevelTime(num_threads);

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		Particle* ptcl;  
		int success = 0;

#pragma omp for
		for (int i=0; i<LevelList[0]->ParticleThisLevel.size(); i++) {
			ptcl = LevelList[0]->ParticleThisLevel[i];
			
			if (ptcl->NextBlockIrr > NextRegTimeBlock)
				continue;

			success = 1;
			for (int j=0; j<local_LevelList[tid].size(); j++) {
				if (local_LevelList[tid][j]->LevelTime == ptcl->NextBlockIrr) {
					local_LevelList[tid][j]->ParticleThisLevel.push_back(ptcl);
					success = 2;
				}
			}

			if (success == 1) {
				Level *level = new Level();
				level->LevelTime = ptcl->NextBlockIrr;
				level->ParticleThisLevel.push_back(ptcl);
				local_LevelList[tid].push_back(level);
				local_LevelTime[tid].push_back(ptcl->NextBlockIrr);
			}
			success = 0;
		}
	}



	// Delete the memory pointed to by the first element
	delete LevelList.front(); 
	LevelList.erase(LevelList.begin());
	LevelTimes.erase(LevelTimes.begin());

	if (LevelList.size() == 0) {
#ifdef time_trace
		_time.irr_sort.markEnd();
		_time.irr_sort.getDuration();
#endif
		local_LevelList.clear();
		local_LevelList.shrink_to_fit();
		local_LevelTime.clear();
		local_LevelTime.shrink_to_fit();
		return false;
	}


	size_t total_size = 0;
	for (const auto &vec :local_LevelTime) {
		total_size += vec.size();
	}


	size_t offset = LevelTimes.size();
	LevelTimes.resize(offset+total_size);

		/*******************************************************
		  Concatenate local_LevelTime vectors into LevelTimes
		*******************************************************/
#pragma omp parallel
	{
#pragma omp for schedule(static) nowait
		for (size_t i = 0; i < local_LevelTime.size(); ++i) {
			size_t local_offset = offset;

			for (size_t j = 0; j < i; ++j) {
				local_offset += local_LevelTime[j].size();
			}
			std::copy(local_LevelTime[i].begin(), local_LevelTime[i].end(), LevelTimes.begin() + local_offset);
		}
	}
	/*******************************************************/


	//Step 1: Sort the vector
	std::sort(LevelTimes.begin(), LevelTimes.end());

	// Step 2: Remove duplicates
	auto last = std::unique(LevelTimes.begin(), LevelTimes.end());

	// Step 3: Erase the removed elements
	LevelTimes.erase(last, LevelTimes.end());
	

	if (this->debug) {
		fprintf(stderr, "LevelTimes: ");
		for (const auto &vec:LevelTimes) {
			fprintf(stderr, "%e, ", vec*time_step*EnzoTimeStep*1e4);
		}
		fprintf(stderr, "\n");


		fprintf(stderr, "offset=%lu, size of LevelList=%lu\n", offset, LevelList.size());
	}
	if (offset != LevelTimes.size()) {
		int j = 0;
		//LevelList.resize(LevelTimes.size());
		for (int i=0; i<LevelTimes.size(); i++) {
			if (LevelTimes[i] == LevelList[j]->LevelTime) {
				j++;
				continue;
			}
			else {
				Level *level = new Level();
				level->LevelTime = LevelTimes[i];
				//LevelList[offset+j] = level;
				LevelList.push_back(level);
			}
		}
	}

		/*******************************************************
		  Concatenate local_LevelList vectors into LevelList 
		*******************************************************/
#pragma omp parallel
	{
		int tid = omp_get_thread_num();

#pragma omp for schedule(static) nowait
		for (size_t i = 0; i < LevelTimes.size(); i++) {
			for (int j=0; j<local_LevelList.size(); j++) {
				for (int k=0; k<local_LevelList[j].size(); k++) {
					if (local_LevelList[j][k]->LevelTime == LevelList[i]->LevelTime) {

						// Reserve space in LevelList[i] to avoid reallocation
						LevelList[i]->ParticleThisLevel.reserve(
								LevelList[i]->ParticleThisLevel.size() + local_LevelList[j][k]->ParticleThisLevel.size());

						// Insert 
						LevelList[i]->ParticleThisLevel.insert(
								LevelList[i]->ParticleThisLevel.end(),
								local_LevelList[j][k]->ParticleThisLevel.begin(),
								local_LevelList[j][k]->ParticleThisLevel.end());
					}
				}
			}
		}
	}
	/*******************************************************/



	//ttb::parallel_sort(LevelList.begin(), LevelList.end(), [](Level *l1, Level *l2){
			//return l1->LevelTime < l2-> LevelTime;});
	std::sort(LevelList.begin(), LevelList.end(), [](Level *l1, Level *l2){
			return l1->LevelTime < l2-> LevelTime;});
	//std::sort(LevelTimes.begin(), LevelTimes.end());

#ifdef time_trace
	_time.irr_sort.markEnd();
	_time.irr_sort.getDuration();
#endif

	if (this->debug) {
		fprintf(stderr, "Update: size of LevelList=%lu\n",
				LevelList.size());
		for (Level* lv: LevelList)
			lv->print_level();
		fprintf(stderr, "This is it.\n\n\n");
	}

	local_LevelList.clear();
	local_LevelList.shrink_to_fit();
	local_LevelTime.clear();
	local_LevelTime.shrink_to_fit();
	return true;
}




