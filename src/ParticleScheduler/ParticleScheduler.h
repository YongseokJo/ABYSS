#ifndef PARTICLE_SCHEDULER_H
#define PARTICLE_SCHEDULER_H
#include <vector>
#include <iostream>
//#include "../Particle/Particle.h"
#include "../defs.h"
#include "../global.h"
#include "SkipList.h"

/*
class Level {
	private:

	public:
		ULL LevelTime;
		Pvector ParticleThisLevel;
		Level* NextLevel;



		Level(void) {
			LevelTime = 0;
			NextLevel = nullptr;
			ParticleThisLevel.clear();
		};

		Level(Particle* ptcl) {
			ParticleThisLevel.clear();
			NextLevel = nullptr;
			LevelTime = ptcl->NextBlockIrr;
			ParticleThisLevel.push_back(ptcl);
		};

		void print_level() {
			fprintf(stderr,"LevelTime=%.4e Myr\n", LevelTime*time_step*EnzoTimeStep*1e4);
			for (Particle* ptcl:ParticleThisLevel) {
				fprintf(stderr,"PID=%d, NextBlockIrr=%.4e Myr\n", ptcl->PID, ptcl->NextBlockIrr*time_step*EnzoTimeStep*1e4);
			}
		};

		~Level(void) {
			ParticleThisLevel.clear();
			ParticleThisLevel.shrink_to_fit();
		};
};


class ParticleScheduler {
	private:
		Pvector *__particle;
		Uvector LevelTimes;
		std::vector<Level*> LevelList;
		bool debug = false;

	public:
		bool isEnd;

		ParticleScheduler(Pvector &particle) {
			this->__particle = &particle;
			isEnd = create_level();
		};

		void run(void);
		bool update_level(void);
		bool create_level(void);


		~ParticleScheduler() {
			LevelList.clear();
			LevelList.shrink_to_fit();
		};
};
*/



class ParticleScheduler {
	private:
		Pvector *__particle;
		bool debug = false;

	public:
		int max_level = 5;
		double prob = 0.5;
		SkipList *skiplist;

		ParticleScheduler(Pvector *particle) {
			__particle = particle;
			skiplist = new SkipList(max_level, prob);
		};

		void run(void);
		bool update_level(void);
		bool create_level(void);


		~ParticleScheduler() {
			//delete skiplist;
		};
};

#endif

