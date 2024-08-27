#ifndef PARTICLE_SCHEDULER_H
#define PARTICLE_SCHEDULER_H
#include <vector>
#include <iostream>
//#include "../Particle/Particle.h"
#include "../defs.h"
#include "../global.h"

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


		~ParticleScheduler() {};
};

#endif
