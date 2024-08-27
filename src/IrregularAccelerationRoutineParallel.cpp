#include <iostream>
#include "global.h"
#include <cassert>
#ifdef OMP
#include <omp.h>
#include <sched.h>
#endif
#include "ParticleScheduler/ParticleScheduler.h"


bool IrregularAccelerationRoutineParallel(std::vector<Particle*> &particle)
{
	fprintf(stdout, "Irr starts!\n");
	ParticleScheduler ptclSchdlr(particle);

	if (ptclSchdlr.isEnd)
		return true;

	do {
	ptclSchdlr.run();
	} while (ptclSchdlr.update_level());
	fprintf(stdout, "Irr ends!\n");

	return true;
}
