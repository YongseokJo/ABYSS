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
	fflush(stdout);
	ParticleScheduler ptclSchdlr(&particle);

	if (ptclSchdlr.create_level())
		return true;

	do {

		ptclSchdlr.run();

	} while (ptclSchdlr.update_level());

	fprintf(stdout, "Irr ends!\n");

	return true;
}
