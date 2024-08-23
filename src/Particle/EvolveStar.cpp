#include <iostream>
#include <vector>
#include <cmath>
#include "../global.h"
//#include "../../macros_and_parameters.h"

//REAL star_feedback7(REAL t1, REAL t2, Particle *self);

REAL Particle::evolveStarMass(REAL t1, REAL t2) {

	if ((isStarEvolution == false) || (t2 == 0)) {
		return 0;
	}
	REAL dm;

	int StarParticleFeedback = 0;

	switch (StarParticleFeedback) {
		case 0:
			dm = 0.;
			isStarEvolution = false;
			break;

		case 128:
			//dm = star_feedback7(t1*EnzoTimeStep, t2*EnzoTimeStep, this);
			break;

		default:
			std::cout << "EvolveStarMass: Invalid choice" << std::endl;
	}

	return dm;
}

/*
REAL star_feedback7(REAL t1, REAL t2, Particle *self) {
	REAL dm, tau1, tau2, m1, m2;

	t1 += EnzoCurrentTime;
	t2 += EnzoCurrentTime;

	tau1 = (t1 - self->CreationTime)/self->DynamicalTime;
	tau2 = (t2 - self->CreationTime)/self->DynamicalTime;

	dm = self->InitialMass*((1+tau1)*std::exp(-tau1)-(1+tau2)*std::exp(-tau2));
	dm = -max(min(dm, self->Mass), 0.)*StarMassEjectionFraction;

	if (tau1 > 12.0) {
		return 0.;
	}

	return dm;
}

*/
