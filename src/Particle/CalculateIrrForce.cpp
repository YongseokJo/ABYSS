#include "Particle.h"
#include <vector>
#include <iostream>
#include <cmath>


void direct_sum(double *x, double *v, double r2, double vx,
	 	        double mass, double a[3], double adot[3]);


/*
 *  Purporse: Calculate the and update the irregular acceleration of particles
 *            using a_0, d(a_0)/dt, a_p, d(a_p)/dt; Nbody6++ manual Eq. 8.9 and 8.10
 * 			  Note that a_p, d(a_p)/dt are calculated from predicted positions and velocities
 *
 *  Date    : 2024.01.09  by Seoyoung Kim
 *  Modified: 2024.01.11  by Seoyoung Kim
 *
 */

void Particle::calculateIrrForce() {

	if (this->NumberOfAC == 0) {
		CurrentTimeIrr += TimeStepIrr;
		return;
	}

	double dt;
	double tirr[2]; // 0 for current and 1 for advanced times

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a0_irr[2][Dim], a0dot_irr[2][Dim]; // 0 for current and 1 for predicted accelerations


	dt      = TimeStepIrr; // interval of time step
	tirr[0] = CurrentTimeIrr; // current time of particle
	tirr[1] = CurrentTimeIrr + TimeStepIrr; // the time to be advanced to



	// initialize irregular force terms for ith particle just in case
	for (int i=0; i<2; i++) {
		for (int dim=0; dim<Dim; dim++){
			a0_irr[i][dim]    = 0.0;
			a0dot_irr[i][dim] = 0.0;
		}
	}

	// scan over the neighbor lists to find the irregular acceleration components
	std::cout <<  "Looping single particles to calculate irregular acceleration...\n" << std::flush;

	for (Particle* ptcl: ACList) {

		// reset temporary variables at the start of a new calculation
		for (int i=0; i<2; i++) {
			r2 = 0.0;
			vx = 0.0;

			ptcl->predictParticleSecondOrder(tirr[i]);
			this->predictParticleSecondOrder(tirr[i]);
			for (int dim=0; dim<Dim; dim++) {
				// calculate position and velocity differences for current time
				x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
				v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];

				// calculate the square of radius and inner product of r and v for each case
				r2 += x[dim]*x[dim];
				vx += v[dim]*x[dim];
			}
			// add the contribution of jth particle to acceleration of current and predicted times
			direct_sum(x ,v, r2, vx, ptcl->Mass, a0_irr[i], a0dot_irr[i]);
		}
	} // end of scanning over neighbor lists


	// correct the force using the 4th order hermite method
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4;

	// set the values of calculation variables
	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;

	// calculate the correction terms and ...
	// calculate the final acceleration using correction terms

	for (int dim=0; dim<Dim; dim++) {
		da_dt2  = (   a0_irr[0][dim] - a0_irr[1][dim]   ) / dt2;
		adot_dt = (a0dot_irr[0][dim] + a0dot_irr[1][dim]) / dt;

		a2 =  12*da_dt2 + 6*adot_dt;
		a3 = (-6*da_dt2 - 2*adot_dt)/dt;

		a_irr[dim][0] = a0_irr[0][dim]; 
		a_irr[dim][1] = a0dot_irr[0][dim]; 
		a_irr[dim][2] = a2; 
		a_irr[dim][3] = a3;

		a_tot[dim][0] = a_reg[dim][0] + a_irr[dim][0];
		a_tot[dim][1] = a_reg[dim][1] + a_irr[dim][1];
		a_tot[dim][2] = a_reg[dim][2] + a_irr[dim][2];
		a_tot[dim][3] = a_reg[dim][3] + a_irr[dim][3];
	}

	// update the current irregular time and irregular time steps
	CurrentTimeIrr = tirr[1];
	this->calculateTimeStepIrr(a_tot); // calculate irregular time step based on total force
	this->updateParticle();
}


