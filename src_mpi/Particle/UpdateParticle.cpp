#include "../global.h"
#include "../def.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>






void Particle::predictParticleSecondOrder(double dt, double pos[], double vel[]) {
	// Doubling check
	// temporary variables for calculation

	dt = dt*EnzoTimeStep;

	// only predict the positions if necessary
	// how about using polynomial correction here?
	for (int dim=0; dim<Dim; dim++) {
		pos[dim] = ((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim];
		vel[dim] =  (a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim];
	}
	return;
}


/*
 *  Purporse: Correct particle positions and velocities up to fourth order
 *            using a_p and d(a_p)/dt; refer to Nbody6++ manual Eq. 8.9 and 8.10
 *
 *  Date    : 2024.01.10  by Yongseok Jo
 *  Modified: 2024.01.11  by Seoyoung Kim
 *
 */
void Particle::correctParticleFourthOrder(double dt, double next_time, double pos[], double vel[], double a[3][4]) {
	double dt3,dt4,dt5;

	dt = dt*EnzoTimeStep;

	dt3 = dt*dt*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;

	// correct the predicted values positions and velocities at next_time
	// and save the corrected values to particle positions and velocities
	// the latest values of a2dot and a3dots (obtained from hermite method) are used
	for (int dim=0; dim<Dim; dim++) {
		NewPosition[dim] = pos[dim]+ a[dim][2]*dt4/24 + a[dim][3]*dt5/120;
		NewVelocity[dim] = vel[dim]+ a[dim][2]*dt3/6  + a[dim][3]*dt4/24;
	}
}


/*
void Particle::polynomialPrediction(double current_time) {

}
*/


void Particle::updateParticle() {
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] = NewPosition[dim];
		Velocity[dim] = NewVelocity[dim];
	}
	//updateTimeStep();
}



void Particle::updateRadius() {

	/* exponential (aggressive) */
	/*
		 const double c = 0.5;
		 const double b = std::log(2) / (NumNeighborMax);  // ln(2) / 40
		 double exp = a * (std::exp(b * NumberOfAC) - 1);
		 */

	/* n=2 polynomial (mild) as n increases it grows mild */

	if (NumberOfNeighbor > FixNumNeighbor) {
		const int n = 2;
		const double c = (NumNeighborMax-FixNumNeighbor);
		const double b = 0.9 / std::pow(c,n);  // ln(2) / 40
		double x = NumberOfNeighbor-FixNumNeighbor;
		double a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
		//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
		//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
		RadiusOfNeighbor *= (1.-a);
	}
	else if (NumberOfNeighbor < FixNumNeighbor) {
		const int n = 3;
		const double c = (NumNeighborMax-FixNumNeighbor);
		const double b = 0.5 / std::pow(c,n);  // ln(2) / 40
		double x = NumberOfNeighbor-FixNumNeighbor;
		double a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
		//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
		//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
		RadiusOfNeighbor *= (1.-a);
	}
}







