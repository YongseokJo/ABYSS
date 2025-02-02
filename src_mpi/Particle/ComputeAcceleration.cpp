#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include "../global.h"
#include "../def.h"
#include <unordered_set>


void calculateSingleAcceleration(Particle *ptcl2, double *pos, double *vel, double (&a)[3], double (&adot)[3], int sign);

void Particle::computeAccelerationIrr() {

	if (this->NumberOfNeighbor == 0) {
		for (int dim=0; dim<Dim; dim++){
			this->NewPosition[dim] = this->Position[dim];
			this->NewVelocity[dim] = this->Velocity[dim];
		}
		return;
	}

	double dt, mdot, epsilon=1e-6;
	double new_time; // 0 for current and 1 for advanced times

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a_tmp[Dim], adot_tmp[Dim]; // 0 for current and 1 for predicted accelerations
	double pos[Dim], vel[Dim];
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double m_r3;
	Particle* ptcl;
	new_time = this->CurrentTimeIrr + this->TimeStepIrr; // the time to be advanced to
	dt       = this->TimeStepIrr*EnzoTimeStep; // interval of time step


	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<Dim; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion Calculation
	 ********************************************************/
	this->predictParticleSecondOrder(this->TimeStepIrr, pos, vel);

	for (int i=0; i<this->NumberOfNeighbor; i++) {

		ptcl = &particles[this->Neighbors[i]];

		if (!ptcl->isActive) {
			if (ptcl->CMPtclIndex != -1) {
				CMPtclsSet.insert(ptcl->CMPtclIndex);
				continue;
			}
			else {
				if (ptcl->Mass == 0) continue;
				else {
					fprintf(stderr, "Why not zero mass? this PID: %d, neighbor PID: %d\n", this->PID, ptcl->PID);
					assert(ptcl->Mass == 0);
				}
			}
		}


		/*
		if (ptcl->isCMptcl) {
			fprintf(stderr, "my = %d , pid of cm = %d\n", this->PID, ptcl->PID);
			fflush(stderr);
		}

	 if (ptcl->Position[0]!=ptcl->Position[0]) {
			fprintf(stderr, "Nan occurs, %lf", ptcl->Position[0]);
			fflush(stderr);
			assert(this->Position[0] ==  this->Position[0]);
			exit(EXIT_FAILURE);
	 }
		if (ptcl->PID == this->PID)  {
			fprintf(stderr, "Myself in neighbor (%d)", PID);
			fflush(stderr);
			exit(EXIT_FAILURE);
			continue;
		}
		*/

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		//mdot = ptcl->evolveStarMass(CurrentTimeIrr,
				//CurrentTimeIrr+TimeStepIrr*1.01)/TimeStepIrr*1e-2; // derivative can be improved
																													 //
																													 // add the contribution of jth particle to acceleration of current and predicted times

		m_r3 = ptcl->Mass/(r2*sqrt(r2));

		for (int dim=0; dim<Dim; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	} // endfor ptcl

	if (!CMPtclsSet.empty()) {
		for (int i: CMPtclsSet) {
			ptcl = &particles[i];
			if (!ptcl->isActive) {
				fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", this->PID, ptcl->PID);
				assert(ptcl->isActive);
			}

			// reset temporary variables at the start of a new calculation
			r2 = 0.0;
			vx = 0.0;

			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

			for (int dim=0; dim<Dim; dim++) {
				// calculate position and velocity differences for current time
				x[dim] = pos_neighbor[dim] - pos[dim];
				v[dim] = vel_neighbor[dim] - vel[dim];

				// calculate the square of radius and inner product of r and v for each case
				r2 += x[dim]*x[dim];
				vx += v[dim]*x[dim];
			}

			//mdot = ptcl->evolveStarMass(CurrentTimeIrr,
					//CurrentTimeIrr+TimeStepIrr*1.01)/TimeStepIrr*1e-2; // derivative can be improved
																														//
																														// add the contribution of jth particle to acceleration of current and predicted times

			m_r3 = ptcl->Mass/(r2*sqrt(r2));

			for (int dim=0; dim<Dim; dim++){
				a_tmp[dim]    += m_r3*x[dim];
				adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
	}


	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;
	double dt_ex = (new_time - this->CurrentTimeReg)*EnzoTimeStep;

	double A, B, C;

	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	//dt5 = dt4*dt;

	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	for (int dim=0; dim<Dim; dim++) {


#define noOptimization // I don't think it works 
#ifdef Optimization
		// do the higher order correcteion


		A = -12   *( this->a_irr[dim][0]   - a_tmp[dim] );
		B = -4*dt *( 2*this->a_irr[dim][1] + adot_tmp[dim] );
		C =  6*dt *( this->a_irr[dim][1]   + adot_tmp[dim] );

		a2 = dt *(A+B)/48;
		a3 = dt *(C-A)/120;

		//fprintf(stderr, "da_dt2=%.2e, adot_dt=%.2e, a2=%.2e, a3=%.2e\n", da_dt2, adot_dt, a2, a3);

		// 4th order correction
		// save the values in the temporary variables
		this->NewPosition[dim] = pos[dim] + a2*dt + a3*dt;
		this->NewVelocity[dim] = vel[dim] + 4*a2  + 5*a3;


		// note that these higher order terms and lowers have different neighbors
		this->a_irr[dim][0] = a_tmp[dim];
		this->a_irr[dim][1] = adot_tmp[dim];
		this->a_irr[dim][2] = a2*24/dt3;
		this->a_irr[dim][3] = a3*120/dt4;
#else
		// do the higher order correcteion
		da_dt2  = (this->a_irr[dim][0] - a_tmp[dim]) / dt2; 
		adot_dt = (this->a_irr[dim][1] + adot_tmp[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*this->a_irr[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		// 4th order correction
		// save the values in the temporary variables
		this->NewPosition[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->NewVelocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;

		// note that these higher order terms and lowers have different neighbors
		this->a_irr[dim][0] = a_tmp[dim];
		this->a_irr[dim][1] = adot_tmp[dim];
		this->a_irr[dim][2] = a2;
		this->a_irr[dim][3] = a3;
#endif
	}



	for (int dim=0; dim<Dim; dim++) {
		this->a_tot[dim][0] = this->a_reg[dim][0] + this->a_irr[dim][0] + this->a_reg[dim][1]*dt_ex; // affect the next
		this->a_tot[dim][1] = this->a_reg[dim][1] + this->a_irr[dim][1];
		this->a_tot[dim][2] = this->a_reg[dim][2] + this->a_irr[dim][2];
		this->a_tot[dim][3] = this->a_reg[dim][3] + this->a_irr[dim][3];
	}


	/*
	if (this->NumberOfNeighbor == 0) {
		//CurrentTimeIrr += TimeStepIrr;
		std::cout << "Error: No neighbor in Irregular force!!" << std::endl;
		return;
	}
	*/
}






void Particle::computeAccelerationReg() {

	double dt, mdot, epsilon=1e-6;
	double new_time; // 0 for current and 1 for advanced times

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a[Dim], adot[Dim], a_new[Dim], adot_new[Dim]; // 0 for current and 1 for predicted accelerations
	double pos[Dim], vel[Dim];
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double m_r3;
	int    j=0;
	Particle* ptcl;
	new_time = this->CurrentTimeReg+this->TimeStepReg; // the time to be advanced to
	dt       = this->TimeStepReg*EnzoTimeStep; // interval of time step
	this->NewNumberOfNeighbor = 0;

	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<Dim; dim++){
		a[dim]        = 0.0;
		adot[dim]     = 0.0;
		a_new[dim]    = 0.0;
		adot_new[dim] = 0.0;
		this->a_irr[dim][0] = 0.;
		this->a_irr[dim][1] = 0.;
	}



	/*******************************************************
	 * Regular Acceleartion Calculation
	 ********************************************************/
	if (this->NumberOfNeighbor == 0)
		this->predictParticleSecondOrder(this->TimeStepReg, pos, vel);
	else
		this->predictParticleSecondOrder(0, pos, vel);


	for (int i=0; i<=global_variable->LastParticleIndex; i++) {
		ptcl = &particles[i];

		if (!ptcl->isActive) continue;

		if (this->PID == ptcl->PID)
			continue;


		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, pos_neighbor, vel_neighbor);
		else
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		//mdot = ptcl->evolveStarMass(CurrentTimeIrr,
		//CurrentTimeIrr+TimeStepIrr*1.01)/TimeStepIrr*1e-2; // derivative can be improved
		//
		// add the contribution of jth particle to acceleration of current and predicted times

		m_r3 = ptcl->Mass/r2/sqrt(r2);


		//std::cout << "PIDs=" <<  this->Neighbors[j] << ', ' << ptcl->PID << NumberOfNeighbor<< std::endl;
		//if (this->Neighbors[j] == ptcl->PID) {
		if (j < this->NumberOfNeighbor && this->Neighbors[j] == ptcl->ParticleIndex) {
			//std::cout << this->PID << ", PIDs=" <<  this->Neighbors[j] << ", " << ptcl->PID << std::endl;
			j++;
		} 
		else {
			for (int dim=0; dim<Dim; dim++){
				a[dim]    += m_r3*x[dim];
				adot[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}


		if (r2 < this->RadiusOfNeighbor) {
			this->NewNeighbors[this->NewNumberOfNeighbor] = ptcl->ParticleIndex;
			this->NewNumberOfNeighbor++;
			for (int dim=0; dim<Dim; dim++){
				this->a_irr[dim][0] += m_r3*x[dim];
				this->a_irr[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++){
				a_new[dim]    += m_r3*x[dim];
				adot_new[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
	} // endfor ptcl

	//std::cout << this->PID << " Number of Neighbor=" << j << ", " << NumberOfNeighbor<< std::endl;
	assert(j == this->NumberOfNeighbor);
	//std::cout << "Number of Neighbor=" << j << ", " << NumberOfNeighbor<< std::endl;
	//std::cout << "New Number of Neighbor=" << NewNumberOfNeighbor<< std::endl;


	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;

	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;
	for (int dim=0; dim<Dim; dim++) {

		// do the higher order correcteion
		da_dt2  = (this->a_reg[dim][0] - a[dim]) / dt2;
		adot_dt = (this->a_reg[dim][1] + adot[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*this->a_reg[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		//fprintf(stderr, "DIM=%d, pid=%d, a=%.2e, da_dt2=%.2e, adot_dt=%.2e, a2=%.2e, a3=%.2e\n", dim,this->PID,a[dim], da_dt2, adot_dt, a2, a3);

		// 4th order correction
		// save the values in the temporary variables
		this->NewPosition[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->NewVelocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;


		this->a_reg[dim][2] = a2;
		this->a_reg[dim][3] = a3;
	}


	for (int dim=0; dim<Dim; dim++) {
		this->a_reg[dim][0] = a_new[dim];
		this->a_reg[dim][1] = adot_new[dim];
		this->a_tot[dim][0] = this->a_reg[dim][0] + this->a_irr[dim][0]; 
		this->a_tot[dim][1] = this->a_reg[dim][1] + this->a_irr[dim][1];
		if (this->NewNumberOfNeighbor == 0) {
			this->a_tot[dim][2] = this->a_reg[dim][2];
			this->a_tot[dim][3] = this->a_reg[dim][3];
		}
	}



	/*
	std::cout << "\ntotal acceleartion\n" << std::flush;
	for (int order=0; order<HERMITE_ORDER; order++) {
		for (int dim=0; dim<Dim; dim++)	 {
			std::cout << a_tot[dim][order] << " ";
		}
		std::cout << std::endl;
	} // endfor dim
		//
	std::cout << "\nreg acceleartion\n" << std::flush;
	for (int order=0; order<HERMITE_ORDER; order++) {
		for (int dim=0; dim<Dim; dim++)	 {
			std::cout << a_reg[dim][order] << " ";
		}
		std::cout << std::endl;
	} // endfor dim
	std::cout << std::endl;

	std::cout << "\nirr acceleartion\n" << std::flush;
	for (int order=0; order<HERMITE_ORDER; order++) {
		for (int dim=0; dim<Dim; dim++)	 {
			std::cout << a_irr[dim][order] << " ";
		}
		std::cout << std::endl;
	} // endfor dim
	std::cout << std::endl;

	*/
	// *position_unit/time_unit/time_unit

		 //std::cout << "\nIrregular Calculation\n" << std::flush;
		 //std::cout <<  "3. a_irr= "<< a_irr[0][0]<< ',' << a_irr[1][0]\
		 //<< ',' << a_irr[2][0] << std::endl;
		 //std::cout <<  "4. a_irr= "<< a_irr[0][0]<< ',' << a_irr[1][0]\
		 << ',' << a_irr[2][0] << std::endl;
		 //std::cout <<  "5. a_irr= "<< a_irr[0][0]<< ',' << a_irr[1][0]\
		 << ',' << a_irr[2][0] << std::endl;

	// update the current irregular time and irregular time steps
	//this->updateParticle((CurrentTimeIrr+TimeStepIrr)*EnzoTimeStep, a_irr);
	//this->updateParticle(CurrentTimeIrr, CurrentTimeIrr+TimeStepIrr, a_tot);
	//this->correctParticleFourthOrder(CurrentTimeIrr, CurrentTimeIrr+TimeStepIrr, a_irr);


	//this->updateParticle();
	//CurrentTimeIrr += TimeStepIrr; // in sorting
	//this->calculateTimeStepIrr(a_tot, a_irr); // calculate irregular time step based on total force

}




// Modified by EW 2025.1.30

void Particle::updateRegularParticleCuda(int *NewNeighborsGPU, int NewNumberOfNeighborGPU, double *new_a, double *new_adot) {
/*
	std::cerr <<  "in here!" << std::endl;
	int NeighborIndex;
	std::cerr <<  "Compute: MyPID=" <<  this->PID;
	std::cerr <<  "(" << NumberOfNeighbor << ") " << std::endl;
	std::cerr <<  "NeighborIndex = ";
	for (int j=0;  j<NumberOfNeighbor; j++) {
		NeighborIndex = Neighbors[j];  // gained neighbor particle (in next time list)
		std::cerr <<  NeighborIndex << "  (" << particles[NeighborIndex].PID << "), ";
		//std::cout <<  particles[NeighborIndex].PID << ", ";
	}
	std::cerr << std::endl;

	std::cerr <<  "(" << NewNumberOfNeighborGPU << ") " << std::endl;
	std::cerr <<  "NewNeighborIndex = ";
	for (int j=0;  j<NewNumberOfNeighborGPU; j++) {
		NeighborIndex = NewNeighborsGPU[j];  // gained neighbor particle (in next time list)
		std::cerr <<  NeighborIndex << "  (" << particles[NeighborIndex].PID << "), ";
		//std::cout <<  particles[NeighborIndex].PID << ", ";
	}
	std::cerr << std::endl;
*/
/*
	if (this->PID == 28150) {
		std::cerr << "NumberOfNeighbor: " << this->NumberOfNeighbor << "\n" << "NewNumberOfNeighborGPU: " << NewNumberOfNeighborGPU << std::endl;
		for (int j=0;  j<NewNumberOfNeighborGPU; j++) {
			int NeighborIndex = NewNeighborsGPU[j];  // gained neighbor particle (in next time list)
			std::cerr <<  NeighborIndex << "  (" << particles[NeighborIndex].PID << "), ";
			//std::cout <<  particles[NeighborIndex].PID << ", ";
		}
		std::cerr << std::endl << std::endl;
	}
*/
	double new_time = this->CurrentTimeReg+this->TimeStepReg;
	double pos[Dim], vel[Dim];
	if (this->NumberOfNeighbor == 0)
		this->predictParticleSecondOrder(this->TimeStepReg, pos, vel);
	else
		this->predictParticleSecondOrder(0, pos, vel);

	double a_tmp[Dim];
	double adot_tmp[Dim];

	for (int dim=0; dim<Dim; dim++) {
		a_tmp[dim]          = 0.;
		adot_tmp[dim]       = 0.;
		this->a_irr[dim][0] = 0.;
		this->a_irr[dim][1] = 0.;
	}


	std::unordered_set<int> hashTableOld;
	std::unordered_set<int> hashTableNew;
	
	int RealNeighbors[MaxNumberOfNeighbor]; // this->Neighbors is containing members, not CM ptcls, but this is containing CM ptcls, not members
	int RealNumberOfNeighbor = 0;

	int size = this->NumberOfNeighbor > NewNumberOfNeighborGPU ? this->NumberOfNeighbor : NewNumberOfNeighborGPU;

	for (int i=0; i<size; i++) {
		if (i < NewNumberOfNeighborGPU) {
			NewNeighborsGPU[i] = ActiveIndexToOriginalIndex[NewNeighborsGPU[i]];
			// assert(particles[NewNeighborsGPU[i]].isActive); // for debugging by EW 2025.1.23
			// /* // for debugging by EW 2025.1.23
			if (!particles[NewNeighborsGPU[i]].isActive) {
				fprintf(stderr, "In GPU, this PID: %d, inActive PID: %d\n", this->PID, particles[NewNeighborsGPU[i]].PID);
				//fprintf(stderr, "ActiveIndexToPID[%d]=%d\n", this->PID, ActiveIndexToPID[this->PID]);
				assert(particles[NewNeighborsGPU[i]].isActive); // for debugging by EW 2025.1.23
			}
			// */
			hashTableNew.insert(NewNeighborsGPU[i]);
		}
		if (i < this->NumberOfNeighbor) {
			if (particles[this->Neighbors[i]].isActive) {
				hashTableOld.insert(this->Neighbors[i]);
				RealNeighbors[RealNumberOfNeighbor++] = this->Neighbors[i];
			}
			else {
				if (particles[this->Neighbors[i]].CMPtclIndex != -1) {
					if (hashTableOld.find(particles[this->Neighbors[i]].CMPtclIndex) == hashTableOld.end()) {
						hashTableOld.insert(particles[this->Neighbors[i]].CMPtclIndex);
						RealNeighbors[RealNumberOfNeighbor++] = particles[this->Neighbors[i]].CMPtclIndex;
					}
				}
				else if (particles[this->Neighbors[i]].Mass != 0) {
					fprintf(stderr, "Why not zero mass? this PID: %d, neighbor PID: %d\n", this->PID, particles[this->Neighbors[i]].PID);
					assert(particles[this->Neighbors[i]].Mass == 0);
				}
			}
		}
	}
	assert(RealNumberOfNeighbor == hashTableOld.size()); // for debugging by EW 2025.1.30

	size = hashTableNew.size() > hashTableOld.size() ? hashTableNew.size() : hashTableOld.size();

	Particle* ptcl;
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double dx[Dim], dv[Dim];
	double dr2;
	double dxdv;
	double m_r3;


	// Aceeleration correction
	for (int i=0; i<size; i++) {
		if ( i < RealNumberOfNeighbor ) {
				
			// neighbor in old but not in new
		 	if ( hashTableNew.find(RealNeighbors[i]) == hashTableNew.end() ) {
				//fprintf(stderr, "in old, not in new = %d\n",this->Neighbors[i]);
				//std::cerr <<  "in old, not in new =" <<  this->Neighbors[i] << std::endl;
				ptcl = &particles[RealNeighbors[i]];
				/* // for debugging by EW 2025.1.23
				if (!ptcl->isActive) {
					fprintf(stderr, "this PID: %d, inActive PID: %d\n", this->PID, ptcl->PID);
					assert(ptcl->isActive); // for debugging by EW 2025.1.23
				}
				*/

				if (ptcl->NumberOfNeighbor == 0)
					ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, pos_neighbor, vel_neighbor);
				else
					ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

				dr2  = 0.0;
				dxdv = 0.0;
				for (int dim=0; dim<Dim; dim++) {
					dx[dim] = pos_neighbor[dim] - pos[dim];
					dv[dim] = vel_neighbor[dim] - vel[dim];
					dr2    += dx[dim]*dx[dim];
					dxdv   += dx[dim]*dv[dim];
				}

				m_r3 = ptcl->Mass/dr2/sqrt(dr2);

				for (int dim=0; dim<Dim; dim++){
					a_tmp[dim]    -= m_r3*dx[dim];
					adot_tmp[dim] -= m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
				}
			}
		}



		if ( i < NewNumberOfNeighborGPU ) {
			ptcl = &particles[NewNeighborsGPU[i]];
			/* // for debugging by EW 2025.1.23
			if (!ptcl->isActive) {
				fprintf(stderr, "this PID: %d, inActive PID: %d\n", this->PID, ptcl->PID);
				assert(ptcl->isActive); // for debugging by EW 2025.1.23
			}
			*/

			if (ptcl->NumberOfNeighbor == 0)
				ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, pos_neighbor, vel_neighbor);
			else
				ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

			dr2  = 0.0;
			dxdv = 0.0;
			for (int dim=0; dim<Dim; dim++) {
				dx[dim] = pos_neighbor[dim] - pos[dim];
				dv[dim] = vel_neighbor[dim] - vel[dim];
				dr2    += dx[dim]*dx[dim];
				dxdv   += dx[dim]*dv[dim];
			}

			m_r3 = ptcl->Mass/dr2/sqrt(dr2);

			for (int dim=0; dim<Dim; dim++){
				a_irr[dim][0] += m_r3*dx[dim];
				a_irr[dim][1] += m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
			}

			// neighbor in new but not in old
			if ( hashTableOld.find(NewNeighborsGPU[i]) == hashTableOld.end() ) {
				//fprintf(stderr, "in new, not in old = %d\n",this->Neighbors[i]);
				//std::cerr <<  "in new, not in old =" <<  NewNeighbors[i] << std::endl;
				for (int dim=0; dim<Dim; dim++){
					a_tmp[dim]    += m_r3*dx[dim];
					adot_tmp[dim] += m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
				}
			}
		}
	}



	/*******************************************************
	 * Acceleartion correction according to past neighbor
	 ********************************************************/


	//fprintf(stderr,"%d Neighbor Correction new=%d, old=%d\n", ptcl->PID, NumNeighborReceive[i], ptcl->NumberOfAC);
	/*
		 if (NumNeighborReceive[i]>NumNeighborMax) {
		 std::cerr <<  "MyPID=" <<  ptcl->PID << ", NN=" << NumNeighborReceive[i] << std::endl;
		 }
		 */

	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	double dt  = this->TimeStepReg*EnzoTimeStep;  // unit conversion
	double dt2 = dt*dt;
	double dt3 = dt2*dt;
	double dt4 = dt3*dt;
	double dt5 = dt4*dt;
	double da_dt2, adot_dt, a2, a3;

	//fprintf(stdout, "PID=%d\n", ptcl->PID);
	for (int dim=0; dim<Dim; dim++) {
		da_dt2  = (this->a_reg[dim][0] - new_a[dim]   - a_tmp[dim]    ) / dt2;
		adot_dt = (this->a_reg[dim][1] + new_adot[dim] + adot_tmp[dim]) / dt;


		a2 =  -6*da_dt2 - 2*adot_dt - 2*this->a_reg[dim][1]/dt;
		a3 = (12*da_dt2 + 6*adot_dt)/dt;

		// note that these higher order terms and lowers have different neighbors

		//fprintf(stdout, "da_dt2 =%.3e, adot_dt =%.3e, dt=%.3e\n", da_dt2, adot_dt, dt);
		//fprintf(stdout, "a2     =%.3e, a3      =%.3e\n", a2, a3);
		/*
			 if (ptcl->PID == 753) {
			 fprintf(stderr, "dim=%d, a2=%.3e, a3=%.3e/a0=%.3e, atot=%.3e, a_tmp=%.3e, adot_tmp=%.3e, dt=%.3e\n", 
			 dim, a2,a3,ptcl->a_reg[dim][0],AccRegReceive[i][dim],a_tmp[dim],adot_tmp[dim],dt*1e10/1e6);
			 fprintf(stderr, "dim=%d, da_dt2=%.3e, adot_dt=%.3e\n", 
			 dim, da_dt2, adot_dt);
			 }
			 */

		// 4th order correction
		// save the values in the temporary variables
		this->NewPosition[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->NewVelocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;

		this->a_reg[dim][2] = a2;
		this->a_reg[dim][3] = a3;

		// reset for future use
		a_tmp[dim]    = 0.;
		adot_tmp[dim] = 0.;
	}

	int _NewNumberOfNeighbor = 0;
	for (int i=0; i<NewNumberOfNeighborGPU; i++) {
		if (particles[NewNeighborsGPU[i]].isCMptcl) {
			for (int j=0; j<particles[NewNeighborsGPU[i]].NumberOfMember; j++) {
				this->NewNeighbors[_NewNumberOfNeighbor++] = particles[NewNeighborsGPU[i]].Members[j];
			}
		}
		else
			this->NewNeighbors[_NewNumberOfNeighbor++] = NewNeighborsGPU[i];
	}
	this->NewNumberOfNeighbor = _NewNumberOfNeighbor;


	for (int dim=0; dim<Dim; dim++) {
		this->a_reg[dim][0] = new_a[dim];
		this->a_reg[dim][1] = new_adot[dim];
		this->a_tot[dim][0] = this->a_reg[dim][0] + this->a_irr[dim][0];
		this->a_tot[dim][1] = this->a_reg[dim][1] + this->a_irr[dim][1];
		if (this->NewNumberOfNeighbor == 0) {
			this->a_tot[dim][2] = this->a_reg[dim][2];
			this->a_tot[dim][3] = this->a_reg[dim][3];
		}
	}
}

