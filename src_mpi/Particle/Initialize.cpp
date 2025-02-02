#include <iostream>
#include "../particle.h"
#include "../global.h"
#include <unordered_set>

void CalculateAcceleration01(Particle* ptcl);
void CalculateAcceleration23(Particle* ptcl);



/*******************************************************
 *  Initialize Accelerations
 *******************************************************/



void CalculateAcceleration01(Particle* ptcl1) {

	double x[Dim], v[Dim];
	double m_r3;
	double v2;
	double r2 = 0;
	double vx = 0;

	std::unordered_set<int> CMPtclsSet;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]    = 0.;
		v[dim]    = 0.;
	}

	ptcl1->NumberOfNeighbor = 0;
	for(int dim=0; dim<Dim; dim++) {
		for (int order=0; order<HERMITE_ORDER; order++) {
			ptcl1->a_reg[dim][order] = 0.0;
			ptcl1->a_irr[dim][order] = 0.0;
			ptcl1->a_tot[dim][order] = 0.0;
		}
	}

	//std::cout << "nbody+: Entering CalculateInitialAcceleration  ..." << std::endl;
	//ptcl1->predictParticleSecondOrder(newTime);
	
	//fprintf(stdout, "pid=%d, nn=%d, numpart=%d\n", ptcl1->PID, ptcl1->NumberOfNeighbor, LastParticleIndex);
	Particle *ptcl2;
	for (int i=0; i<global_variable->LastParticleIndex; i++) {
		ptcl2 = &particles[i];

		if (ptcl1->PID == ptcl2->PID) {
			continue;
		}

		if (!ptcl2->isActive) {
			if (ptcl2->CMPtclIndex != -1) {
				CMPtclsSet.insert(ptcl2->CMPtclIndex);
				continue;
			}
			else {
				if (ptcl2->Mass == 0) continue;
				else {
					fprintf(stderr, "Why not zero mass? this PID: %d, neighbor PID: %d\n", ptcl1->PID, ptcl2->PID);
					assert(ptcl2->Mass == 0);
				}
			}
		}

		r2 = 0;
		vx = 0;
		v2 = 0;

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		//ptcl2->predictParticleSecondOrder(newTime);
		for (int dim=0; dim<Dim; dim++) {
			x[dim] = ptcl2->Position[dim] - ptcl1->Position[dim];
			v[dim] = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
			r2    += x[dim]*x[dim];
			vx    += v[dim]*x[dim];
			v2    += v[dim]*v[dim];
		}

		m_r3 = ptcl2->Mass/r2/sqrt(r2); 

		if (r2 > ptcl1->RadiusOfNeighbor) {
			for (int dim=0; dim<Dim; dim++) {
				// Calculate 0th and 1st derivatives of acceleration
				ptcl1->a_reg[dim][0] += m_r3*x[dim];
				ptcl1->a_reg[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				ptcl1->a_irr[dim][0] += m_r3*x[dim];
				ptcl1->a_irr[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
			ptcl1->Neighbors[ptcl1->NumberOfNeighbor]=ptcl2->ParticleIndex; // Eunwoo: PID -> ParticleIndex
			ptcl1->NumberOfNeighbor++;
			//fprintf(stdout, "pid=%d, nn=%d\n", ptcl1->PID, ptcl1->NumberOfNeighbor);
		} // endfor dim
	} // endfor ptcl2

	if (!CMPtclsSet.empty()) {
		for (int i: CMPtclsSet) {
			ptcl2 = &particles[i];
			if (!ptcl2->isActive) {
				fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", ptcl1->PID, ptcl2->PID);
				assert(ptcl2->isActive);
			}

			r2 = 0;
			vx = 0;
			v2 = 0;

			// updated the predicted positions and velocities just in case
			// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
			//ptcl2->predictParticleSecondOrder(newTime);
			for (int dim=0; dim<Dim; dim++) {
				x[dim] = ptcl2->Position[dim] - ptcl1->Position[dim];
				v[dim] = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
				r2    += x[dim]*x[dim];
				vx    += v[dim]*x[dim];
				v2    += v[dim]*v[dim];
			}

			m_r3 = ptcl2->Mass/r2/sqrt(r2); 

			if (r2 > ptcl1->RadiusOfNeighbor) {
				for (int dim=0; dim<Dim; dim++) {
					// Calculate 0th and 1st derivatives of acceleration
					ptcl1->a_reg[dim][0] += m_r3*x[dim];
					ptcl1->a_reg[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
				}
			}
			else {
				for (int dim=0; dim<Dim; dim++) {
					ptcl1->a_irr[dim][0] += m_r3*x[dim];
					ptcl1->a_irr[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
				}
				for (int j; j < ptcl2->NumberOfMember; j++) {
					ptcl1->Neighbors[ptcl1->NumberOfNeighbor]=ptcl2->Members[j]; // Eunwoo: PID -> ParticleIndex
					ptcl1->NumberOfNeighbor++;
				}
				//fprintf(stdout, "pid=%d, nn=%d\n", ptcl1->PID, ptcl1->NumberOfNeighbor);
			} // endfor dim
		} // endfor ptcl2
	}
		//
	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<2; order++) {
			ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order]; 
		}
	}
	return;
}


	/*
		// Calculate 2nd and 3rd derivatives of acceleration
		if (restart) {
			;
		}
	*/

void CalculateAcceleration23(Particle* ptcl1) {

	double x[Dim], v[Dim], a21[Dim], a21dot[Dim], a1[Dim], a2[Dim], a1dot[Dim], a2dot[Dim];
	double a, b, c;
	double rdf_r2, vdf_r2, rdfdot_r2, v2, r2, r3, vr, m_r3;
	double adot2, adot3;

	std::unordered_set<int> CMPtclsSet;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = ptcl1->a_tot[dim][0];
		a1dot[dim]  = ptcl1->a_tot[dim][1];
	}

	Particle *ptcl2;
	for (int i=0; i<global_variable->LastParticleIndex; i++) {
		ptcl2 = &particles[i];

		if (ptcl1->PID == ptcl2->PID) {
			continue;
		}

		if (!ptcl2->isActive) {
			if (ptcl2->CMPtclIndex != -1) {
				CMPtclsSet.insert(ptcl2->CMPtclIndex);
				continue;
			}
			else {
				if (ptcl2->Mass == 0) continue;
				else {
					fprintf(stderr, "Why not zero mass? this PID: %d, neighbor PID: %d\n", ptcl1->PID, ptcl2->PID);
					assert(ptcl2->Mass == 0);
				}
			}
		}
		
		r2 = 0;
		r3 = 0;
		v2 = 0;
		vr = 0;
		rdf_r2 = 0;
		vdf_r2 = 0;
		rdfdot_r2 = 0;

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		for (int dim=0; dim<Dim; dim++) {
			a2[dim]    = ptcl2->a_tot[dim][0];
			a2dot[dim] = ptcl2->a_tot[dim][1];
			x[dim]     = ptcl2->Position[dim] - ptcl1->Position[dim];
			v[dim]     = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
			r2        += x[dim]*x[dim];
			vr        += v[dim]*x[dim];
			v2        += v[dim]*v[dim];
		}

		r3   = r2*sqrt(r2);
		m_r3 = ptcl2->Mass/r3; 

		for (int dim=0; dim<Dim; dim++) {
			a21[dim]    = m_r3*x[dim];
			a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vr/r2);
			rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
			vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
			rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
		}

		a = vr/r2;
		b = v2/r2 + rdf_r2 + a*a;
		c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);


		if (r2 > ptcl1->RadiusOfNeighbor) {
			for (int dim=0; dim<Dim; dim++) {
				adot2 = -ptcl2->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
				adot3 = -ptcl2->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
				ptcl1->a_reg[dim][2] += adot2;
				ptcl1->a_reg[dim][3] += adot3;
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				adot2 = -ptcl2->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
				adot3 = -ptcl2->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
				ptcl1->a_irr[dim][2] += adot2;
				ptcl1->a_irr[dim][3] += adot3;
			}
		} // endfor if
	} //endfor ptcl2

	if (!CMPtclsSet.empty()) {
		for (int i: CMPtclsSet) {
			ptcl2 = &particles[i];
			if (!ptcl2->isActive) {
				fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", ptcl1->PID, ptcl2->PID);
				assert(ptcl2->isActive);
			}

			r2 = 0;
			r3 = 0;
			v2 = 0;
			vr = 0;
			rdf_r2 = 0;
			vdf_r2 = 0;
			rdfdot_r2 = 0;

			// updated the predicted positions and velocities just in case
			// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
			for (int dim=0; dim<Dim; dim++) {
				a2[dim]    = ptcl2->a_tot[dim][0];
				a2dot[dim] = ptcl2->a_tot[dim][1];
				x[dim]     = ptcl2->Position[dim] - ptcl1->Position[dim];
				v[dim]     = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
				r2        += x[dim]*x[dim];
				vr        += v[dim]*x[dim];
				v2        += v[dim]*v[dim];
			}

			r3   = r2*sqrt(r2);
			m_r3 = ptcl2->Mass/r3; 

			for (int dim=0; dim<Dim; dim++) {
				a21[dim]    = m_r3*x[dim];
				a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vr/r2);
				rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
				vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
				rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
			}

			a = vr/r2;
			b = v2/r2 + rdf_r2 + a*a;
			c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);


			if (r2 > ptcl1->RadiusOfNeighbor) {
				for (int dim=0; dim<Dim; dim++) {
					adot2 = -ptcl2->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
					adot3 = -ptcl2->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
					ptcl1->a_reg[dim][2] += adot2;
					ptcl1->a_reg[dim][3] += adot3;
				}
			}
			else {
				for (int dim=0; dim<Dim; dim++) {
					adot2 = -ptcl2->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
					adot3 = -ptcl2->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
					ptcl1->a_irr[dim][2] += adot2;
					ptcl1->a_irr[dim][3] += adot3;
				}
			} // endfor if
		} // endfor ptcl2
	}

	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=2; order<4; order++) {
			ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order];
		}
	}
}




/*******************************************************
 *  Initialize Time Steps 
 *******************************************************/


double getNewTimeStep(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int& TimeLevel, ULL &TimeBlock, double &TimeStep);


void Particle::initializeTimeStep() {

	//std::cout << "Initializing timesteps ..." << std::endl;

	double dtIrr, dtReg;

	dtReg = getNewTimeStep(this->a_reg, this->a_reg);
	//std::cout << "dtReg=" << dtReg << std::endl;
	getBlockTimeStep(dtReg, this->TimeLevelReg, this->TimeBlockReg, this->TimeStepReg);
	//std::cout << "TimeStepReg=" << this->TimeStepReg*EnzoTimeStep*1e10/1e6 << std::endl;


	this->TimeStepReg  = MIN(1,this->TimeStepReg);
	this->TimeBlockReg = std::min(block_max, this->TimeBlockReg);
	this->TimeLevelReg = std::min(0, this->TimeLevelReg);

	if (this->NumberOfNeighbor != 0) {
		dtIrr = getNewTimeStep(this->a_tot, this->a_irr);
		getBlockTimeStep(dtIrr, this->TimeLevelIrr, this->TimeBlockIrr, this->TimeStepIrr);
	}
	else {
		this->TimeBlockIrr = this->TimeBlockReg;
		this->TimeLevelIrr = this->TimeLevelReg;
		this->TimeStepIrr  = this->TimeStepReg;
	}

	this->CurrentTimeIrr  = 0;
	this->CurrentTimeReg  = 0;
	this->CurrentBlockIrr = 0;
	this->CurrentBlockReg = 0;

}


/* this is for few body initialization
 *
int InitializeTimeStep(Particle* particle, int size) {
	std::cout << "Initializing timesteps ..." << std::endl;
	double dtIrr, dtReg;
	Particle *ptcl;

	for (int i=0; i<size; i++){
		ptcl = &particle[i];
		dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
		getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeBlockReg, ptcl->TimeStepReg);

		if (ptcl->NumberOfAC != 0) {
			dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
			getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeBlockIrr, ptcl->TimeStepIrr);
		}
		else {
			ptcl->TimeBlockIrr = ptcl->TimeBlockReg;
			ptcl->TimeLevelIrr = ptcl->TimeLevelReg;
			ptcl->TimeStepIrr  = ptcl->TimeStepReg;
		}

		ptcl->TimeStepReg  = MIN(1,ptcl->TimeStepReg);
		ptcl->TimeBlockReg = std::min(block_max, ptcl->TimeBlockReg);
		ptcl->TimeLevelReg = std::min(0, ptcl->TimeLevelReg);

		ptcl->CurrentTimeIrr  = 0;
		ptcl->CurrentTimeReg  = 0;
		ptcl->CurrentBlockIrr = 0;
		ptcl->CurrentBlockReg = 0;
	} // endfor size

	for (int i=0; i<size; i++){
		ptcl = &particle[i];

		if (ptcl->NumberOfAC != 0) {
			while (ptcl->TimeLevelIrr >= ptcl->TimeLevelReg) {
				ptcl->TimeLevelIrr--;
			}
		}
		ptcl->TimeBlockIrr = static_cast<ULL>(pow(2, ptcl->TimeLevelIrr-time_block));
		ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg-time_block));
	} //endfor size
	return true;
}



*/
