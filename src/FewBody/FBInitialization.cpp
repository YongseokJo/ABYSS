#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"

void direct_sum(REAL *x, REAL *v, REAL r2, REAL vx,
		REAL mass, REAL (&a)[3], REAL (&adot)[3]);

REAL getNewTimeStepIrr(REAL f[3][4], REAL df[3][4]);
void getBlockTimeStep(REAL dt, int& TimeLevel, ULL &TimeBlock, REAL &TimeStep);
bool UpdateComputationChain(Particle* ptcl);
void UpdateNextRegTime(std::vector<Particle*> &particle);
bool CreateComputationList(Particle* ptcl);
bool CreateComputationChain(std::vector<Particle*> &particle);


//       /////////////////////        //
//       /////////////////////        //
//       CALCULATEACCELERATION        //
//       /////////////////////        //
//       /////////////////////        //


void CalculateFBAcceleration(Particle* ptclI, Particle* ptclCM, std::vector<Particle*> &particle) {

	int j=0;
	Particle* ptcl1;
	size_t length = ptclI->GroupParticles.size();
	REAL x[Dim], v[Dim];
	REAL a21[Dim], a21dot[Dim], a1[Dim], a2[Dim], a1dot[Dim], a2dot[Dim];
	REAL rdf_r2, vdf_r2, rdfdot_r2, v2, r2, r3, vr, m_r3;
	REAL adot2, adot3;
	REAL a,b,c;
	REAL current_time = ptclCM->CurrentTimeIrr;

	// Initialize relevant variables

	for (int dim=0; dim<Dim; dim++) {
		x[dim]    = 0.;
		v[dim]    = 0.;
	}

	//std::sort(particle.begin(),particle.end(),
				//[](Particle* p1, Particle* p2) { return p1->PID > p2->PID; });

	for (int i=0; i<length+1; i++) {

		if (i==0) 
			ptcl1 = ptclI;
		else		
			ptcl1 = ptclI->GroupParticles[i-1];

		// updated the predicted positions and velocities just in case
		ptcl1->predictParticleSecondOrderIrr(current_time);
		//std::sort(ptcl1->ACList.begin(),ptcl1->ACList.end());
		std::sort(ptcl1->ACList.begin(),ptcl1->ACList.end(),
				[](Particle* p1, Particle* p2) { return p1->ParticleOrder < p2->ParticleOrder;});


		// initialization of relevant variables 
		j = 0;
		for(int dim=0; dim<Dim; dim++) {
			for (int order=0; order<HERMITE_ORDER; order++) {
				ptcl1->a_reg[dim][order] = 0.0;
				ptcl1->a_irr[dim][order] = 0.0;
				ptcl1->a_tot[dim][order] = 0.0;

				ptclCM->a_reg[dim][order] = 0.0;
				ptclCM->a_irr[dim][order] = 0.0;
				ptclCM->a_tot[dim][order] = 0.0;
			}
		}

		for (Particle *ptcl2: particle) {
			r2 = 0;
			vr = 0;
			v2 = 0;

			if (ptcl1->PID == ptcl2->PID) {
				continue;
			}

			ptcl2->predictParticleSecondOrderIrr(current_time);
			for (int dim=0; dim<Dim; dim++) {
				x[dim] = ptcl2->PredPosition[dim] - ptcl1->PredPosition[dim];
				v[dim] = ptcl2->PredVelocity[dim] - ptcl1->PredVelocity[dim];
				r2    += x[dim]*x[dim];
				vr    += v[dim]*x[dim];
				v2    += v[dim]*v[dim];
			}

			m_r3 = ptcl2->Mass/r2/sqrt(r2); 

			if ((ptcl1->NumberOfAC==0) || j >= ptcl1->NumberOfAC || (ptcl2 != ptcl1->ACList[j])) {
				for (int dim=0; dim<Dim; dim++) {
					// Calculate 0th and 1st derivatives of acceleration
					ptcl1->a_reg[dim][0] += m_r3*x[dim];
					ptcl1->a_reg[dim][1] += m_r3*(v[dim] - 3*x[dim]*vr/r2);
				}
			}
			else {
				for (int dim=0; dim<Dim; dim++) {
					ptcl1->a_irr[dim][0] += m_r3*x[dim];
					ptcl1->a_irr[dim][1] += m_r3*(v[dim] - 3*x[dim]*vr/r2);
				}
				j++;
			} // endfor dim
		} // end of loop for ptcl2 (full particle)


		// update total acceleration as well
		for (int dim=0; dim<Dim; dim++)	 {
			for (int order=0; order<HERMITE_ORDER; order++) {
				ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order]; 
			}
		}
	}// end of loop for pair particle, ptclI and ptclJ


	// copy the calculated values to CM particle
	// and initialize the 3rd and 4th order derivatives just in case

	for (Particle* ptclJ : ptclI->GroupParticles) {
		for (int dim = 0; dim < Dim; dim++) {
			for (int order = 0; order < 2; order++) {
				ptclCM->a_reg[dim][order] += ptclJ->a_reg[dim][order]*ptclJ->Mass;
				ptclCM->a_irr[dim][order] += ptclJ->a_irr[dim][order]*ptclJ->Mass;
				ptclCM->a_tot[dim][order] += ptclJ->a_tot[dim][order]*ptclJ->Mass;
			}
		}
	}

	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<2; order++) { // Eunwoo: Why order<2?? -> order<HERMITE_ORDER
			ptclCM->a_reg[dim][order] += ptclI->a_reg[dim][order]*ptclI->Mass;
			ptclCM->a_reg[dim][order] /= ptclCM->Mass; 
			ptclCM->a_irr[dim][order] += ptclI->a_irr[dim][order]*ptclI->Mass;
			ptclCM->a_irr[dim][order] /= ptclCM->Mass; 
			ptclCM->a_tot[dim][order] += ptclI->a_tot[dim][order]*ptclI->Mass;
			ptclCM->a_tot[dim][order] /= ptclCM->Mass; 
		}
	}


	// updated the predicted positions and velocities just in case

	ptcl1 = ptclCM;
	ptcl1->predictParticleSecondOrderIrr(current_time);
	std::cout << "CM particle neighbor = " << ptcl1->NumberOfAC << std::endl;
	j = 0;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = ptcl1->a_tot[dim][0];
		a1dot[dim]  = ptcl1->a_tot[dim][1];
	}


	for (Particle *ptcl2: particle) {

		r2 = 0;
		r3 = 0;
		v2 = 0;
		vr = 0;
		rdf_r2 = 0;
		vdf_r2 = 0;
		rdfdot_r2 = 0;

		if ((ptcl2->PID == ptclI->PID) || 
			(std::find_if(ptclI->GroupParticles.begin(), ptclI->GroupParticles.end(), 
						[ptcl2](Particle* ptcl) { return ptcl->PID == ptcl2->PID; }) != ptclI->GroupParticles.end())) {
			continue;
		} // Eunwoo: fixed this parts

		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		ptcl2->predictParticleSecondOrderIrr(current_time);

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


		if ((ptcl1->NumberOfAC==0) || j >= ptcl1->NumberOfAC ||(ptcl2 != ptcl1->ACList[j])) {
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
			j++;
		} // endfor if
	} // end of loop for ptcl2 (full particle)

	for (int dim=0; dim<Dim; dim++)	 {
		for (int order=2; order<HERMITE_ORDER; order++) { // Eunwoo: Why order = 2 ??
			ptclCM->a_tot[dim][order] = ptclCM->a_reg[dim][order] + ptclCM->a_irr[dim][order]; 
		}
	}

	// save the regular and irregular time steps as well

}






//        ///////////////////        //
//        ///////////////////        //
//           ISFBCANDIDATE           //
//        ///////////////////        //
//        ///////////////////        //



// made 2024.08.08 by Eunwoo Chung

// from particles with shortest time steps...
// Group detection - by distance for the simplest case
// reference - not yet but PeTar/search_group_candidate.hpp might be refered laler

void Particle::isFBCandidate() {

	// temporary calculation variables
	REAL x[Dim];
	REAL r2;
	REAL current_time;
	REAL r_group = 1e-4/position_unit; // group detecting radius: 1e-4 pc (tentative) for FewBody physics
	// Eunwoo: idea - search bound/unbound particles within neighbor as BIFROST do
	std::vector<Particle*> groupParticles; // Vector to store pointers of particles in the group not including mother

	// predict the particle position to obtain more information
	// particle regularlized if the conditions are satisfied at a future time

	// need to consider case when c.m particles are the cause of the small steps
	// need to be added later - CHECK

	// fprintf(binout, "Function isFBCandidate started\n"); // Eunwoo debug

	for (Particle* ptcl: ACList) {

		// if particle time step is too large, skip
		// if the neighbor step is larger then 8 times of the candidate particle, then skip
		r2 = 0.0;

		if ((ptcl->TimeLevelIrr) > (this->TimeLevelIrr+3) || ptcl->isGroup || ptcl->isCMptcl) // Eunwoo didn't fully understand this
			continue;


		// find out what the paired particle is
		current_time = this->CurrentTimeIrr > ptcl->CurrentTimeIrr ? \
									 this->CurrentTimeIrr : ptcl->CurrentTimeIrr;
		this->predictParticleSecondOrderIrr(current_time);
		ptcl->predictParticleSecondOrderIrr(current_time);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences
			x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
		}

		// find out particles which can be detected as a group

		// if (r2<r_group*r_group) {
		if ((r2<r_group*r_group) && (this->PID != ptcl->PID)) { // Eunwoo for test
			groupParticles.push_back(ptcl);
			fprintf(binout, "group member added!\n"); // Eunwoo debug
			fprintf(stdout, "group member added!\n"); // Eunwoo check
		}
	}

	if (!groupParticles.empty()) {
		this->isGroup = true; // Assuming isGroup is a member of Particle to mark group status
		for (Particle* ptcl : groupParticles) {
			ptcl->isGroup = true;
		}
	}

	this->GroupParticles = groupParticles; // 'this' is the mother of 'groupParticles'

}




	//        ///////////////////        //
	//        ///////////////////        //
	//           STARTNEWKSREG           //
	//        ///////////////////        //
	//        ///////////////////        //



	// start new KS regularlization



	//        ///////////////////        //
	//        ///////////////////        //
	//        NEWKSINITIALIZATION        //
	//        ///////////////////        //
	//        ///////////////////        //



	// initialize conditions for new Few body group

void NewFBInitialization(Particle* ptclI, std::vector<Particle*> &particle, std::vector<Particle*> &ComputationList) {

	// basic variables for calculation
	Particle *ptclCM;
	Group *ptclGroup;
	// std::vector<Particle*> KSNeighborList;

	// int ptclIIndex;
	int numMembers = 0;

	REAL dtReg, dtIrr;

	std::cout <<"Starting Routine NewFBInitialization" << std::endl;

	// BASIC DEFINITIONS


	// define the pair particle for particle I
	// ptclJ = ptclI->BinaryPairParticle;
	// ptclJ->BinaryPairParticle = ptclI;

	//std::cout << "Predicting Particle Positions" << std::endl;
	fprintf(binout, "\n-------------------------NEW-GROUP------------------------\n");
	fprintf(stdout, "\n-------------------------NEW-GROUP------------------------\n"); // Eunwoo check

	// fprintf(binout, "Radius = %e, \n", dist(ptclI->Position, ptclJ->Position));
	fprintf(binout, "I. Total Acceleration - ax:%e, ay:%e, az:%e \n", ptclI->a_tot[0][0], ptclI->a_tot[1][0], ptclI->a_tot[2][0]);
	fprintf(binout, "I. Time Steps - irregular:%e, regular:%e \n", ptclI->TimeStepIrr*EnzoTimeStep*1e4, ptclI->TimeStepReg*EnzoTimeStep*1e4);
	for (Particle* ptclJ : ptclI->GroupParticles) {
		numMembers += 1;
        ptclJ->isGroup = true; // Mark particle as part of a group // Eunwoo: redundant - this is already done on isFBCandidate
        fprintf(binout, "I_%d. Distance from I: %e\n", numMembers, dist(ptclI->Position, ptclJ->Position));
        fprintf(binout, "I_%d. Total Acceleration - ax:%e, ay:%e, az:%e \n", numMembers, ptclJ->a_tot[0][0], ptclJ->a_tot[1][0], ptclJ->a_tot[2][0]);
        fprintf(binout, "I_%d. Time Steps - irregular:%e, regular:%e\n", numMembers, ptclJ->TimeStepIrr*EnzoTimeStep*1e4, ptclJ->TimeStepReg*EnzoTimeStep*1e4);
    }

	//fprintf(binout, "\nPosition: ptclI - x:%e, y:%e, z:%e, \n", ptclI->Position[0], ptclI->Position[1], ptclI->Position[2]);
	//fprintf(binout, "Velocity: ptclI - vx:%e, vy:%e, vz:%e, \n", ptclI->Velocity[0], ptclI->Velocity[1], ptclI->Velocity[2]);

	/*
	fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclI->a_tot[0][1], ptclI->a_tot[1][1], ptclI->a_tot[2][1]);
	fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclI->a_tot[0][2], ptclI->a_tot[1][2], ptclI->a_tot[2][2]);
	fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclI->a_tot[0][3], ptclI->a_tot[1][3], ptclI->a_tot[2][3]);
	*/

	//fprintf(binout, "\nPosition: ptclJ - x:%e, y:%e, z:%e, \n", ptclJ->Position[0], ptclJ->Position[1], ptclJ->Position[2]);
	//fprintf(binout, "Velocity: ptclJ - vx:%e, vy:%e, vz:%e, \n", ptclJ->Velocity[0], ptclJ->Velocity[1], ptclJ->Velocity[2]);


	// need to put option if there aren't any close neighbors
	// define the new center of mass particle
	//std::cout << "Make new particle and binary information" << std::endl;
	ptclCM  = new Particle();
	ptclGroup = new Group();

	// fprintf(stderr, "time error max: %e", ptclGroup->manager.time_error_max); // Eunwoo debug


	Particle* ptcl = ptclI;
	for (Particle* ptclJ : ptclI->GroupParticles) {
    	if (ptclJ->CurrentTimeIrr > ptcl->CurrentTimeIrr) {
        	ptcl = ptclJ;
    	}
	}

	ptclCM->CurrentTimeIrr  = ptcl->CurrentTimeIrr;
	ptclCM->CurrentTimeReg  = ptcl->CurrentTimeReg;
	ptclCM->CurrentBlockIrr = ptcl->CurrentBlockIrr; 
	ptclCM->CurrentBlockReg = ptcl->CurrentBlockReg;

	ptclCM->TimeStepIrr     = ptcl->TimeStepIrr;
	ptclCM->TimeBlockIrr    = ptcl->TimeBlockIrr;
	ptclCM->TimeLevelIrr    = ptcl->TimeLevelIrr;

	ptclCM->TimeStepReg     = ptcl->TimeStepReg;
	ptclCM->TimeBlockReg    = ptcl->TimeBlockReg;
	ptclCM->TimeLevelReg    = ptcl->TimeLevelReg;

	ptclCM->PID             = -(ptcl->PID + NNB);
	ptclCM->GroupMother		= ptclI;
	ptclCM->GroupParticles	= ptclI->GroupParticles;
	ptclCM->GroupInfo		= ptclGroup;
	ptclCM->isCMptcl        = true;

	fprintf(binout, "The ID of CM is %d.\n",ptclCM->PID);
	fprintf(stdout, "The ID of CM is %d.\n",ptclCM->PID); // Eunwoo check
	fprintf(binout, "The ID of members are");
	fprintf(stdout, "The ID of members are"); // Eunwoo check
	for (Particle* ptclJ : ptclI->GroupParticles) {
		fprintf(binout, " %d, ", ptclJ->PID);
		fprintf(stdout, " %d, ", ptclJ->PID); // Eunwoo check
	}
	fprintf(binout, "The ID of mother is %d.\n", ptclI->PID);
	fprintf(stdout, "The ID of mother is %d.\n", ptclI->PID); // Eunwoo check
	//fflush(binout);

	ptclCM->Mass = ptclI->Mass; // Start with the mass of the mother particle
	for (Particle* ptclJ : ptclI->GroupParticles) {
		ptclCM->Mass += ptclJ->Mass; // Add the mass of each group member
	}

	ptclI->predictParticleSecondOrderIrr(ptclCM->CurrentTimeIrr);
	for (Particle* ptclJ : ptclI->GroupParticles){
		ptclJ->predictParticleSecondOrderIrr(ptclCM->CurrentTimeIrr);
	}

	for (Particle* ptclJ : ptclI->GroupParticles){
		for (int dim=0; dim<Dim; dim++){
			ptclCM->Position[dim]	+= ptclJ->PredPosition[dim]*ptclJ->Mass;
			ptclCM->Velocity[dim]	+= ptclJ->PredVelocity[dim]*ptclJ->Mass;
		}
	}

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->Position[dim]	= (ptclCM->Position[dim] + ptclI->PredPosition[dim]*ptclI->Mass)/ptclCM->Mass;
		ptclCM->Velocity[dim]	= (ptclCM->Velocity[dim] + ptclI->PredVelocity[dim]*ptclI->Mass)/ptclCM->Mass;
		ptclCM->PredPosition[dim] = ptclCM->Position[dim];
		ptclCM->PredVelocity[dim] = ptclCM->Velocity[dim];
	}

	ptclGroup->groupCM       = ptclCM;
	ptclGroup->CurrentTime  = ptcl->CurrentTimeIrr;
	// ptclBin->CurrentTau = 0.0; // already initialized at formation, but added for safety.
	// ptclBin->PredTau = 0.0;



	// copy the neighbor list for c.m particle
	ptclCM->NumberOfAC = 0;
	ptclCM->RadiusOfAC = ptcl->RadiusOfAC;

	for (Particle* neighbor: ptcl->ACList) { // Eunwoo: fixed this part
		if (neighbor->PID == ptclI->PID || 
		std::find_if(ptclI->GroupParticles.begin(), ptclI->GroupParticles.end(), 
                 [neighbor](Particle* ptcl) { return ptcl->PID == neighbor->PID; }) != ptclI->GroupParticles.end()) {
			continue;
		}
		ptclCM->ACList.push_back(neighbor);
		ptclCM->NumberOfAC++;
	}

	// Eunwoo check
	std::cout << ptclCM->PID << "("<< ptclCM->NumberOfAC <<")" << "=" <<std::flush;
	for (Particle * nn:ptclCM->ACList) {
		std::cout << nn->PID << ", ";
	}
	std::cout << std::endl; 
	// Eunwoo check


	// calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately for the binary pair particle and the cm particle
	//std::cout << "CalculateKSAcceleration starts" << std::endl;
	CalculateFBAcceleration(ptclI,ptclCM,particle);

	std::cout << "Calculating Time steps for the CM particle" << std::endl;
	ptclCM->calculateTimeStepReg();
	if (ptclCM->TimeLevelReg <= ptcl->TimeLevelReg-1 
			&& ptcl->TimeBlockReg/2+ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr+ptcl->TimeBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg-1;
	}
	else if  (ptclCM->TimeLevelReg >= ptcl->TimeLevelReg+1) {
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg+1;
	}
	else 
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg;


	ptclCM->TimeStepReg  = static_cast<REAL>(pow(2, ptclCM->TimeLevelReg));
	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));

	ptclCM->calculateTimeStepIrr(ptclCM->a_tot, ptclCM->a_irr);
	while (ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr <= global_time_irr 
			&& ptclCM->TimeLevelIrr <= ptcl->TimeLevelIrr) { //first condition guarantees that ptclcm is small than ptcl
		ptclCM->TimeLevelIrr++;
		ptclCM->TimeStepIrr  = static_cast<REAL>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}







	//update particle vector, nextregtime, computation chain and list
	ptclI->isErase = true;
	for (Particle* ptclJ : ptclI->GroupParticles) {
		ptclJ->isErase = true;
	}
	

	// update ParticleOrder
	int j = 0;
	for (Particle* ptcl : particle) {
		ptcl->ParticleOrder -= j;
		
		// Check if the particle is the mother particle or one of the group members
		if (ptcl == ptclI || std::find(ptclI->GroupParticles.begin(), ptclI->GroupParticles.end(), ptcl) != ptclI->GroupParticles.end()) {
			j++;
		}
	}



	/*
	fprintf(stderr, "particle:");
	for (Particle* ptcl:particle) {
		fprintf(stderr, "%d, ", ptcl->PID);
	}
	fprintf(stderr, "\n");
	*/

	// erase particle I and its group particles
	particle.erase(
			std::remove_if(particle.begin(), particle.end(),
				[](Particle* p) {
				bool to_remove = p->isErase;
				return to_remove;
				}),
			particle.end());
	ptclI->isErase = false;
	for (Particle* ptclJ : ptclI->GroupParticles) {
		ptclJ->isErase = false;
	}

	/*
	fprintf(stderr, "particle (%d, %d):", ptclI->PID, ptclJ->PID);
	for (Particle* ptcl:particle) {
		fprintf(stderr, "%d, ", ptcl->PID);
	}
	fprintf(stderr, "\n");
	*/


	// update particle order and add it to particle vector
	ptclCM->ParticleOrder=particle.size();
	particle.push_back(ptclCM); // Eunwoo: This part is not the problem
	// if (std::find(particle.begin(), particle.end(), ptclCM) == particle.end()) {
	// 	particle.push_back(ptclCM);
	// }
	UpdateNextRegTime(particle);

  CreateComputationChain(particle);
	CreateComputationList(FirstComputation);

	// Add the binary to binary integration list
	//std::cout << "Add the ptclBin to Binary List" << std::endl;
	GroupList.push_back(ptclGroup);


	fprintf(binout, "\nFBInitialization.cpp: result of CM particle value calculation\n");
	fprintf(binout, "from function NewFBInitialization\n");
	//fflush(binout);

	fprintf(binout, "The ID of ith particle is %d of %d",ptclCM->PID, ptclI->PID);
	fprintf(stdout, "The ID of ith particle is %d of %d",ptclCM->PID, ptclI->PID); // Eunwoo check
	for (Particle* ptclJ : ptclI->GroupParticles) {
		fprintf(binout, ", %d", ptclJ->PID);
		fprintf(stdout, ", %d", ptclJ->PID); // Eunwoo check
	}
	fprintf(binout, "\n");
	fprintf(stdout, "\n"); // Eunwoo check

	// fprintf(binout, "Position - x:%e, y:%e, z:%e, \n", ptclCM->Position[0], ptclCM->Position[1], ptclCM->Position[2]);
	// fprintf(binout, "Velocity - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0], ptclCM->Velocity[1], ptclCM->Velocity[2]);
	// //fflush(binout);

	// fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	// fprintf(binout, "Time Steps - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr, ptclCM->TimeStepReg);




	// calculate the initial values of relevant variables
	std::cout << "Initializing Group Information" << std::endl;
	// std::cout << "Initializing Binary Information" << std::endl;

	// Eunwoo debug

	ptclGroup->Members = ptclI->GroupParticles;
	ptclGroup->Members.push_back(ptclI);

	// Eunwoo debug

	// Eunwoo: Where is initialization of REAL PredTime, REAL TimeStep, int TimeLevel ?????

	// ptclGroup->InitializeGroup(ptclCM->CurrentTimeIrr);

	// fprintf(binout, "\nKS coordinates - u1:%e, u2:%e, u3:%e, u4:%e \n", ptclBin->u[0], ptclBin->u[1], ptclBin->u[2], ptclBin->u[3]);
	// fprintf(binout, "KS coordinates - udot1:%e, udot2:%e, udot3:%e, udot4:%e \n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
	// fprintf(binout, "KS coordinates - udot1:%e, udot2:%e, udot3:%e, udot4:%e \n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
	// fprintf(binout, "Other important KS variables - r:%e, h:%e, gamma: %e, tau:%e, step:%e \n", ptclBin->r, ptclBin->h, ptclBin->gamma, ptclBin->dTau, ptclBin->TimeStep);
	//fflush(binout);


	// we also need to change the neighbor list of Particles
	// assuming that all the neighbors are bidirectional
	// may need to update later if the radius for neighbor differs depending on the particle
	//
	// first for particle I
	std::cout << "Changing the Neighbor List of Particles" << std::endl;

	ptclI->isErase = true;
	for (Particle* ptclJ : ptclI->GroupParticles) {
		ptclJ->isErase = true;
	}

	int size = 0;
	for (Particle* ptcl: particle) {
		size = ptcl->ACList.size();

		ptcl->ACList.erase(
				std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					return to_remove; }),
				ptcl->ACList.end());
		ptcl->NumberOfAC = ptcl->ACList.size();

		if (size != ptcl->NumberOfAC) {
			ptcl->ACList.push_back(ptclCM);
			ptcl->NumberOfAC++;
		}
		//fflush(stderr);
	}

	// this makes binary particles have their CM particle as a neighbor
	ptclI->isErase = false;
	for (Particle* ptclJ : ptclI->GroupParticles) {
		ptclJ->isErase = false;
	}

	fprintf(binout, "\n---------------------END-OF-NEW-GROUP---------------------\n\n");
	fprintf(stdout, "\n---------------------END-OF-NEW-GROUP---------------------\n\n"); // Eunwoo check
	fflush(binout);
	fflush(stdout); // Eunwoo check
}


