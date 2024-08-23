#include <vector>
#include <iostream>
#include "../global.h"
#include <cmath>
#include "defs.h"
#include "cuda_functions.h"

/*
 *  Purporse: calculate acceleration and neighbors of all the particles on GPU.
 *  -> calculates based on their current positions
 *
 *  Date    : 2024.01.17  by Seoyoung Kim
 *
 */

void SendAllParticlesToGPU(std::vector <Particle*> &particle);

void CalculateAllAccelerationOnGPU(std::vector<Particle*> &particle){

	// variables for opening GPU
	int numGpuOpen = NNB+10;
	int numGpuSend = 1024;
	int numGpuCal;
	int mpi_rank = 0;

	// variables for saving variables to send to GPU
	CUDA_REAL* MassSend;
	CUDA_REAL* MdotSend;
	CUDA_REAL(*PositionSend)[Dim];
	CUDA_REAL(*VelocitySend)[Dim];

	CUDA_REAL* r2OfACSend;
	CUDA_REAL  TimeStepRegTmp, ri2;
	CUDA_REAL* TimeStepRegSend;

	CUDA_REAL(*AccSend)[Dim];
	CUDA_REAL(*AccDotSend)[Dim];
	//CUDA_REAL *PotSend;

	int **ACListReceive;
	int *NumNeighborReceive;

	// temporary variables for calculating the irregular force
	CUDA_REAL dx[Dim];
	CUDA_REAL dv[Dim];
	CUDA_REAL rij2,dr2i,dr3i,drdv;

	// extra variables
	bool neighborOK;
	int ACnumi2;
	int ACjid;
	CUDA_REAL dt = particle[0]->TimeStepReg*EnzoTimeStep;



	// first let's open the GPU
	//OpenDevice(&mpi_rank);


	// need to make array to send to GPU
	// allocate memory to the temporary variables
	MassSend        = new CUDA_REAL[NNB];
	MdotSend        = new CUDA_REAL[NNB];
	PositionSend    = new CUDA_REAL[NNB][Dim];
	VelocitySend    = new CUDA_REAL[NNB][Dim];

	r2OfACSend      = new CUDA_REAL[NNB];
	TimeStepRegSend = new CUDA_REAL[NNB];

	AccSend         = new CUDA_REAL[NNB][Dim];
	AccDotSend      = new CUDA_REAL[NNB][Dim];

	NumNeighborReceive = new int[NNB];
	ACListReceive      = new int*[NNB];
	for (int i=0; i<NNB; i++) {
		ACListReceive[i] = new int[NumNeighborMax];
	}


	// copy the data of particles to the arrays to be sent
	for (int i=0; i<NNB; i++) {
		ri2           = 0;
		MassSend[i]   = particle[i]->Mass;
		r2OfACSend[i] = (particle[i]->RadiusOfAC)*(particle[i]->RadiusOfAC);
		MdotSend[i]   = particle[i]->Mass;

		for (int dim=0; dim<Dim; dim++) {
			PositionSend[i][dim] = particle[i]->Position[dim];
			VelocitySend[i][dim] = particle[i]->Velocity[dim];
			ri2                 += particle[i]->Position[dim]*particle[i]->Position[dim];
		}

		TimeStepRegTmp     = 1.0/8.0*sqrt(1.0 + ri2);
		TimeStepRegSend[i] = MIN(TimeStepRegTmp,1.0);
	}

	// send the arrays to GPU
	//SendToDevice(&NNB, MassSend, PositionSend, VelocitySend, MdotSend, &NumNeighborMax);


	// calculate the force by sending the particles to GPU in multiples of 1024
	for (int i=0; i<NNB; i+=numGpuSend) {
		numGpuCal = std::min(1024,(NNB-i));

		//CalculateAccelerationOnDevice(&NNB, PositionSend, VelocitySend,
					//AccSend, AccDotSend, MdotSend, r2OfACSend, NumNeighborReceive, ACListReceive, dt);

		// copy the values of regular forces and neighbors obtained in GPU to particles
		// and also calculate the irregular forces in the process
		for (int i2=i; i2<(i+1024); i2++){
			for (int dim=0; dim<Dim; dim++) {
				particle[i2]->a_reg[dim][0] += AccSend[i2][dim];
				particle[i2]->a_reg[dim][1] += AccDotSend[i2][dim];
			}
			//ACnumi2 = AClistGpu[i][0];

			// if there are neighbors, save them to neighbor list and calculate the irregular force
			if (ACnumi2 > 0){
				particle[i2]->NumberOfAC = ACnumi2;

				for (int j=1; j<(ACnumi2+1); j++) {
					//ACjid = AClistGpu[i][j];
					particle[i2]->ACList.push_back(particle[ACjid]);

					rij2 = 0.0;
					drdv = 0.0;

					for (int dim=0; dim<Dim; dim++) {
						dx[dim] = particle[ACjid]->Position - particle[i2]->Position;
						dv[dim] = particle[ACjid]->Velocity - particle[i2]->Velocity;
						rij2   += dx[dim]*dx[dim];
						drdv   += dx[dim]*dv[dim];
					}

					dr2i = 1.0/rij2;
					dr3i = particle[ACjid]->Mass*dr2i*sqrt(dr2i);

					for (int dim=0; dim<Dim; dim++) {
						particle[i2]->a_irr[dim][0] += dx[dim]*dr3i;
						particle[i2]->a_irr[dim][1] += (dv[dim]-dx[dim]*drdv)*dr3i;
					}
				} // endfor neighbors 
			} else {  // in case of no neighbors, just set the neighbor number to 0 just in case. 
				particle[i2]->NumberOfAC = 0;
			} // if statement on neighbor number ends


			// copy the values to other values as well
			for (int dim=0; dim<3; dim++){
				particle[i2]->a_tot[dim][0] = particle[i2]->a_reg[dim][0] + particle[i2]->a_irr[dim][0];
				particle[i2]->a_tot[dim][1] = particle[i2]->a_reg[dim][1] + particle[i2]->a_irr[dim][1];
			}
		} // saving the 1024 results from GPU to local class ends
	} // endfor total calculation 


	// free all temporary variables
	delete[] MassSend;
	delete[] PositionSend;
	delete[] VelocitySend;

	delete[] r2OfACSend;
	delete[] TimeStepRegSend;

	delete[] AccSend;
	delete[] AccDotSend;

	// close GPU
	//
	//CloseDevice();

} // calculate 0th, 1st derivative of force + neighbors on GPU ends



/*
 *  Purporse: calculate acceleration and neighbors of some particles on GPU.
 *  the particles subject for calculation of regular force are given by std::vector<Particle*> list
 *  predicted positions and velocities are used for calculation on GPU for now
 * 
 *  -> calculates based on their current positions
 *
 *  Date    : 2024.01.17  by Seoyoung Kim
 *
 */

void CalculateListAccelerationOnGPU(std::vector<int> &IndexList, std::vector<Particle*> &particle){

	// variables for opening GPU
	int ListSize   = IndexList.size();
	int numGpuSend = 1024;
	int mpi_rank   = 0;
	int numGpuCal;

	// variables for GPU
	CUDA_REAL *MassSend;
	CUDA_REAL *MdotSend;
	CUDA_REAL(*PositionSend)[Dim];
	CUDA_REAL(*VelocitySend)[Dim];

	CUDA_REAL *r2OfACSend;
	CUDA_REAL *TimeStepRegSend;

	CUDA_REAL(*AccSend)[Dim];
	CUDA_REAL(*AccDotSend)[Dim];
	CUDA_REAL *PotSend;
	int **AClistGpu;
	int massFlag;

	// temporary variables for calculating the irregular force

	CUDA_REAL dx[Dim], dv[Dim];
	CUDA_REAL rij2, dr2i, dr3i, drdv;

	// extra variables
	int ACnumi2, ACjid, i2reg;


	// need to make array to send to GPU
	// allocate memory to the temporary variables
	MassSend     = new CUDA_REAL[ListSize];
	MdotSend     = new CUDA_REAL[ListSize];
	PositionSend = new CUDA_REAL[ListSize][Dim];
	VelocitySend = new CUDA_REAL[ListSize][Dim];

	r2OfACSend   = new CUDA_REAL[ListSize];
	TimeStepRegSend  = new CUDA_REAL[ListSize];

	AccSend       = new CUDA_REAL[ListSize][Dim];
	AccDotSend    = new CUDA_REAL[ListSize][Dim];
	PotSend       = new CUDA_REAL[ListSize];

	AClistGpu = new int*[ListSize];
	for (int i=0; i<(NNB); i++) {
		AClistGpu[i] = new int[NumNeighborMax];
	}


	// copy the data of particles to the arrays to be sent
	// predicted positions and velocities are sent

	for (int i=0; i<ListSize; i++) {
		rij2          = 0;
		MassSend[i]   = particle[(IndexList[i])]->Mass;
		MdotSend[i]   = particle[(IndexList[i])]->Mass;
		r2OfACSend[i] = particle[(IndexList[i])]->RadiusOfAC * particle[(IndexList[i])]->RadiusOfAC;

		for (int dim=0; dim<Dim; dim++) {
			PositionSend[i][dim] = particle[(IndexList[i])]->PredPosition[dim];
			VelocitySend[i][dim] = particle[(IndexList[i])]->PredVelocity[dim];
		}
		TimeStepRegSend[i] = particle[(IndexList[i])]->TimeStepReg;
	}

	// send the arrays to GPU
	SendAllParticlesToGPU(particle);

	// calculate the force by sending the particles to GPU in multiples of 1024
	for (int i=0; i<ListSize; i+=numGpuSend) {
		numGpuCal = std::min(1024,(NNB-i));

		//CalculateAccelerationOnDevice(&NNB,&i, PositionSend, VelocitySend,
		//			AccSend, AccDotSend, MdotSend, r2OfACSend);

		// copy the values of regular forces and neighbors obtained in GPU to particles
		// and also calculate the irregular forces in the process
		for (int i2=i; i2<(i+1024); i2++){
			i2reg = (IndexList[i2]);

			// need to clear the AClist and the number of neighbors before copying the new values
			particle[i2reg]->NumberOfAC = 0;
			particle[i2reg]->ACList.clear();

			for (int dim=0; dim<Dim; dim++) {
				particle[i2reg]->a_reg[dim][0] += AccSend[i2][dim];
				particle[i2reg]->a_reg[dim][1] += AccDotSend[i2][dim];
			}

			ACnumi2 = AClistGpu[i][0];

			// if there are neighbors, save them to neighbor list and calculate the irregular force
			if (ACnumi2 > 0){
				particle[i2reg]->NumberOfAC = ACnumi2;

				for (int j=1; j<(ACnumi2+1); j++) {
					ACjid = AClistGpu[i][j];
					particle[i2reg]->ACList.push_back(particle[ACjid]);
					rij2 = 0.0;
					drdv = 0.0;

					for (int dim=0; dim<Dim; dim++) {
						dx[dim] = particle[ACjid]->PredPosition - particle[i2reg]->PredPosition;
						dv[dim] = particle[ACjid]->PredVelocity - particle[i2reg]->PredVelocity;
						rij2   += dx[dim]*dx[dim];
						drdv   += dx[dim]*dv[dim];
					}
					dr2i = 1.0/rij2;
					dr3i = particle[ACjid]->Mass*dr2i*sqrt(dr2i);

					for (int dim=0; dim<Dim; dim++) {
						particle[i2reg]->a_irr[dim][0] += dx[dim]*dr3i;
						particle[i2reg]->a_irr[dim][1] += (dv[dim]-dx[dim]*drdv)*dr3i;
					}
				} // loop on neighbors end
			} else { // in case of no neighbors, just set the neighbor number to 0 just in case. 
				particle[i2]->NumberOfAC = 0;
			} // if statement on neighbor number ends

			for (int dim=0; dim<3; dim++){ // copy the values to other values as well
				particle[i2reg]->a_tot[dim][0] = particle[i2reg]->a_reg[dim][0] + particle[i2reg]->a_irr[dim][0];
				particle[i2reg]->a_tot[dim][1] = particle[i2reg]->a_reg[dim][1] + particle[i2reg]->a_irr[dim][1];
			}
		} // saving the 1024 results from GPU to local class ends
	} // loop of total calculation ends


	// free all temporary variables
	delete[] MassSend;
	delete[] PositionSend;
	delete[] VelocitySend;

	delete[] r2OfACSend;
	delete[] TimeStepRegSend;

	delete[] AccSend;
	delete[] AccDotSend;
	delete[] PotSend;
	delete[] AClistGpu;

} // calculate 0th, 1st derivative of force + neighbors on GPU ends



/*
 *  Purporse: send the information of all particles to GPU in regular integration steps
 *  send the predicted positions and velocities (consistent prediction must be performed before sending)
 *
 *  -> calculates based on their current positions
 *
 *  Date    : 2024.01.17  by Seoyoung Kim
 *
 */

void SendAllParticlesToGPU(std::vector <Particle*> &particle) {

	// variables for saving variables to send to GPU
	CUDA_REAL * Mass;
	CUDA_REAL * Mdot;
	CUDA_REAL(*Position)[Dim];
	CUDA_REAL(*Velocity)[Dim];

	// allocate memory to the temporary variables
	Mass     = new CUDA_REAL[NNB];
	Mdot     = new CUDA_REAL[NNB];
	Position = new CUDA_REAL[NNB][Dim];
	Velocity = new CUDA_REAL[NNB][Dim];

	// copy the data of particles to the arrays to be sent
	for (int i=0; i<NNB; i++) {
		Mass[i] = particle[i]->Mass;
		Mdot[i] = 0; //particle[i]->Mass;

		for (int dim=0; dim<Dim; dim++) {
			Position[i][dim] = particle[i]->PredPosition[dim];
			Velocity[i][dim] = particle[i]->PredVelocity[dim];
		}
	}

	// send the arrays to GPU
	//SendToDevice(&NNB, Mass, Position, Velocity, Mdot, &NumNeighborMax);

	// free the temporary variables
	delete[] Mass;
	delete[] Mdot;
	delete[] Position;
	delete[] Velocity;
}


void SendAllParticlesToGPU(CUDA_REAL time, std::vector <Particle*> &particle) {

	// variables for saving variables to send to GPU
	CUDA_REAL * Mass;
	CUDA_REAL * Mdot;
	CUDA_REAL * Radius2;
	CUDA_REAL(*Position)[Dim];
	CUDA_REAL(*Velocity)[Dim];
	int size = (int) particle.size();

	// allocate memory to the temporary variables
	Mass     = new CUDA_REAL[size];
	Mdot     = new CUDA_REAL[size];
	Radius2  = new CUDA_REAL[size];
	Position = new CUDA_REAL[size][Dim];
	Velocity = new CUDA_REAL[size][Dim];


	// copy the data of particles to the arrays to be sent
	for (int i=0; i<size; i++) {
		Mass[i]    = particle[i]->Mass;
		Mdot[i]    = 0; //particle[i]->Mass;
		Radius2[i] = particle[i]->RadiusOfAC*particle[i]->RadiusOfAC; // mass wieght?
		if (particle[i]->NumberOfAC == 0)
			particle[i]->predictParticleSecondOrder(time);
		else
			particle[i]->predictParticleSecondOrderIrr(time);

		for (int dim=0; dim<Dim; dim++) {
			Position[i][dim] = particle[i]->PredPosition[dim];
			Velocity[i][dim] = particle[i]->PredVelocity[dim];
		}
	}

	//fprintf(stdout, "Sending particles to GPU...\n");
	//fflush(stdout);
	// send the arrays to GPU
	SendToDevice(&size, Mass, Position, Velocity, Radius2, Mdot);

	// fprintf(stdout, "Done.\n"); // Eunwoo debug
	fflush(stdout);
	// free the temporary variables
	delete[] Mass;
	delete[] Mdot;
	delete[] Radius2;
	delete[] Position;
	delete[] Velocity;
}

