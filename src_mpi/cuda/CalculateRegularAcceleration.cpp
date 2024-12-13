#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <mpi.h>
#include "../def.h"
#include "../global.h"
#include "cuda_functions.h"

void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG);
void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int* data, int NumTask, int TAG);
void sendAllParticlesToGPU(double new_time);
void CalculateAccelerationOnDevice(int *NumTargetTotal, int *h_target_list, double acc[][3], double adot[][3], int NumNeighbor[], int *NeighborList);

/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU
 *
 *  Date    : 2024.01.18  by Seoyoung Kim
 *
 */
void calculateRegAccelerationOnGPU(std::vector<int> RegularList){



	// regIds are the list of positions of particles subject to regular force calculation in std::vector list particle

	// variables for opening GPU
	// const int buffer = 10;
	// int numGpuOpen = NNB+buffer;
	//const int NumPtclPerEachCalMax = 2048; // this also caps the number of particles computed each iteration
	int NeighborIndex; // this size should coincide with number of threads
	int ListSize = RegularList.size();
	int *IndexList = new int[ListSize];

	//int NumGpuCal;

	// variables for saving variables to send to GPU
	// only regular particle informations are stored here
	double (*AccRegReceive)[Dim];
	double (*AccRegDotReceive)[Dim];
	double (*AccIrr)[Dim];
	double (*AccIrrDot)[Dim];
	//int (*ACListReceive)[NumNeighborMax];

	//double* PotSend;
	// int **ACListReceive;
	int *ACListReceive;
	int *NumNeighborReceive;
	int MassFlag;


	double a_tmp[Dim]{0}, adot_tmp[Dim]{0};
	double da, dadot;
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;


	double DFR, FRD, SUM, AT3, BT2;
	double DTR, DTSQ, DT2, DT6,DTSQ12, DTR13;

	Particle *ptcl;


	double new_time = NextRegTimeBlock*time_step;  // next regular time


	// need to make array to send to GPU
	// allocate memory to the temporary variables
	//PotSend         = new double[ListSize];
	AccRegReceive    = new double[ListSize][Dim];
	AccRegDotReceive = new double[ListSize][Dim];
	AccIrr           = new double[ListSize][Dim];
	AccIrrDot        = new double[ListSize][Dim];

	NumNeighborReceive  = new int[ListSize];

	// ACListReceive      = new int*[ListSize];
	ACListReceive = new int[ListSize * NumNeighborMax];

	for (int i=0; i<ListSize; i++) {
		// ACListReceive[i] = new int[NumNeighborMax];
		for (int dim=0; dim<Dim; dim++) {
			AccRegReceive[i][dim]    = 0;
			AccRegDotReceive[i][dim] = 0;
			AccIrr[i][dim]           = 0;
			AccIrrDot[i][dim]        = 0;
		}
	}
	sendAllParticlesToGPU(new_time);  // needs to be updated
	
	
	for (int i=0; i<ListSize; i++) {
		IndexList[i] = RegularList[i];
	} // endfor copy info


	//std::cout <<  "Starting Calculation On Device ..." << std::endl;
	// send information of all the particles to GPU
	// includes prediction
/*
#ifdef time_trace
	_time.reg_sendall.markStart();
#endif

	// Particles have been already at T_new through irregular time step

#ifdef time_trace
	_time.reg_sendall.markEnd();
	_time.reg_sendall.getDuration();
#endif

#ifdef time_trace
	_time.reg_gpu.markStart();
#endif
*/

	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive, AccRegDotReceive, NumNeighborReceive, ACListReceive);

	/*
#ifdef time_trace
	_time.reg_gpu.markEnd();
	_time.reg_gpu.getDuration();

	_time.reg_cpu1.markStart();
#endif
*/



	/*
#ifdef time_trace
	_time.reg_cpu2.markEnd();
	_time.reg_cpu2.getDuration();

	_time.reg_cpu3.markStart();
#endif
*/

	//UpdateNextRegTime(particle);

	//for (Particle* ptcl: RegularList)
		//ptcl->calculateTimeStepIrr(ptcl->a_tot,ptcl->a_irr);
	//std::cout <<  "Calculation On Device Done ..." << std::endl;



	// Regular Gravity
	int task = 4;
	int completed_tasks = 0;
	int total_tasks = RegularList.size();
	double next_time = NextRegTimeBlock*time_step;
	int remaining_tasks, ptcl_id, completed_rank;


	MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object
	
	//std::cout << "TotalTask=" << total_tasks << std::endl;

				/*
				std::cout << "RegularList, PID= ";
				for (int i=0; i<total_tasks; i++) {
					std::cout << RegularList[i]<< ", ";
				}*/
				//std::cout << std::endl;

	InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
	InitialAssignmentOfTasks(RegularList, next_time, total_tasks, PTCL_TAG);
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= total_tasks) break;
		//std::cout << "InitialAssignmentOfTasks out of" << NumTask<< ": " << i << std::endl;
		MPI_Isend(&NumNeighborReceive[i] ,                     1, MPI_INT,
			 	i+1, 10, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Isend(&ACListReceive[i*NumNeighborMax], NumNeighborReceive[i], MPI_INT,
			 	i+1, 11, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Isend(&AccRegReceive[i][0]   ,                     3, MPI_DOUBLE,
			 	i+1, 12, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Isend(&AccRegDotReceive[i][0],                     3, MPI_DOUBLE,
			 	i+1, 13, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		/*
		for (int dim=0; dim<Dim; dim++)
			fprintf(stderr, "Root: new_a=%.3e", AccRegReceive[i][dim]);
		fprintf(stderr, "\n");
		*/
	}

	MPI_Waitall(NumberOfCommunication, requests, statuses);
	NumberOfCommunication = 0;


	// further assignments
	int j=0;
	remaining_tasks = total_tasks-NumberOfWorker;
	while (completed_tasks < total_tasks) {
		// Check which worker is done
		MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
		// Poll until send completes
		/*
			 flag=0;
			 while (!flag) {
			 MPI_Test(&request, &flag, &status);
		// Perform other work while waiting
		}
		*/
		MPI_Wait(&request, &status);
		completed_rank = status.MPI_SOURCE;
		//printf("Rank %d: Send operation completed (%d).\n",completed_rank, ptcl_id_return);

		if (remaining_tasks > 0) {
			j = NumberOfWorker + completed_tasks;
			ptcl_id = RegularList[j];
			MPI_Isend(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD,
				 	&requests[NumberOfCommunication++]);
			MPI_Isend(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD,
				 	&requests[NumberOfCommunication++]);

			MPI_Isend(&NumNeighborReceive[j] ,                     1, MPI_INT,
					completed_rank, 10, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
			MPI_Isend(&ACListReceive[j*NumNeighborMax], NumNeighborReceive[j], MPI_INT,
					completed_rank, 11, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
			MPI_Isend(&AccRegReceive[j][0]   ,                     3, MPI_DOUBLE,
					completed_rank, 12, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
			MPI_Isend(&AccRegDotReceive[j][0],                     3, MPI_DOUBLE,
					completed_rank, 13, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);

			MPI_Waitall(NumberOfCommunication, requests, statuses);
			NumberOfCommunication = 0;

			remaining_tasks--;
		}
		completed_tasks++;
	}



	// Update Regular Gravity 
	task = 5;
	completed_tasks = 0;


	//std::cout << "TotalTask=" << total_tasks << std::endl;

				/*
				std::cout << "RegularList, PID= ";
				for (int i=0; i<total_tasks; i++) {
					std::cout << RegularList[i]<< ", ";
				}*/
				//std::cout << std::endl;

	InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
	InitialAssignmentOfTasks(RegularList, next_time, total_tasks, PTCL_TAG);
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= total_tasks) break;
		MPI_Isend(&NumNeighborReceive[i] ,                     1, MPI_INT,
				i+1, 10, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Isend(&ACListReceive[i*NumNeighborMax], NumNeighborReceive[i], MPI_INT,
				i+1, 11, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
	}
	MPI_Waitall(NumberOfCommunication, requests, statuses);
	NumberOfCommunication = 0;


	// further assignments
	j=0;
	remaining_tasks = total_tasks-NumberOfWorker;
	while (completed_tasks < total_tasks) {
		// Check which worker is done
		MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
		// Poll until send completes
		/*
			 flag=0;
			 while (!flag) {
			 MPI_Test(&request, &flag, &status);
		// Perform other work while waiting
		}
		*/
		MPI_Wait(&request, &status);
		completed_rank = status.MPI_SOURCE;
		//printf("Rank %d: Send operation completed (%d).\n",completed_rank, ptcl_id_return);

		if (remaining_tasks > 0) {
			j = NumberOfWorker + completed_tasks;
			ptcl_id = RegularList[j];
			MPI_Isend(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD,
				 	&requests[NumberOfCommunication++]);
			MPI_Isend(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD,
				 	&requests[NumberOfCommunication++]);
			MPI_Isend(&NumNeighborReceive[j] ,                     1, MPI_INT,
					completed_rank, 10, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
			MPI_Isend(&ACListReceive[j*NumNeighborMax], NumNeighborReceive[j], MPI_INT,
					completed_rank, 11, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
			MPI_Waitall(NumberOfCommunication, requests, statuses);
			NumberOfCommunication = 0;
			remaining_tasks--;
		}
		completed_tasks++;
	}





	delete[] NumNeighborReceive;
	delete[] AccRegReceive;
	delete[] AccRegDotReceive;
	delete[] AccIrr;
	delete[] AccIrrDot;
	// for (int i=0; i<ListSize; i++) {
	//	delete[] ACListReceive[i];
	// }
	delete[] ACListReceive;


#ifdef time_trace
	_time.reg_cpu3.markEnd();
	_time.reg_cpu3.getDuration();
#endif
	//CloseDevice();
} // calculate 0th, 1st derivative of force + neighbors on GPU ends





void sendAllParticlesToGPU(double new_time) {

	// variables for saving variables to send to GPU
	double * Mass;
	double * Mdot;
	double * Radius2;
	double(*Position)[Dim];
	double(*Velocity)[Dim];
	int size = NumberOfParticle;

	// allocate memory to the temporary variables
	Mass     = new double[size];
	Mdot     = new double[size];
	Radius2  = new double[size];
	Position = new double[size][Dim];
	Velocity = new double[size][Dim];

	Particle *ptcl;

	// copy the data of particles to the arrays to be sent
	for (int i=0; i<size; i++) {
		ptcl       = &particles[i];
		Mass[i]    = ptcl->Mass;
		Mdot[i]    = 0; //particle[i]->Mass;
		Radius2[i] = ptcl->RadiusOfNeighbor; // mass wieght?

		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position[i], Velocity[i]);
		else
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, Position[i], Velocity[i]);
	}

	//fprintf(stdout, "Sending particles to GPU...\n");
	//fflush(stdout);
	// send the arrays to GPU
	SendToDevice(&size, Mass, Position, Velocity, Radius2, Mdot);

	//fprintf(stdout, "Done.\n");
	//fflush(stdout);
	// free the temporary variables
	delete[] Mass;
	delete[] Mdot;
	delete[] Radius2;
	delete[] Position;
	delete[] Velocity;
}

