#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <mpi.h>
#include "../global.h"
#include "../QueueScheduler.h"
#include "cuda_functions.h"

void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG);
void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int* data, int NumTask, int TAG);
void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList);
void CalculateAccelerationOnDevice(int *NumTargetTotal, int *h_target_list, double acc[][3], double adot[][3], int NumNeighbor[], int *NeighborList);

/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU
 *
 *  Date    : 2024.01.18  by Seoyoung Kim
 *
 */
void calculateRegAccelerationOnGPU(std::unordered_set<int> RegularList, QueueScheduler &queue_scheduler){



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
#ifdef DEBUG
	std::cout << "sendAllParticlesToGPU starts" << std::endl;
#endif
	sendAllParticlesToGPU(new_time, RegularList, IndexList);  // needs to be updated
#ifdef DEBUG
	std::cout << "sendAllParticlesToGPU ended" << std::endl;
#endif
	
	
	/*
	for (int i=0; i<ListSize; i++) {
		IndexList[i] = RegularList[i];
	} // endfor copy info
	*/


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
#ifdef DEBUG
	std::cout << "CalculateAccelerationOnDevice starts" << std::endl;
#endif
	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive, AccRegDotReceive, NumNeighborReceive, ACListReceive);
#ifdef DEBUG
	std::cout << "CalculateAccelerationOnDevice ended" << std::endl;
#endif

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



/*
	std::cout << "(REG_CUDA) RegularList, PID= ";
	for (int i=0; i<RegularList.size(); i++) {
		std::cout << RegularList[i]<< ", ";
	}
	std::cout << std::endl;
	*/

#ifdef DEBUG
	std::cout << "Adjust Regular Gravity starts" << std::endl;
#endif

	// Adjust Regular Gravity
	int i=0, task=REG_CUDA;
	queue_scheduler.initialize(REG_CUDA);
	queue_scheduler.takeQueueRegularList(RegularList);
	do
	{
		queue_scheduler.assignQueueRegularList();

        for (auto worker = queue_scheduler.WorkersToGo.begin(); worker != queue_scheduler.WorkersToGo.end();)
        {
            if ((*worker)->NumberOfQueues > 0) // original
            {
				//std::cout << "(REG_CUDA) My Rank =" << (*worker)->MyRank << std::endl;
				MPI_Send(&task, 1, MPI_INT, (*worker)->MyRank, TASK_TAG, MPI_COMM_WORLD);
				MPI_Send(&ActiveIndexToOriginalIndex[IndexList[i]], 1, MPI_INT, (*worker)->MyRank, PTCL_TAG, MPI_COMM_WORLD);
				MPI_Send(&NumNeighborReceive[i], 1, MPI_INT, (*worker)->MyRank, 10, MPI_COMM_WORLD);
				MPI_Send(&ACListReceive[i * NumNeighborMax], NumNeighborReceive[i], MPI_INT, (*worker)->MyRank, 11, MPI_COMM_WORLD);
				MPI_Send(&AccRegReceive[i][0], 3, MPI_DOUBLE, (*worker)->MyRank, 12, MPI_COMM_WORLD);
				MPI_Send(&AccRegDotReceive[i][0], 3, MPI_DOUBLE, (*worker)->MyRank, 13, MPI_COMM_WORLD);
				((*worker))->onDuty = true;
				/*
				(*worker)->CurrentQueue++;
				(*worker)->CurrentQueue %= MAX_QUEUE;
				(*worker)->NumberOfQueues--;
				*/
                worker = queue_scheduler.WorkersToGo.erase(worker);
				i++;
#ifdef DEBUG
				std::cout << "i: " << i << std::endl;
#endif
            }
			else
			{
                ++worker;
#ifdef DEBUG
				std::cout << "worker MyRank: " << (*worker)->MyRank << std::endl;
#endif
			}
        }
#ifdef DEBUG
		std::cout << "queue_scheduler.waitQueue(0) starts" << std::endl;
#endif
		queue_scheduler.waitQueue(0); // blocking wait
#ifdef DEBUG
		std::cout << "queue_scheduler.waitQueue(0) ended" << std::endl;
#endif
	} while (queue_scheduler.isComplete());

#ifdef DEBUG
	std::cout << "Adjust Regular Gravity ended" << std::endl;
#endif


#ifdef nouse
	ws.initialize();
	int return_value, i = 0;
	ws._setTask(4);
	ws._total_tasks = RegularList.size();
	ws._completed_tasks = 0;

	do
	{
		if (ws._FreeWorkers.size() == 0 || ws._assigned_tasks == ws._total_tasks)
		{
			// have to add check all the sends are recved.
			// MPI_Waitall(NumberOfCommunication, requests, statuses);
			// NumberOfCommunication = 0;
			ws._checkCompletion(return_value);
			/* we can do something here */
			ws._Callback();
			ws._completed_tasks++;
		}
		if (ws._FreeWorkers.size() != 0 && ws._assigned_tasks < ws._total_tasks)
		{
			ws._WorkerTmp = ws._FreeWorkers.back();
			ws._FreeWorkers.pop_back();
			MPI_Send(&ws._WorkerTmp->task, 1, MPI_INT, ws._WorkerTmp->MyRank, TASK_TAG, MPI_COMM_WORLD);
			MPI_Send(&RegularList[i], 1, MPI_INT, ws._WorkerTmp->MyRank, PTCL_TAG, MPI_COMM_WORLD);
			MPI_Send(&NumNeighborReceive[i], 1, MPI_INT, ws._WorkerTmp->MyRank, 10, MPI_COMM_WORLD);
			MPI_Send(&ACListReceive[i * NumNeighborMax], NumNeighborReceive[i], MPI_INT, ws._WorkerTmp->MyRank, 11, MPI_COMM_WORLD);
			MPI_Send(&AccRegReceive[i][0], 3, MPI_DOUBLE, ws._WorkerTmp->MyRank, 12, MPI_COMM_WORLD);
			MPI_Send(&AccRegDotReceive[i][0], 3, MPI_DOUBLE, ws._WorkerTmp->MyRank, 13, MPI_COMM_WORLD);
			ws._WorkerTmp->onDuty = true;
			ws._assigned_tasks++;
			//fprintf(stdout, "assigned_tasks = %d/%d, number of free worker = %d pid = %d rank = %d\n",
					//ws._assigned_tasks, ws._total_tasks, ws._FreeWorkers.size(), RegularList[i], ws._WorkerTmp->MyRank);
		}
		i++;
	} while (ws._completed_tasks < ws._total_tasks);
#endif

	delete[] IndexList;

	delete[] AccRegReceive;
	delete[] AccRegDotReceive;
	delete[] AccIrr;
	delete[] AccIrrDot;

	delete[] NumNeighborReceive;
	delete[] ACListReceive;


#ifdef time_trace
	_time.reg_cpu3.markEnd();
	_time.reg_cpu3.getDuration();
#endif
	//CloseDevice();
} // calculate 0th, 1st derivative of force + neighbors on GPU ends





void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList) {

	// variables for saving variables to send to GPU
	double * Mass;
	double * Mdot;
	double * Radius2;
	double(*Position)[Dim];
	double(*Velocity)[Dim];
	//int size = NumberOfParticle;
	int size=0, j=0;


	// allocate memory to the temporary variables
	Mass     = new double[NumberOfParticle];
	Mdot     = new double[NumberOfParticle];
	Radius2  = new double[NumberOfParticle];
	Position = new double[NumberOfParticle][Dim];
	Velocity = new double[NumberOfParticle][Dim];

	Particle *ptcl;

	// copy the data of particles to the arrays to be sent
		
	for (int i=0; i<=LastParticleIndex; i++) {
		ptcl = &particles[i];

		if (!ptcl->isActive) {
			// fprintf(stdout, "Skipping inactive particle (%d)\n", ptcl->PID);
			continue;
		}

		if (RegularList.find(i) != RegularList.end()) {
			IndexList[j] = size;
			j++;
		}


		Mass[size]    = ptcl->Mass;
		Mdot[size]    = 0; //particle[i]->Mass;
		Radius2[size] = ptcl->RadiusOfNeighbor; // mass wieght?

		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position[size], Velocity[size]);
		else
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, Position[size], Velocity[size]);

		ActiveIndexToOriginalIndex[size] = i;
		// std::cout << "(size , i) = "  << size << " " << i << std::endl;
		size++;
	} 

	// fprintf(stdout, "in sendAllParticlesToGPU, NumberOfParticle = %d, size=%d, TotalNumberOfParticle=%d\n", NumberOfParticle, size, LastParticleIndex+1);

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

