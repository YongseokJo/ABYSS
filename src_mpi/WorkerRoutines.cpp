#include <iostream>
#include <vector>
#include <errno.h>
#include "global.h"
#include "Queue.h"

void ComputeAcceleration(int ptcl_id, double next_time);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void makePrimordialGroup(Particle* ptclCM);
void NewFBInitialization(Particle* ptclCM);
void deleteGroup(Particle* ptclCM);
void NewFBInitialization3(Group* group);

void WorkerRoutines() {

	//std::cout << "Processor " << MyRank << " is ready." << std::endl;

	TaskName task = Error;
	MPI_Status status;
	MPI_Request request;
	int ptcl_id;
	double next_time;
	int NewNumberOfNeighbor;
	int NewNeighbors[NumNeighborMax];
	int size=0;
	double new_a[Dim];
	double new_adot[Dim];
	Particle *ptcl;
	std::chrono::high_resolution_clock::time_point start_point;
	std::chrono::high_resolution_clock::time_point end_point;

	while (true) {
		MPI_Recv(&task, 1, MPI_INT, ROOT, TASK_TAG, MPI_COMM_WORLD, &status);
		//MPI_Irecv(&task, 1, MPI_INT, ROOT, TASK_TAG, MPI_COMM_WORLD, &request);
		//MPI_Wait(&request, &status);
		//if (status.MPI_TAG == TERMINATE_TAG) break;
		//std::cerr << "Processor " << MyRank << " received task " << task << std::endl;

		switch (task) {
			case IrrForce: // Irregular Acceleration
				MPI_Recv(&ptcl_id,   1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				//std::cout << "(IRR_FORCE) Processor " << MyRank<< ": PID= "<<ptcl_id << std::endl;
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status); // (Query to myself) it seems like it's not needed.
#ifdef PerformanceTrace
				ptcl = &particles[ptcl_id];
				start_point = std::chrono::high_resolution_clock::now();
				ptcl->computeAccelerationIrr();
				end_point = std::chrono::high_resolution_clock::now();
				performance.IrregularForce +=
					std::chrono::duration_cast<std::chrono::nanoseconds>(end_point - start_point).count();
#else
				ptcl->computeAccelerationIrr();
#endif

				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->calculateTimeStepIrr();
				ptcl->NextBlockIrr = ptcl->NewCurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->isUpdateToDate = true;
				//std::cout << "IrrCal done " << MyRank << std::endl;
				break;

			case RegForce: // Regular Acceleration
				//std::cout << "RegCal start " << MyRank << std::endl;
				MPI_Recv(&ptcl_id,   1, MPI_INT,    ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);

				particles[ptcl_id].computeAccelerationReg();
				//ComputeAcceleration(ptcl_id, next_time);
				//std::cout << "RegCal end" << MyRank << std::endl;
				break;

			case IrrUpdate: // Irregular Update Particle
				//std::cout << "IrrUp Processor " << MyRank << std::endl;
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl = &particles[ptcl_id];

				if (ptcl->NumberOfNeighbor != 0) // IAR modified
					ptcl->updateParticle();
				ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
				ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;
				//std::cout << "pid=" << ptcl_id << ", CurrentBlockIrr=" << particles[ptcl_id].CurrentBlockIrr << std::endl;
				//std::cout << "IrrUp end " << MyRank << std::endl;
				break;

			case RegUpdate: // Regular Update Particle
				//std::cout << "RegUp start " << MyRank << std::endl;
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				//std::cout << "ptcl " << ptcl_id << std::endl;

				ptcl = &particles[ptcl_id];
				ptcl->updateParticle();

				for (int i=0; i<ptcl->NewNumberOfNeighbor; i++)
					ptcl->Neighbors[i] = ptcl->NewNeighbors[i];
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

				ptcl->CurrentBlockReg += ptcl->TimeBlockReg;
				ptcl->CurrentTimeReg   = ptcl->CurrentBlockReg*time_step;
				ptcl->calculateTimeStepReg();
				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockReg;
				ptcl->calculateTimeStepIrr();
				ptcl->updateRadius();
				if (ptcl->NumberOfNeighbor == 0) {
					ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
					ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg*time_step;
				}
				ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				break;

			case RegCuda: // Update Regular Particle CUDA
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//std::cout << "(REG_CUDA) Processor " << MyRank<< ": PID= "<<ptcl_id << std::endl;
				MPI_Recv(&NewNumberOfNeighbor, 1, MPI_INT, ROOT, 10, MPI_COMM_WORLD, &status);
				MPI_Recv(NewNeighbors, NewNumberOfNeighbor, MPI_INT, ROOT, 11, MPI_COMM_WORLD, &status);
				MPI_Recv(new_a, 3, MPI_DOUBLE, ROOT, 12, MPI_COMM_WORLD, &status);
				MPI_Recv(new_adot, 3, MPI_DOUBLE, ROOT, 13, MPI_COMM_WORLD, &status);
				particles[ptcl_id].updateRegularParticleCuda(NewNeighbors, NewNumberOfNeighbor, new_a, new_adot);
				break;

			case RegCudaUpdate: // Update Regular Particle CUDA II
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				ptcl = &particles[ptcl_id];

				for (int j = 0; j < ptcl->NewNumberOfNeighbor; j++)
					ptcl->Neighbors[j] = ptcl->NewNeighbors[j];
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

				ptcl->updateParticle();
				ptcl->CurrentBlockReg = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;
				ptcl->CurrentTimeReg = ptcl->CurrentBlockReg * time_step;
				ptcl->calculateTimeStepReg();
				ptcl->calculateTimeStepIrr();
				if (ptcl->NumberOfNeighbor == 0) {
					ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
					ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg*time_step;
				}
				ptcl->updateRadius();
				ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of ptcl particle
				break;

			case InitAcc1: // Initialize Acceleration(01)
				//std::cout << "Processor " << MyRank<< " initialization starts." << std::endl;
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				//std::cout << "Processor " << MyRank<< ": PID= "<<ptcl_id << std::endl;
				ptcl = &particles[ptcl_id];
				//std::cerr << ptcl_id+i << std::endl;
				CalculateAcceleration01(ptcl);
				//std::cout << "Processor " << MyRank<< " done." << std::endl;
				break;

			case InitAcc2: // Initialize Acceleration(23)
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				//std::cout << "Processor " << MyRank<< ": PID= "<<ptcl_id << std::endl;
				ptcl = &particles[ptcl_id];
				CalculateAcceleration23(ptcl);
				break;

			case InitTime: // Initialize Time Step
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				//std::cout << "Processor " << MyRank<< ": PID= "<<ptcl_id << std::endl;
				ptcl = &particles[ptcl_id];
				if (ptcl->isActive)
					ptcl->initializeTimeStep();
				break;

			case TimeSync: // Initialize Timestep variables
				broadcastFromRoot(time_block);
				broadcastFromRoot(block_max);
				broadcastFromRoot(time_step);
				//MPI_Win_sync(win);  // Synchronize memory
				//MPI_Barrier(shared_comm);
				//MPI_Win_fence(0, win);
				fprintf(stderr, "(%d) nbody+:time_block = %d, EnzoTimeStep=%e\n", MyRank, time_block, EnzoTimeStep);
				fflush(stderr);
				break;


#ifdef FEWBODY
			case SearchPrimordialGroup: // Primordial binary search
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl = &particles[ptcl_id];

				ptcl->NewNumberOfNeighbor = 0;
				ptcl->checkNewGroup2();

				break;

			case SearchGroup: // Few-body group search
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl = &particles[ptcl_id];
				// std::cerr << "FB search of particle  " << ptcl_id << " is initiated on rank " << MyRank << "." <<std::endl;
				ptcl->NewNumberOfNeighbor = 0;

				if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody) {
					ptcl->checkNewGroup2();
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
				}
				else {
					if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 < TSEARCH)
						ptcl->checkNewGroup();
				}
				/*
				if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody) {
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
					std::cout << "ptcl PID: " << ptcl->PID << ", ptcl NewNumberOfNeighbor: " << ptcl->NewNumberOfNeighbor << std::endl;
				}
				else {
					ptcl->NewNumberOfNeighbor = 0;
					if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 < TSEARCH)
						ptcl->checkNewGroup();
				}
				*/
				// std::cerr << "FB search of particle  " << ptcl_id << " is successfully finished on rank " << MyRank << "." <<std::endl;

				break;

			case MakePrimordialGroup: // Make a primordial group
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl = &particles[ptcl_id];

				makePrimordialGroup(ptcl);

				break;

			case MakeGroup: // Make a group
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl = &particles[ptcl_id];

				NewFBInitialization(ptcl);
				std::cout << "FewBody object of particle " << ptcl->PID
						  << " is successfully initialized on rank " << MyRank << "." <<std::endl;
				break;

			case DeleteGroup: // Delete a Group struct
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl = &particles[ptcl_id];
				deleteGroup(ptcl);
				break;

			case ARIntegration: // SDAR for few body encounters
				MPI_Recv(&ptcl_id,   1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);
				
				ptcl = &particles[ptcl_id];
				// std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << std::endl;

				/* (Query) this will be done already. 
				ptcl->computeAccelerationIrr();
				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->calculateTimeStepIrr();
				ptcl->NextBlockIrr = ptcl->NewCurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				*/

				if (!ptcl->isCMptcl || ptcl->GroupInfo == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->isCMptcl=%d ptcl->GroupInfo=%p\n", ptcl->isCMptcl, ptcl->GroupInfo);
					exit(EXIT_FAILURE);
				}
				
				ptcl->GroupInfo->ARIntegration(next_time);
				if (!ptcl->GroupInfo->isMerger && !ptcl->GroupInfo->isTerminate)
					ptcl->GroupInfo->isTerminate = ptcl->GroupInfo->CheckBreak();

				if (ptcl->GroupInfo->isTerminate) {
					if (ptcl->getBinaryInterruptState() == BinaryInterruptState::none)
						ptcl->setBinaryInterruptState(BinaryInterruptState::terminated);

					delete ptcl->GroupInfo;

					std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << " deleted!" <<std::endl;
				}
				else
					std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << " done!" <<std::endl;
				break;
			
			case MergeManyBody: // Merger insided many-body (>2) group

				MPI_Recv(&ptcl_id,   1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				
				ptcl = &particles[ptcl_id];
				std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << std::endl;

				if (!ptcl->isCMptcl || ptcl->GroupInfo == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->isCMptcl=%d ptcl->GroupInfo=%p\n", ptcl->isCMptcl, ptcl->GroupInfo);
					exit(EXIT_FAILURE);
				}

				NewFBInitialization3(ptcl->GroupInfo);

				ptcl->GroupInfo->isMerger = false;
				ptcl->setBinaryInterruptState(BinaryInterruptState::none);

				std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << " NewFBInitialization3 done!" <<std::endl;
				break;
#endif 

			case Synchronize: // Synchronize
				MPI_Win_sync(win);  // Synchronize memory
				MPI_Barrier(shared_comm);
				break;

			case Ends: // Simualtion ends
				std::cout << "Processor " << MyRank<< " returns." << std::endl;
				//return;
				break;

			case Error:
				perror("Error task assignments");
				exit(EXIT_FAILURE);
				break;
			default:
				break;
		}

		// return that it's over
		//task = -1;
		if (task == IrrForce || task == RegForce || task == IrrUpdate || task == RegUpdate)
			MPI_Isend(&ptcl_id, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD,&request);
		else
			MPI_Isend(&task, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD,&request);

		MPI_Wait(&request, &status);
		//std::cerr << "Processor " << MyRank << " done." << std::endl;
	}
}


void ComputeAcceleration(int ptcl_id, double next_time) {
	Particle *ptcl = &particles[ptcl_id];
	Particle *neighbor;
	double neighbor_pos[Dim], neighbor_vel[Dim];
	for (int i=0; i<ptcl->NumberOfNeighbor; i++) {
		neighbor = &particles[ptcl->Neighbors[i]];
		//neighbor->predictParticle(next_time, neighbor_pos, neighbor_vel);
		//ptcl->getAcceleration(neighbor_pos, neighbor_vel);
	}
}

