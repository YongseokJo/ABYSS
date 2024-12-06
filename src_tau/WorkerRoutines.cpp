#include <iostream>
#include <vector>
#include <errno.h>
#include "global.h"
#include "particle.h"
#include <vector>

void ComputeAcceleration(int ptcl_id, double next_time);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);

void WorkerRoutines() {

	std::cout << "Processor " << MyRank << " is ready." << std::endl;

	int task=-1;
	MPI_Status status;
	MPI_Request request;
	int ptcl_id;
	double next_time;

	while (true) {
		//MPI_Recv(&task, 1, MPI_INT, ROOT, TASK_TAG, MPI_COMM_WORLD, &status);
		MPI_Irecv(&task, 1, MPI_INT, ROOT, TASK_TAG, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
		//if (status.MPI_TAG == TERMINATE_TAG) break;
		//std::cerr << "Processor " << MyRank << " received task " << task << std::endl;

		switch (task) {
			case 0: // Irregular Acceleration
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);
				particles[ptcl_id].computeAccelerationIrr();
				particles[ptcl_id].NewCurrentBlockIrr = particles[ptcl_id].CurrentBlockIrr + particles[ptcl_id].TimeBlockIrr; // of this particle
				//std::cout << "pid=" << ptcl_id << ", NextBlockIrr=" << particles[ptcl_id].NextBlockIrr \
				//	<< "TimeBlockIrr="<< particles[ptcl_id].TimeBlockIrr << "CurrentBlockIrr="<<particles[ptcl_id].CurrentBlockIrr \
				//	<< "NewCurrentBlockIrr="<< particles[ptcl_id].NewCurrentBlockIrr << std::endl;
			  particles[ptcl_id].calculateTimeStepIrr();
				particles[ptcl_id].NextBlockIrr = particles[ptcl_id].NewCurrentBlockIrr + particles[ptcl_id].TimeBlockIrr; // of this particle
				//std::cout << "pid=" << ptcl_id << ", NextBlockIrr=" << particles[ptcl_id].NextBlockIrr \
				/<< "TimeBlockIrr="<< particles[ptcl_id].TimeBlockIrr 	<< std::endl;
				//std::cout << "IrrCal done " << MyRank << std::endl;
				break;

			case 1: // Regular Acceleration
				//std::cout << "RegCal start " << MyRank << std::endl;
				MPI_Recv(&ptcl_id,   1, MPI_INT,    ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);
				particles[ptcl_id].computeAccelerationReg();
				//ComputeAcceleration(ptcl_id, next_time);
				//std::cout << "RegCal end" << MyRank << std::endl;
				break;

			case 2: // Irregular Update Particle
				//std::cout << "IrrUp Processor " << MyRank << std::endl;
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
			  particles[ptcl_id].updateParticle();
				particles[ptcl_id].CurrentBlockIrr = particles[ptcl_id].NewCurrentBlockIrr;
				particles[ptcl_id].CurrentTimeIrr  = particles[ptcl_id].CurrentBlockIrr*time_step;
				//std::cout << "pid=" << ptcl_id << ", CurrentBlockIrr=" << particles[ptcl_id].CurrentBlockIrr << std::endl;
				//std::cout << "IrrUp end " << MyRank << std::endl;
				break;

			case 3: // Regular Update Particle
				//std::cout << "RegUp start " << MyRank << std::endl;
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				//std::cout << "ptcl " << ptcl_id << std::endl;
			  particles[ptcl_id].updateParticle();

				for (int i=0; i<particles[ptcl_id].NewNumberOfNeighbor; i++)
					particles[ptcl_id].Neighbors[i] = particles[ptcl_id].NewNeighbors[i];
				particles[ptcl_id].NumberOfNeighbor = particles[ptcl_id].NewNumberOfNeighbor;

				particles[ptcl_id].CurrentBlockReg += particles[ptcl_id].TimeBlockReg;
				particles[ptcl_id].CurrentTimeReg   = particles[ptcl_id].CurrentBlockReg*time_step;
			  particles[ptcl_id].calculateTimeStepReg();
			  particles[ptcl_id].calculateTimeStepIrr();
			  particles[ptcl_id].updateRadius();
				if (particles[ptcl_id].NumberOfNeighbor == 0) {
					particles[ptcl_id].CurrentBlockIrr = particles[ptcl_id].CurrentBlockReg;
					particles[ptcl_id].CurrentTimeIrr = particles[ptcl_id].CurrentBlockReg*time_step;
				}
				particles[ptcl_id].NextBlockIrr = particles[ptcl_id].CurrentBlockIrr + particles[ptcl_id].TimeBlockIrr; // of this particle
				//std::cout << "RegUp end " << MyRank << std::endl;
				break;

			case 8: // Initialize Acceleration
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				CalculateAcceleration01(&particles[ptcl_id]);
				break;

			case 9: // Initialize Acceleration
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				CalculateAcceleration23(&particles[ptcl_id]);
				particles[ptcl_id].initializeTimeStep();
				break;

			case 10: // Initialize Timestep variables
				broadcastFromRoot(time_block);
				broadcastFromRoot(block_max);
				broadcastFromRoot(time_step);
				MPI_Win_sync(win);  // Synchronize memory
				MPI_Barrier(shared_comm);
				//MPI_Win_fence(0, win);
				//fprintf(stderr, "nbody+:time_block = %d, EnzoTimeStep=%e\n", time_block, EnzoTimeStep);
				//fflush(stderr);
				break;


			case 100: // Synchronize
				MPI_Win_sync(win);  // Synchronize memory
				MPI_Barrier(shared_comm);
				break;

			case -100: // Simualtion ends
				std::cout << "Processor " << MyRank<< " returns." << std::endl;
				return;
				break;

			case -1:
				perror("Error task assignments");
				exit(EXIT_FAILURE);
				break;
			default:
				break;
		}

		// return that it's over
		//task = -1;
		if (task == 0 || task == 1 || task == 2 || task == 3)
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

