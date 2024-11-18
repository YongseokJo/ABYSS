#include <iostream>
#include <vector>
#include <errno.h>
#include "global.h"
#include "particle.h"
#include <vector>

void ComputeAcceleration(int ptcl_id, double next_time);

void WorkerRoutines() {
	int task;
	MPI_Status status;
	MPI_Request request;
	int ptcl_id;
	double next_time;

	while (true) {
		MPI_Recv(&task, 1, MPI_INT, ROOT, TASK_TAG, MPI_COMM_WORLD, &status);
		//if (status.MPI_TAG == TERMINATE_TAG) break;
		std::cout << "Processor " << MyRank << " received task " << task << std::endl;

		switch (task) {
			case 0: // Acceleration
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);
				//ComputeAcceleration(ptcl_id, next_time);
				break;

			case 1: // Update Particle
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);
			  //particles[ptcl_id].update(next_time);
				break;


			case 9: // Initialize Acceleration
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				particles[ptcl_id].initializeAcceleration();
				particles[ptcl_id].initializeTimeStep();
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
		MPI_Isend(&task, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD,&request);
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

