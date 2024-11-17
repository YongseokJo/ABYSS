#include <mpi.h>
#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "global.h"


void DefaultGlobal();
void WorkerRoutines();
void RootRoutines(int num_workers, std::vector<int> tasks);


int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);



	if (size < 2) {
		std::cerr << "This program requires at least 2 processes.\n";
		MPI_Finalize();
		return 1;
	}

	/*
	// comm for each node
	MPI_Comm shmcomm;
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
	int local_rank, local_size;
	MPI_Comm_rank(shmcomm, &local_rank);
	MPI_Comm_size(shmcomm, &local_size);
	*/


	DefaultGlobal();

	Particle *particles;
	MPI_Win win;

	// Declare shared memory pointer and window object
	int NumberOfParticle = 1000;

	// Allocate shared memory window
	MPI_Win_allocate_shared(sizeof(Particle) * MaxNumberOfParticle, sizeof(Particle), MPI_INFO_NULL, MPI_COMM_WORLD, &particles, &win);

	// Write Particles

	// Synchronize all processes before reading from the shared memory
	MPI_Win_fence(0, win);



	if (MyRank == ROOT) {
		/* Initial configuration */
		int num_tasks = 10;
		std::vector<int> tasks(num_tasks);
		for (int i = 0; i < num_tasks; ++i) tasks[i] = i + 1;
		int num_workers = size - 1;

		RootRoutines(num_workers, tasks)



	} else {
		WorkerRoutines();
	}

	// Finalize the window and MPI environment
	MPI_Win_free(&win);
	MPI_Finalize();
	return 0;
}

