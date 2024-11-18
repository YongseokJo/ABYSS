#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "global.h"
#include <mpi.h>



void initializeMPI(int argc, char *argv[]) {
	/* MPI Initialization */
	//MPI_Win win;
	//int MyRank;
	//int NumberOfProcessor;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
	MPI_Comm_size(MPI_COMM_WORLD, &NumberOfProcessor);


	if (NumberOfProcessor < 2) {
		std::cerr << "This program requires at least 2 processes.\n";
		MPI_Finalize();
	}
	/*
	// comm for each node
	MPI_Comm shmcomm;
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
	int local_rank, local_size;
	MPI_Comm_rank(shmcomm, &local_rank);
	MPI_Comm_size(shmcomm, &local_size);
	*/

	// Declare shared memory pointer and window object
	// Allocate shared memory window
	MPI_Win_allocate_shared(sizeof(Particle) * MaxNumberOfParticle, sizeof(Particle), MPI_INFO_NULL, MPI_COMM_WORLD, &particles, &win);

}
