#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "global.h"
#include <mpi.h>


void DefaultGlobal();
void initializeMPI(int argc, char *argv[]);
void WorkerRoutines();
void RootRoutines();
int Parser(int argc, char *argv[]);
int readData();


int main(int argc, char *argv[]) {

	std::cout << "Starting ABYSS ..." <<  std::endl;
	//binout = fopen("binary_output.txt", "w");
	//fprintf(binout, "Starting nbody - Binary OUTPUT\n"); 

	/* Initialize global variables */
	DefaultGlobal();

	/* Input options */
	Parser(argc, argv);

	/* MPI Initialization */
	initializeMPI(argc, argv);

	// Write Particles
	if (MyRank == ROOT && readData() == FAIL)
		fprintf(stderr, "Read Data Failed!\n");

	// Synchronize all processes before reading from the shared memory
	MPI_Win_fence(0, win);

	if (MyRank == ROOT) {
		RootRoutines();
	} else {
		WorkerRoutines();
	}

	// Finalize the window and MPI environment
	MPI_Win_free(&win);
	MPI_Finalize();
	return 0;
}

