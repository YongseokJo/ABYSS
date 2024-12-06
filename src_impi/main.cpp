#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "global.h"
#include <mpi.h>
#include <unistd.h>


void broadcastFromRoot(int &data);
void DefaultGlobal();
void initializeMPI(int argc, char *argv[]);
void WorkerRoutines();
void RootRoutines();
int Parser(int argc, char *argv[]);
int readData();
int readParameterFile();


int main(int argc, char *argv[]) {

	//binout = fopen("binary_output.txt", "w");
	//fprintf(binout, "Starting nbody - Binary OUTPUT\n"); 




	/* Initialize global variables */
	DefaultGlobal();

	/* MPI Initialization */
	initializeMPI(argc, argv);


	/*
	// Insert this function definition at the top of your code after the include directives.
	char hostname[256];
	gethostname(hostname, sizeof(hostname));

	// Insert this code right after the  MPI initialization routines (though not a mandatory requirement 
	// to add there only). Please make a judgement based on your code.
	// Retrieve process ID and hostname
	pid_t pid = getpid();

	volatile int i = 0;
	while (0 == i)
	{
		std::cout << "My rank = " << MyRank << " PID = " << pid << " running on Host = " << hostname << " in sleep " << std::endl;
		sleep(5);
	}
	*/

	/* Input options */
	//Parser(argc, argv);
	readParameterFile();

	// Write Particles
	if (MyRank == ROOT && readData() == _FAIL)
		fprintf(stderr, "Read Data Failed!\n");

	broadcastFromRoot(NumberOfParticle);
	//MPI_Win_fence(0, win);

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

