#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "global.h"
#include <mpi.h>


template <typename T>
void InitialAssignmentOfTasks(T data, int NumTask);




void initializeMPI(int argc, char *argv[]) {
	/* MPI Initialization */
	//MPI_Win win;
	//int MyRank;
	//int NumberOfProcessor;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
	MPI_Comm_size(MPI_COMM_WORLD, &NumberOfProcessor);
	NumberOfWorker = NumberOfProcessor - 1;

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

	int *shared_mem;

	// Create a shared memory communicator
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, MyRank, MPI_INFO_NULL, &shared_comm);

	// Get rank and size in the shared communicator
	int shared_rank, shared_size;
	MPI_Comm_rank(shared_comm, &shared_rank);
	MPI_Comm_size(shared_comm, &shared_size);
	fprintf(stderr,"My Rank =%d : Shared Rank = %d, Shared size = %d\n", MyRank, shared_rank, shared_size);

	// Allocate shared memory
	if (shared_rank == 0) {
		//MPI_Win_allocate_shared(sizeof(int), sizeof(int), MPI_INFO_NULL, shared_comm, &shared_mem, &win);
		MPI_Win_allocate_shared(sizeof(Particle) * MaxNumberOfParticle, sizeof(Particle), MPI_INFO_NULL, shared_comm, &particles_original, &win);
	} else {
		MPI_Win_allocate_shared(0, sizeof(Particle), MPI_INFO_NULL, shared_comm, &particles_original, &win);
		//MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, shared_comm, &shared_mem, &win);
	}
	    // Query shared memory of rank 0
	
	MPI_Aint size_bytes;
	int disp_unit;

	MPI_Win_shared_query(win, 0, &size_bytes, &disp_unit, &particles);
}

void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG) {
	MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
	MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
	int i;
	for (i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		MPI_Isend(&data[i],   1, MPI_INT,    i+1, TAG, MPI_COMM_WORLD, &requests[i]);
	}
	MPI_Waitall(i, requests, statuses);
}





void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG) {
	MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
	MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
	int i;
	for (i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		MPI_Isend(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, MPI_COMM_WORLD, &requests[i]);
		MPI_Isend(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, MPI_COMM_WORLD, &requests[i]);
	}
	MPI_Waitall(i, requests, statuses);
}


void InitialAssignmentOfTasks(int* data, int NumTask, int TAG) {
	MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
	MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
	int i;
	for (i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		MPI_Isend(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[i]);
	}
	MPI_Waitall(i, requests, statuses);
}

void InitialAssignmentOfTasks(int data, int NumTask, int TAG) {
	MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
	MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
	int i;
	for (i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//std::cout << "InitialAssignmentOfTasks out of" << NumTask<< ": " << i << std::endl;
		MPI_Isend(&data, 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[i]);
	}
	//fprintf(stderr, "Number of tasks assigned = %d\n", i);
	//fflush(stderr);
	MPI_Waitall(i, requests, statuses);
}


void broadcastFromRoot(int &data) {
	MPI_Bcast(&data, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
}
void broadcastFromRoot(double &data) {
	MPI_Bcast(&data, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
}
void broadcastFromRoot(ULL &data) {
	MPI_Bcast(&data, 1, MPI_UNSIGNED_LONG_LONG, ROOT, MPI_COMM_WORLD);
}

void ParticleSynchronization() {
		std::cout << "Particle synchronization." << std::endl;
		MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
		MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
		int task=100;
		int completed_tasks = 0;
		//int ptcl_id=104;
		InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
		//std::cerr << "before, Rank=" << MyRank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].CurrentTimeIrr << std::endl;
		//MPI_Win_sync(win);  // Synchronize memory
		//MPI_Win_sync(win);  // Synchronize memory
		//MPI_Win_fence(0, win);
		//MPI_Barrier(shared_comm);
		//std::cerr << "after, Rank=" << MyRank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].CurrentTimeIrr << std::endl;
		for (int i=0; i<NumberOfWorker; i++) {

			MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &requests[i]);
		}
		MPI_Waitall(NumberOfWorker, requests, statuses);
		//std::cerr << "Done" << std::endl;
}

void FurtherAssignmentOfTasks() {

}

