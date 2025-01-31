#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
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
		MPI_Win_allocate_shared(sizeof(Particle) * MaxNumberOfParticle, sizeof(Particle),
								MPI_INFO_NULL, shared_comm, &particles_original, &win);
		MPI_Win_allocate_shared(sizeof(GlobalVariable), sizeof(GlobalVariable),
								MPI_INFO_NULL, shared_comm, &global_variable_original, &win2);
		MPI_Win_allocate_shared(sizeof(int) * MaxNumberOfParticle, sizeof(int),
								MPI_INFO_NULL, shared_comm, &ActiveIndexToOriginalIndex_original, &win3);
		MPI_Win_allocate_shared(sizeof(int) * NumberOfProcessor, sizeof(int),
								MPI_INFO_NULL, shared_comm, &tasks_original, &win4);
		MPI_Win_allocate_shared(sizeof(int) * NumberOfProcessor*(MAX_QUEUE+1), sizeof(int),
								MPI_INFO_NULL, shared_comm, &queues_original, &win5);
	} else {
		MPI_Win_allocate_shared(0, sizeof(Particle), MPI_INFO_NULL, shared_comm, &particles_original, &win);
		MPI_Win_allocate_shared(0, sizeof(GlobalVariable), MPI_INFO_NULL, shared_comm, &global_variable_original, &win2);
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, shared_comm, &ActiveIndexToOriginalIndex_original, &win3);
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, shared_comm, &tasks_original, &win4);
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, shared_comm, &queues_original, &win5);
	}
	// Query shared memory of rank 0
	
	MPI_Aint size_bytes;
	int disp_unit;

	MPI_Win_shared_query(win, 0, &size_bytes, &disp_unit, &particles);
	MPI_Win_shared_query(win2, 0, &size_bytes, &disp_unit, &global_variable);
	MPI_Win_shared_query(win3, 0, &size_bytes, &disp_unit, &ActiveIndexToOriginalIndex);
	MPI_Win_shared_query(win4, 0, &size_bytes, &disp_unit, &tasks);
	MPI_Win_shared_query(win5, 0, &size_bytes, &disp_unit, &queues);
}

void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i],   1, MPI_INT,    i+1, TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Send(&data[i],   1, MPI_INT,    i+1, TAG, MPI_COMM_WORLD);
	}
}


void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		//MPI_Isend(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);

		MPI_Send(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, MPI_COMM_WORLD);
		MPI_Send(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, MPI_COMM_WORLD);
	}
}


void InitialAssignmentOfTasks(int* data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Send(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD);
	}
}

void InitialAssignmentOfTasks(int data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		//std::cerr << "InitialAssignmentOfTasks out of" << NumTask<< ": " << i << std::endl;
		if (i >= NumTask)
			break;
		//MPI_Isend(&data, 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Send(&data, 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD);
	}
	//fprintf(stderr, "Number of tasks assigned = %d\n", i);
	//fflush(stderr);
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
		//std::cout << "Particle synchronization." << std::endl;
		MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
		MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
		int task=100;
		int completed_tasks = 0;
		//int ptcl_id=104;
		InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
		//std::cerr << "before, Rank=" << MyRank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].CurrentTimeIrr << std::endl;
		MPI_Win_sync(win);  // Synchronize memory
		MPI_Barrier(shared_comm);
		//std::cerr << "after, Rank=" << MyRank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].CurrentTimeIrr << std::endl;
		for (int i=0; i<NumberOfWorker; i++) {

			MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &requests[i]);
		}
		MPI_Waitall(NumberOfWorker, requests, statuses);
		//std::cerr << "Done" << std::endl;
}

void FurtherAssignmentOfTasks() {

}

