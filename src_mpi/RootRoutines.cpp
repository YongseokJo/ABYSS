#include "global.h"
#include <vector>
#include <iostream>
#include <mpi.h>


void RootRoutines() {


	Particle* ptcl;
	int min_time_level=0;
	int workers=0, worker_rank;
	int task = 0, task_recv, total_tasks;
	int remaining_tasks=0, completed_tasks=0, completed_rank;
	
	MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object

	int sender_rank, sender_tag;
	int ptcl_id;


	/* Initialization */
	{
		task = 9;
		total_tasks = NumberOfParticle;

		// Initial assignments
		worker_rank = 1;
		for (ptcl_id=0; ptcl_id<NumberOfWorker; ptcl_id++) {
			if (total_tasks < ptcl_id)
				break;
			MPI_Send(&task,    1, MPI_INT   , worker_rank, TASK_TAG, MPI_COMM_WORLD);
			MPI_Send(&ptcl_id, 1, MPI_INT   , worker_rank, PTCL_TAG, MPI_COMM_WORLD);
			worker_rank++;
		}

		// further assignments
		remaining_tasks = total_tasks-NumberOfWorker;
		while (completed_tasks < total_tasks) {
			// Check which worker is done
			MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			completed_rank = status.MPI_SOURCE;
			completed_tasks++;

			if (remaining_tasks > 0) {
				ptcl_id = NumberOfWorker + completed_tasks;
				MPI_Send(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD);
				MPI_Send(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD);
				remaining_tasks--;
			} else {
				printf("Rank %d: No more tasks to assign\n", completed_rank);
			}
		}
	}

	/* timestep correction */
	{
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			if (ptcl->NumberOfNeighbor != 0) {
				while (ptcl->TimeLevelIrr >= ptcl->TimeLevelReg) {
					ptcl->TimeStepIrr *= 0.5;
					ptcl->TimeBlockIrr *= 0.5;
					ptcl->TimeLevelIrr--;
				}
			}
			if (ptcl->TimeLevelIrr < min_time_level) {
				min_time_level = ptcl->TimeLevelIrr;
			}
		}

		// resetting time_block based on the system
		time_block = std::max(-60, min_time_level-MIN_LEVEL_BUFFER);
		block_max = static_cast<ULL>(pow(2, -time_block));
		time_step = std::pow(2,time_block);

		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			ptcl->TimeBlockIrr = static_cast<ULL>(pow(2, ptcl->TimeLevelIrr-time_block));
			ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg-time_block));
#ifdef IRR_TEST
			ptcl->TimeStepReg = 1;
			ptcl->TimeLevelReg = 0;
			ptcl->TimeBlockReg = block_max;
#endif
			ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
		}
	}




	fprintf(stdout, "nbody+:time_block = %d, EnzoTimeStep=%e\n", time_block, EnzoTimeStep);


	/* Assign tasks 
	for (int worker = 1; worker < size && next_task < total_tasks; ++worker) {
		MPI_Send(&tasks[next_task], 1, MPI_INT, worker, TASK_TAG, MPI_COMM_WORLD);
		++next_task;
	}*/

	/* Receive the restuls */
	/*
	while (next_task < total_tasks || active_workers > 0) {
		int worker_rank;
		MPI_Status status;
		MPI_Recv(&worker_rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);


		// we have to synchronize after one time step you got it?
		MPI_Win_fence(0, win);

		if (next_task < total_tasks) {
			MPI_Send(&tasks[next_task], 1, MPI_INT, worker_rank, TASK_TAG, MPI_COMM_WORLD);
			++next_task;
		} else {
			MPI_Send(nullptr, 0, MPI_INT, worker_rank, TERMINATE_TAG, MPI_COMM_WORLD);
			--active_workers;
		}
	}
	*/


}




