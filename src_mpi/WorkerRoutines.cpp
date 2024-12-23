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
	std::vector<int> ptcl_id_vector;
	double next_time;
	int *NewNumberOfNeighbor;
	int *NewNeighbors;
	int size=0;
	double *new_a;
	double *new_adot;
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
			case 0: // Irregular Acceleration
#ifdef LoadBalance
				MPI_Probe(ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_INT, &size);
				ptcl_id_vector.resize(size);
				MPI_Recv(ptcl_id_vector.data(), size, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time,      1   , MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);

				for (int i=0; i<size; i++) {
					ptcl = &particles[ptcl_id_vector[i]];

#ifdef PerformanceTrace
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
				}
#else
				MPI_Recv(&ptcl_id,   1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);

				ptcl = &particles[ptcl_id];
				ptcl->computeAccelerationIrr();
				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->calculateTimeStepIrr();
				ptcl->NextBlockIrr = ptcl->NewCurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
#endif
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
#ifdef LoadBalance
				MPI_Probe(ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_INT, &size);
				ptcl_id_vector.resize(size);
				MPI_Recv(ptcl_id_vector.data(), size, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);

				for (int i=0; i<size; i++) {
					ptcl = &particles[ptcl_id_vector[i]];

					if (ptcl->GroupInfo) {
						ptcl->GroupInfo->sym_int.particles.cm.NumberOfNeighbor = ptcl->NumberOfNeighbor;
						for (int i=0; i<ptcl->NumberOfNeighbor; i++)
							ptcl->GroupInfo->sym_int.particles.cm.Neighbors[i] = ptcl->Neighbors[i];
						bool int_normal = ptcl->GroupInfo->ARIntegration(ptcl->NewCurrentBlockIrr*time_step);
						if (!int_normal) break;
					}
					if (ptcl->NumberOfNeighbor != 0)
						ptcl->updateParticle();
					ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
					ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;
					//std::cout << "pid=" << ptcl_id << ", CurrentBlockIrr=" << particles[ptcl_id].CurrentBlockIrr << std::endl;
				}
#else
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);

				ptcl = &particles[ptcl_id];

				if (ptcl->GroupInfo) {
					ptcl->GroupInfo->sym_int.particles.cm.NumberOfNeighbor = ptcl->NumberOfNeighbor;
					for (int i=0; i<ptcl->NumberOfNeighbor; i++)
						ptcl->GroupInfo->sym_int.particles.cm.Neighbors[i] = ptcl->Neighbors[i];
					bool int_normal = ptcl->GroupInfo->ARIntegration(ptcl->NewCurrentBlockIrr*time_step);
					if (!int_normal) break;
				}
				if (ptcl->NumberOfNeighbor != 0) // IAR modified
					ptcl->updateParticle();
				ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
				ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;
#endif
				//std::cout << "pid=" << ptcl_id << ", CurrentBlockIrr=" << particles[ptcl_id].CurrentBlockIrr << std::endl;
				//std::cout << "IrrUp end " << MyRank << std::endl;
				break;

			case 3: // Regular Update Particle
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
																																				 //std::cout << "RegUp end " << MyRank << std::endl;
				break;

			case 4: // Update Regular Particle CUDA
				//MPI_Recv(&size               , 1                  , MPI_INT, ROOT, ANY_TAG, MPI_COMM_WORLD, &status);
				// Probe to find the size of the incoming message
				MPI_Probe(ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_INT, &size);
				ptcl_id_vector.resize(size);
				MPI_Recv(ptcl_id_vector.data(), size                    , MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				NewNeighbors        = new int[size*NumNeighborMax];
				NewNumberOfNeighbor = new int[size];
				new_a               = new double[size*3];
				new_adot            = new double[size*3];
				MPI_Recv(NewNumberOfNeighbor , size                    , MPI_INT, ROOT, 10, MPI_COMM_WORLD, &status);
				MPI_Recv(NewNeighbors         , size*NumNeighborMax     , MPI_INT, ROOT, 11, MPI_COMM_WORLD, &status);
				MPI_Recv(new_a                , size*3                  , MPI_DOUBLE, ROOT, 12, MPI_COMM_WORLD, &status);
				MPI_Recv(new_adot             , size*3                  , MPI_DOUBLE, ROOT, 13, MPI_COMM_WORLD, &status);
				//MPI_Waitall(NumberOfCommunication, requests, statuses);
				//NumberOfCommunication = 0;
				for (int i=0; i<size; i++) {
					ptcl = &particles[ptcl_id_vector[i]];
					ptcl->updateRegularParticleCuda(NewNeighbors, NewNumberOfNeighbor[i], new_a, new_adot, i);
				}
				delete NewNeighbors;
				delete NewNumberOfNeighbor;
				break;

			case 5: // Update Regular Particle CUDA II
				MPI_Probe(ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_INT, &size);
				ptcl_id_vector.resize(size);
				MPI_Recv(ptcl_id_vector.data(), size                    , MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


				for (int i=0; i<size; i++) {
					ptcl = &particles[ptcl_id_vector[i]];

					for (int j=0; j<ptcl->NewNumberOfNeighbor; j++)
						ptcl->Neighbors[j] = ptcl->NewNeighbors[j];
					ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

					ptcl->updateParticle();
					ptcl->CurrentBlockReg = ptcl->CurrentBlockReg+ptcl->TimeBlockReg;
					ptcl->CurrentTimeReg  = ptcl->CurrentBlockReg*time_step;
					ptcl->calculateTimeStepReg();
					ptcl->calculateTimeStepIrr();
					/* // IAR original
					if (ptcl->NumberOfNeighbor == 0) {
						ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
						ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg*time_step;
					}
					*/ // IAR original
					ptcl->updateRadius();
					ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of ptcl particle
				}
				break;

			case 7: // Initialize Acceleration(01)
				//std::cout << "Processor " << MyRank<< " initialization starts." << std::endl;
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl_id = ptcl_id*LoadBalanceParticle;

				for (int i=0; i<LoadBalanceParticle; i++) {
					//std::cerr << ptcl_id+i << std::endl;
					CalculateAcceleration01(&particles[ptcl_id+i]);
				}
				//std::cout << "Processor " << MyRank<< " done." << std::endl;
				break;

			case 8: // Initialize Acceleration(23)
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl_id = ptcl_id*LoadBalanceParticle;

				for (int i=0; i<LoadBalanceParticle; i++) {
					CalculateAcceleration23(&particles[ptcl_id+i]);
				}
				break;

			case 9: // Initialize Time Step
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl_id = ptcl_id*LoadBalanceParticle;

				for (int i=0; i<LoadBalanceParticle; i++) {
					if (!particles[ptcl_id+i].isActive)
						continue;
					particles[ptcl_id+i].initializeTimeStep();
				}
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

			case 20: // Primordial binary search
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);

				ptcl = &particles[ptcl_id];
				ptcl->checkNewGroup2();
				break;

			case 21: // CheckBreak
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);

				ptcl = &particles[ptcl_id];
				if (ptcl->GroupInfo) {
					if (!ptcl->GroupInfo->isMerger)
						ptcl->GroupInfo->isTerminate = ptcl->GroupInfo->CheckBreak();
					else if (ptcl->GroupInfo->isMerger && !ptcl->GroupInfo->isTerminate)
						ptcl->GroupInfo->isMerger = false;
				}
				break;

			case 22: // Few-body group search
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);

				ptcl = &particles[ptcl_id];
				if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 > tbin) break;
				ptcl->checkNewGroup();
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

