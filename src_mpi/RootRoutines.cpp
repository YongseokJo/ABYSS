#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cassert>
#include <mpi.h>
#include "global.h"
#include "SkipList.h"
#include "Queue.h"
#include "Worker.h"
#include "QueueScheduler.h"


void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG);
void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int* data, int NumTask, int TAG);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void ParticleSynchronization();
void updateNextRegTime(std::vector<int>& RegularList);
bool createSkipList(SkipList *skiplist);
bool updateSkipList(SkipList *skiplist, int ptcl_id);
int writeParticle(double current_time, int outputNum);
void calculateRegAccelerationOnGPU(std::vector<int> RegularList, QueueScheduler &queue_scheduler);

void formPrimordialBinaries(int beforeNumberOfParticle);
void insertNeighbors(int Order);
void formBinaries(std::vector<int>& ParticleList, std::vector<int>& newCMptcls, std::unordered_map<int, int>& existing, std::unordered_map<int, int>& terminated);

Worker* workers;

void RootRoutines() {

	std::cout << "Root processor is ready." << std::endl;

	Particle* ptcl;
	int min_time_level=0;
	//int worker_rank;
	int task, task_recv, total_tasks;
	int remaining_tasks=0, completed_tasks=0, completed_rank;
	int flag=0;
	// std::map<int, int> CMPtcl_Worker_Map; // commented out by EW 2025.1.5
	// std::map<int, int> Worker_CMPtcl_Map; // commented out by EW 2025.1.5
	// std::vector<int> CMPtclList; // commented out by EW 2025.1.5
	// std::vector<int> CMWorkerList; // commented out by EW 2025.1.5

	std::unordered_map<int, int> CMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11
	std::unordered_map<int, int> PrevCMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11
	std::vector<int> newCMptcls; // by EW 2025.1.6 // unordered_set? by EW 2025.1.11
	std::vector<int> EmptyIndex; // by EW 2025.1.7  empty slots in particles e.g., due to mergers
	// unordered_set? by EW 2025.1.11
	// merged particles & PISN will be contained here
	// new single Particle formed in Enzo can be formed in ParticleIndex of these ptcls
	// if empty, NumberOfParticle++
	


	std::vector<int> RegularList;
	//MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
	//MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
	MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object

	//int sender_rank, sender_tag;
	int ptcl_id;

	workers = new Worker[NumberOfWorker+1];

	for (int i=0; i<=NumberOfWorker; i++) {
		workers[i].initialize(i);
	}

	QueueScheduler queue_scheduler;

	/* Particle loading Check */
	/*
	{
		//, NextRegTime= %.3e Myr(%llu),
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			fprintf(stdout, "PID=%d, pos=(%lf, %lf, %lf), vel=(%lf, %lf, %lf)\n",
					ptcl->PID,
					ptcl->Position[0],
					ptcl->Position[1],
					ptcl->Position[2],
					ptcl->Velocity[0],
					ptcl->Velocity[1],
					ptcl->Velocity[2]
					);
		}
		fflush(stdout);
	}*/

	/* Initialization */
	{

		std::vector<int> PIDs;
		PIDs.reserve(NumberOfParticle);
		PIDs.resize(NumberOfParticle);
		for (int i = 0; i < NumberOfParticle; i++)
		{
			PIDs[i] = i;
		}

		std::cout << "Initialization of particles starts." << std::endl;
		queue_scheduler.initialize(INIT_ACC1);
		queue_scheduler.takeQueue(PIDs);
		do
		{
			//queue_scheduler.printFreeWorker();
			//queue_scheduler.printWorkerToGo();
			queue_scheduler.assignQueueAuto();
			//queue_scheduler.printFreeWorker();
			//queue_scheduler.printWorkerToGo();
			queue_scheduler.runQueueAuto();
			queue_scheduler.waitQueue(0); //blocking wait
		} while(queue_scheduler.isComplete());

		std::cout << "Init 01 done" << std::endl;
		queue_scheduler.initialize(INIT_ACC2);
		queue_scheduler.takeQueue(PIDs);
		do
		{
			queue_scheduler.assignQueueAuto();
			queue_scheduler.runQueueAuto();
			queue_scheduler.waitQueue(0); //blocking wait
		} while(queue_scheduler.isComplete());

		std::cout << "Init 02 done" << std::endl;

#ifdef FEWBODY
		// Primordial binary search

		queue_scheduler.initialize(20);
		queue_scheduler.takeQueue(PIDs);
		do
		{
			queue_scheduler.assignQueueAuto();
			queue_scheduler.runQueueAuto();
			queue_scheduler.waitQueue(0); //blocking wait
		} while(queue_scheduler.isComplete());

		std::cout << "Primordial binary search done" << std::endl;

		// example code by EW 2025.1.7
		Queue queue;
		int rank;
		int OriginalNumberOfParticle = NumberOfParticle;
		formPrimordialBinaries(OriginalNumberOfParticle);
		assert(OriginalNumberOfParticle <= NumberOfParticle); // for debugging by EW 2025.1.4
		assert(CMPtclWorker.empty()); // for debugging by EW 2025.1.4
		if (OriginalNumberOfParticle != NumberOfParticle) {
			std::cout << "In total, " << OriginalNumberOfParticle - NumberOfParticle
					  << " primordial binaries are created." << std::endl;
			queue_scheduler.initialize(23);
			for (int i=OriginalNumberOfParticle; i<NumberOfParticle; i++) {
				std::cout << "New Primordial Binary of PID="
						  << i << " is created with being assigned to a worker of rank "
						  << rank << "." << std::endl;
				ptcl = &particles[i];
				CMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker.size() % NumberOfWorker + 1});
				PIDs.push_back(ptcl->ParticleIndex);
				rank = CMPtclWorker[ptcl->ParticleIndex];

				queue.task = 23;
				queue.pid = ptcl->ParticleIndex;
				workers[rank].addQueue(queue);
				queue_scheduler.WorkersToGo.insert(&workers[rank]);
			}
			queue_scheduler.setTotalQueue(CMPtclWorker.size());
			do {
				queue_scheduler.runQueueAuto();
				queue_scheduler.waitQueue(0);
			} while(queue_scheduler.isComplete());
		}
		else {
			std::cout << "There is no primordial binary." << std::endl;
		}
		fprintf(stdout, "PrimordialBinariesRoutine has ended...\nThe total number of particles is  %d\n", NumberOfParticle);
		fflush(stdout);
#endif

		// Initialize Time Step
		queue_scheduler.initialize(INIT_TIME);
		queue_scheduler.takeQueue(PIDs);
		do
		{
			queue_scheduler.assignQueueAuto();
			queue_scheduler.runQueueAuto();
			queue_scheduler.waitQueue(0); //blocking wait
		} while(queue_scheduler.isComplete());

		/*
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			if (ptcl->isActive)
				fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
								"NumNeighbor= %d\n",
						ptcl->PID,
						ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockIrr,
						ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockReg,
						// NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
						// NextRegTimeBlock,
						ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeBlockIrr,
						ptcl->TimeLevelIrr,
						ptcl->TimeBlockReg,
						ptcl->TimeLevelReg,
						ptcl->NumberOfNeighbor);
		}
		*/
	}



	/* synchronization */
	//ParticleSynchronization();

	/* timestep correction */
	{
		std::cout << "Time Step correction." << std::endl;
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];

			if (!ptcl->isActive)
				continue;

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
		time_step = pow(2,time_block);

		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];

			if (!ptcl->isActive)
				continue;

			ptcl->TimeBlockIrr = static_cast<ULL>(pow(2, ptcl->TimeLevelIrr-time_block));
			ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg-time_block));
#ifdef IRR_TEST
			ptcl->TimeStepReg = 1;
			ptcl->TimeLevelReg = 0;
			ptcl->TimeBlockReg = block_max;
#endif
			ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
		}
		std::cout << "Time Step done." << std::endl;
	}

	/*
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			if (ptcl->isActive)
				fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
								"NumNeighbor= %d\n",
						ptcl->PID,
						ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockIrr,
						ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockReg,
						// NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
						// NextRegTimeBlock,
						ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeBlockIrr,
						ptcl->TimeLevelIrr,
						ptcl->TimeBlockReg,
						ptcl->TimeLevelReg,
						ptcl->NumberOfNeighbor);
		}
		*/

	/* timestep variable synchronization */
	{
		std::cout << "Time Step synchronization." << std::endl;
		task=10;
		completed_tasks = 0; total_tasks = NumberOfWorker;
		InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
		//MPI_Waitall(NumberOfCommunication, requests, statuses);
		//NumberOfCommunication = 0;
		broadcastFromRoot(time_block);
		broadcastFromRoot(block_max);
		broadcastFromRoot(time_step);
		//MPI_Win_sync(win);  // Synchronize memory
		//MPI_Barrier(shared_comm);
		while (completed_tasks < total_tasks) {
			MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			completed_tasks++;
		}
		fprintf(stderr, "nbody+:time_block = %d, EnzoTimeStep=%e\n", time_block, EnzoTimeStep);
		fflush(stderr);
	}

	/* Particle Initialization Check */
	/*
	{
		//, NextRegTime= %.3e Myr(%llu),
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			fprintf(stdout, "%d(%d)=",ptcl->PID,ptcl->NumberOfNeighbor);
			for (int j=0;j<ptcl->NumberOfNeighbor;j++) {
				fprintf(stdout, "%d, ",ptcl->Neighbors[j]);
			}
			fprintf(stdout, "\n");
		}
	}
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"\
					"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"\
					"NumNeighbor= %d\n",
					ptcl->PID,
					ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
					ptcl->CurrentBlockIrr,
					ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
					ptcl->CurrentBlockReg,
					//NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
					//NextRegTimeBlock,
					ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
					ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
					ptcl->TimeBlockIrr,
					ptcl->TimeLevelIrr,
					ptcl->TimeBlockReg,
					ptcl->TimeLevelReg,
					ptcl->NumberOfNeighbor
					);

			fprintf(stdout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
					a1_reg = (%.4e,%.4e,%.4e), a2_reg = (%.4e,%.4e,%.4e), a3_reg = (%.4e,%.4e,%.4e)\n\
					a1_irr = (%.4e,%.4e,%.4e), a2_irr = (%.4e,%.4e,%.4e), a3_irr = (%.4e,%.4e,%.4e)\n", 
					ptcl->a_tot[0][0],
					ptcl->a_tot[1][0],
					ptcl->a_tot[2][0],
					ptcl->a_reg[0][0],
					ptcl->a_reg[1][0],
					ptcl->a_reg[2][0],
					ptcl->a_irr[0][0],
					ptcl->a_irr[1][0],
					ptcl->a_irr[2][0],
					ptcl->NumberOfNeighbor,
					ptcl->RadiusOfNeighbor,
					ptcl->a_reg[0][1],
					ptcl->a_reg[1][1],
					ptcl->a_reg[2][1],
					ptcl->a_reg[0][2],
					ptcl->a_reg[1][2],
					ptcl->a_reg[2][2],
					ptcl->a_reg[0][3],
					ptcl->a_reg[1][3],
					ptcl->a_reg[2][3],
					ptcl->a_irr[0][1],
					ptcl->a_irr[1][1],
					ptcl->a_irr[2][1],
					ptcl->a_irr[0][2],
					ptcl->a_irr[1][2],
					ptcl->a_irr[2][2],
					ptcl->a_irr[0][3],
					ptcl->a_irr[1][3],
					ptcl->a_irr[2][3]
						);

		}
		fflush(stdout);
	}
	*/
	/* Particle Initialization Check */
	// /*
	{
		//, NextRegTime= %.3e Myr(%llu),
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			if (ptcl->isActive)
				fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
								"NumNeighbor= %d\n",
						ptcl->PID,
						ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockIrr,
						ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockReg,
						// NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
						// NextRegTimeBlock,
						ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeBlockIrr,
						ptcl->TimeLevelIrr,
						ptcl->TimeBlockReg,
						ptcl->TimeLevelReg,
						ptcl->NumberOfNeighbor);
		}
	}
	// */


	/* Actual Loop */
	{
		int max_level = 5;
		double prob = 0.5;
		SkipList *skiplist;
		Node* ThisLevelNode;
		double current_time_irr=0;
		double next_time=0;
		Worker* worker;


		//ParticleSynchronization();
		while (1) {
			updateNextRegTime(RegularList);
			/*
			std::cout << "NextRegTimeBlock=" << NextRegTimeBlock << std::endl;
			std::cout << "PID= ";
			for (int i=0; i<RegularList.size(); i++)
				std::cout << RegularList[i]<< ", ";
			std::cout << std::endl;
			std::cout << "size of regularlist= " << RegularList.size() << std::endl;
			*/

			skiplist = new SkipList(max_level, prob);
			if (createSkipList(skiplist) == FAIL)
				fprintf(stderr, "There are no irregular particles!\nBut is it really happening? check skiplist->display()\n");

			bool bin_termination = false;
			bool new_binaries = false;


			// Irregular
			while ( skiplist->getFirstNode() != nullptr) {
				// update_idx=0; // commented out by EW 2025.1.6
				ThisLevelNode = skiplist->getFirstNode();
				next_time     = particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr\
									 	    + particles[ThisLevelNode->ParticleList[0]].TimeStepIrr;

				// print out particlelist
				/*
				fprintf(stdout, "(IRR_FORCE) PID (%d) =", ThisLevelNode->ParticleList.size());
				for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
					ptcl = &particles[ThisLevelNode->ParticleList[i]];
					fprintf(stdout, "%d, ", ptcl->PID);
				}
				fprintf(stdout, "\n");
				fflush(stdout);
				*/

				// Irregular Force
#ifdef FEWBODY
				// std::cout << "Irr force starts" << std::endl;
				int cm_pid;
				Queue queue;
				queue_scheduler.initializeIrr(IRR_FORCE, next_time, ThisLevelNode->ParticleList);
				auto iter = queue_scheduler.CMPtcls.begin();
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					do { 
						worker = queue_scheduler.waitQueue(1); // non-blocking wait
						// if there's any CMPtcl
						if (queue_scheduler.CMPtcls.size() > 0)
						{
							/* check if there's any CM ptcl ready to go for SDAR*/
							if (iter == queue_scheduler.CMPtcls.end())
								iter = queue_scheduler.CMPtcls.begin();
							cm_pid = *(iter);
							ptcl = &particles[cm_pid];
							for (int j = 0; j < ptcl->NumberOfNeighbor; j++)
							{
								if (particles[ptcl->Neighbors[j]].isUpdateToDate == false)
								{
									iter++;
									goto skip_to_next;
								}
							}
							queue.task = FB_SDAR;
							queue.pid = cm_pid;
							queue.next_time = next_time;
							workers[CMPtclWorker[cm_pid]].addQueue(queue);
							queue_scheduler.assignWorker(&workers[CMPtclWorker[cm_pid]]);
							iter = queue_scheduler.CMPtcls.erase(iter);
						skip_to_next:;
						}
						if (worker != nullptr) 
						{
							//fprintf(stdout, "Worker rank: %d\n", worker->MyRank);
						}
					} while (worker == nullptr);
					queue_scheduler.callback(worker);
				} while (queue_scheduler.isComplete() || queue_scheduler.CMPtcls.size() > 0);

				// std::cout << "Irregular Force done" << std::endl;
#else
				queue_scheduler.initialize(IRR_FORCE, next_time);
				queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());
#endif

				//ParticleSynchronization();


				// if remaining, further sorting
				/*
				if (update_idx < total_tasks) {
						updateSkipList(skiplist, ThisLevelNode->ParticleList[update_idx]);
						update_idx++;
					}				updateSkipList(skiplist, ptcl_id_return);
					*/

				// Irregular Update
				queue_scheduler.initialize(IRR_UPDATE);
				queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());

#ifdef FEWBODY
				int OriginalSize = ThisLevelNode->ParticleList.size();
				for (int i=0; i<OriginalSize; i++ ){
					ptcl = &particles[ThisLevelNode->ParticleList[i]];
					if (ptcl->isCMptcl && !ptcl->isActive) {

						PrevCMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker[ptcl->ParticleIndex]});
						CMPtclWorker.erase(ptcl->ParticleIndex);

						insertNeighbors(ptcl->ParticleIndex);

						for (int j=0; j < ptcl->NewNumberOfNeighbor; j++) {
							ThisLevelNode->ParticleList.push_back(ptcl->NewNeighbors[j]);
							particles[ptcl->NewNeighbors[j]].isActive = true;
						}

						bin_termination = true;
					}
				}

				// Erase terminated CM particles by EW 2025.1.6
				ThisLevelNode->ParticleList.erase(
					std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
						[](int i) {
						return !particles[i].isActive;
						}
					),
					ThisLevelNode->ParticleList.end()
				);

				// Few-body group search
				queue_scheduler.initialize(22);
				queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());

				int OriginalParticleListSize = ThisLevelNode->ParticleList.size();
				int rank_delete, rank_new;
				formBinaries(ThisLevelNode->ParticleList, newCMptcls, CMPtclWorker, PrevCMPtclWorker);
				if (OriginalParticleListSize != ThisLevelNode->ParticleList.size()) {
					new_binaries = true;
					for (int i=OriginalParticleListSize; i<ThisLevelNode->ParticleList.size(); i++) {
						ptcl = &particles[ThisLevelNode->ParticleList[i]];
						fprintf(stdout, "New CM Particle PID: %d\n", ptcl->PID);
						fflush(stdout);
					}
					//work_scheduler.doSimpleTask(24, newCMptcls); 
					// This should be fixed by EW 2025.1.6
					// First, delete CM particle inside NewNeighbors using void deleteGroup(Particle* ptclCM) // task 25?
					// Should adjust CMPtclWorker & PrevCMPtclWorker
					// Then, use task 24

					// example code by EW 2025.1.7
					Particle* ptclCM;
					Particle* mem_ptclCM;
					int total_queues = 0;
					for (int i=0; i<newCMptcls.size(); i++) {
						ptclCM = &particles[newCMptcls[i]]; // 2025.01.10 edited to newCMptcls[i] by YS
						queue_scheduler.initialize(25);
						total_queues = 0;
						for (int j=0; j<ptclCM->NewNumberOfNeighbor; j++) {
							mem_ptclCM = &particles[ptclCM->NewNeighbors[j]];
							if (mem_ptclCM->isCMptcl) {
								rank_delete = CMPtclWorker[mem_ptclCM->ParticleIndex];
								queue.task = 25;
								queue.pid = mem_ptclCM->ParticleIndex;
								workers[rank_delete].addQueue(queue);
								queue_scheduler.WorkersToGo.insert(&workers[rank_delete]);
								total_queues++;
							}
						}
						if (total_queues > 0)
						{
							queue_scheduler.setTotalQueue(total_queues);
							do
							{
								queue_scheduler.runQueueAuto();
								queue_scheduler.waitQueue(0);
							} while (queue_scheduler.isComplete());
						}

						queue_scheduler.initialize(24);
						rank_new = CMPtclWorker[ptclCM->ParticleIndex];
						fprintf(stdout, "New CM ptcl is %d\n", newCMptcls[0]);
						fprintf(stdout, "Rank of CM ptcl %d: %d\n", ptclCM->ParticleIndex, rank_new);
						queue.task = 24;
						queue.pid = ptclCM->ParticleIndex;
						workers[rank_new].addQueue(queue);
						workers[rank_new].runQueue();
						workers[rank_new].callback();
					}
					std::cout << "All new fewbody objects are initialized." << std::endl;

					ThisLevelNode->ParticleList.erase(
						std::remove_if(
								ThisLevelNode->ParticleList.begin(), 
								ThisLevelNode->ParticleList.end(),
								[](int i) { return !particles[i].isActive; }
						),
						ThisLevelNode->ParticleList.end()
					);
				}
				newCMptcls.clear();

				// std::cout << "erase success" << std::endl;
#endif

				for (int i=0; i<ThisLevelNode->ParticleList.size(); i++)
					updateSkipList(skiplist, ThisLevelNode->ParticleList[i]);

				//std::cout << "update success" << std::endl;
				//skiplist->display();
				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
						ptcl = &particles_original[ThisLevelNode->ParticleList[i]];
						fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
								"NumNeighbor= %d\n",
								ptcl->PID,
								ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockIrr,
								ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockReg,
								NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
								NextRegTimeBlock,
								ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
								ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
								ptcl->TimeBlockIrr,
								ptcl->TimeLevelIrr,
								ptcl->TimeBlockReg,
								ptcl->TimeLevelReg,
								ptcl->NextBlockIrr*time_step*EnzoTimeStep*1e10/1e6,
								ptcl->NextBlockIrr,
								ptcl->NumberOfNeighbor
								);

						fprintf(stdout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
								a1_reg = (%.4e,%.4e,%.4e), a2_reg = (%.4e,%.4e,%.4e), a3_reg = (%.4e,%.4e,%.4e)\n\
								a1_irr = (%.4e,%.4e,%.4e), a2_irr = (%.4e,%.4e,%.4e), a3_irr = (%.4e,%.4e,%.4e)\n", 
								ptcl->a_tot[0][0],
								ptcl->a_tot[1][0],
								ptcl->a_tot[2][0],
								ptcl->a_reg[0][0],
								ptcl->a_reg[1][0],
								ptcl->a_reg[2][0],
								ptcl->a_irr[0][0],
								ptcl->a_irr[1][0],
								ptcl->a_irr[2][0],
								ptcl->NumberOfNeighbor,
								ptcl->RadiusOfNeighbor,
								ptcl->a_reg[0][1],
								ptcl->a_reg[1][1],
								ptcl->a_reg[2][1],
								ptcl->a_reg[0][2],
								ptcl->a_reg[1][2],
								ptcl->a_reg[2][2],
								ptcl->a_reg[0][3],
								ptcl->a_reg[1][3],
								ptcl->a_reg[2][3],
								ptcl->a_irr[0][1],
								ptcl->a_irr[1][1],
								ptcl->a_irr[2][1],
								ptcl->a_irr[0][2],
								ptcl->a_irr[1][2],
								ptcl->a_irr[2][2],
								ptcl->a_irr[0][3],
								ptcl->a_irr[1][3],
								ptcl->a_irr[2][3]
									);
					}
					fflush(stdout);
				}
				*/
				current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr*time_step;
				skiplist->deleteFirstNode();

				//std::cout << "deleteFirstNode success" << std::endl;
				//ThisLevelNode = skiplist->getFirstNode();

				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
						ptcl = &particles[ThisLevelNode->ParticleList[i]];
						fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"\
								"NumNeighbor= %d\n",
								ptcl->PID,
								ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockIrr,
								ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockReg,
								NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
								NextRegTimeBlock,
								ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
								ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
								ptcl->TimeBlockIrr,
								ptcl->TimeLevelIrr,
								ptcl->TimeBlockReg,
								ptcl->TimeLevelReg,
								ptcl->NumberOfNeighbor
								);
					}
				}
				*/
				//exit(SUCCESS);
#ifdef IRR_TEST
				std::cout << "current_time_irr=" << current_time_irr<< std::endl;
				// create output at appropriate time intervals
				if (current_time_irr >= outputTime) {
					writeParticle(current_time_irr, outNum++);
					outputTime += outputTimeStep;
				}

				// end if the global time exceeds the end time
				if (current_time_irr >= 1) {
					task=-100;
					InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
					MPI_Waitall(NumberOfCommunication, requests, statuses);
					NumberOfCommunication = 0;
					std::cout << EnzoTimeStep << std::endl;
					std::cout << "Simulation Done!" << std::endl;
					return;
				}
#endif
			} // Irr
			delete skiplist;
			skiplist = nullptr;
			//exit(SUCCESS);

#ifdef FEWBODY
			if (bin_termination || new_binaries)
				updateNextRegTime(RegularList);
			// std::cout << "updateNextRegTime done" << std::endl;
#endif

#ifdef CUDA
			{
				//total_tasks = RegularList.size();
				next_time = NextRegTimeBlock*time_step;

				//fprintf(stdout, "Regular starts\n");
				calculateRegAccelerationOnGPU(RegularList, queue_scheduler);

				// Update Regular
				queue_scheduler.initialize(REG_CUDA_UPDATE);
				queue_scheduler.takeQueue(RegularList);
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());

				//fprintf(stdout, "Regular done\n");
#ifdef unused

				// Few-body group search
				//work_scheduler.initialize();
				//work_scheduler.doSimpleTask(22, RegularList);

				int beforeRegularListSize = RegularList.size();

				formBinaries(RegularList, newCMptcls, CMPtclWorker, PrevCMPtclWorker);
				if (beforeRegularListSize != RegularList.size()) {
					for (int i=beforeRegularListSize; i<RegularList.size(); i++) {
						ptcl = &particles[RegularList[i]];
						fprintf(stdout, "New CM Particle PID: %d\n", ptcl->PID);
						fflush(stdout);
					}
					//work_scheduler.doSimpleTask(24, newCMptcls); 
					// This should be fixed by EW 2025.1.6
					// First, delete CM particle inside NewNeighbors using void deleteGroup(Particle* ptclCM) // task 25?
					// Should adjust CMPtclWorker & PrevCMPtclWorker
					// Then, use task 24

					// example code by EW 2025.1.7
					Particle* ptclCM;
					Particle* mem_ptclCM;
					for (int i=0; i<newCMptcls.size(); i++) {
						ptclCM = &particles[i];
						for (int j=0; j<ptclCM->NewNumberOfNeighbor; j++) {
							mem_ptclCM = &particles[ptclCM->NewNeighbors[j]];
							if (mem_ptclCM->isCMptcl) {
								int rank1 = CMPtclWorker[mem_ptclCM->ParticleIndex];
								// Send and do task 25 in worker number: rank1
							}
						}
						int rank2 = CMPtclWorker[ptclCM->ParticleIndex]; 
						// Send and do task 24 in worker number: rank2
					}

					RegularList.erase(
						std::remove_if(
								RegularList.begin(), 
								RegularList.end(),
								[](int i) { return !particles[i].isActive; }
						),
						RegularList.end()
					);
				}
				newCMptcls.clear();

				std::cout << "erase success" << std::endl;
#endif
			}
				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<RegularList.size(); i++) {
						ptcl = &particles[RegularList[i]];
						fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
								"NumNeighbor= %d\n",
								ptcl->PID,
								ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockIrr,
								ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockReg,
								NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
								NextRegTimeBlock,
								ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
								ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
								ptcl->TimeBlockIrr,
								ptcl->TimeLevelIrr,
								ptcl->TimeBlockReg,
								ptcl->TimeLevelReg,
								ptcl->NextBlockIrr*time_step*EnzoTimeStep*1e10/1e6,
								ptcl->NextBlockIrr,
								ptcl->NumberOfNeighbor
								);

						fprintf(stdout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
								a1_reg = (%.4e,%.4e,%.4e), a2_reg = (%.4e,%.4e,%.4e), a3_reg = (%.4e,%.4e,%.4e)\n\
								a1_irr = (%.4e,%.4e,%.4e), a2_irr = (%.4e,%.4e,%.4e), a3_irr = (%.4e,%.4e,%.4e)\n", 
								ptcl->a_tot[0][0],
								ptcl->a_tot[1][0],
								ptcl->a_tot[2][0],
								ptcl->a_reg[0][0],
								ptcl->a_reg[1][0],
								ptcl->a_reg[2][0],
								ptcl->a_irr[0][0],
								ptcl->a_irr[1][0],
								ptcl->a_irr[2][0],
								ptcl->NumberOfNeighbor,
								ptcl->RadiusOfNeighbor,
								ptcl->a_reg[0][1],
								ptcl->a_reg[1][1],
								ptcl->a_reg[2][1],
								ptcl->a_reg[0][2],
								ptcl->a_reg[1][2],
								ptcl->a_reg[2][2],
								ptcl->a_reg[0][3],
								ptcl->a_reg[1][3],
								ptcl->a_reg[2][3],
								ptcl->a_irr[0][1],
								ptcl->a_irr[1][1],
								ptcl->a_irr[2][1],
								ptcl->a_irr[0][2],
								ptcl->a_irr[1][2],
								ptcl->a_irr[2][2],
								ptcl->a_irr[0][3],
								ptcl->a_irr[1][3],
								ptcl->a_irr[2][3]
									);
					}
					//fflush(stdout); 
				}
				*/
#else
			//std::cout << "Regular Routine Starts." << std::endl;
			// Regular
			{
				// Regular Gravity
				task = 1;
				completed_tasks = 0;
				total_tasks = RegularList.size();
				next_time = NextRegTimeBlock*time_step;

				//std::cout << "TotalTask=" << total_tasks << std::endl;

				/*
				std::cout << "RegularList, PID= ";
				for (int i=0; i<total_tasks; i++) {
					std::cout << RegularList[i]<< ", ";
				}*/
				//std::cout << std::endl;

				InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
				InitialAssignmentOfTasks(RegularList, next_time, total_tasks, PTCL_TAG);
				MPI_Waitall(NumberOfCommunication, requests, statuses);
				NumberOfCommunication = 0;


				// further assignments
				remaining_tasks = total_tasks-NumberOfWorker;
				while (completed_tasks < total_tasks) {
					// Check which worker is done
					MPI_Irecv(&ptcl_id_return, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
					// Poll until send completes
					/*
						 flag=0;
						 while (!flag) {
						 MPI_Test(&request, &flag, &status);
					// Perform other work while waiting
					}
					*/
					MPI_Wait(&request, &status);
					completed_rank = status.MPI_SOURCE;
					//printf("Rank %d: Send operation completed (%d).\n",completed_rank, ptcl_id_return);

					if (remaining_tasks > 0) {
						ptcl_id = RegularList[NumberOfWorker + completed_tasks];
						MPI_Send(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD);
						MPI_Send(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD);
						MPI_Send(&next_time, 1, MPI_DOUBLE, completed_rank, TIME_TAG, MPI_COMM_WORLD);
						remaining_tasks--;
					} else {
						//printf("Rank %d: No more tasks to assign\n", completed_rank);
					}
					//updateSkipList(skiplist, ptcl_id_return);
					completed_tasks++;
				}


				//ParticleSynchronization();


				// Regular Update
				//std::cout<< "Reg Acc Done." <<std::endl;
				task = 3;
				completed_tasks = 0;

				InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
				InitialAssignmentOfTasks(RegularList, total_tasks, PTCL_TAG);
				MPI_Waitall(NumberOfCommunication, requests, statuses);
				NumberOfCommunication = 0;


				// further assignments
				remaining_tasks = total_tasks-NumberOfWorker;
				while (completed_tasks < total_tasks) {
					// Check which worker is done
					MPI_Irecv(&ptcl_id_return, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
					MPI_Wait(&request, &status);
					completed_rank = status.MPI_SOURCE;

					if (remaining_tasks > 0) {
						ptcl_id = RegularList[NumberOfWorker + completed_tasks];
						//MPI_Isend(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD, &request);
						//MPI_Isend(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD, &request);
						MPI_Send(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD);
						MPI_Send(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD);
						remaining_tasks--;
					} else {
						//printf("Rank %d: No more tasks to assign\n", completed_rank);
					}
					completed_tasks++;
				}

				//ParticleSynchronization();
				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<total_tasks; i++) {
						ptcl = &particles[RegularList[i]];
						fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
								"NumNeighbor= %d\n",
								ptcl->PID,
								ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockIrr,
								ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockReg,
								NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
								NextRegTimeBlock,
								ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
								ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
								ptcl->TimeBlockIrr,
								ptcl->TimeLevelIrr,
								ptcl->TimeBlockReg,
								ptcl->TimeLevelReg,
								ptcl->NextBlockIrr*time_step*EnzoTimeStep*1e10/1e6,
								ptcl->NextBlockIrr,
								ptcl->NumberOfNeighbor
								);

						fprintf(stdout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
								a1_reg = (%.4e,%.4e,%.4e), a2_reg = (%.4e,%.4e,%.4e), a3_reg = (%.4e,%.4e,%.4e)\n\
								a1_irr = (%.4e,%.4e,%.4e), a2_irr = (%.4e,%.4e,%.4e), a3_irr = (%.4e,%.4e,%.4e)\n", 
								ptcl->a_tot[0][0],
								ptcl->a_tot[1][0],
								ptcl->a_tot[2][0],
								ptcl->a_reg[0][0],
								ptcl->a_reg[1][0],
								ptcl->a_reg[2][0],
								ptcl->a_irr[0][0],
								ptcl->a_irr[1][0],
								ptcl->a_irr[2][0],
								ptcl->NewNumberOfNeighbor,
								ptcl->RadiusOfNeighbor,
								ptcl->a_reg[0][1],
								ptcl->a_reg[1][1],
								ptcl->a_reg[2][1],
								ptcl->a_reg[0][2],
								ptcl->a_reg[1][2],
								ptcl->a_reg[2][2],
								ptcl->a_reg[0][3],
								ptcl->a_reg[1][3],
								ptcl->a_reg[2][3],
								ptcl->a_irr[0][1],
								ptcl->a_irr[1][1],
								ptcl->a_irr[2][1],
								ptcl->a_irr[0][2],
								ptcl->a_irr[1][2],
								ptcl->a_irr[2][2],
								ptcl->a_irr[0][3],
								ptcl->a_irr[1][3],
								ptcl->a_irr[2][3]
									);
					}
					//fflush(stdout); 
				}
				*/
				//current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr;
			} // Regular Done.
#endif
			global_time = NextRegTimeBlock*time_step;

			// create output at appropriate time intervals
			if (global_time >= outputTime) {
				writeParticle(global_time, outNum++);
				outputTime += outputTimeStep;
			}

			// end if the global time exceeds the end time
			if (global_time >= 1) {
				task=-100;
				InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
				//MPI_Waitall(NumberOfCommunication, requests, statuses);
				//NumberOfCommunication = 0;
				std::cout << EnzoTimeStep << std::endl;
				std::cout << "Simulation Done!" << std::endl;
				return;
			}
			//exit(SUCCESS);
		} // While(1)
	} // Actual Loop
}




void updateNextRegTime(std::vector<int>& RegularList) {

	ULL time_tmp=0, time=block_max;
	Particle *ptcl;

	RegularList.clear();

	for (int i=0; i<NumberOfParticle; i++)
	{
		//std::cout << i << std::endl;
		ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;
		// Next regular time step
		time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;

		// Find the minum regular time step
		if (time_tmp <= time) {
			//fprintf(stderr, "PID=%d, time_tme=%llu\n", ptcl->PID, time_tmp);
			if (time_tmp < time) {
				RegularList.clear();
				time = time_tmp;
			}
			RegularList.push_back(ptcl->ParticleIndex);
		}
	}
	NextRegTimeBlock = time;
}







bool createSkipList(SkipList *skiplist) {

	bool debug = false;
	//fprintf(stdout, "create level starts!\n");
	//fflush(stdout);

	if (debug) {
		fprintf(stdout, "create level starts!\n");
		fflush(stdout);
	}

	/*
#ifdef time_trace
	_time.irr_chain.markStart();
#endif
*/

	Particle* ptcl;

	for (int i=0; i<NumberOfParticle; i++) {
		ptcl =  &particles[i];

		// if ((ptcl->NumberOfNeighbor != 0) && (ptcl->NextBlockIrr <= NextRegTimeBlock)) { // IAR original
		if (ptcl->isActive && ptcl->NextBlockIrr <= NextRegTimeBlock) {	// IAR modified
			//fprintf(stdout, "PID=%d, NBI=%llu\n", ptcl->PID, ptcl->NextBlockIrr);
			if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleIndex))
				skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleIndex);
		}
	}

	/*
#ifdef time_trace
	_time.irr_chain.markEnd();
	_time.irr_chain.getDuration();
#endif
*/

	/*
	if (debug) {
		skiplist->display();
	}
	*/

	//fprintf(stdout, "create level ends!\n");
	//fflush(stdout);

	if (skiplist->getFirstNode() == nullptr)
		return FAIL;
	else
		return SUCCESS;
}



bool updateSkipList(SkipList *skiplist, int ptcl_id) {
	bool debug = false;
	if (debug) {
		fprintf(stdout, "update level starts!\n");
		fflush(stdout);
	}

	/*
#ifdef time_trace
	_time.irr_sort.markStart();
#endif
*/

	/* Update New Time Steps */
	//Node* ThisLevelNode = skiplist->getFirstNode();

	/*
		 if (this->debug) {
		 fprintf(stdout, "PID=%d, NBI=%llu, size=%lu\n", ptcl->PID, ptcl->NextBlockIrr, ThisLevelNode->particle_list.size());
		 fprintf(stdout, "NextBlockIrr=%llu\n",ptcl->NextBlockIrr);
		 fflush(stdout);
		 }
		 */

	Particle * ptcl = &particles[ptcl_id];

	//std::cout << "NextBlockIrr of "<< ptcl_id<<" = " << ptcl->NextBlockIrr << std::endl;
	if (ptcl->NextBlockIrr > NextRegTimeBlock)
		return true;

	if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleIndex))
		skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleIndex);

	if (debug) {
	}

	/*
#ifdef time_trace
	_time.irr_sort.markEnd();
	_time.irr_sort.getDuration();
#endif
*/

	if (debug) {
		//skiplist->display();
		//fprintf(stderr, "This is it.\n\n\n");
	}

	return true;
}




