#include <iostream>
#include <vector>
#include <map>
#include <mpi.h>
#include "global.h"
#include "SkipList.h"
#include "WorkScheduler.h"


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
void calculateRegAccelerationOnGPU(std::vector<int> RegularList, WorkScheduler &ws);

void SetPrimordialBinaries();
void FBTermination(int Order);
void FBTermination2(int Order);
void SetBinaries(std::vector<int>& ParticleList);

void formPrimordialBinaries(int& beforeNumberOfParticle);

Worker* workers;

void RootRoutines() {

	std::cout << "Root processor is ready." << std::endl;

	Particle* ptcl;
	int min_time_level=0;
	//int worker_rank;
	int task, task_recv, total_tasks;
	int remaining_tasks=0, completed_tasks=0, completed_rank;
	int flag=0;
	std::map<int, int> CMPtcl_Worker_Map;
	std::map<int, int> Worker_CMPtcl_Map;
	std::vector<int> CMPtclList;
	std::vector<int> CMWorkerList;

	std::vector<int> existingCMPtclIndexList; // by EW 2025.1.4
	std::vector<int> terminatedCMPtclIndexList; // by EW 2025.1.4
	

	std::vector<int> RegularList;
	//MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
	//MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
	MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object

	int sender_rank, sender_tag;
	int ptcl_id;

	workers = new Worker[NumberOfWorker];

	for (int i=0; i<NumberOfWorker; i++) {
		workers[i].MyRank = i+1;
	}

	WorkScheduler work_scheduler;

	std::vector<int> PIDs;
	PIDs.reserve(NumberOfParticle);
	PIDs.resize(NumberOfParticle);
	for (int i = 0; i < NumberOfParticle; i++)
	{
		PIDs[i] = i;
	}

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
		std::cout << "Initialization of particles starts." << std::endl;
		work_scheduler.initialize();
		work_scheduler.doSimpleTask(7, PIDs);

		work_scheduler.initialize();
		work_scheduler.doSimpleTask(8, PIDs);

#ifdef FEWBODY
		// Primordial binary search
		std::vector<int> BinaryList;

		work_scheduler.initialize();
		work_scheduler.doSimpleTask(20, NumberOfParticle);

		// Example codes by EW 2025.1.4
		int beforeNumberOfParticle = NumberOfParticle;
		formPrimordialBinaries(beforeNumberOfParticle);
		assert(beforeNumberOfParticle >= NumberOfParticle); // for debugging by EW 2025.1.4
		assert(existingCMPtclIndexList.empty()); // for debugging by EW 2025.1.4
		if (beforeNumberOfParticle != NumberOfParticle) {
			Particle* ptcl;
			for (int i=beforeNumberOfParticle; i<NumberOfParticle; i++) {
				// i is ParticleOrder
				ptcl = &particles[i];
				existingCMPtclIndexList.push_back(i); // this contains ParticleOrders by EW 2025.1.4
				CMPtclList.push_back(ptcl->PID); // this contains PIDs by EW 2025.1.4
				CMWorkerList.push_back(rank); // this contains ranks by EW 2025.1.4
				CMPtcl_Worker_Map.insert(ptcl->PID, rank); // this maps between PIDs and ranks by EW 2025.1.4
			}

			work_scheduler.doSimpleTask(23, NumberOfParticle - beforeNumberOfParticle);
			// rank and CMPtcl might not be matched correctly inside work_scheduler by this simple code by EW 2025.1.4
		}
		// Example codes by EW 2025.1.4

		if (BinaryList.size() >	0) {
			CMPtclList.push_back(PID);
			CMWorkerList.push_back(rank);
			CMPtcl_Worker_Map.insert(PID, rank);
			SetPrimordialBinaries();

		}
		fprintf(stdout, "SetPrimordialBinaries ends... NumberOfParticle: %d\n", NumberOfParticle);
		fflush(stdout);
#endif

		// Initialize Time Step
		work_scheduler.initialize();
		work_scheduler.doSimpleTask(9, PIDs);
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

			if (!ptcl->isActive) // Eunwoo: If I'm correct, there should be no isActive=false particle here!
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

		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			if (ptcl->isActive) // Eunwoo: If I'm correct, there should be no isActive=false particle here!
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
	{
		//, NextRegTime= %.3e Myr(%llu),
		for (int i=0; i<NumberOfParticle; i++) {
			ptcl = &particles[i];
			if (ptcl->isActive) // Eunwoo: If I'm correct, there should be no isActive=false particle here!
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



	/* Actual Loop */
	{
		int max_level = 5;
		double prob = 0.5;
		SkipList *skiplist;
		int update_idx;
		Node* ThisLevelNode;
		flag = 0;
		double current_time_irr=0;
		double next_time=0;
		int ptcl_id_return;
		int AdaptiveLoadBalancing;
		int size=0;

		//ParticleSynchronization();
		while (1) {
			updateNextRegTime(RegularList);
			std::cout << "NextRegTimeBlock=" << NextRegTimeBlock << std::endl;
			std::cout << "PID= ";
			for (int i=0; i<RegularList.size(); i++)
				std::cout << RegularList[i]<< ", ";
			std::cout << std::endl;
			std::cout << "size of regularlist= " << RegularList.size() << std::endl;

			skiplist = new SkipList(max_level, prob);
			if (createSkipList(skiplist) == FAIL)
				fprintf(stderr, "There are no irregular particles!\nBut is it really happening? check skiplist->display()\n");

			bool bin_termination = false;
			bool new_binaries = false;


			// Irregular
			while ( skiplist->getFirstNode() != nullptr) {
				update_idx=0;
				ThisLevelNode = skiplist->getFirstNode();
				next_time     = particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr\
									 	    + particles[ThisLevelNode->ParticleList[0]].TimeStepIrr;

				// Irregular Force
				work_scheduler.initialize();
				work_scheduler.doIrregularForce(ThisLevelNode->ParticleList, next_time);

				//ParticleSynchronization();


				// if remaining, further sorting
				/*
				if (update_idx < total_tasks) {
						updateSkipList(skiplist, ThisLevelNode->ParticleList[update_idx]);
						update_idx++;
					}				updateSkipList(skiplist, ptcl_id_return);
					*/

				// Irregular Update
				work_scheduler.initialize();
				work_scheduler.doSimpleTask(2, ThisLevelNode->ParticleList);


#ifdef FEWBODY
				// Check few-body groups should be broken or not
				work_scheduler.initialize();
				work_scheduler.doSimpleTask(21, ThisLevelNode->ParticleList);

                // Check fewbody termination
				int OriginalSize = ThisLevelNode->ParticleList.size();
				for (int i=0; i<OriginalSize; i++) {
					Particle* ptclCM = &particles[ThisLevelNode->ParticleList[i]];
					if (ptclCM->GroupOrder < 0) continue;

					if (groups[ptclCM->GroupOrder].isTerminate) {
						if (groups[ptclCM->GroupOrder].isMerger) {
							for (int i=0; i<groups[ptclCM->GroupOrder].sym_int.particles.getSize(); i++) {
								Particle* ptcl = &particles[groups[ptclCM->GroupOrder].sym_int.particles[i].ParticleOrder];
								if (ptcl->Mass != 0.0)
									ThisLevelNode->ParticleList[i] = ptcl->ParticleOrder;
							}
							FBTermination2(ThisLevelNode->ParticleList[i]);
							bin_termination = true;
						}
						else {
							if (groups[ptclCM->GroupOrder].sym_int.particles.getSize() > 2) { // Terminate many-body (> 2) group

								ThisLevelNode->ParticleList[i] = groups[ptclCM->GroupOrder].sym_int.particles[0].ParticleOrder; // CM particle -> 1st group member
								for (int i=1; i<groups[ptclCM->GroupOrder].sym_int.particles.getSize(); i++) {
									Particle* ptcl = &particles[groups[ptclCM->GroupOrder].sym_int.particles[i].ParticleOrder]; 
									ThisLevelNode->ParticleList.push_back(ptcl->ParticleOrder); // push_back other group members
								}
								/* Eunwoo: This might be needed later!
								std::vector<Particle*> members = ptcl->GroupInfo->Members;
								FBTermination(ThisLevelNode->ParticleList[i]);
								AddNewGroupsToList2(members, particle);
								members.clear();
								*/
								FBTermination(ThisLevelNode->ParticleList[i]);
								bin_termination = true;
							}
							else { // Terminate binary

								ThisLevelNode->ParticleList[i] = groups[ptclCM->GroupOrder].sym_int.particles[0].ParticleOrder; // CM particle -> 1st binary member
								ThisLevelNode->ParticleList.push_back(groups[ptclCM->GroupOrder].sym_int.particles[1].ParticleOrder); // push_back 2nd binary member

								FBTermination(ThisLevelNode->ParticleList[i]);
								bin_termination = true;
							}
						}
					}
				}
				if (bin_termination)
					global_variable->NumberOfParticle = NumberOfParticle;

				// Few-body group search
				work_scheduler.initialize();
				work_scheduler.doSimpleTask(22, ThisLevelNode->ParticleList);

				int beforeNumberOfParticle = NumberOfParticle;
				SetBinaries(ThisLevelNode->ParticleList);
				if (beforeNumberOfParticle != NumberOfParticle) {
					new_binaries = true;
					for (int i=beforeNumberOfParticle; i<NumberOfParticle; i++) {
						Particle* ptcl = &particles[i];
						ThisLevelNode->ParticleList.push_back(ptcl->ParticleOrder);
						fprintf(stdout, "PID(%d) is pushed back to the ParticleList\n", ptcl->PID);
						fflush(stdout);
					}
				}

				ThisLevelNode->ParticleList.erase(
					std::remove_if(
							ThisLevelNode->ParticleList.begin(), 
							ThisLevelNode->ParticleList.end(),
							[&particles](int i) { return !particles[i].isActive; }
					),
					ThisLevelNode->ParticleList.end()
				);

				std::cout << "erase success" << std::endl;
#endif

				for (int i=0; i<ThisLevelNode->ParticleList.size(); i++)
					updateSkipList(skiplist, ThisLevelNode->ParticleList[i]);

				std::cout << "update success" << std::endl;
				//skiplist->display();
				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<total_tasks; i++) {
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

				std::cout << "deleteFirstNode success" << std::endl;
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
			std::cout << "updateNextRegTime done" << std::endl;
#endif

#ifdef CUDA
			{
				//total_tasks = RegularList.size();
				next_time = NextRegTimeBlock*time_step;

				//fprintf(stdout, "Regular starts\n");
				calculateRegAccelerationOnGPU(RegularList, work_scheduler);

				// Update Regular
				work_scheduler.initialize();
				work_scheduler.doSimpleTask(5, RegularList);
				//fprintf(stdout, "Regular done\n");
#ifdef FEWBODY

				fprintf(stdout, "Regular list: ");
				fflush(stdout);
				Particle* ptcl;
				for (int i=0; i<RegularList.size(); i++) {
					ptcl = &particles[RegularList[i]];
					fprintf(stdout, "%d ", ptcl->PID);
				}
				fprintf(stdout, "\n");

				fprintf(stdout, "Regular group search\n");
				fflush(stdout);
				// Few-body group search
				work_scheduler.initialize();
				work_scheduler.doSimpleTask(22, RegularList);

				SetBinaries(RegularList);
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
			RegularList.push_back(ptcl->ParticleOrder);
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
			if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleOrder))
				skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleOrder);
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

	if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleOrder))
		skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleOrder);

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




