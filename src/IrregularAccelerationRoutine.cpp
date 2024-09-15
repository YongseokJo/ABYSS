#include <iostream>
#include "global.h"
#include <cassert>

#define no_chain_debug

std::vector<Particle*> ComputationList{};
int writeParticle(std::vector<Particle*> &particle, REAL MinRegTime, int outputNum);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool UpdateComputationChain(Particle* ptcl);
//bool UpdateComputationChainParallel(std::vector<Particle*> &particle);
bool CreateComputationListParallel(std::vector<Particle*> &particle);
bool CreateComputationList(Particle* ptcl);
bool AddNewGroupsToList(std::vector<Particle*> &particle);
void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle);



Particle *FirstComputation;
// 1. sort particles according to next irregular time
// 2. the sorted list should be update regularly
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle)
{

	//fprintf(binout, "Starting irregular force\n");
	//std::cout << "Creating a computation chain ...\n";

	bool first = true;
	Particle *ParticleForComputation;
	//AddNewBinariesToList(particle, particle);
	//std::cout << "Create Chain\n" << std::flush;
#ifdef time_trace
	_time.irr_chain.markStart();
#endif

	if (CreateComputationChain(particle) == false) {
		std::cout << "No irregular particle to update ...\n";
		return true;
	}
#ifdef time_trace
	_time.irr_chain.markEnd();
	_time.irr_chain.getDuration();
#endif



	//std::cout << "Calculating irregular force ...\n" << std::endl;

	// Caculating irregular acceleration


	while (CreateComputationList(FirstComputation) && ComputationList.size() != 0) {
			
		//std::cout << "List size=" << ComputationList.size() << std::endl;

#define binary
#ifdef binary
#ifdef time_trace
	_time.irr_bin.markStart();
#endif
		//if (AddNewBinariesToList(ComputationList, particle) && ComputationList.size() == 0) {
		if (AddNewGroupsToList(particle) && ComputationList.size() == 0) {
			//std::cout << "No irregular particle to update afte binary formation." << std::endl;
			fprintf(stdout, "No irregular particle to update after binary formation.\n");
			break;
		}




		//fflush(stdout);
		for (Particle* ptcl:particle) {
				if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
						|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
					fprintf(stdout, "before, myself = %d\n", ptcl->PID);
					fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
					fflush(stdout);
					//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}

		for (Particle* ptcl:particle) {
			if (ptcl->CurrentTimeIrr > 1) 
				fprintf(stderr, "outside before bin, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
		}

		for (Particle* ptcl:particle) {
			if (ptcl->CurrentTimeIrr > 1)  {
				fprintf(stderr, "outside after bin, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
				fflush(stderr);
			}
		}

		for (Particle* ptcl:particle) {
				if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
						|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
					fprintf(stdout, "after, myself = %d\n", ptcl->PID);
					fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
					fflush(stdout);
					//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}

#ifdef time_trace
	_time.irr_bin.markEnd();
	_time.irr_bin.getDuration();
#endif
#endif

		//std::cout << "Start IRR\n" << std::flush;
		//while (ParticleForComputation != nullptr) {

#ifdef chain_debug
		std::cerr << "BU, ComputationList=" << std::flush;
		for (Particle* ptcl:ComputationList) {
			fprintf(stderr,"%.3e (%d), ", ptcl->CurrentTimeIrr+ptcl->TimeStepIrr, ptcl->PID);
		}
		std::cerr << std::endl;

		fprintf(stderr,"BU,IrrChain=");
		Particle* next=FirstComputation;
		while (next!= nullptr) {
			fprintf(stderr,"%.3e (%d), ", next->CurrentTimeIrr+next->TimeStepIrr, next->PID);
			next = next->NextParticleForComputation;
		}
		fprintf(stderr,"\n");

#endif



		for (Particle* ptcl:ComputationList) {


#ifdef time_trace
			_time.irr_force.markStart();
#endif

			ptcl->calculateIrrForce(); // this includes particle position
#ifdef time_trace
			_time.irr_force.markEnd();
			_time.irr_force.getDuration();
#endif


			//ParticleForComputation = ParticleForComputation->NextParticleForComputation;// this includes time evolution.
			//ParticleForComputation = SortComputationChain(ParticleForComputation);
		}


		for (Particle* ptcl:particle) {
			if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
					|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
				fprintf(stdout, "before term, myself = %d\n", ptcl->PID);
				fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
				fflush(stdout);
				//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}

		// update particles and chain
		// The next particle of the particle calculated lastly should be the start of the next iteration.
#ifdef binary
		binary_time_prev = binary_time;
		binary_block = ComputationList[0]->CurrentBlockIrr+ComputationList[0]->TimeBlockIrr;
		binary_time  = binary_block*time_step;
		global_time_irr = ComputationList[0]->CurrentBlockIrr+ComputationList[0]->TimeBlockIrr;
		//std::cout << "ComputationList of " << ComputationList.size() << " : " ;
		// bool bin_termination=false;	
#endif




		for (Particle* ptcl:ComputationList) {
			if (ptcl->CurrentTimeIrr > 1) 
				fprintf(stderr, "outside, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
			//std::cout << ptcl->PID << " " ;
			ptcl->updateParticle();
			ptcl->CurrentBlockIrr += ptcl->TimeBlockIrr;
			ptcl->CurrentTimeIrr   = ptcl->CurrentBlockIrr*time_step;
			//std::cout << "before TimeStepCal\n" << std::flush;
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
			//std::cout << "after TimeStepCal\n" << std::flush;
			ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
#ifdef binary
			if (ptcl->isCMptcl) {
				ptcl->GroupInfo->ARIntegration(binary_time, particle);

/*  // Terminate a group if there is a particle is far away from the every other particles.
	// This version needs too many integration step for more than two body case (Too high stab: unstable)
	
				for (size_t i = 0; i < ptcl->GroupInfo->Members.size(); ++i) {

					Particle* ptclI = ptcl->GroupInfo->Members[i];
					bool isFarFromAll = true;

					for (size_t j = 0; j < ptcl->GroupInfo->Members.size(); ++j) {

						if (i == j) continue;

						Particle* ptclJ = ptcl->GroupInfo->Members[j];
						if (dist(ptclI->Position, ptclJ->Position) < 1e-2/position_unit) {
							isFarFromAll =false;
							break;

						}
					}
					if (isFarFromAll) {
						fprintf(binout, "Terminating Binary at time : %e \n", binary_time);
						fprintf(stdout, "Terminating Binary at time : %e \n", binary_time);
						FBTermination(ptcl, particle);
						bin_termination=true;
						break;	
					}
				}
*/

				for (size_t i = 0; i < ptcl->GroupInfo->Members.size(); ++i) {
					Particle* ptclI = ptcl->GroupInfo->Members[i];
					for (size_t j = i + 1; j < ptcl->GroupInfo->Members.size(); ++j) {
						Particle* ptclJ = ptcl->GroupInfo->Members[j];
						REAL rv = (ptclI->Position[0]-ptclJ->Position[0])*(ptclI->Velocity[0]-ptclJ->Velocity[0])
									+ (ptclI->Position[1]-ptclJ->Position[1])*(ptclI->Velocity[1]-ptclJ->Velocity[1])
									+ (ptclI->Position[2]-ptclJ->Position[2])*(ptclI->Velocity[2]-ptclJ->Velocity[2]);

						// if ((dist(ptclI->Position, ptclJ->Position) > 5e-3/position_unit) && (rv > 0)) {
						if (dist(ptclI->Position, ptclJ->Position) > 5e-3/position_unit) {

							fprintf(binout, "Terminating Binary at time : %e \n", binary_time);
							fprintf(stdout, "Terminating Binary at time : %e \n", binary_time);
							FBTermination(ptcl, particle);
							bin_termination=true;
							break;							
						}
					}
					if (bin_termination=true) 
						break;
				}

			}

			if (bin_termination=true) continue;

#endif


#ifdef time_trace
		_time.irr_sort.markStart();
#endif
			UpdateComputationChain(ptcl);
#ifdef time_trace
		_time.irr_sort.markEnd();
		_time.irr_sort.getDuration();
#endif
			if (ptcl->CurrentTimeIrr > 1) {
				fprintf(stderr, "outside after, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
				fflush(stderr);
			}
		}

		for (Particle* ptcl:particle) {
			if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
					|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
				fprintf(stdout, "after term, myself = %d\n", ptcl->PID);
				fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
				fflush(stdout);
				//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}
#ifdef binary
		if (bin_termination)  {
			//std::cerr << "After termination, NextRegTimeBlock=" << NextRegTimeBlock << std::endl; 
			CreateComputationChain(particle); // because NextRegTime might have been changed.
			/* 
				 fprintf(stderr, "in irr, particle: ");
				 for (Particle* ptcl: particle) {
				 fprintf(stderr,"%d, ",ptcl->PID);	
				 }
				 fprintf(stderr,"\n");	
				 */
			bin_termination = false; // Change its value to false again.
		}
#endif



		//std::cout << "\ndone!" << std::endl ;

#ifdef chain_debug
		std::cerr << "AU, ComputationList=" << std::flush;
		for (Particle* ptcl:ComputationList) {
			fprintf(stderr,"%.3e (%d), ", ptcl->CurrentTimeIrr+ptcl->TimeStepIrr, ptcl->PID);
		}
		fprintf(stderr,"\n");

		fprintf(stderr,"AC,AU,IrrChain=");
		next=FirstComputation;
		while (next!= nullptr) {
			fprintf(stderr,"%.3e (%d), ", next->CurrentTimeIrr+next->TimeStepIrr, next->PID);
			next = next->NextParticleForComputation;
		}
		fprintf(stderr,"\n");
#endif


#define no_IRR_TEST
#ifdef IRR_TEST
		// create output at appropriate time intervals just for irr
		if (outputTime <= ComputationList[0]->CurrentTimeIrr ) {
			writeParticle(particle, ComputationList[0]->CurrentTimeIrr, outNum++);
			outputTime += outputTimeStep;
		}
#endif

		//fflush(stdout);
	}
	//std::cout << "Finishing irregular force ..." << std::endl;
	//kstd::cerr << "Finishing irregular force ..." << std::endl;
	//fflush(stderr);



	//std::cout << "Finishing irregular force ...\n" << std::endl;
	return true;
}



