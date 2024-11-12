#include <iostream>
#include "global.h"
#include <cassert>

#define no_chain_debug
#define no_Eunwoo_debug // Eunwoo

std::vector<Particle*> ComputationList{};
int writeParticle(std::vector<Particle*> &particle, REAL MinRegTime, int outputNum);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool UpdateComputationChain(Particle* ptcl);
//bool UpdateComputationChainParallel(std::vector<Particle*> &particle);
bool CreateComputationListParallel(std::vector<Particle*> &particle);
bool CreateComputationList(Particle* ptcl);
bool AddNewGroupsToList(std::vector<Particle*> &particle);
bool AddNewGroupsToList2(std::vector<Particle*> &members, std::vector<Particle*> &particle);
void NewFBInitialization2(std::vector<Particle*> &stillGroup, std::vector<Particle*> &particle);
void FBTermination(Particle* ptclCM, std::vector<Particle*> &particle);
void FBTermination2(Particle* ptclCM, std::vector<Particle*> &stillGroup, std::vector<Particle*> &noMoreGroup, std::vector<Particle*> &particle);



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
		_time.irr_bin_search.markStart();
#endif

		bool binSearch = AddNewGroupsToList(particle);

		// if (AddNewGroupsToList(particle) && ComputationList.size() == 0) { // original
		if (binSearch && ComputationList.size() == 0) {
			//std::cout << "No irregular particle to update afte binary formation." << std::endl;
			fprintf(stdout, "No irregular particle to update after binary formation.\n");
			break;
		}

#ifdef time_trace
		_time.irr_bin_search.markEnd();
		_time.irr_bin_search.getDuration();
#endif

#ifdef Eunoo_debug // Eunwoo optimization test
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
				if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
						|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
					fprintf(stdout, "after, myself = %d\n", ptcl->PID);
					fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
					fflush(stdout);
					//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}
#endif // Eunwoo optimization test
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


#ifdef time_trace
		_time.irr_force.markStart();
#endif
		for (Particle* ptcl:ComputationList) {

			ptcl->calculateIrrForce(); // this includes particle position

			//ParticleForComputation = ParticleForComputation->NextParticleForComputation;// this includes time evolution.
			//ParticleForComputation = SortComputationChain(ParticleForComputation);
		}
#ifdef time_trace
		_time.irr_force.markEnd();
		_time.irr_force.getDuration();
#endif

#ifdef Eunwoo_debug // Eunwoo optimization debug
		for (Particle* ptcl:particle) {
			if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
					|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
				fprintf(stdout, "before term, myself = %d\n", ptcl->PID);
				fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
				fflush(stdout);
				//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}
#endif // Eunwoo optimization debug

		// update particles and chain
		// The next particle of the particle calculated lastly should be the start of the next iteration.
#ifdef binary
		binary_time_prev = binary_time;
		binary_block = ComputationList[0]->CurrentBlockIrr+ComputationList[0]->TimeBlockIrr;
		binary_time  = binary_block*time_step;
		// global_time_irr = ComputationList[0]->CurrentBlockIrr+ComputationList[0]->TimeBlockIrr;
		//std::cout << "ComputationList of " << ComputationList.size() << " : " ;
#endif


		for (Particle* ptcl:ComputationList) {
#ifdef Eunwoo_debug // Eunwoo optimization debug
			if (ptcl->CurrentTimeIrr > 1) 
				fprintf(stderr, "outside, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
#endif // Eunwoo optimizatino debug
			// std::cout << ptcl->PID << " " ;
#ifdef binary
			// if (ptcl->isCMptcl) {
			// 	if (ptcl->GroupInfo->Members[0]->PID == 716 && ptcl->GroupInfo->Members[1]->PID == 44) {
			// 		fprintf(mergerout, "TimeStepIrr: %e Myr\n", ptcl->TimeStepIrr*EnzoTimeStep*1e4);
			// 		fprintf(mergerout, "TimeStepReg: %e Myr\n", ptcl->TimeStepReg*EnzoTimeStep*1e4);
			// 		fprintf(mergerout, "ACnum: %d\n", ptcl->NumberOfAC);
			// 	}
			// }
			if (ptcl->isCMptcl) {
				// if (ptcl->NumberOfAC == 0)
				// 	fprintf(stderr, "No neighbor! PID: %d\n", ptcl->PID);
#ifdef time_trace
				_time.irr_bin_int.markStart();
#endif
				bool int_normal = ptcl->GroupInfo->ARIntegration(ptcl->CurrentTimeIrr + ptcl->TimeStepIrr, particle);
#ifdef time_trace
				_time.irr_bin_int.markEnd();
				_time.irr_bin_int.getDuration();
#endif
				if (int_normal) {
				// if (ptcl->GroupInfo->ARIntegration(ptcl->CurrentTimeIrr + ptcl->TimeStepIrr, particle)){ // true: integrated normally
					if (ptcl->GroupInfo->CheckBreak()) {
						if (ptcl->GroupInfo->Members.size() > 2) {
							std::vector<Particle*> members = ptcl->GroupInfo->Members;
							FBTermination(ptcl, particle);
							AddNewGroupsToList2(members, particle);
							members.clear();
							bin_termination = true;
							continue;
						}
						else {
							FBTermination(ptcl, particle);
							bin_termination = true;
							continue; // Do not UpdateComputationChain!
						}
					}
				}
				else { // false: terminated by stellar merger, TDE, GW merger, etc.
					bin_termination = true;
					continue; // Do not UpdateComputationChain!
				}	
			}
#endif
			if (ptcl->NumberOfAC != 0) // IAR modified
				ptcl->updateParticle();
			ptcl->CurrentBlockIrr += ptcl->TimeBlockIrr;
			ptcl->CurrentTimeIrr   = ptcl->CurrentBlockIrr*time_step;
// #ifdef time_trace
// 			_time.irr_calctime.markStart();
// #endif
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
// #ifdef time_trace
// 			_time.irr_calctime.markEnd();
// 			_time.irr_calctime.getDuration();
// #endif
			ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle


#ifdef time_trace
		_time.irr_sort.markStart();
#endif
			UpdateComputationChain(ptcl);
#ifdef time_trace
		_time.irr_sort.markEnd();
		_time.irr_sort.getDuration();
#endif

#ifdef Eunwoo_debug // Eunwoo optimization debug
			if (ptcl->CurrentTimeIrr > 1) {
				fprintf(stderr, "outside after, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
				fflush(stderr);
			}
#endif // Eunwo optimization debug
		}

#ifdef Eunwoo_debug // Eunwoo optimization debug
		for (Particle* ptcl:particle) {
			if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
					|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
				fprintf(stdout, "after term, myself = %d\n", ptcl->PID);
				fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
				fflush(stdout);
				//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}
#endif // Eunwo optimization debug

#ifdef binary
		if (bin_termination)  {
			//std::cerr << "After termination, NextRegTimeBlock=" << NextRegTimeBlock << std::endl;
#ifdef time_trace
			_time.irr_chain.markStart();
#endif
			CreateComputationChain(particle); // because NextRegTime might have been changed.
			/* 
			fprintf(stderr, "in irr, particle: ");
			for (Particle* ptcl: particle) {
			fprintf(stderr,"%d, ",ptcl->PID);	
			}
			fprintf(stderr,"\n");	
			*/
			bin_termination = false; // Change its value to false again.
#ifdef time_trace
			_time.irr_chain.markEnd();
			_time.irr_chain.getDuration();
#endif
		}
#endif


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

	}


	//std::cout << "Finishing irregular force ...\n" << std::endl;
	return true;
}



