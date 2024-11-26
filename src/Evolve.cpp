#include <cmath>
#include <iostream>
#include "defs.h"
#include "global.h"
#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

// #define SEVN
#ifdef SEVN
void UpdateEvolution(Particle* ptcl);
#endif

int writeParticle(std::vector<Particle*> &particle, REAL MinRegTime, int outputNum);
// int ReceiveFromEzno(std::vector<Particle*> &particle);
// int SendToEzno(std::vector<Particle*> &particle);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool RegularAccelerationRoutine(std::vector<Particle*> &particle);
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle);
bool IrregularAccelerationRoutineParallel(std::vector<Particle*> &particle);


bool IsOutput           = false;
REAL binary_time      = 0;
REAL binary_time_prev = 0;
ULL binary_block        = 0;
REAL outputTime       = 0;
ULL NextRegTimeBlock    = 0;
int outNum              = 0;
// REAL global_time_irr  = 0;
ULL global_time_irr  = 0;
std::vector<Particle*> ComputationChain{};
#ifdef time_trace
TimeTracer _time;
#endif



void Evolve(std::vector<Particle*> &particle) {

	std::cout << "Evolve starting ..." << std::endl;
	int freq   = 0;


//	if (NNB == 0) {
//		std::cout << "No particle to be calculated ..." << std::endl;
//		goto Communication;
//	}

	//CreateComputationChain(particle);

	writeParticle(particle, global_time, outNum++);
	outputTime = outputTimeStep;

	while (true) {

		// It's time to compute regular force.
#ifdef time_trace
		_time.reg.markStart();
#endif
#ifdef NSIGHT
nvtxRangePushA(("RegularAccelerationRoutine " + std::to_string(RegularList.size())).c_str());
#endif
		RegularAccelerationRoutine(particle); // do not update particles unless NNB=0
		
#ifdef NSIGHT
nvtxRangePop();
#endif


#ifdef time_trace
		_time.reg.markEnd();
		_time.reg.getDuration();
		_time.irr.markStart();
#endif
		
#ifdef NSIGHT
nvtxRangePushA("IrregularAccelerationRoutine");
#endif


#ifdef OMP
		IrregularAccelerationRoutineParallel(particle);
#else
		IrregularAccelerationRoutine(particle);
#endif

#ifdef NSIGHT
nvtxRangePop();
#endif

#ifdef time_trace
		_time.irr.markEnd();
		_time.irr.getDuration();
#ifndef SEVN
		_time.output();
#endif
#ifdef SEVN
		_time.sevn.markStart();
		_time.output();
#endif
#endif

		global_time = NextRegTimeBlock*time_step;

#ifdef SEVN
		bool evolved;

		for(Particle* ptcl: particle) {

			evolved = false;
		
			if (ptcl->star == nullptr || ptcl->star->amiremnant()) // CM particle & Stars with mass < 2.2 Msol don't have star class.
				continue;

			while (ptcl->EvolutionTime + ptcl->star->getp(Timestep::ID) < global_time*EnzoTimeStep*1e4) {
				// ptcl->star->sync_with(global_time*EnzoTimeStep*1e4 - ptcl->EvolutionTime); // Eunwoo: sometimes, this call some errors!
				// fprintf(SEVNout, "Before. PID: %d, GT: %e, ET: %e, dt: %e\n", ptcl->PID, global_time*EnzoTimeStep*1e4, ptcl->EvolutionTime, ptcl->star->getp(Timestep::ID));
				// fprintf(SEVNout, "dt_a: %e, dt_b: %e\n", ptcl->star->getp(Timestep::ID), global_time*EnzoTimeStep*1e4 - ptcl->EvolutionTime);
				// fflush(SEVNout);
				// assert(ptcl->star->getp(Timestep::ID) == global_time*EnzoTimeStep*1e4 - ptcl->EvolutionTime);
				ptcl->EvolutionTime += ptcl->star->getp(Timestep::ID);
				ptcl->star->evolve();
				// if (ptcl->Mass*mass_unit > 1000) {
				// 	fprintf(SEVNout, "PID: %d, Mass: %e, Radius: %e, EvolutionTime: %e\n", ptcl->PID, ptcl->Mass*mass_unit, ptcl->radius*position_unit, ptcl->EvolutionTime);
				// 	fprintf(SEVNout, "EvolutionTimeStep: %e\n", ptcl->star->getp(Timestep::ID));
				// 	fflush(SEVNout);
				// }
				// fprintf(SEVNout, "After. PID: %d, ET: %e, WT: %e\n", ptcl->PID, ptcl->EvolutionTime, ptcl->star->getp(Worldtime::ID));
				// assert(ptcl->EvolutionTime == ptcl->star->getp(Worldtime::ID));
				// fflush(SEVNout);
				evolved = true;
			}

			if (evolved) {
				UpdateEvolution(ptcl);
			}
							
		}
		if (MasslessList.size() != 0) {
			particle.erase(
				std::remove_if(particle.begin(), particle.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					//if (to_remove) delete p;
					return to_remove;
					}),
				particle.end());

			for (int i=0; i<particle.size(); i++) {
				particle[i]->ParticleOrder = i;
			}
			// This might be not necessary because we're moving to the regular acceleration routine, and re-set neighbors.
			// Should be checked again later.
			for (Particle* massless: MasslessList) {
				int index = 0;
				for (Particle* ptcl: particle) {
					index = 0;
					for (Particle* neighbor: ptcl->ACList) {
						if (neighbor->PID == massless->PID) {
							ptcl->ACList.erase(ptcl->ACList.begin() + index);
							ptcl->NumberOfAC--;
							break; // Eunwoo check
						}
						index++;
					}
				}
				// fprintf(SEVNout, "3\n");
				// fflush(SEVNout);
				// delete massless->star; // I don't know why but there is an issue with Star::default_destructor(). This must be fixed!
				// delete massless;
				// fprintf(SEVNout, "4\n");
				// fflush(SEVNout);
			}
		}
#endif

#ifdef SEVN
#ifdef time_trace
		_time.sevn.markEnd();
		_time.sevn.getDuration();
		_time.output();
#endif
#endif

		// create output at appropriate time intervals
		if (outputTime <= global_time ) {
			writeParticle(particle, global_time, outNum++);
			outputTime += outputTimeStep;
		}

		// end if the global time exceeds the end time
		if (global_time >= 1) {
			std::cout << EnzoTimeStep << std::endl;
			return;
			//exit(EXIT_FAILURE);
		}


	//Communication:
	//	do
	//	{
	//		SendToEzno(particle);
	//		ReceiveFromEzno(particle);
	//	} while (NNB == 0);
	//	global_time = 0.;
	//	NextRegTime = 0.;
	}
}



