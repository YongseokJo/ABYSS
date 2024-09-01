#include <cmath>
#include <iostream>
#include "defs.h"
#include "global.h"
#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

int writeParticle(std::vector<Particle*> &particle, REAL MinRegTime, int outputNum);
// int ReceiveFromEzno(std::vector<Particle*> &particle);
// int SendToEzno(std::vector<Particle*> &particle);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool RegularAccelerationRoutine(std::vector<Particle*> &particle);
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle);

bool IsOutput           = false;
REAL binary_time      = 0;
REAL binary_time_prev = 0;
ULL binary_block        = 0;
REAL outputTime       = 0;
ULL NextRegTimeBlock    = 0;
int outNum              = 0;
REAL global_time_irr  = 0;
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

	fprintf(stdout, "output time step: %e\n", outputTimeStep); // Eunwoo debug

	while (true) {

		// It's time to compute regular force.
#ifdef time_trace
		_time.reg.markStart();
#endif
#ifdef NSIGHT
nvtxRangePushA("RegularAccelerationRoutine");
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

		IrregularAccelerationRoutine(particle);

#ifdef NSIGHT
nvtxRangePop();
#endif

#ifdef time_trace
		_time.irr.markEnd();

		_time.irr.getDuration();
		_time.output();
#endif

		global_time = NextRegTimeBlock*time_step;

		// fprintf(stdout, "global time: %e\n", global_time); // Eunwoo debug

		// create output at appropriate time intervals
		if (outputTime <= global_time ) {
			writeParticle(particle, global_time, outNum++);
			outputTime += outputTimeStep;
/*
			for (Particle* pp : particle) {
				if (pp->PID == 718) {
					fprintf(binout, "outputTime: %e\n", global_time);
					fprintf(binout, "PID: %d. Position - x:%e, y:%e, z:%e, \n", pp->PID, pp->Position[0], pp->Position[1], pp->Position[2]);
					fprintf(binout, "PID: %d. Velocity - vx:%e, vy:%e, vz:%e, \n", pp->PID, pp->Velocity[0], pp->Velocity[1], pp->Velocity[2]);
					fprintf(binout, "PID: %d. Mass - %e, \n", pp->PID, pp->Mass);
					// fprintf(binout, "PID: %d. Distance from mother particle (pc): %e\n", pp->PID, dist(ptclI->Position, pp->Position)*position_unit);
					fprintf(binout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", pp->PID, pp->a_tot[0][0], pp->a_tot[1][0], pp->a_tot[2][0]);
					fprintf(binout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", pp->PID, pp->a_tot[0][1], pp->a_tot[1][1], pp->a_tot[2][1]);
					fprintf(binout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", pp->PID, pp->a_tot[0][2], pp->a_tot[1][2], pp->a_tot[2][2]);
					fprintf(binout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", pp->PID, pp->a_tot[0][3], pp->a_tot[1][3], pp->a_tot[2][3]);
					fprintf(binout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", pp->PID, pp->a_irr[0][0], pp->a_irr[1][0], pp->a_irr[2][0]);
					fprintf(binout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", pp->PID, pp->a_irr[0][1], pp->a_irr[1][1], pp->a_irr[2][1]);
					fprintf(binout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", pp->PID, pp->a_irr[0][2], pp->a_irr[1][2], pp->a_irr[2][2]);
					fprintf(binout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", pp->PID, pp->a_irr[0][3], pp->a_irr[1][3], pp->a_irr[2][3]);
					fprintf(binout, "PID: %d. Time Steps (Myr) - irregular:%e, regular:%e \n", pp->PID, pp->TimeStepIrr*EnzoTimeStep*1e4, pp->TimeStepReg*EnzoTimeStep*1e4);
					fprintf(binout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", pp->PID, pp->TimeBlockIrr, pp->TimeBlockReg);
					fprintf(binout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", pp->PID, pp->CurrentBlockIrr, pp->CurrentBlockReg);
				}
			}
*/
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



