#include <iostream>
#include <cmath>
#include "../global.h"


REAL getNewTimeStepReg(REAL v[3], REAL f[3][4]);
REAL getNewTimeStepIrr(REAL f[3][4], REAL df[3][4]);
void getBlockTimeStep(REAL dt, int& TimeLevel, ULL &TimeBlock, REAL &TimeStep);

// Update TimeStepIrr
void Particle::calculateTimeStepIrr(REAL f[3][4],REAL df[3][4]) {
	REAL TimeStepTmp;
	int TimeLevelTmp, TimeLevelTmp0;
	ULL TimeBlockTmp;

	if (this->NumberOfAC == 0) {
		TimeLevelIrr = TimeLevelReg;
		TimeStepIrr = static_cast<REAL>(pow(2, TimeLevelIrr));
		TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-time_block));
		return;
	}

	getBlockTimeStep(getNewTimeStepIrr(a_tot, a_irr), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	/*
	if (PID == 965) {
		fprintf(stderr, "stepirr=%.3e, levelirr=%d, blockirr=%llu\n", TimeStepTmp, TimeLevelTmp, TimeBlockTmp);
		fprintf(stderr, "currentblokirr=%llu, nexttimereg=%llu \n", CurrentBlockIrr, NextRegTimeBlock);
		fprintf(stderr, "currentblockreg=%llu, levelreg=%d, blockreg=%llu\n", CurrentBlockReg, TimeLevelReg, TimeBlockReg);
	}
	*/

	/*
	if (NextRegTimeBlock < TimeBlockReg && isRegular) {
		fprintf(stderr, "PID=%d, NextRegTimeBlock=%llu, TimeBlockReg=%llu\n", PID, NextRegTimeBlock, TimeBlockReg);
	}
	*/

	/*
	if (PID == 430) {
		std::cerr << "Before TimeLevelIrr=" << TimeLevelIrr << ", TimeLevelReg="<< TimeLevelReg << std::endl;
		std::cerr << "       TimeLevelTmp=" << TimeLevelTmp << std::endl;
	}
	*/

	if (TimeLevelTmp > TimeLevelIrr) {
		if (fmod(CurrentBlockIrr, 2*TimeBlockIrr)==0) {
			TimeLevelTmp = TimeLevelIrr+1;
			TimeBlockTmp = 2*TimeBlockIrr;
		}
		else {
			TimeLevelTmp = TimeLevelIrr;
			TimeBlockTmp = TimeBlockIrr;
		}
	}
	else if (TimeLevelTmp < TimeLevelIrr) {
		if (TimeLevelTmp < TimeLevelIrr-1) {
			TimeLevelTmp = TimeLevelIrr - 2;
			TimeBlockTmp = TimeBlockIrr/4;
		}
		else {
			TimeLevelTmp = TimeLevelIrr - 1;
			TimeBlockTmp = TimeBlockIrr/2;
		}
	} else {
		TimeLevelTmp = TimeLevelIrr;
		TimeBlockTmp = TimeBlockIrr;
	}

	/*
	if (PID == 430) {
		std::cerr << "Middle TimeLevelTmp=" << TimeLevelTmp << std::endl;
		std::cerr << std::endl;
	}
	*/

	//std::cout << "TimeStepIrrTmp=" << TimeStepIrrTmp << std::endl;
	while (((CurrentBlockIrr < CurrentBlockReg+TimeBlockReg) && (CurrentBlockIrr+TimeBlockTmp > CurrentBlockReg+TimeBlockReg)) || (TimeLevelTmp >= TimeLevelReg)) {
		/*
		fprintf(stderr,"CurrentBlockIrr = %llu\n",
				CurrentBlockIrr);
		fprintf(stderr,"CurrentBlockReg = %llu\n",
				CurrentBlockReg);
		fprintf(stderr,"TimeBlockReg    = %llu\n",
				TimeBlockReg);
		fprintf(stderr,"TimeBlockIrr    = %llu\n",
				TimeBlockIrr);
				*/
		//if (TimeLevelTmp > TimeLevelReg)
			//fprintf(stderr, "PID=%d, Irr=%d, Reg=%d\n", PID, TimeLevelTmp0, TimeLevelReg);

		TimeLevelTmp--;
		TimeBlockTmp *= 0.5;
	}

	TimeLevelIrr = TimeLevelTmp;

	if (TimeLevelIrr < time_block) {
		//std::cerr << "TimeLevelIrr is too small" << std::endl;
		TimeLevelIrr = std::max(time_block, TimeLevelIrr);
	}


	TimeStepIrr = static_cast<REAL>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-time_block));

	if (TimeStepIrr > 1) {
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",TimeStepIrr, TimeLevelIrr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}


	if (TimeStepIrr*EnzoTimeStep*1e4 < KSTime && (this->isCMptcl == false)) // Eunwoo: this should be fixed.
	// Eunwoo: I don't want to use KSTime but let's think about it.
	// if (this->isCMptcl == false) // Eunwoo: just for debug
		GroupCandidateList.push_back(this); // Eunwoo edited
		// fprintf(binout, "GroupCandidateList added!\n");
	/*
	if (PID == 430) {
		std::cerr << "After TimeLevelIrr=" << TimeLevelIrr << std::endl;
		std::cerr << std::endl;
	}
	*/
}



// Update TimeStepReg // need to review
void Particle::calculateTimeStepReg() {
	//fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	REAL TimeStepTmp;
	ULL TimeBlockTmp;
	int TimeLevelTmp, TimeLevelTmp0;

	getBlockTimeStep(getNewTimeStepReg(Velocity, a_reg), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);

	//fprintf(stderr, "in CalReg, raw time step=%.2eMyr, ", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);


	//std::cout << "NBODY+: TimeStepRegTmp = " << TimeStepRegTmp << std::endl;

	TimeLevelTmp0 = TimeLevelTmp;

	if (TimeLevelTmp >= TimeLevelReg+1) {
		if (fmod(CurrentBlockReg, 2*TimeBlockReg)==0 \
				&& CurrentTimeReg != 0) {
			TimeLevelTmp   = TimeLevelReg+1;
			TimeBlockTmp   = TimeBlockReg*2;

			while ((TimeLevelTmp0 > TimeLevelTmp) \
					&& (fmod(CurrentBlockReg, 2*TimeBlockTmp) == 0)) {
				TimeBlockTmp = 2*TimeBlockTmp;
				TimeLevelTmp   = TimeLevelTmp + 1;
			}
		}
		else {
			TimeLevelTmp   = TimeLevelReg;
		}
	}
	else if (TimeLevelTmp < TimeLevelReg) {
		TimeLevelTmp = TimeLevelReg - 1;
		if (TimeLevelTmp0 < TimeLevelTmp)
			TimeLevelTmp--;
	}
	else {
		TimeLevelTmp = TimeLevelReg;
	}


	// update needed. regcor_gpu.for:725 (Makino, ApJ, 369)
	/*
	if (TimeStepRegTmp > 0.1 && TimeStepRegTmp > TimeStepReg) {
		REAL v2 = 0., a2=0., dt;
		for (int dim=0; dim<Dim; dim++) {
			v2 += (PredVelocity[dim]-NewVelocity[dim])*(PredVelocity[dim]-NewVelocity[dim]);
			a2 += a_reg[dim][0]*a_reg[dim][0];
		}
		dt = TimeStepReg*std::pow((1e-4*TimeStepReg*TimeStepReg*a2/v2),0.1);
		if (dt < TimeStepRegTmp) {
			TimeStepRegTmp = TimeStepReg;
		}	
	}
	*/

	//fprintf(stderr, " final time step=%.2eMyr\n", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);

	TimeLevelReg = std::max(time_block,TimeLevelTmp);
	//TimeLevelReg = std::max(time_block, TimeLevelReg);

	if (this->NumberOfAC == 0) {
		TimeLevelIrr = TimeLevelReg;
	}

	TimeStepReg  = static_cast<REAL>(pow(2, TimeLevelReg));
	TimeBlockReg = static_cast<ULL>(pow(2, TimeLevelReg-time_block));

	if (TimeStepReg*EnzoTimeStep*1e4 < 1e-7) {
		fprintf(stderr, "TimeStep = %.3e, TimeStepTmp0 = %.3e\n",
			 	TimeStepReg*EnzoTimeStep*1e4, static_cast<REAL>(pow(2, TimeLevelTmp0))*EnzoTimeStep*1e4);
		fflush(stderr);
	}

	if (CurrentTimeReg+TimeStepReg > 1 && CurrentTimeReg != 1.0) {
		TimeStepReg = 1 - CurrentTimeReg;
		TimeBlockReg = block_max-CurrentBlockReg;
	}

	if (TimeStepReg*EnzoTimeStep*1e4<1e-9) {
		throw std::runtime_error("TimeStepReg is too small.");
	}
	if (TimeStepReg > 1) {
		fprintf(stderr, "TimeStepReg=%e, TimeLevelReg=%d, TimeLevelTmp0=%d\n",TimeStepReg, TimeLevelReg, TimeLevelTmp0);
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",TimeStepIrr, TimeLevelIrr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}

	//std::cout << "NBODY+: TimeStepReg = " << TimeStepReg << std::endl;
}
