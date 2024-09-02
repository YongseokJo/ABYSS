#include "global.h"
#include <iostream>
#include <cmath>


int time_block = -30;
ULL block_max = static_cast<ULL>(pow(2, -time_block));
REAL time_step = std::pow(2,time_block);

REAL getNewTimeStep(REAL f[3][4], REAL df[3][4]);
void getBlockTimeStep(REAL dt, int& TimeLevel, ULL &TimeBlock, REAL &TimeStep);

/*
 *  Purporse: Initialize timsteps 
 *
 *  Modified: 2024.01.16  by Yongseok Jo
 *  Modified: 2024.05.29  by Yongseok Jo
 *
 */
int InitializeTimeStep(std::vector<Particle*> &particle) {

	std::cout << "Initializing timesteps ..." << std::endl;

	int min_time_level=0;
	REAL dtIrr, dtReg;

	for (Particle* ptcl: particle) {
		dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
		//std::cout << "dtReg=" << dtReg << std::endl;
		getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeBlockReg, ptcl->TimeStepReg);

		if (ptcl->NumberOfAC != 0) {
			dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
			getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeBlockIrr, ptcl->TimeStepIrr);
		}
		else {
			ptcl->TimeBlockIrr = ptcl->TimeBlockReg;
			ptcl->TimeLevelIrr = ptcl->TimeLevelReg;
			ptcl->TimeStepIrr  = ptcl->TimeStepReg;
		}

		ptcl->TimeStepReg  = MIN(1,ptcl->TimeStepReg);
		ptcl->TimeBlockReg = std::min(block_max, ptcl->TimeBlockReg);
		ptcl->TimeLevelReg = std::min(0, ptcl->TimeLevelReg);

		ptcl->CurrentTimeIrr  = 0;
		ptcl->CurrentTimeReg  = 0;
		ptcl->CurrentBlockIrr = 0;
		ptcl->CurrentBlockReg = 0;

#define no_IRR_TEST
#ifdef IRR_TEST
		ptcl->TimeStepReg = 1;
		ptcl->TimeLevelReg = 0;
		ptcl->TimeBlockReg = block_max;
#endif
	}



	// Irregular Time Step Correction
	for (Particle* ptcl: particle) {
		if (ptcl->NumberOfAC != 0) {
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

	for (Particle* ptcl: particle) {
		ptcl->TimeBlockIrr = static_cast<ULL>(pow(2, ptcl->TimeLevelIrr-time_block));
		ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg-time_block));
#ifdef IRR_TEST
		ptcl->TimeStepReg = 1;
		ptcl->TimeLevelReg = 0;
		ptcl->TimeBlockReg = block_max;
#endif
		ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
	}




	fprintf(stdout, "nbody+:time_block = %d, EnzoTimeStep=%e\n", time_block, EnzoTimeStep);
	return true;
}


int InitializeTimeStep(Particle* particle, int size) {
	std::cout << "Initializing timesteps ..." << std::endl;
	REAL dtIrr, dtReg;
	Particle *ptcl;

	for (int i=0; i<size; i++){
		ptcl = &particle[i];
		dtReg = getNewTimeStep(ptcl->a_reg, ptcl->a_reg);
		getBlockTimeStep(dtReg, ptcl->TimeLevelReg, ptcl->TimeBlockReg, ptcl->TimeStepReg);

		if (ptcl->NumberOfAC != 0) {
			dtIrr = getNewTimeStep(ptcl->a_tot, ptcl->a_irr);
			getBlockTimeStep(dtIrr, ptcl->TimeLevelIrr, ptcl->TimeBlockIrr, ptcl->TimeStepIrr);
		}
		else {
			ptcl->TimeBlockIrr = ptcl->TimeBlockReg;
			ptcl->TimeLevelIrr = ptcl->TimeLevelReg;
			ptcl->TimeStepIrr  = ptcl->TimeStepReg;
		}

		ptcl->TimeStepReg  = MIN(1,ptcl->TimeStepReg);
		ptcl->TimeBlockReg = std::min(block_max, ptcl->TimeBlockReg);
		ptcl->TimeLevelReg = std::min(0, ptcl->TimeLevelReg);

		ptcl->CurrentTimeIrr  = 0;
		ptcl->CurrentTimeReg  = 0;
		ptcl->CurrentBlockIrr = 0;
		ptcl->CurrentBlockReg = 0;
	} // endfor size

	for (int i=0; i<size; i++){
		ptcl = &particle[i];

		if (ptcl->NumberOfAC != 0) {
			while (ptcl->TimeLevelIrr >= ptcl->TimeLevelReg) {
				ptcl->TimeLevelIrr--;
			}
		}
		ptcl->TimeBlockIrr = static_cast<ULL>(pow(2, ptcl->TimeLevelIrr-time_block));
		ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg-time_block));
	} //endfor size
	ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
	return true;
}



