#include <iostream>
#include <stdio.h>
#include "global.h"

Particle *particles_original;
Particle *particles;

int *queues_original;
int *queues;

int *tasks_original;
int *tasks;

int *ActiveIndexToOriginalIndex;
int *ActiveIndexToOriginalIndex_original;

MPI_Win win;
MPI_Win win2;
MPI_Win win3;
MPI_Win win4;
MPI_Win win5;

MPI_Comm shared_comm;
int MyRank;
int NumberOfProcessor;
int NumberOfWorker;
int NumberOfCommunication;


GlobalVariable *global_variable;
GlobalVariable *global_variable_original;
int LastParticleIndex; // The last index of particle array
int NumberOfParticle; // The number of active particles
int NewPID;

// Task
int Task[NumberOfTask];

// Time
double global_time;
double global_time_irr;
ULL NextRegTimeBlock;
int time_block;
double time_step;
ULL block_max;
double outputTimeStep;
double endTime;

double binary_time;
double binary_time_prev;
ULL binary_block;

// Enzo to Nbody
Particle* FirstEnzoParticle;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
double EnzoTimeStep;



// i/o
char* fname;
double inputTime;
bool restart;
char* foutput;
bool IsOutput;
double outputTime;
int outNum;

FILE* binout;
FILE* mergerout;
#ifdef SEVN
FILE* SEVNout;
#endif
FILE* workerout;

#ifdef PerformanceTrace
Performance performance;
#endif

void DefaultGlobal() {

	/* Task initialization */
	//int Task[NumberOfTask];
	for (int i=0;i<NumberOfTask; i++) {
		Task[i] = i;
	}

	NumberOfCommunication = 0;

	/* Timesteps */
	endTime = 1;
	EnzoTimeStep   = endTime/1e10; // endTime should be Myr
	outputTimeStep = outputTimeStep/endTime; // endTime should be Myr

	time_block = -30;
	block_max = static_cast<ULL>(pow(2, -time_block));
	time_step = std::pow(2,time_block);

	inputTime = 0.0;
	endTime = 0.0;
	outputTimeStep = 0.;

	global_time = 0.;
	outputTime = 0.;

}
