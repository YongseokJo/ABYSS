#ifndef GLOBAL_H
#define GLOBAL_H
#include "def.h"
#include "particle.h"
#include <mpi.h>

extern Particle *particles;
extern Particle *particles_original;

/* Communicators */
extern MPI_Win win;
extern MPI_Comm shared_comm;
extern int MyRank;
extern int NumberOfProcessor;
extern int NumberOfWorker;
const int ROOT = 0;



extern int NumberOfParticle;

// Task
const int TASK_TAG = 1;
const int PTCL_TAG = 2;
const int TIME_TAG = 3;
const int ANY_TAG = 100;
const int TERMINATE_TAG = 666;
extern int Task[NumberOfTask];

// Time
extern double global_time;
extern double global_time_irr;
extern ULL NextRegTimeBlock;
extern int time_block;
extern double time_step;
extern ULL block_max;


extern double binary_time;
extern double binary_time_prev;
extern ULL binary_block;

// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTimeStep;



// i/o
extern char* fname;
extern double inputTime;
extern double endTime;
extern bool restart;
extern char* foutput;
extern bool IsOutput;
extern double outputTime;
extern int outNum;
extern double outputTimeStep;


#endif
