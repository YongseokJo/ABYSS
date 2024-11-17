#include "def.h"
#include "particle.h"
#include <mpi.h>

extern int TaskType[NumberOfTaskTypes];
extern Particle *particles;


/* Communicators */
extern MPI_Win win;
extern int MyRank;
extern int NumberOfParticle;

// Task
const int ROOT = 0;
const int TASK_TAG = 1;
const int PTCL_TAG = 2;
const int TIME_TAG = 3;
const int ANY_TAG = 100;
const int TERMINATE_TAG = 666;

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
extern double EnzoTime;
extern double EnzoTimeStep;



