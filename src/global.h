#include <vector>
// #include <mpi.h>
#include <iostream>
#include "Particle/Particle.h"


#ifdef time_trace
#include "TimeTrace.h"
extern TimeTracer _time;
#endif





extern std::vector<int> LevelList;
extern int NNB;
extern int newNNB;
extern double global_time;
extern double NextRegTime;
extern const double dt_min;
extern const int dt_level_min;
extern double dt_block;  // this stores the minimum time step
extern int dt_block_level;
extern std::vector<Particle*> ComputationChain;
extern std::vector<int> RegIndexList; // how about changing this to particle list
extern Particle* FirstComputation;
extern int NumNeighborMax;


//extern bool debug;
extern char* fname;
extern double inputTime;
extern double endTime;
extern bool restart;


// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTime;
extern double EnzoTimeStep;


// i/o
extern char* foutput;
extern bool IsOutput;
extern double outputTime;
extern double outputTimeStep;
extern int outNum;
//
//


// MPI variables
//extern MPI_Comm inter_comm;
//extern MPI_Comm nbody_comm;

// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTime;
extern double EnzoTimeStep;



