#ifndef GLOBAL_H
#define GLOBAL_H
#include <vector>
// #include <mpi.h>
#include <iostream>
#include "Particle/Particle.h"
#include "FewBody/Group.h" // Eunwoo added
//#include "ParticleScheduler/ParticleScheduler.h"
#include <stdio.h>
#include <stdexcept>

#ifdef time_trace
#include "TimeTrace.h"
extern TimeTracer _time;
#endif



extern std::vector<int> LevelList;
extern int NNB;
extern int newNNB;
//extern int NumNeighborMax;

// Time
extern REAL global_time;
// extern REAL global_time_irr;
extern ULL global_time_irr;
extern ULL NextRegTimeBlock;
extern int time_block;
extern REAL time_step;
extern ULL block_max;

extern REAL min_timestep; // Eunwoo added


extern REAL binary_time;
extern REAL binary_time_prev;
extern ULL binary_block;
extern bool bin_termination; // Eunwoo added

// ComputationChain
extern std::vector<Particle*> ComputationChain;
extern Particle* FirstComputation;
extern std::vector<Particle*> ComputationList;
extern int ComputationTimeMarker;
extern std::vector<Particle*> RegularList;
extern std::vector<Group*> GroupCandidateList; // Newly detected group candidates // Eunwoo added
// extern std::vector<Group*> GroupList; // List of groups to calculate // Eunwoo added

//extern bool debug;
extern char* fname;
extern REAL inputTime;
extern REAL endTime;
extern bool restart;


// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern REAL EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern REAL EnzoTime;
extern REAL EnzoTimeStep;


// i/o
extern char* foutput;
extern bool IsOutput;
extern REAL outputTime;
extern REAL outputTimeStep;
extern int outNum;
//
//

extern FILE* binout;
extern FILE* mergerout;
// #define SEVN
#ifdef SEVN
extern FILE* SEVNout;
// extern IO* sevnio;
extern std::vector<Particle*> MasslessList;
#endif

typedef std::vector<int> Ivector;
typedef std::vector<ULL> Uvector;
typedef std::vector<double> Dvector;
#endif
