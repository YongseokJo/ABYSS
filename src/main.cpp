#include <iostream>
//#include "defs.h"
#include "global.h"
#include "nbody.h"
#include "cuda/cuda_functions.h"


using namespace std;

//Global Variables
int NNB; REAL global_time; //bool debug;
bool bin_termination = false; // Eunwoo added
std::vector<int> LevelList;
REAL EnzoTimeStep;
int newNNB = 0;
std::vector<Particle*> RegularList;
std::vector<Group*> GroupCandidateList; // Eunwoo edited
std::vector<Group*> GroupList; // Eunwoo edited
FILE* binout;

//int NumNeighborMax = 100;



int main(int argc, char *argv[]) {
	cout << "Staring Nbody+ ..." << endl;
	binout = fopen("binary_output.txt", "w");                                                                          
	fprintf(binout, "Starting nbody - Binary OUTPUT\n"); 
	std::vector<Particle*> particle{};
	int irank=0;
	std::ios::sync_with_stdio(false);
	//comm        = com;
	//inter_comm  = inter_com;
	//nbody_comm  = nbody_com;

	global_time = 0.;
	//debug = true;

	//InitialCommunication(particle);
	Parser(argc, argv);

	EnzoTimeStep   = endTime/1e10; // endTime should be Myr
	outputTimeStep = outputTimeStep/endTime; // endTime should be Myr

	cout << "EnzoTimeStep = "   << EnzoTimeStep   << endl;
	cout << "outputTimeStep = " << outputTimeStep << endl;

	if (readData(particle) == FAIL)
		fprintf(stderr, "Read Data Failed!\n");

	/***
		for (Particle* elem: particle) {
		std::cout << elem->Position[0] <<" ";
		}
		std::cout << std::endl;
	 ***/

	fprintf(stderr, "Initializing Device!\n");
	InitializeDevice(&irank);
	fprintf(stderr, "Initializing Particles!\n");
	InitializeParticle(particle);


	//createComputationChain(particle);

	/*
	for (Particle* elem: particle) {
		std::cout << elem->TimeStepIrr <<" ";
	}
	*/

	/*
	for (Particle* elem: particle) {
		fprintf(stdout, "PID=%d, TReg=%e, TIrr=%e\n", elem->getPID(),elem->TimeStepReg, elem->TimeStepIrr);
		fprintf(stdout, "%e, %e, %e, %e, %e, %e\n",
				elem->Force[0], elem->Force[1], elem->Force[2], elem->ForceDot[0], elem->ForceDot[1], elem->ForceDot[2]);
		fprintf(stdout, "%e, %e, %e, %e, %e, %e\n",
				elem->dFReg[0][0], elem->dFReg[0][1], elem->dFReg[0][2], elem->dFIrr[0][0], elem->dFIrr[0][1], elem->dFIrr[0][3]);
		fprintf(stdout, "%e, %lf, %lf, %lf, %e, %e, %e\n\n",
				elem->Mass, elem->Position[0], elem->Position[1], elem->Position[2], elem->Velocity[0], elem->Velocity[1], elem->Velocity[2]);

	}
	std::cout << std::endl;
	*/



	Evolve(particle);

	// Particle should be deleted at some point
	fclose(binout);

	return 0;
}
