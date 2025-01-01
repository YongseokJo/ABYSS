#ifndef GROUP_H
#define GROUP_H

#include <vector>
#include <iostream>
#include "../def.h"
#include <cassert>

#include "Common/binary_tree.h"
#include "Common/Float.h"
#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "ar_interaction.hpp"
#include "ar_perturber.hpp"
// #define SEVN

struct Particle;
struct Group
{

	Particle* groupCM;
	bool isTerminate; // For later use: I will use this when Binary Interrupt State is being used
	bool isMerger; // True if merger happened

	double CurrentTime;  // this show how much the binary system has evolved

	AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, AR::Information<Particle,Particle>> sym_int;
	AR::TimeTransformedSymplecticManager<Interaction> manager;

#ifdef SEVN
	bool useSEVN;		// true if one of the group members has Star class as its member.
	double EvolutionTime; // Myr
	double EvolutionTimeStep; // Myr
#endif

	// Constructor
#ifdef SEVN
	Group(void) 
		: groupCM(nullptr),
		isTerminate(false),
		isMerger(false),
		CurrentTime(0.0),
		sym_int(),
		manager(),
		useSEVN(false),
		EvolutionTime(0.0),
		EvolutionTimeStep(0.0)
	{}
#else
	Group(void) 
		: groupCM(nullptr),
		isTerminate(false),
		isMerger(false),
		CurrentTime(0.0),
		sym_int(),
		manager()
	{}
#endif

	bool ARIntegration(double next_time);
	bool CheckBreak();
	void initialManager();
	void initialIntegrator(int NumMembers);

	~Group() {

		sym_int.clear(); // Delocate memory
		sym_int.particles.clear(); // Delocate memory
		manager.step.clear(); // Delocate memory
		groupCM = nullptr;
	};

};
#endif
