#pragma once

#include <vector>
#include <iostream>
#include "../defs.h"
#include <cassert>

#include "FB_defs.h"
#include "Common/binary_tree.h"
#include "Common/Float.h"
#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "interaction.h"
#include "perturber.h"

class Particle;
class Group
{
	private:
	public:

		std::vector<Particle*> Members;     // list of Group members

		Particle* groupCM;
		bool isTerminate;
		bool isErase;

		REAL PredTime;
		REAL CurrentTime;  // this show how much the binary system has evolved
		REAL TimeStep;
		int TimeLevel;

		AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, AR::Information<Particle,Particle>> sym_int;
		AR::TimeTransformedSymplecticManager<Interaction> manager;

		// Constructor
		Group(void) 
			: groupCM(nullptr),
			isTerminate(false),
			isErase(false),
			PredTime(0),
			CurrentTime(0),
			TimeStep(0),
			TimeLevel(0),
			sym_int(),
			manager()

		{
			Members.clear();
		}

		// void InitializeGroup(REAL current_time);
		// void getStumpffCoefficients(REAL z);
		void ARIntegration(REAL next_time);
		void predictGroup(REAL next_time);
		void initialManager();
		void initialIntegrator();

		~Group() {
            groupCM = nullptr;
            Members.clear();
            Members.shrink_to_fit();

			sym_int.clear(); // Delocate memory
			sym_int.particles.clear(); // Delocate memory
			manager.step.clear(); // Delocate memory
        };

};