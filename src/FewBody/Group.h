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
#include "ar_interaction.hpp"
#include "ar_perturber.hpp"

class Particle;
class Group
{
	private:
	public:

		std::vector<Particle*> Members;     // list of Group members

		Particle* groupCM;
		bool isTerminate; // For later use: I will use this when Binary Interrupt State is being used
		bool isErase;

		REAL CurrentTime;  // this show how much the binary system has evolved

		AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, AR::Information<Particle,Particle>> sym_int;
		AR::TimeTransformedSymplecticManager<Interaction> manager;

		bool PNon;

		// Constructor
		Group(void) 
			: groupCM(nullptr),
			isTerminate(false),
			isErase(false),
			CurrentTime(0),
			sym_int(),
			manager(),
			PNon(false)

		{
			Members.clear();
		}

		bool ARIntegration(REAL next_time, std::vector<Particle*> &particle);
		bool CheckBreak();
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