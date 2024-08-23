// #define ASSERT(x) assert(x)

#include <vector>
#include <iostream>
#include "../defs.h"
#include <cassert>

// #include "AR/symplectic_integrator.h"
// #include "interaction.h"
// #include "perturber.h"

class Particle;
class Group
{
	private:


	public:

		// corresponding single particles and neighbor particles
		// overlapping variables with ptclCM are deleted. 
		// int NumberOfAC;     // number of neighbors - not currently in use
		// REAL RadiusOfAC;
		std::vector<Particle*> Members;     // list of Group members
		// std::vector<Particle> Members;

		Particle* groupCM;
		bool isTerminate;
		bool isErase;

		// information of binary particles in cartesian coordinates
		// REAL Mass;

		// REAL Position[Dim];
		// REAL Velocity[Dim];
		// REAL PredPosition[Dim];
		// REAL PredVelocity[Dim];
		REAL PredTime;
		REAL CurrentTime;  // this show how much the binary system has evolved
		REAL TimeStep;
		int TimeLevel;

		// REAL a_tot[Dim][HERMITE_ORDER];
		// REAL a_reg[Dim][HERMITE_ORDER];
		// REAL a_irr[Dim][HERMITE_ORDER];
		// REAL BackgroundAcceleration[Dim];



		// Constructor
		Group(void) {

			//std::cout << "Constructor called" << std::endl;

            groupCM     = nullptr;

			isTerminate = false;
			isErase     = false;
			// NumberOfAC     = 0; // number of neighbors
			// RadiusOfAC     = InitialRadiusOfAC;
			// Mass           = 0;
			PredTime       = 0;
			CurrentTime    = 0;
			TimeStep       = 0;
			TimeLevel      = 0;

            Members.clear();


			// for (int i=0; i<Dim; i++) {

			// 	Position[i]     = 0;
			// 	Velocity[i]     = 0;
			// 	PredPosition[i] = 0;
			// 	PredVelocity[i] = 0;

			// 	BackgroundAcceleration[i] = 0;

			// 	for (int j=0; j<HERMITE_ORDER; j++) {
			// 		a_reg[i][j] = 0;
			// 		a_irr[i][j] = 0;
			// 		a_tot[i][j] = 0;
			// 	}
			// }  -> it is all saved in the CM particle information

		}

		// void InitializeGroup(REAL current_time);
		// void getStumpffCoefficients(REAL z);
		void ARIntegration(REAL next_time);
		void predictGroup(REAL next_time);
		void IntegrateGroup(REAL next_time);

		~Group() {
            groupCM = nullptr;
            Members.clear();
            Members.shrink_to_fit();
        };

};
