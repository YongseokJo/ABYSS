/*
*  Purporse: Particle Class
*
*  Date    : 2024.01.02  by Yongseok Jo
*
*/
#ifndef PARTICLE_H
#define PARTICLE_H

#pragma once


#include <vector>
#include <iostream>
#include "../defs.h"
#include <iomanip> // Eunwoo
#include "Common/Float.h" // Eunwoo

enum class BinaryInterruptState:int {none = 0, form = 1, exchange = 2, collision = 3}; // Eunwoo
#define BINARY_STATE_ID_SHIFT 4 // Eunwoo
#define BINARY_INTERRUPT_STATE_MASKER 0xF // Eunwoo

class Group;
class Particle
{
	private:

	public:
		int ParticleOrder;
		int PID;
		int ParticleType;
		REAL Mass;
		REAL InitialMass;
		REAL CreationTime;
		REAL DynamicalTime;
		REAL Velocity[Dim];
		REAL Position[Dim];
		REAL NewVelocity[Dim];
		REAL NewPosition[Dim];
		REAL PredTime;
		REAL PredTimeIrr;
		REAL PredTimeReg;
		REAL CurrentTimeIrr;
		REAL CurrentTimeReg;
		ULL CurrentBlockIrr;
		ULL NextBlockIrr;
		ULL CurrentBlockReg;
		REAL TimeStepIrr;
		REAL TimeStepReg;
		ULL TimeBlockIrr;
		ULL TimeBlockReg;
		int TimeLevelIrr;
		int TimeLevelReg;
		REAL PredPosition[Dim];
		REAL PredVelocity[Dim];
		REAL a_tot[Dim][HERMITE_ORDER];
		REAL a_reg[Dim][HERMITE_ORDER];
		REAL a_irr[Dim][HERMITE_ORDER];
		REAL BackgroundAcceleration[Dim];
		REAL LocalDensity;
		Particle* NextParticleInEnzo;
		Particle* NextParticleForComputation;
		std::vector<Particle*> ACList;     // list of AC neighbor 
		int NumberOfAC; // number of neighbors
		REAL RadiusOfAC;
		bool isStarEvolution;
		bool isCMptcl; // check if this particle is center-of-mass particle
		bool isErase;

		// Eunwoo added for SDAR

		bool isGroup; // Eunwoo, check whether this is a member of the group
		REAL radius;
		REAL dm;
		REAL time_check; // time to check next interrupt
		long long int binary_state; // contain two parts, low bits (first BINARY_STATE_ID_SHIFT bits) is binary interrupt state and high bits are pair ID
		Group* GroupInfo;
		REAL a_spin[3]; // dimensionless spin parameter a

		// Constructor
		Particle(void) {__initializer__();};
		Particle(REAL *data, int PID);
		void __initializer__(void) {
			//std::cout << "Constructor called" << std::endl;
			PID             = -1;
			ParticleOrder   = -1;
			Mass            = 0;
			InitialMass     = 0;
			NumberOfAC      = 0; // number of neighbors
			RadiusOfAC      = -1;
			ParticleType    = -9999;
			CurrentTimeIrr  = 0.; // consistent with actual current time
			CurrentTimeReg  = 0.;
			CurrentBlockIrr = 0; // consistent with actual current time
			CurrentBlockReg = 0;	
			NextBlockIrr    = 0;
			PredTimeIrr     = 0;
			PredTimeReg     = 0;
			TimeStepIrr     = 0;
			TimeStepReg     = 0;
			TimeLevelIrr    = 0;
			TimeLevelReg    = 0;
			TimeBlockIrr    = 0;
			TimeBlockReg    = 0;
			LocalDensity    = 0;
			isStarEvolution = true;
			// isBinary        = false; // Eunwoo deleted
			isCMptcl        = false;
			isErase         = false;
			for (int i=0; i<Dim; i++) {
				Velocity[i]     = 0.;
				Position[i]     = 0.;
				PredPosition[i] = 0.;
				PredVelocity[i] = 0.;
				BackgroundAcceleration[i] = 0.;
				for (int j=0; j<HERMITE_ORDER; j++) {
					a_tot[i][j] = 0.;
					a_reg[i][j] = 0.;
					a_irr[i][j] = 0.;
				}
			}
			NextParticleInEnzo         = nullptr;
			NextParticleForComputation = nullptr;

			ACList.clear();

			// Eunwoo added for SDAR

			radius			= 0.0;
			dm				= 0.0;
			binary_state	= 0;
			time_check		= NUMERIC_FLOAT_MAX;
			isGroup			= false;
			GroupInfo		= nullptr;
			for (int i=0; i<3; i++)
				a_spin[i]	= 0.;

		};

		void updateParticle(REAL mass, REAL *vel, REAL pos[], int particletype) {

			Mass = mass;
			ParticleType = particletype;

			for (int i=0; i<Dim; i++) {
				Velocity[i] = vel[i];
				Position[i] = pos[i];
			}
		};
		REAL getDistanceTo(Particle *particle);
		void setParticleInfo(REAL *data, int PID);
		void setParticleInfo(REAL *data, int PID, Particle* NextParticleInEnzo);
		void setParticleInfo(int *PID, REAL *Mass, REAL *CreationTime, REAL *DynamicalTime,
		REAL *Position[Dim], REAL *Velocity[Dim],
		REAL *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i);
		void setParticleInfo(int *PID, REAL *Mass, REAL *CreationTime, REAL *DynamicalTime,
		REAL *Position[Dim], REAL *Velocity[Dim],
		REAL *BackgroundAcceleration[Dim],  int ParticleType, Particle* NextParticleInEnzo, int i);
		void setParticleInfo(int *PID, REAL *BackgroundAcceleration[Dim],
		Particle* NextParticleInEnzo, int i);
		void setParticleInfo(REAL *Mass, REAL *BackgroundAcceleration[Dim],
		Particle* NextParticleInEnzo, int i);
		void initializeTimeStep();
		int getPID() {return PID;};
		void calculateIrrForce();
		void calculateRegAccelerationSecondOrder(std::vector<Particle*> &particle);
		void calculateRegAccelerationFourthOrder(std::vector<Particle*> &particle);

		void predictParticleSecondOrder(REAL time);
		void predictParticleSecondOrderIrr(REAL time);
		void correctParticleFourthOrder(REAL current_time, REAL next_time, REAL a[3][4]);

		void normalizeParticle();
		void calculateTimeStepIrr(REAL f[3][4], REAL df[3][4]);
		void calculateTimeStepIrr2(REAL f[3][4], REAL df[3][4]); // Eunwoo add
		void calculateTimeStepReg();
		void calculateTimeStepReg2();
		bool checkNeighborForEvolution();
		void updateEvolveParticle(std::vector<Particle*> &particle);
		void updateParticle();
		REAL evolveStarMass(REAL t1, REAL t2);
		// void isKSCandidate(); // Eunwoo deleted
		// void convertBinaryCoordinatesToCartesian(); // Eunwoo deleted
		void polynomialPrediction(REAL current_time);
		void UpdateRadius();
		void UpdateNeighbor(std::vector<Particle*> &particle);
		void isFBCandidate(); // Eunwoo
		void checkNewGroup(); // Eunwoo
		void checkNewGroup2(); // Eunwoo

		// Eunwoo added for SDAR

		//! save pair id in binary_state with shift bit size of BINARY_STATE_ID_SHIFT
		void setBinaryPairID(const int _id) {
			binary_state = (binary_state&BINARY_INTERRUPT_STATE_MASKER) | (_id<<BINARY_STATE_ID_SHIFT);
		}

		//! save binary interrupt state in the first  BINARY_STATE_ID_SHIFT bit in binary_state
		void setBinaryInterruptState(const BinaryInterruptState _state) {
			binary_state = ((binary_state>>BINARY_STATE_ID_SHIFT)<<BINARY_STATE_ID_SHIFT) | int(_state);
		}

		//! get binary interrupt state from binary_state
		BinaryInterruptState getBinaryInterruptState() const {
			return static_cast<BinaryInterruptState>(binary_state&BINARY_INTERRUPT_STATE_MASKER);
		}

		//! get pair ID from binary_state 
		int getBinaryPairID() const {
			return (binary_state>>BINARY_STATE_ID_SHIFT);
		}

		REAL* getPos() {
			return Position;
		}

		REAL* getVel() {
			return Velocity;
		}

		static void printColumnTitle(std::ostream & _fout, const int _width=20) {
			_fout<<std::setw(_width)<<"mass"
				<<std::setw(_width)<<"pos.x"
				<<std::setw(_width)<<"pos.y"
				<<std::setw(_width)<<"pos.z"
				<<std::setw(_width)<<"vel.x"
				<<std::setw(_width)<<"vel.y"
				<<std::setw(_width)<<"vel.z"
				<<std::setw(_width)<<"radius"
				<<std::setw(_width)<<"id";
		}

		//! print data of class members using column style (required)
		/*! print data of class members in one line for column style. Notice no newline is printed at the end
		@param[out] _fout: std::ostream output object
		@param[in] _width: print width (defaulted 20)
		*/

		void printColumn(std::ostream & _fout, const int _width=20){
			_fout<<std::setw(_width)<<Mass
				<<std::setw(_width)<<Position[0]
				<<std::setw(_width)<<Position[1]
				<<std::setw(_width)<<Position[2]
				<<std::setw(_width)<<Velocity[0]
				<<std::setw(_width)<<Velocity[1]
				<<std::setw(_width)<<Velocity[2]
				<<std::setw(_width)<<radius
				<<std::setw(_width)<<PID;
		}


		~Particle() {
			ACList.clear();
			// Deallocate memory
			ACList.shrink_to_fit();
			NextParticleInEnzo         = nullptr;
			NextParticleForComputation = nullptr;
			GroupInfo					= nullptr;
			// fprintf(stderr, "deleting particle, pid=%d\n", PID); // Eunwoo deleted for debug
		};
};


typedef std::vector<Particle*> Pvector;
#endif
