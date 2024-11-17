#ifndef PARTICLE_H
#define PARTICLE_H

#include "def.h"
#include <cmath>

struct Particle {
	int PID;
	int ParticleType;
	double Position[3];      // Position in 3D space (x, y, z)
	double Velocity[3];      // Velocity in 3D space (vx, vy, vz)
	double Mass;             // Mass of the particle
	double a_tot[Dim][HERMITE_ORDER];
	double a_reg[Dim][HERMITE_ORDER];
	double a_irr[Dim][HERMITE_ORDER];
	int    Neighbors[MaxNumberOfNeighbor];    // Velocity in 3D space (vx, vy, vz)
	int    NumberOfNeighbor;    // Velocity in 3D space (vx, vy, vz)

	double CurrentTimeIrr;
	double CurrentTimeReg;
	ULL CurrentBlockIrr;
	ULL CurrentBlockReg;
	ULL NextBlockIrr;
	double TimeStepIrr;
	double TimeStepReg;
	ULL TimeBlockIrr;
	ULL TimeBlockReg;
	int TimeLevelIrr;
	int TimeLevelReg;
	double NewPosition[Dim];
	double NewVelocity[Dim];
	double BackgroundAcceleration[Dim];
	double RadiusOfNeighbor; //this should be squared of it.



	Particle() : Mass(1.0) {
		Position[0] = Position[1] = Position[2] = 0.0;
		Velocity[0] = Velocity[1] = Velocity[2] = 0.0;
		PID             = -1;
		Mass            = 0;
		RadiusOfNeighbor= -1;
		ParticleType    = -9999;
		CurrentTimeIrr  = 0.; // consistent with actual current time
		CurrentTimeReg  = 0.;
		CurrentBlockIrr = 0; // consistent with actual current time
		CurrentBlockReg = 0;
		NextBlockIrr    = 0;
		TimeStepIrr     = 0;
		TimeStepReg     = 0;
		TimeLevelIrr    = 0;
		TimeLevelReg    = 0;
		TimeBlockIrr    = 0;
		TimeBlockReg    = 0;
		for (int i=0; i<Dim; i++) {
			Velocity[i]     = 0.;
			Position[i]     = 0.;
			NewPosition[i] = 0.;
			NewVelocity[i] = 0.;
			BackgroundAcceleration[i] = 0.;
			for (int j=0; j<HERMITE_ORDER; j++) {
				a_tot[i][j] = 0.;
				a_reg[i][j] = 0.;
				a_irr[i][j] = 0.;
			}
		}
	}

	Particle(double x, double y, double z, double vx, double vy, double vz, double m, double q)
		: Mass(m) {
			Position[0] = x; Position[1] = y; Position[2] = z;
			Velocity[0] = vx; Velocity[1] = vy; Velocity[2] = vz;
		}

	void updateParticle(); 

	void getAcceleration(const double pos[], const double vel[]) {
		double dx[Dim], dr2;
		for (int dim=0; dim<Dim; dim++) {
			dx[dim] = pos[dim] - this->Position[dim];
			dr2 += dx[dim]*dx[dim];
		}

		for (int dim=0; dim<Dim; dim++) {
			a_[dim] += dx[dim]/dr2/sqrt(dr2);
		}
	}

	void predictParticleSecondOrder(double dt, double pos[], double vel[]);
	void correctParticleFourthOrder(double dt, double next_time, double pos[], double vel[], double a[3][4]);


	void update_timestep() {
		double acc = mag(acceleration);
		double vel = mag(Velocity);

		time_step = eta*sqrt(std::abs(vel/acc));
	}

	void updateRadius();

	//void initializeNeighbor();
	void initializeAcceleration();
	void initializeTimeStep();

};




#endif
