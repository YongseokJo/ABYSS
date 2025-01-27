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
	int    Neighbors[MaxNumberOfNeighbor];    // 
	int    NumberOfNeighbor;    // 
	int    NewNeighbors[MaxNumberOfNeighbor];    // 
	int    NewNumberOfNeighbor;    // 

	double CurrentTimeIrr;
	double CurrentTimeReg;
	ULL CurrentBlockIrr;
	ULL NewCurrentBlockIrr;
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
	double RadiusOfNeighbor; //this should be squared of it.
	double BackgroundAcceleration[Dim];



	Particle() {
		Position[0] = Position[1] = Position[2] = 0.0;
		Velocity[0] = Velocity[1] = Velocity[2] = 0.0;
		PID             = -1;
		Mass            = 0;
		RadiusOfNeighbor= -1;
		NumberOfNeighbor= 0;
		NewNumberOfNeighbor= 0;
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

	/*
	Particle(double x, double y, double z, double vx, double vy, double vz, double m, double q)
		: Mass(m) {
			Position[0] = x; Position[1] = y; Position[2] = z;
			Velocity[0] = vx; Velocity[1] = vy; Velocity[2] = vz;
		}
		*/

	void initialize(double *data, int PID) {
		this->PID          = PID;
		this->Position[0]  = data[0];
		this->Position[1]  = data[1];
		this->Position[2]  = data[2];
		this->Velocity[0]  = data[3];
		this->Velocity[1]  = data[4];
		this->Velocity[2]  = data[5];
		this->Mass         = data[6];
		this->ParticleType = 0; //NormalStar+SingleParticle;
		//this->NextParticleInEnzo = NextParticleInEnzo;
		this->CurrentTimeReg             = 0;
		this->CurrentTimeIrr             = 0;
		this->RadiusOfNeighbor           = 0.1*0.1;
		//this->RadiusOfNeighbor         = 1;
		//this->RadiusOfNeighbor         = 0;
		this->NumberOfNeighbor           = 0;
		this->NewNumberOfNeighbor        = 0;
	}

	void normalizeParticle() {
		// pc to computing unit, km/s to computing unit
		this->Mass *= 1e9;
		this->Mass /= mass_unit;
		for (int dim=0; dim<Dim; dim++) {
			this->Position[dim] *= 1000; // kpc to pc
			this->Position[dim] /= position_unit;
			this->Velocity[dim] *= 1e5*yr/pc; // km/s to pc/yr
			this->Velocity[dim] /= velocity_unit;
		}
	}



	void updateParticle(); 

	/*
	void getAcceleration(const double pos[], const double vel[]) {
		double dx[Dim], dr2;
		for (int dim=0; dim<Dim; dim++) {
			dx[dim] = pos[dim] - this->Position[dim];
			dr2 += dx[dim]*dx[dim];
		}

		for (int dim=0; dim<Dim; dim++) {
			a_reg[dim] += dx[dim]/dr2/sqrt(dr2);
		}
	}
	*/

	void predictParticleSecondOrder(double dt, double pos[], double vel[]);
	void correctParticleFourthOrder(double dt, double pos[], double vel[], double a[3][4]);


	/*
	void update_timestep() {
		double acc = mag(acceleration);
		double vel = mag(Velocity);

		time_step = eta*sqrt(std::abs(vel/acc));
	}
	*/

	void updateRadius();

	//void initializeNeighbor();
	//void initializeAcceleration();
	void initializeTimeStep();

	void computeAccelerationIrr();
	void computeAccelerationReg();


	void calculateTimeStepIrr();
	void calculateTimeStepReg();

};




#endif
