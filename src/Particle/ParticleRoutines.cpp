#include "../global.h"
#include <cmath>

void generate_Matrix(REAL a[3], REAL (&A)[3][4]);

REAL Particle::getDistanceTo(Particle *particle) {
	REAL d=0; 
	for (int i=0; i<Dim; i++)
		d += std::pow(this->Position[i]-particle->Position[i],2);
	return std::sqrt(d);
}

// Eunwoo: Set initial spin for BH. (Geneva model, MESA model, Gaussian distribution, etc.)
void initialBHspin(Particle *particle) {

	// 1. Aligned spin with the same magnitude
	// /*
	for (int i=0; i<2; i++)
		particle->a_spin[i] = 0.0;
	particle->a_spin[2] = 0.5;
	// */

	//2. Random spin generator
	/*
	std::random_device rd; // Obtain a random number from hardware
	std::mt19937 mt(rd()); // Seed the generator
	std::uniform_real_distribution<> distr(0.0, 1.0); // Define the range (0 to 1)
	while (1) {
		for (int i=0; i<2; i++)
			particle->a_spin[i] = (2 * distr(mt) - 1) * 0.8;
		if (mag(particle->a_spin) < 0.8*0.8)
			break;
	}
	*/
}

Particle::Particle(REAL *data, int PID) {
	__initializer__();
	this->PID          = PID;
	this->Position[0]  = data[0];
	this->Position[1]  = data[1];
	this->Position[2]  = data[2];
	this->Velocity[0]  = data[3];
	this->Velocity[1]  = data[4];
	this->Velocity[2]  = data[5];
	this->Mass         = data[6];
	// this->ParticleType = NormalStar+SingleParticle;
	// /* // Eunwoo added
	if (this->Mass*1e9 > 8) {
		this->ParticleType = Blackhole+SingleParticle;
		this->radius = 6*this->Mass*1e9/mass_unit/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // innermost stable circular orbit around a Schwartzshild BH = 3 * R_sch
		initialBHspin(this);
	}
	else {
		this->ParticleType = NormalStar+SingleParticle;
		this->radius = 2.25461e-8/position_unit*pow(this->Mass*1e9, 0.8); // stellar radius in code unit
	}
	// */ // Eunwoo added
	//this->NextParticleInEnzo = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
}

void Particle::setParticleInfo(REAL *data, int PID) {
	this->PID          = PID;
	this->Position[0]  = data[0];
	this->Position[1]  = data[1];
	this->Position[2]  = data[2];
	this->Velocity[0]  = data[3];
	this->Velocity[1]  = data[4];
	this->Velocity[2]  = data[5];
	this->Mass         = data[6];
	this->ParticleType = NormalStar+SingleParticle;
	//this->NextParticleInEnzo = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
}

/*
void Particle::setParticleInfo(REAL *data, int PID, Particle* NextParticleInEnzo) {
	this->PID          = PID;
	this->Position[0]  = data[0];
	this->Position[1]  = data[1];
	this->Position[2]  = data[2];
	this->Velocity[0]  = data[3];
	this->Velocity[1]  = data[4];
	this->Velocity[2]  = data[5];
	this->Mass         = data[6];
	this->ParticleType = NormalStar+SingleParticle;
	this->NextParticleInEnzo = NextParticleInEnzo;
}

void Particle::setParticleInfo(int *PID, REAL *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i) {
	this->PID                        = PID[i];
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->NextParticleInEnzo = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
}

void Particle::setParticleInfo(REAL *Mass, REAL *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i) {
	this->Mass                       = Mass[i]*EnzoMass;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->NextParticleInEnzo = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
}


void Particle::setParticleInfo(int *PID, REAL *Mass, REAL *CreationTime, REAL *DynamicalTime,
	 	REAL *Position[Dim], REAL *Velocity[Dim], REAL *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i) {
	this->PID                        = PID[i];
	this->Mass                       = Mass[i]*EnzoMass;
	this->InitialMass                = this->Mass;
	this->CreationTime               = CreationTime[i]*EnzoTime;
	this->DynamicalTime              = DynamicalTime[i]*EnzoTime;
	this->Position[0]                = Position[0][i]*EnzoLength;
	this->Position[1]                = Position[1][i]*EnzoLength;
	this->Position[2]                = Position[2][i]*EnzoLength;
	this->Velocity[0]                = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]                = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]                = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->ParticleType               = NormalStar+SingleParticle;
	this->NextParticleInEnzo         = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
}


void Particle::setParticleInfo(int *PID, REAL *Mass, REAL *CreationTime, REAL *DynamicalTime, REAL *Position[Dim],
		REAL *Velocity[Dim], REAL *BackgroundAcceleration[Dim], int ParticleType, Particle* NextParticleInEnzo, int i) {
	this->PID          = PID[i];
	this->Mass         = Mass[i]*EnzoMass;
	this->InitialMass  = this->Mass;
	this->Position[0]  = Position[0][i]*EnzoLength;
	this->Position[1]  = Position[1][i]*EnzoLength;
	this->Position[2]  = Position[2][i]*EnzoLength;
	this->Velocity[0]  = Velocity[0][i]*EnzoVelocity;
	this->Velocity[1]  = Velocity[1][i]*EnzoVelocity;
	this->Velocity[2]  = Velocity[2][i]*EnzoVelocity;
	this->BackgroundAcceleration[0]  = BackgroundAcceleration[0][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[1]  = BackgroundAcceleration[1][i]\
																		 *EnzoAcceleration;
	this->BackgroundAcceleration[2]  = BackgroundAcceleration[2][i]\
																		 *EnzoAcceleration;
	this->ParticleType       = ParticleType;
	this->NextParticleInEnzo = NextParticleInEnzo;
	this->CurrentTimeReg             = 0;
	this->CurrentTimeIrr             = 0;
}
*/


void Particle::normalizeParticle() {
	// pc to computing unit, km/s to computing unit
	Mass *= 1e9;
	Mass /= mass_unit;
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] *= 1000; // kpc to pc
		Position[dim] /= position_unit;
		Velocity[dim] *= 1e5*yr/pc; // km/s to pc/yr
		Velocity[dim] /= velocity_unit;
	}
}

// Eunwoo deleted

/*
void Particle::convertBinaryCoordinatesToCartesian() {
	if (!this->isCMptcl)  {
		fprintf(stderr,"This is NOT a CM particle!\n");
		return;
	}
	fprintf(stdout,"Converting the KS coordinates to physical coordinates of ptclI and ptclJ\n");
	Binary* ptclBin = BinaryInfo;
	Particle* ptclI = BinaryParticleI; 
	Particle* ptclJ = BinaryParticleJ;

	REAL R[Dim], Rdot[Dim];
	REAL Rinv;
	REAL ratioM;
	REAL L[3][4];

	// update the values of positions of ptclI and ptcl J
	R[0]   = ptclBin->u[0]*ptclBin->u[0] - ptclBin->u[1]*ptclBin->u[1] - ptclBin->u[2]*ptclBin->u[2] + ptclBin->u[3]*ptclBin->u[3];
	R[1]   = 2*(ptclBin->u[0]*ptclBin->u[1] - ptclBin->u[2]*ptclBin->u[3]);
	R[2]   = 2*(ptclBin->u[0]*ptclBin->u[2] + ptclBin->u[1]*ptclBin->u[3]);
	ratioM = ptclJ->Mass/this->Mass;

	for (int dim=0; dim<Dim; dim++) {
		ptclI->Position[dim] = this->Position[dim] + ratioM*R[dim];
		ptclJ->Position[dim] = this->Position[dim] - R[dim];
	}


	// do the same thing for velocity components
	generate_Matrix(ptclBin->u,L);

	Rinv = 1/(ptclBin->u[0]*ptclBin->u[0] + ptclBin->u[1]*ptclBin->u[1] + ptclBin->u[2]*ptclBin->u[2] + ptclBin->u[3]*ptclBin->u[3]) ;


	for (int dim=0; dim<Dim; dim++) {
		Rdot[dim] = 0.0;
		for (int dimu=0; dimu<4; dimu++) {
			Rdot[dim] += 2*L[dim][dimu]*ptclBin->udot[dim]*Rinv;
		}
	}


	for (int dim=0; dim<Dim; dim++) {
		ptclI->Velocity[dim] = this->Velocity[dim] + ratioM*Rdot[dim];
		ptclJ->Velocity[dim] = ptclI->Velocity[dim] - Rdot[dim];
	}

	ptclI->CurrentBlockIrr = this->CurrentBlockIrr;
	ptclI->CurrentBlockReg = this->CurrentBlockReg;
	ptclI->CurrentTimeIrr = this->CurrentBlockIrr*time_step;
	ptclI->CurrentTimeReg = this->CurrentBlockReg*time_step;

	ptclI->TimeStepIrr     = this->TimeStepIrr;
	ptclI->TimeBlockIrr    = this->TimeBlockIrr;
	ptclI->TimeLevelIrr    = this->TimeLevelIrr;

	ptclI->TimeStepReg     = this->TimeStepReg;
	ptclI->TimeBlockReg    = this->TimeBlockReg;
	ptclI->TimeLevelReg    = this->TimeLevelReg;

	ptclJ->CurrentBlockIrr = ptclI->CurrentBlockIrr;
	ptclJ->CurrentBlockReg = ptclI->CurrentBlockReg;
	ptclJ->CurrentTimeIrr = ptclI->CurrentTimeIrr;
	ptclJ->CurrentTimeReg = ptclI->CurrentTimeReg;

	ptclJ->TimeStepIrr     = this->TimeStepIrr;
	ptclJ->TimeBlockIrr    = this->TimeBlockIrr;
	ptclJ->TimeLevelIrr    = this->TimeLevelIrr;

	ptclJ->TimeStepReg     = this->TimeStepReg;
	ptclJ->TimeBlockReg    = this->TimeBlockReg;
	ptclJ->TimeLevelReg    = this->TimeLevelReg;

	fprintf(stdout,"END CONVERTING THE COORDINATES\n \n");
}

*/

