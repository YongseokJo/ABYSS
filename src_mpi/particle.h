#ifndef PARTICLE_H
#define PARTICLE_H

#include "def.h"
#include <cmath>

// SDAR
#include "Common/Float.h"
#include <iostream>
#include <iomanip>
#ifdef SEVN
#include "star.h" // Eunwoo added for SEVN
#endif

enum class BinaryInterruptState:int {
	none = 0, 
	form = 1, 
	exchange = 2, 
	collisioncandidate = 3,
	collision = 4,
	manybody = 5,	// many-body (>2) group terminated // added by EW 2025.1.19
	kicked = 6,		// kicked or exploded as PISN during stellar evolution // added by EW 2025.1.19
	terminated = 7,	// CM particle only; terminated // added by EW 2025.1.20
	merger = 8		// CM particle only; merger inside // added by EW 2025.1.20
};
#define BINARY_STATE_ID_SHIFT 4 // Eunwoo
#define BINARY_INTERRUPT_STATE_MASKER 0xF // Eunwoo

struct Group;
struct Particle {
	int PID;
	int ParticleIndex; // (Query) please change it to index when you get a chance // (Answer) done!
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
	double NewPosition[Dim]; // this might be initialized in NewFBInitialization by EW 2025.1.7
	double NewVelocity[Dim]; // this might be initialized in NewFBInitialization by EW 2025.1.7
	double RadiusOfNeighbor; //this should be squared of it. // this might be initialized in NewFBInitialization by EW 2025.1.7
	double BackgroundAcceleration[Dim]; // this might be initialized in NewFBInitialization by EW 2025.1.7

	// For SDAR
	bool isActive;
	bool isUpdateToDate;
	double radius; // used in SEVN too
	double dm; // Stellar mass which will be distributed to nearby gas cells // used in SEVN too
	double time_check; // time to check next interrupt
	long long int binary_state; // contain two parts, low bits (first BINARY_STATE_ID_SHIFT bits) is binary interrupt state and high bits are pair ID
	double a_spin[3]; // dimensionless spin parameter a
	Group* GroupInfo;
	bool isCMptcl; // do we need this? (Query) we can simply use if GroupInfo == nullptr right?
	int CMPtclIndex; // added for write_out_group function by EW 2025.1.6
	int Members[10]; // ParticleIndex of group members; only used for cm ptcls by EW 2025.1.30
	int NumberOfMember; // Number of group members; only used for cm ptcls by EW 2025.1.30

#ifdef SEVN
	// For SEVN
	Star* StellarEvolution;
	double FormationTime; // Myr // for restart
	double WorldTime; // Myr // FormationTime + EvolutionTime
#endif

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
			a_spin[i] = 0.;
		}
		isActive = true;
		ParticleIndex = -1;
		radius = 0.;
		dm = 0.0;
		time_check = NUMERIC_FLOAT_MAX;
		setBinaryInterruptState(BinaryInterruptState::none);
		GroupInfo = nullptr;
		isCMptcl = false; //(Query)
		isUpdateToDate = true;
		CMPtclIndex = -1;
		NumberOfMember = 0;
#ifdef SEVN
		StellarEvolution = nullptr;
		FormationTime = 0.0; // Myr
		WorldTime = 0.0; // Myr
#endif
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
		// this->ParticleType = 0; //NormalStar+SingleParticle;
		//this->NextParticleInEnzo = NextParticleInEnzo;
		this->CurrentTimeReg		= 0;
		this->CurrentTimeIrr		= 0;
		this->RadiusOfNeighbor		= ACRadius*ACRadius;
		//this->RadiusOfNeighbor	= 1;
		//this->RadiusOfNeighbor	= 0;
		this->NumberOfNeighbor		= 0;
		this->NewNumberOfNeighbor	= 0;

		this->isActive				= true;
		this->ParticleIndex			= PID;
		this->dm = 0.0;
		this->time_check = NUMERIC_FLOAT_MAX;
		this->setBinaryInterruptState(BinaryInterruptState::none);
		this->GroupInfo = nullptr;
		this->isCMptcl = false;
		this->a_spin[0] = 0.;
		this->a_spin[1] = 0.;
		this->a_spin[2] = 0.;
		this->CMPtclIndex = -1;
		this->isUpdateToDate = true;
		this->NumberOfMember = 0;

#ifndef SEVN
		this->ParticleType = NormalStar+SingleStar;
		this->radius = 2.25461e-8/position_unit*pow(this->Mass*1e9, 1./3); // stellar radius in code unit
		/*
		if (this->Mass*1e9 > 8) {
			this->ParticleType = Blackhole+SingleStar;
			this->radius = 2*this->Mass*1e9/mass_unit/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzshild radius in code unit
			// initialBHspin(this);
		}
		else {
			this->ParticleType = NormalStar+SingleStar;
			this->radius = 2.25461e-8/position_unit*pow(this->Mass*1e9, 1./3); // stellar radius in code unit
		}
		*/
#endif
	}

	// Clear function // Simplified by EW 2025.1.7
    void clear() {
        PID = -1;
		ParticleIndex = -1;

        NumberOfNeighbor = 0;
        NewNumberOfNeighbor = 0;

		isUpdateToDate = true;
        isActive = false;
		GroupInfo = nullptr;
		isCMptcl = false;
		CMPtclIndex = -1;
		NumberOfMember = 0;
		setBinaryInterruptState(BinaryInterruptState::none);
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

	void predictParticleSecondOrder(double dt, CUDA_REAL pos[], CUDA_REAL vel[]);
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
	void calculateTimeStepIrr2();
	void calculateTimeStepReg();

	void updateRegularParticleCuda(int *NewNeighbors, int NewNumberOfNeighbor, double *new_a, double *new_adot);

	// SDAR
	void checkNewGroup();
	void checkNewGroup2();

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

	double* getPos() {
		return Position;
	}

	double* getVel() {
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

	// made by EW 2025.1.6
	void copyNewNeighbor(Particle* ptcl) {
		this->NewNumberOfNeighbor = ptcl->NewNumberOfNeighbor;
		for (int i = 0; i < ptcl->NewNumberOfNeighbor; i++)
			this->NewNeighbors[i] = ptcl->NewNeighbors[i];
	}
	
};




#endif
