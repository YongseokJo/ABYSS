

#define FixNumNeighbor 30
#define HERMITE_ORDER 4
#define Dim 3

#define UpdateFrequency 10



// Particle Type Definition
// Particle type defined by combination of binary switches
// fictious (2^8), new (2^7), single (2^6), black hole (2^5), dead star (2^4)
// TypeIII (2^3), TypeII (2^2), TypeI (2^1), normal star or particle (2^0) 
// for example, normal signle particle is 001000001 = 2^6 + 2^0 = 65
// new black hole single particle is 011100000 = 2^7+2^6+2^5 = 
// fictious particles are usually COM particles for regularization
#define NormalStar 1
#define TypeIStar 2
#define TypeIIStar 4
#define TypeIIIStar 8
#define DeadStar 16
#define Blackhole 32
#define SingleParticle 64
#define NewParticle 128
#define FicticiousParticle 256


#define DONE    1
#define FAIL   -1


// numerical values
#define mytolerance 5.4210109e-20

// Physics related parameters
#define eta 0.01 // by YS Jo subject to modifty
#define EPS2 1e-4 // smoothing length
#define InitialRadiusOfAC 0.04 // 0.04 pc
//#define InitialRadiusOfAC 1000. // 0.04 pc
#define MIN_LEVEL_BUFFER 30

// Physical constants
#define G_cgs 6.67430e-8
#define G // pc, yr, Msun

// KS regularlization related variables

#define NumNeighborMax 100
#define stumpffN 12 // the order of approximation for stumpff functions
#define KSDistance 0.000002  // the distance of binary particles from each other
//#define KSDistance 1e-10  // the distance of binary particles from each other
#define KSTime 2e-7  // short timestep criteria for binary search
//#define KSTime 1e-10  // short timestep criteria for binary search
#define PerturberDistance 0.01  // the distance of binary perturbers from the pair

// Physical units in cgs
#define pc 3.08567758149137e18
#define yr 3.1536e7
#define Msun 1.98847e33

// Code unit in pc, yr, Msun
#define time_unit 1e10 // in 1e10 yr
#define position_unit 4. // in 4 pc
#define velocity_unit 4e-10 // in 4e-10 pc/yr
//#define mass_unit 256e-20  // 256e-20 Msun in the unit that G = 1.
#define mass_unit 0.0001424198  // Msun in the unit that G = 1.

#define no_time_trace

#define mag(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define mag0(a) (a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0])
#define dist(a,b) std::sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))


typedef unsigned long long ULL;

//#define FLOAT
#define DOUBLE
#ifdef FLOAT
typedef float REAL;
#define RN f
#define RCAST(a)  static_cast<float>(a)
#define ABS(a) static_cast<float>(std::abs(a))
#define MIN(a,b) std::min(RCAST(a),RCAST(b))
#else 
typedef double REAL;
#define RN d
#define RCAST(a)  static_cast<double>(a)
#define ABS(a) static_cast<double>(std::abs(a))
#define MIN(a,b) std::min(RCAST(a),RCAST(b))
#endif

//typedef float CUDA_REAL;
typedef double CUDA_REAL;



#define NAN_CHECK(val) assert((val) == (val));

#define NSIGHT // for nsight debugging
