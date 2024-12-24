//#define CUDA

// #define LoadBalance
#define PerformanceTrace

#define FAIL    -1
#define SUCCESS  1


#define NumberOfTask 20
#define MaxNumberOfParticle 1000000
#define MaxNumberOfNeighbor 10000



#define MaxNumberOfCommunication 10000


//#define FixNumNeighbor 20
#define FixNumNeighbor 500
//#define NumNeighborMax 100
#define NumNeighborMax 5000


// SDAR
#define rbin 0.25e-3 // pc
#define tbin 1e-6 // Myr
#define NormalStar 1
#define Blackhole 32
#define SingleParticle 64


#define MIN_LEVEL_BUFFER 30

#define Dim 3
#define eta 0.01
#define HERMITE_ORDER 4


typedef unsigned long long ULL;



#define mag(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define mag0(a) (a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0])
#define dist(a,b) std::sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))



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
// Physical constants
// #define G_cgs 6.67430e-8 // Eunwoo: crash with SEVN
// #define G // pc, yr, Msun // Eunwoo: crash with SEVN

#define RCAST(a)  static_cast<double>(a)
#define ABS(a) static_cast<double>(std::abs(a))
#define MIN(a,b) std::min(RCAST(a),RCAST(b))

#define CUDA_REAL double
#define nbodymax 100000000 //100000000 for node14
#define NSIGHT // for nsight debugging
#define BatchSize 64 // each thread calculates BatchSize particles
#define GridDimY 32 // each block calcuates NNB/GridDimY particles
#define NNB_per_block 256
//#define BatchSize 32 // each thread calculates BatchSize particles
//#define GridDimY 16 // each block calcuates NNB/GridDimY particles
//#define NNB_per_block 128
