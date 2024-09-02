# Compiler
#CXX = icc
CXX = icc
#CXX = g++
NVCC = nvcc

# Compiler flags
<<<<<<< HEAD
CXXFLAGS  = -O3 -std=c++11  -Wall -Wextra -Wuninitialized -I./Particle -I./ -I./../SDAR/src -g # Eunwoo edited -fopenmp for sympletic_integrator.h?
NVCCFLAGS = -O2 --std=c++11 -G -Xcompiler -rdynamic -g -lineinfo -Xptxas=-v -lcudart -lcompute-sanitizer
#NVCC_FLAGS_DEBUG = -g -G --device-debug
=======
CXXFLAGS  = -O3 -march=native -ftree-vectorize -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -qopenmp  #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -qopenmp  #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  -fsanitize=thread
#CXXFLAGS  = -O2 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  #-fsanitize=thread
NVCCFLAGS = -O2 --std=c++11 -G -Xcompiler -rdynamic -g -lineinfo -Xptxas=-v -lcudart -lcompute-sanitizer 
NVCC_FLAGS_DEBUG = -g -G --device-debug -lineinfo
>>>>>>> f688615c336c02f68febfc5c6d41376999e08ae3
#NVCC_FLAGS_PROFILE = -g -lineinfo

# SDARFLAGS = -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE # Eunwoo added for SDAR

# Libraries
LDFLAGS = -L/usr/local/cuda-11.8/lib64/ -lcudart -lcublas -lnvToolsExt #-ltbb
#add -lnvToolsExt for <nvToolsExt.h>
#LDFLAGS = -L/mnt/sw/nix/store/zi2wc26znf75csf5hhz77p0d2bbz53ih-cuda-11.8.0/lib64/ -lcudart

<<<<<<< HEAD
# Source files // Eunwoo: Should change this line
CXX_SRCS = main.cpp ReadWrite.cpp Particle/ParticleRoutines.cpp ComputationChainRoutine.cpp SortComputationChain.cpp InitializeParticles.cpp TimeStepRoutines.cpp Particle/CalculateIrrForce.cpp IrregularAccelerationRoutine.cpp Particle/CalculateRegForce.cpp RegularAccelerationRoutine.cpp Particle/UpdateParticle.cpp Particle/CalculateTimeStep.cpp InitializeTimeStep.cpp Evolve.cpp Particle/CheckNeighborForEvolution.cpp Particle/UpdateEvolveParticle.cpp parser.cpp cuda/CalculateRegularAcceleration.cpp cuda/CalculateAccelerationForAll.cpp Particle/EvolveStar.cpp FewBody/FBInitialization.cpp FewBody/FBIntegration.cpp FewBody/FBTermination.cpp GroupAccelerationRoutine.cpp 
=======
# Source files
CXX_SRCS = main.cpp ReadWrite.cpp Particle/ParticleRoutines.cpp ComputationChainRoutine.cpp SortComputationChain.cpp InitializeParticles.cpp TimeStepRoutines.cpp Particle/CalculateIrrForce.cpp IrregularAccelerationRoutine.cpp Particle/CalculateRegForce.cpp RegularAccelerationRoutine.cpp Particle/UpdateParticle.cpp Particle/CalculateTimeStep.cpp InitializeTimeStep.cpp Evolve.cpp Particle/CheckNeighborForEvolution.cpp Particle/UpdateEvolveParticle.cpp parser.cpp cuda/CalculateRegularAcceleration.cpp cuda/CalculateAccelerationForAll.cpp Particle/EvolveStar.cpp Binary/BinaryRoutines.cpp Binary/KSRegularlizationInitialization.cpp Binary/KSRegularlizationIntegration.cpp Binary/KSRegularlizationTermination.cpp BinaryAccelerationRoutine.cpp  IrregularAccelerationRoutineParallel.cpp ParticleScheduler/ParticleSchedulerRoutines.cpp
>>>>>>> f688615c336c02f68febfc5c6d41376999e08ae3
#cuda/cuda_routines.cpp


CU_SRCS = cuda/cuda_my_acceleration.cu cuda/cuda_kernels.cu 
# cuda/cuda_routines.cu

#
OBJS = $(CU_SRCS:.cu=.o) $(CXX_SRCS:.cpp=.o)

# Executable name
EXEC = nbodyplus_binary.exe

# Main target
all: $(EXEC)

# Link object files into the executable
$(EXEC): $(OBJS)
	#$(NVCC) $(NVCCFLAGS) $(OBJS) -o $@ $(LDFLAGS)
	#$(NVCC) $(CXXFLAGS) $(NVCCFLAGS) $(OBJS) -o $(EXEC) -v
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ -v $(LDFLAGS)
	#$(NVCC) $(CXXFLAGS) $(NVCCFLAGS) $(OBJS) -o $@ -v

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@


# Clean the project directory
clean:
	rm -f $(OBJS) $(EXEC)


