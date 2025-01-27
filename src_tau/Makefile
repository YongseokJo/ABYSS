# Compiler
#CXX = icc
#CXX = mpic++
#CXX = /mnt/home/yjo10/packages/tau-2.34/build/x86_64/bin/tau_cxx.sh
CXX = /mnt/home/yjo10/packages/tau-2.34/build/x86_64/bin/tau_cc.sh
#CXX = g++
NVCC = nvcc

#exprot TAU_MAKEFILE=/mnt/home/yjo10/packages/tau-2.34/build/x86_64/lib/Makefile.tau-mpi
#TAU_MAKEFILE_INCLUDE =  /mnt/home/yjo10/packages/tau-2.34/build/include
#TAU_MAKEFILE_LIBS =  /mnt/home/yjo10/packages/tau-2.34/build/x86_64/lib

# Compiler flags
#CXXFLAGS  = -O1 -ftree-vectorize -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g  #-fopenmp  -fsanitize=thread
CXXFLAGS =  -O1 -I./Particle -I./ -g  -W -export-dynamic
#CXXFLAGS = -O2 
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -qopenmp  #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  -fsanitize=thread
#CXXFLAGS  = -O2 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  #-fsanitize=thread
NVCCFLAGS = -O2 --std=c++11 -G -Xcompiler -rdynamic -g -lineinfo -Xptxas=-v -lcudart -lcompute-sanitizer 
NVCC_FLAGS_DEBUG = -g -G --device-debug -lineinfo
#NVCC_FLAGS_PROFILE = -g -lineinfo

# Libraries
#LDFLAGS = -L/usr/local/cuda-11.8/lib64/ -lcudart -lcublas -lnvToolsExt #-ltbb
LDFLAGS = -lm 
#add -lnvToolsExt for <nvToolsExt.h>
#LDFLAGS = -L/mnt/sw/nix/store/zi2wc26znf75csf5hhz77p0d2bbz53ih-cuda-11.8.0/lib64/ -lcudart

# Source files
CXX_SRCS = main.cpp DefaultGlobal.cpp MPIRoutines.cpp readParameterFile.cpp WorkerRoutines.cpp RootRoutines.cpp TimeStepRoutines.cpp parser.cpp ReadWrite.cpp Particle/Initialize.cpp Particle/UpdateParticle.cpp Particle/ComputeAcceleration.cpp
#cuda/cuda_routines.cpp


#CU_SRCS = cuda/cuda_my_acceleration.cu cuda/cuda_kernels.cu 
# cuda/cuda_routines.cu

#
#OBJS = $(CU_SRCS:.cu=.o) $(CXX_SRCS:.cpp=.o)
OBJS = $(CXX_SRCS:.cpp=.o)

# Executable name
EXEC = abyss_mpi_tau.exe

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


#%.o: %.cu
#	$(NVCC) $(NVCCFLAGS) $(NVCC_FLAGS_DEBUG) -c $< -o $@


# Clean the project directory
clean:
	rm -f $(OBJS) $(EXEC)


