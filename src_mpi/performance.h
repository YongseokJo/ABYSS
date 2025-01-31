#ifndef PERFORMANCE_H
#define PERFORMANCE_H

#include <iostream>
#include <chrono>

#define MAX_PERFORMANCE_ENTITY 100

#define InitAcc1_Root 0
#define InitAcc1_Queueing 1
#define InitAcc1_Running 2
#define IrrForce_Root 10
#define IrrForce_Worker 11

struct Performance {
	//std::chrono::milliseconds IrregularForce{0};

	std::chrono::high_resolution_clock::time_point start_point[MAX_PERFORMANCE_ENTITY];
	std::chrono::high_resolution_clock::time_point end_point[MAX_PERFORMANCE_ENTITY];

	long duration[MAX_PERFORMANCE_ENTITY];

    void start(int tag) {
		start_point[tag] = std::chrono::high_resolution_clock::now();
	}

	void end(int tag) {
		end_point[tag] = std::chrono::high_resolution_clock::now();
		duration[tag] +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(end_point[tag] - start_point[tag]).count();
	}

	long get(int tag) {
		return duration[tag];
	}
};

#endif

