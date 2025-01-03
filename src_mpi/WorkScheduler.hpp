#ifndef WORK_SCHEDULER_H
#define WORK_SCHEDULER_H
#include <iostream>
#include "def.h"
#include "global.h"

class WorkScheduler {
public:
bool isWorking[NumberOfWorker];

private:
std::vector<int> _ParticleList;


void initialize(std::vector<int>& ParticleList) {
    _ParticleList = ParticleList;
    for (int i = 0; i < NumberOfWorker; i++) {
        isWorking[i] = false;
    }
}

void assignJobs() {

}

};

#endif