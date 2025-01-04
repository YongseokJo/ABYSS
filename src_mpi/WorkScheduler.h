#ifndef WORK_SCHEDULER_H
#define WORK_SCHEDULER_H
#include <iostream>
#include <vector>
#include <map>
#include "def.h"
#include "particle.h"
#include "global.h"

struct Worker {
    int MyRank; 
    int task;
    int PID;
    double next_time;
    bool onDuty;
    std::vector<int> CMPtclIDs;
    bool isCMWorker;

    void initialize() {
        task = -1;
        onDuty = false;
        next_time = -1;
        PID = -1;
        if (CMPtclIDs.size() > 0) {
           isCMWorker = true; 
        }
    }
};


class WorkScheduler
{
public:
    int NumberOfOnGoingWork;
    int NumberOfFreeWorker;


    WorkScheduler() {
        _FreeWorkers.reserve(NumberOfWorker);
    }
    void initialize()
    {
        for (int i = 0; i < NumberOfWorker; i++) {
            workers[i].initialize();
            _FreeWorkers.push_back(&workers[i]);
        }
        NumberOfOnGoingWork = 0;
        NumberOfFreeWorker = NumberOfWorker;
    }


    void doSimpleTask(int task, int total_tasks) {
        int i = 0;
        _task = task;
        _total_tasks = total_tasks;
        _completed_tasks = 0;
        int pid_tmp;

        _setTask(_task);

        do {
            if (_FreeWorkers.size() == 0) {
                // have to add check all the sends are recved.
                // MPI_Waitall(NumberOfCommunication, requests, statuses);
                // NumberOfCommunication = 0;
                _checkCompletion();
                /* we can do something here */
                _Callback();
                _completed_tasks++;
            }

            if (_FreeWorkers.size() != 0) {
                _WorkerTmp = _FreeWorkers.back();
                _FreeWorkers.pop_back();
                _WorkerTmp->PID = pid_tmp;
                _assignJobs(_WorkerTmp);
            }
            pid_tmp++;
        } while(_completed_tasks < _total_tasks);
    }




    void doIrregularForce(std::vector<int> &ParticleList,
                          std::vector<int> &CMPtclList, std::vector<int> &CMWorkerList,
                          std::map<int, int> CMPtcl_Worker_Map, std::map<int, int> Worker_CMPtcl_Map,
                          double &next_time)

    {
        _ParticleList = ParticleList;
        _total_tasks = _ParticleList.size();
        _completed_tasks = 0;
        _assinged_tasks = 0;
        std::vector<int> SingleParticleList;
        std::unordered_map<int, int> CMPtclMap;
        int pid_tmp;
        SingleParticleList.reserve(_total_tasks);

        _setTask(_task);

        for (int i=0; i<_ParticleList.size(); i++) {
            pid_tmp = _ParticleList[i];
            if (particles[pid_tmp].GroupInfo == nullptr)
            {
                SingleParticleList.push_back(pid_tmp);
            }
            else
            {
                CMPtclMap.insert({pid_tmp, i});
            }
        }

        do {
            if (_FreeWorkers.size() == 0) {
                // have to add check all the sends are recved.
                // MPI_Waitall(NumberOfCommunication, requests, statuses);
                // NumberOfCommunication = 0;
                _checkCompletion();
                /* we can do something here */
                _Callback();
                _completed_tasks++;
            }

            if (_FreeWorkers.size() != 0 && _assinged_tasks < _total_tasks)
            {
                _WorkerTmp = _FreeWorkers.back();
                _FreeWorkers.pop_back();

                if (_WorkerTmp->isCMWorker)
                {
                    for (int cm_pid : _WorkerTmp->CMPtclIDs)
                    {
                        _WorkerTmp->isCMWorker = false;
                        if (CMPtclMap.find(cm_pid) != CMPtclMap.end())
                        {
                            _WorkerTmp->PID = cm_pid;
                            _WorkerTmp->isCMWorker = true;
                            CMPtclMap.erase(cm_pid);
                            break;
                        }
                    }
                }

                if (!_WorkerTmp->isCMWorker)
                {
                    pid_tmp = SingleParticleList.back();
                    SingleParticleList.pop_back();
                    _WorkerTmp->PID = pid_tmp;
                }

                _WorkerTmp->next_time = next_time;
                _assignJobs(_WorkerTmp);
                _assinged_tasks++;
            }
        } while (_completed_tasks < _total_tasks);

        if (SingleParticleList.size() != 0 || CMPtclMap.size() != 0) {
            fprintf(stderr, "Something's wrong. The particle list is not empty.\n");
            exit(1);
        }

        SingleParticleList.clear();
        SingleParticleList.shrink_to_fit();
    }

private:
    std::vector<int> _ParticleList;
    int _total_tasks, _completed_tasks, _assinged_tasks, _remaining_tasks;
    int _WorkerRank, _pid;
    int _task;
    std::vector<Worker*> _FreeWorkers;
    MPI_Request _request;  // Pointer to the request handle
    MPI_Status _status;    // Pointer to the status object
    Particle* _ptcl;
    Worker* _WorkerTmp, _CompletedWorker;



    void _assignJobs(Worker *worker) {
        if ((worker->task == 0) || (worker->task == 1)) {
            MPI_Send(&worker->task,      1, MPI_INT,    worker->MyRank, TASK_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->PID,       1, MPI_INT,    worker->MyRank, PTCL_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->next_time, 1, MPI_DOUBLE, worker->MyRank, TIME_TAG, MPI_COMM_WORLD);
        }
        else if ((worker->task == 7) || (worker->task == 8)) {
            MPI_Send(&worker->task, 1, MPI_INT, worker->MyRank, TASK_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->PID,  1, MPI_INT, worker->MyRank, PTCL_TAG, MPI_COMM_WORLD);
        }
        worker->onDuty = true;
    }


    void _checkCompletion() {
        int task;
		MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &_request);
    }


    void _Callback() {
        MPI_Wait(&_request, &_status);
        int rank = _status.MPI_SOURCE;
        if (!workers[rank-1].onDuty) {
            fprintf(stderr, "Something's worng! the worker was not on duty.");
            exit(1);
        }
        workers[rank-1].onDuty = false;
        _FreeWorkers.push_back(&workers[rank-1]);
    }


    void _setTask(int task) {
        for (int i = 0; i < NumberOfWorker; i++) {
            workers[i].task = task;
        }
    }
};




#endif