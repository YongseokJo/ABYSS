#ifndef WORK_SCHEDULER_H
#define WORK_SCHEDULER_H
#include <iostream>
#include <vector>
#include <unordered_map>
#include "mpi.h"
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

    Worker() {
        MyRank = -1;
        task = -1;
        onDuty = false;
        next_time = -1;
        PID = -1;
        CMPtclIDs.clear();
        isCMWorker = false;
    }

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

extern Worker* workers;

class WorkScheduler
{
public:

    WorkScheduler() {
        _FreeWorkers.reserve(NumberOfWorker);
    }
    void initialize()
    {
        _FreeWorkers.clear();
        _FreeWorkers.resize(NumberOfWorker);
        for (int i = 0; i < NumberOfWorker; i++) {
            workers[i].initialize();
            _FreeWorkers[i] = &workers[i];
        }
        _completed_tasks = 0;
        _assigned_tasks = 0;
    }


    void doSimpleTask(int task, std::vector<int> &list) {
        int return_value, i = 0;
        _task = task;
        _total_tasks = list.size();
        _setTask(_task);
        //fprintf(stdout, "task=%d, the size of FreeWorker is %d\n", _task, _FreeWorkers.size());
        do {
            if (_FreeWorkers.size() == 0 || _assigned_tasks == _total_tasks) {
                // have to add check all the sends are recved.
                // MPI_Waitall(NumberOfCommunication, requests, statuses);
                // NumberOfCommunication = 0;
                _checkCompletion(return_value);
                /* we can do something here */
                _Callback();
                _completed_tasks++;
            }

            if (_FreeWorkers.size() != 0 && _assigned_tasks < _total_tasks) {
                _WorkerTmp = _FreeWorkers.back();
                _FreeWorkers.pop_back();
                _WorkerTmp->PID = list[i++];
                _assignJobs(_WorkerTmp);
                _assigned_tasks++;
                //fprintf(stdout, "assigned_tasks = %d, number of free worker = %d pid = %d rank = %d\n",
                //_assigned_tasks, _FreeWorkers.size(), _WorkerTmp->PID, _WorkerTmp->MyRank);
                fflush(stdout);
            }
        } while(_completed_tasks < _total_tasks);
    }




    void doIrregularForce(std::vector<int> &ParticleList,
                          double &next_time)
    {
        _task=0;
        _ParticleList = ParticleList;
        _total_tasks = _ParticleList.size();
        std::vector<int> SingleParticleList;
        std::unordered_map<int, int> CMPtclMap;
        int pid_tmp, return_value;
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
            if (_FreeWorkers.size() == 0 || _assigned_tasks == _total_tasks) {
                // have to add check all the sends are recved.
                // MPI_Waitall(NumberOfCommunication, requests, statuses);
                // NumberOfCommunication = 0;
                _checkCompletion(return_value);
                /* we can do something here */
                _Callback();
                _completed_tasks++;
                // I have to add that if the worker is a CMworker, I should run SDAR integration
            }

            if (_FreeWorkers.size() != 0 && _assigned_tasks < _total_tasks)
            {
                _WorkerTmp = _FreeWorkers.back();
                _FreeWorkers.pop_back();

                // CM ptcl assignment
                // if this worker contains CM ptcl, see if any relevant CM ptcls in ParticleList.
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

                // Single particle assignment
                if (!_WorkerTmp->isCMWorker)
                {
                    pid_tmp = SingleParticleList.back();
                    SingleParticleList.pop_back();
                    _WorkerTmp->PID = pid_tmp;
                }

                _WorkerTmp->next_time = next_time;
                //fprintf(stdout, "assigned_tasks = %d/%d, number of free worker = %d pid = %d rank = %d\n",
                //_assigned_tasks, _total_tasks, _FreeWorkers.size(), _WorkerTmp->PID, _WorkerTmp->MyRank);
                _assignJobs(_WorkerTmp);
                _assigned_tasks++;
            }
        } while (_completed_tasks < _total_tasks);

        if (SingleParticleList.size() != 0 || CMPtclMap.size() != 0) {
            fprintf(stderr, "Something's wrong. The particle list is not empty.\n");
            exit(1);
        }

        SingleParticleList.clear();
        SingleParticleList.shrink_to_fit();
    }

//private:
    std::vector<int> _ParticleList;
    int _total_tasks, _completed_tasks, _assigned_tasks, _remaining_tasks;
    int  _pid;
    int _task;
    std::vector<Worker*> _FreeWorkers;
    MPI_Request _request;  // Pointer to the request handle
    MPI_Status _status;    // Pointer to the status object
    Particle* _ptcl;
    Worker* _WorkerTmp, _CompletedWorker;



    void _assignJobs(Worker *worker) {
        if ((worker->task == 0) || (worker->task == 1))
        {
            MPI_Send(&worker->task,      1, MPI_INT,    worker->MyRank, TASK_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->PID,       1, MPI_INT,    worker->MyRank, PTCL_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->next_time, 1, MPI_DOUBLE, worker->MyRank, TIME_TAG, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(&worker->task, 1, MPI_INT, worker->MyRank, TASK_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->PID,  1, MPI_INT, worker->MyRank, PTCL_TAG, MPI_COMM_WORLD);
        }
        worker->onDuty = true;
    }


    void _checkCompletion(int &return_value) {
        MPI_Irecv(&return_value, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &_request);
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