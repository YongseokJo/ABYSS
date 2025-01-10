#ifndef QUEUE_SCHEDULER_H
#define QUEUE_SCHEDULER_H
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "mpi.h"
#include "def.h"
#include "particle.h"
#include "global.h"
#include "Queue.h"
#include "Worker.h"





class QueueScheduler
{
public:
    std::unordered_set<Worker*> WorkersToGo;


    QueueScheduler()
    {
        _FreeWorkers.reserve(NumberOfWorker);
        WorkersToGo.reserve(NumberOfWorker);
    }

    void initialize(int task, double next_time)
    {
        _initialize();
        _task = task;
        _next_time = next_time;
    }

    void initialize(int task)
    {
        _initialize();
        _task = task;
        _next_time = -1.;
    }

    void takeQueue(std::vector<int> &queue_list) {
        _queue_list = queue_list;
        _total_queues = _queue_list.size();
    }

    void assignQueueAuto() {
        for (auto worker = _FreeWorkers.begin(); worker != _FreeWorkers.end();)
        {
            if (_queue_list.size() > 0)
            {
                _queue.task = _task;
                _queue.pid = _queue_list.back();
                _queue_list.pop_back();
                _queue.next_time = _next_time;
                (*worker)->addQueue(_queue);
                WorkersToGo.insert(*worker);
                _assigned_queues++;
                worker = _FreeWorkers.erase(worker);
                //_queue.print();
            }
            else {
               ++worker; 
            }
        }
    }

    void runQueueAuto() {
        for (auto worker = WorkersToGo.begin(); worker != WorkersToGo.end();)
        {
            if ((*worker)->NumberOfQueues > 0)
            {
                (*worker)->runQueue();
                worker = WorkersToGo.erase(worker);
            }
            else
            {
                ++worker;
            }
        }
    }

    void assignQueue(Worker *worker, Queue *queue) {
        worker->addQueue(*queue);
        _FreeWorkers.erase(worker);
        _assigned_queues++;
        WorkersToGo.insert(worker);
    }

    Worker* waitQueue(int type) {
        if (type == 0) // blocking
        {
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &_status);
            _completed_queues++;
            // Retrieve the rank of the source processor
            _rank = _status.MPI_SOURCE;
            //fprintf(stdout, "returned rank = %d\n",_rank);
            workers[_rank].callback();
            _FreeWorkers.insert(&workers[_rank]);
            return &workers[_rank];
        }
        if (type == 1) // non-blocking
        {
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &_flag, &_status);
            _rank = _status.MPI_SOURCE;
            return &workers[_rank];
        }
    }


    void callback() {

    }




    bool isComplete() {
        if (_completed_queues == _total_queues && _assigned_queues == _total_queues)
            return _complete;
        else
            return _not_complete; 
    }

    void printFreeWorker() {
        std::cout << "<FreeWorkers List>" << std::endl;
        for (Worker* worker: _FreeWorkers) {
            std::cout << "MyRank=" << worker->MyRank << std::endl;
        }
    }

    void printWorkerToGo() {
        std::cout << "<WorkersToGo List>" << std::endl;
        for (Worker* worker: WorkersToGo) {
            std::cout << "MyRank=" << worker->MyRank << std::endl;
        }
    }


private:
    std::unordered_set<Worker*> _FreeWorkers;
    std::vector<int> _queue_list;
    Queue _queue;
    int _task, _rank, _flag;
    double _next_time;
    int _total_queues, _assigned_queues, _completed_queues;
    const bool _complete = false;
    const bool _not_complete = true;
    MPI_Request _request;  // Pointer to the request handle
    MPI_Status _status;    // Pointer to the status object



    void _initialize() {
        WorkersToGo.clear();
        _FreeWorkers.clear();
        for (int i = 1; i <= NumberOfWorker; i++) {
            workers[i].initialize();
            _FreeWorkers.insert(&workers[i]);
        }
        _total_queues=0;
        _assigned_queues=0;
        _completed_queues=0;;
    }

#ifdef unuse


        void doSimpleTask(int task, std::vector<int> &list)
    {
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
        std::unordered_set<int> CMPtcls;
        std::vector<int> UpdatedCMPtcl;
        std::unordered_set<int> ReadyToGoCMPtcl;
        int pid_tmp, return_value;
        SingleParticleList.reserve(_total_tasks);

        _setTask(_task);

        for (int i=0; i<_ParticleList.size(); i++) {
            pid_tmp = _ParticleList[i];
            particles[pid_tmp].isUpdateToDate = false;
            if (particles[pid_tmp].isCMptcl)
            {
                SingleParticleList.push_back(pid_tmp);
            }
            else
            {
                CMPtcls.insert(pid_tmp);
            }
        }

        int fb_total_tasks = CMPtcls.size(); // this is for the SDAR computations by YS 2025.01.06
        int fb_assigned_tasks = 0; // this is for the SDAR computations by YS 2025.01.06
        UpdatedCMPtcl.resize(fb_total_tasks);

        do {
            if (_FreeWorkers.size() == 0 || _assigned_tasks == _total_tasks) 
            {
                // have to add check all the sends are recved.
                // MPI_Waitall(NumberOfCommunication, requests, statuses);
                // NumberOfCommunication = 0;
                _checkCompletion(return_value);
                // this ensures that cmptcls are ready to go (neighbors are update to date)
                _checkCMPtclNeighbor(UpdatedCMPtcl, ReadyToGoCMPtcl); 
                _WorkerTmp = _Callback();
                _completed_tasks++;

            }

            // SDAR assignment
            if (_FreeWokrers.size() != 0 && ReadyToGoCMPtcl.size() != 0) {
                // I have to add that if the worker is a CMworker, I should run SDAR integration (Query to myself)
                if (particles[_WorkerTmp->PID].isCMptcl && _WorkerTmp->isCMWorker) 
                {
                    // all neighbors are up to date.
                
                
                        fprintf(stderr, "Worker does not match!\n");
                        exit(EXIT_FAILURE);
                    }
                    _FreeWorkers.pop_back();
                    _WorkerTmp->task = 26;
                    _assignJobs(_WorkerTmp);
                    fb_assigned_tasks++;
                }
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
                        if (CMPtcls.find(cm_pid) != CMPtcls.end())
                        {
                            _WorkerTmp->PID = cm_pid;
                            _WorkerTmp->isCMWorker = true;
                            CMPtcls.erase(cm_pid);
                            UpdatedCMPtcl.push_back(cm_pid);
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

                _WorkerTmp->task = 0;
                _WorkerTmp->next_time = next_time;
                //fprintf(stdout, "assigned_tasks = %d/%d, number of free worker = %d pid = %d rank = %d\n",
                //_assigned_tasks, _total_tasks, _FreeWorkers.size(), _WorkerTmp->PID, _WorkerTmp->MyRank);
                _assignJobs(_WorkerTmp);
                _assigned_tasks++;
            }
        } while (_completed_tasks < _total_tasks && fb_completed_tasks < fb_total_tasks);

        if (SingleParticleList.size() != 0 || CMPtcls.size() != 0) {
            fprintf(stderr, "Something's wrong. The particle list is not empty.\n");
            exit(1);
        }

        UpdatedCMPtcl.clear();
        UpdatedCMPtcl.shrink_to_fit();
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
    int _rank;
    Particle* _ptcl;
    Worker* _WorkerTmp, _CompletedWorker;



    void _assignJobs(Worker *worker) {
        if ((worker->task == 0) || (worker->task == 1) || (worker->task == 26))
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


    Worker* _Callback() {
        MPI_Wait(&_request, &_status);
        int rank = _status.MPI_SOURCE;
        if (!workers[rank-1].onDuty) {
            fprintf(stderr, "Something's worng! the worker was not on duty.");
            exit(1);
        }
        workers[rank-1].onDuty = false;
        _FreeWorkers.push_back(&workers[rank-1]);
        return &workers[rank-1];
    }


    void _setTask(int task) {
        for (int i = 0; i < NumberOfWorker; i++) {
            workers[i].task = task;
        }
    }

    void _checkCMPtclNeighbor(const std::vector<int>& UpdatedCMPtcl, std::unordered_set<int>& ReadyToGoCMPtcl) {
        Particle *ptcl;
        int i = 0;
        for (int pid:UpdatedCMPtcl) {
           ptcl = &particles[pid]; 
           if (!ptcl->isCMPtcl) {
            fprintf(stderr, "this is not CM ptcl!\n");
            exit(EXIT_FAILURE);
           }
           for (int j = 0; j < _ptcl->NumberOfNeighbor; j++)
           {
               if (particles[_ptcl->Neighbors[j]].isUpdateToDate == false)
                   break;
           }
           ReadyToGoCMPtcl.insert(pid); // (Query to myself) should I add worker rank?
        }
    }
#endif
};




#endif
