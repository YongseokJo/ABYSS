#ifndef WORKER_H
#define WORKER_H
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "mpi.h"
#include "global.h"
#include "Queue.h"



#define MAX_QUEUE 1000
#define MAX_COMMUNICATION 10

struct Worker {
    int MyRank; 
    bool onDuty;
    std::unordered_set<int> CMPtclIDs;
    bool isCMWorker;
    Queue queues[MAX_QUEUE];
    short NumberOfQueues;
    short CurrentQueue;
    short NumberOfCommunications;

    MPI_Request _request[MAX_COMMUNICATION];  // Pointer to the request handle
    MPI_Status _status[MAX_COMMUNICATION];    // Pointer to the status object

    Worker() {
        _initialize();
    }

    void initialize() {
        onDuty = false;
        if (CMPtclIDs.size() > 0) {
           isCMWorker = true; 
        }
        NumberOfQueues = 0;
        CurrentQueue = 0;
        NumberOfCommunications = 0;
    }

    void initialize(int _MyRank) {
        MyRank = _MyRank;
        onDuty = false;
        if (CMPtclIDs.size() > 0) {
           isCMWorker = true; 
        }
        NumberOfQueues = 0;
        NumberOfCommunications = 0;
    }

    void addQueue(Queue queue) {
        int index = CurrentQueue+NumberOfQueues++;
        if (NumberOfQueues == MAX_QUEUE) {
           fprintf(stderr, "NumberOfQueues exceeds MAX_QUEUE!"); 
           exit(EXIT_FAILURE);
        }
        index %= MAX_QUEUE;
        queues[index] = queue;
    }

    void runQueue() {
        if (onDuty) {
            std::cout << "Worker " << MyRank << " is already on duty" << std::endl;
            std::cerr << "Worker " << MyRank << " is already on duty" << std::endl;
            exit(1);
        }
        if (NumberOfQueues == 0) {
            fprintf(stderr, "There is no queue in this worker!\n");
            exit(EXIT_FAILURE);
        }
        //Queue *q = &queues[CurrentQueue];
        //q->print();
        //sendTask(*q);
        sendTask(queues[CurrentQueue]);
        //std::cout << "Worker " << MyRank << " is on duty" << std::endl;
    }

    void removeQueue() {
    }

    void callback() {
        MPI_Waitall(NumberOfCommunications, _request, _status);
        NumberOfCommunications = 0;
        int return_value;
        MPI_Recv(&return_value, 1, MPI_INT, this->MyRank, TERMINATE_TAG, MPI_COMM_WORLD, &_status[0]);
        if (!onDuty) {
            fprintf(stderr, "Something's worng! the worker %d was not on duty.\n", this->MyRank);
            fprintf(stdout, "Something's worng! the worker %d was not on duty.\n", this->MyRank);
            exit(1);
        }
        onDuty = false;
        CurrentQueue++;
        CurrentQueue %= MAX_QUEUE;
        NumberOfQueues--;
        //std::cout << "Worker " << MyRank << " is off duty" << std::endl;
    }

/*
    void callback(int &return_value) {
        MPI_Recv(&return_value, 1, MPI_INT, this->MyRank, TERMINATE_TAG, MPI_COMM_WORLD, &_status);
        if (!onDuty) {
            fprintf(stderr, "Something's worng! the worker was not on duty.");
            exit(1);
        }
        onDuty = false;
    }
    */


    void sendTask(Queue &_queue) {
        if ((_queue.task == IrrForce) || (_queue.task == RegForce) || (_queue.task == ARIntegration))
        {
            MPI_Isend(&_queue.task,      1, MPI_INT,    this->MyRank, TASK_TAG, MPI_COMM_WORLD, &_request[NumberOfCommunications++]);
            MPI_Isend(&_queue.pid,       1, MPI_INT,    this->MyRank, PTCL_TAG, MPI_COMM_WORLD, &_request[NumberOfCommunications++]);
            MPI_Isend(&_queue.next_time, 1, MPI_DOUBLE, this->MyRank, TIME_TAG, MPI_COMM_WORLD, &_request[NumberOfCommunications++]);
        }
        else if (_queue.task == CommunicationSpeedBenchmark1) {
            MPI_Isend(&_queue.task, 1, MPI_INT, this->MyRank, TASK_TAG, MPI_COMM_WORLD, &_request[NumberOfCommunications++]);
        }
        else
        {
            MPI_Isend(&_queue.task, 1, MPI_INT, this->MyRank, TASK_TAG, MPI_COMM_WORLD, &_request[NumberOfCommunications++]);
            MPI_Isend(&_queue.pid,  1, MPI_INT, this->MyRank, PTCL_TAG, MPI_COMM_WORLD, &_request[NumberOfCommunications++] );
        }
        onDuty = true;
    }

    Queue* getCurrentQueue() {return &queues[CurrentQueue];}

private:
    //Queue current_queue;


    void _initialize() {
        MyRank = -1;
        onDuty = false;
        CMPtclIDs.clear();
        isCMWorker = false;
        NumberOfQueues = 0;
        CurrentQueue   = 0;
        NumberOfCommunications = 0;
        //queues.reserve(MAX_QUEUE);
    }
};

extern Worker* workers;
#endif