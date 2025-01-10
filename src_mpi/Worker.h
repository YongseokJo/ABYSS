#ifndef WORKER_H
#define WORKER_H
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "mpi.h"
#include "def.h"
#include "particle.h"
#include "global.h"
#include "Queue.h"



#define MAX_QUEUE 10

struct Worker {
    int MyRank; 
    bool onDuty;
    std::unordered_set<int> CMPtclIDs;
    bool isCMWorker;
    Queue queues[MAX_QUEUE];
    short NumberOfQueues;
    short CurrentQueue;


    Worker() {
        MyRank = -1;
        onDuty = false;
        CMPtclIDs.clear();
        isCMWorker = false;
        NumberOfQueues = 0;
        CurrentQueue   = 0;
        //queues.reserve(MAX_QUEUE);
    }

    void initialize() {
        onDuty = false;
        if (CMPtclIDs.size() > 0) {
           isCMWorker = true; 
        }
        NumberOfQueues = 0;
    }

    void addQueue(Queue queue) {
        int index = CurrentQueue+NumberOfQueues++;
        index %= MAX_QUEUE;
        queues[index] = queue;
    }

    void runQueue() {
        if (NumberOfQueues == 0) {
            fprintf(stderr, "There is no queue in this worker!\n");
            exit(EXIT_FAILURE);
        }
        Queue *q = &queues[CurrentQueue++];
        //q->print();
        sendTask(*q);
        //sendTask(queues[CurrentQueue++]);
        CurrentQueue %= MAX_QUEUE;
        NumberOfQueues--;
    }

    void removeQueue() {
    }

    void callback() {
        int return_value;
        MPI_Recv(&return_value, 1, MPI_INT, this->MyRank, TERMINATE_TAG, MPI_COMM_WORLD, &_status);
        if (!onDuty) {
            fprintf(stderr, "Something's worng! the worker was not on duty.");
            exit(1);
        }
        onDuty = false;
    }

    void callback(int &return_value) {
        MPI_Recv(&return_value, 1, MPI_INT, this->MyRank, TERMINATE_TAG, MPI_COMM_WORLD, &_status);
        if (!onDuty) {
            fprintf(stderr, "Something's worng! the worker was not on duty.");
            exit(1);
        }
        onDuty = false;
    }


    void sendTask(Queue &_queue) {
        if ((_queue.task == 0) || (_queue.task == 1) || (_queue.task == 26))
        {
            MPI_Send(&_queue.task,      1, MPI_INT,    this->MyRank, TASK_TAG, MPI_COMM_WORLD);
            MPI_Send(&_queue.pid,       1, MPI_INT,    this->MyRank, PTCL_TAG, MPI_COMM_WORLD);
            MPI_Send(&_queue.next_time, 1, MPI_DOUBLE, this->MyRank, TIME_TAG, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(&_queue.task, 1, MPI_INT, this->MyRank, TASK_TAG, MPI_COMM_WORLD);
            MPI_Send(&_queue.pid,  1, MPI_INT, this->MyRank, PTCL_TAG, MPI_COMM_WORLD);
        }
        onDuty = true;
    }



private:
    //Queue current_queue;
    MPI_Request _request;  // Pointer to the request handle
    MPI_Status _status;    // Pointer to the status object
};

extern Worker* workers;
#endif