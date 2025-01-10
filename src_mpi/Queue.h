#ifndef QUEUE_H
#define QUEUE_H
/* This should contain all the information that should be transmitted to workers */
struct Queue {
    int task;
    int pid;
    double next_time;
    
    void print() {
        std::cout << "Task: " << task << ", PID: " << pid << ", Next Time: " << next_time << std::endl;
    }
};
#endif