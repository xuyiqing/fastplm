#ifndef FASTPLM_CRUSH_QUEUE_H
#define FASTPLM_CRUSH_QUEUE_H

// [[Rcpp::plugins(cpp17)]]
#include <iostream>
#include <atomic>
#include <vector>
#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <random>
#include <algorithm>
#include <tuple>
#include <queue>

class CrushQueue {
public:
    typedef void* PayloadType;
    typedef std::function<void(PayloadType)> FunctionType;

private:
    std::vector<std::thread> workers;
    std::mutex mutex;
    std::size_t aliveWorkers;
    std::condition_variable hasJobs;
    std::condition_variable noJobs;

    std::queue<PayloadType> payloads;
    FunctionType function;

    bool running = true;

public:
    // Turn off Intel HT.
    CrushQueue(std::size_t threadCount = std::thread::hardware_concurrency() / 2);
    ~CrushQueue();

    void join() {
        std::unique_lock<std::mutex> lock(mutex);
        noJobs.wait(lock, [this](){ return aliveWorkers == 0; });
    }

    void crush();
    void commit(FunctionType&& function, std::queue<PayloadType>&& payloads);
    void commit(PayloadType payload);
};

extern CrushQueue* mainQueue;

#endif
