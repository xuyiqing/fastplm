#ifndef FASTPLM_CRASH_QUEUE_H
#define FASTPLM_CRASH_QUEUE_H

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
    std::mutex taskFetchingMutex;
    std::atomic<long> aliveWorkers;
    std::condition_variable alarmClock;
    std::condition_variable sleepClock;

    std::queue<PayloadType> payloads;
    FunctionType function;

    bool running = true;

public:
    // Turn off Intel HT.
    CrushQueue(std::size_t threadCount = std::thread::hardware_concurrency() / 2);
    ~CrushQueue();
    void crash();
    void commit(FunctionType&& function, std::queue<PayloadType>&& payloads);
};

extern CrushQueue* mainQueue;

#endif
