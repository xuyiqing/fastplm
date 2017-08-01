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

class CrashQueue {
private:
    std::vector<std::thread> workers;
    std::queue<void*> payloads;
    std::function<void(void*)> function;
    std::mutex taskFetchingMutex;
    std::atomic<long> aliveWorkers;

    std::condition_variable alarmClock;
    std::condition_variable sleepClock;

    bool running = true;

public:
    CrashQueue(std::size_t threadCount = std::thread::hardware_concurrency())
    : aliveWorkers(threadCount) {
        for (std::size_t i = 0; i < threadCount; i ++) {
            workers.emplace_back([this]() -> void {
                while (running) {
                    void* payload;
                    {
                        std::unique_lock<std::mutex> lock(taskFetchingMutex);
                        if (payloads.empty()) {
                            aliveWorkers.fetch_sub(1);
                            sleepClock.notify_one();

                            alarmClock.wait(lock);
                            continue;
                        }
                        payload = payloads.front();
                        payloads.pop();
                    }

                    function(payload);
                }
            });
        }

        // Make sure all workers finished running.
        while (aliveWorkers.load() > 0);
        std::unique_lock<std::mutex> lock(taskFetchingMutex);
    }

    ~CrashQueue() {
        running = false;
        alarmClock.notify_all();
        for (auto& worker : workers)
            worker.join();
    }

    void run() {
        this->aliveWorkers = workers.size();
        alarmClock.notify_all();

        while (true) {
            std::unique_lock<std::mutex> lock(taskFetchingMutex);
            if (aliveWorkers.load() == 0)
                break;
            sleepClock.wait(lock);
        }
    }

    void commit(std::function<void(void*)>&& function, std::queue<void*>&& payloads) {
        this->function = std::move(function);
        this->payloads = std::move(payloads);
    }
};

#endif
