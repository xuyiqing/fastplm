#include "CrushQueue.h"

CrushQueue* mainQueue = nullptr;

CrushQueue::CrushQueue(std::size_t threadCount): aliveWorkers(threadCount) {
    for (std::size_t i = 0; i < threadCount; i ++) {
        workers.emplace_back([this]() -> void {
            while (running) {
                void* payload;
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    if (payloads.empty()) {
                        aliveWorkers -= 1;
                        noJobs.notify_one();
                        hasJobs.wait(lock);
                        aliveWorkers += 1;
                        continue;
                    }
                    payload = payloads.front();
                    payloads.pop();
                }

                function(payload);
            }
        });
    }

    // Make sure all workers finished initialization.
    this->join();
}

CrushQueue::~CrushQueue() {
    running = false;
    hasJobs.notify_all();
    for (auto& worker : workers)
        worker.join();
}

void CrushQueue::crush() {
    std::unique_lock<std::mutex> lock(mutex);
    hasJobs.notify_all();
    do {
        noJobs.wait(lock);
    } while (aliveWorkers > 0);
}

void CrushQueue::commit(CrushQueue::FunctionType&& function,
                        std::queue<CrushQueue::PayloadType>&& payloads) {
    this->join();
    this->function = std::move(function);
    this->payloads = std::move(payloads);
}

void CrushQueue::commit(PayloadType paylod) {
    std::unique_lock<std::mutex> lock(mutex);
    payloads.push(paylod);
    hasJobs.notify_one();
}
