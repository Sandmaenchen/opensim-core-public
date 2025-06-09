
#include "UKFThreadPool.h"

// UKFThreadPool methods

// OpenSim::UKFThreadPool::UKFThreadPool() : numTasksPending(0), stop(false) {
//     workers.emplace_back(std::bind(&UKFThreadPool::workerThread, this));
// }

OpenSim::UKFThreadPool::UKFThreadPool(size_t num_threads) : numTasksPending(0), stop(false) {
    for (size_t i = 0; i < num_threads; ++i) {
        workers.emplace_back(std::bind(&UKFThreadPool::workerThread, this));
    }
}

/*
template<class F>
void OpenSim::UKFThreadPool::enqueue(F f) {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        tasks.emplace(std::function<void()>(f));
    }
    numTasksPending++;
    condition.notify_one();
}
*/

void OpenSim::UKFThreadPool::waitUntilCompleted() {
    std::unique_lock<std::mutex> lock(main_mutex);
    if (numTasksPending != 0) {
        main_condition.wait(lock);
    }
    else {
        lock.unlock();
    }
}

OpenSim::UKFThreadPool::~UKFThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (std::thread& worker : workers) {
        worker.join();
    }
}

void OpenSim::UKFThreadPool::workerThread() {
    while (true) {
        std::function<void()> task;
        {
            std::unique_lock<std::mutex> queueLock(queue_mutex);
            condition.wait(queueLock, [this] { return stop || !tasks.empty(); });
            if (stop && tasks.empty()) {
                return;
            }
            task = tasks.front();
            tasks.pop();
        }
        task();
        {
            std::lock_guard<std::mutex> mainLock(main_mutex);
            numTasksPending--;
            if (numTasksPending == 0) {
                main_condition.notify_one();
            }
        }
    }
}
