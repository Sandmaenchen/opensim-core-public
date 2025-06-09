#ifndef OPENSIM_UKFTHREADPOOL_H_
#define OPENSIM_UKFTHREADPOOL_H_

/* -------------------------------------------------------------------------- *
 *                           OpenSim:  UKFThreadPool.h                        *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2024 Stanford University and the Authors                *
 * Author(s): Matti Kortelainen                                               *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include <thread>
#include <mutex>
#include <functional>
#include <vector>
#include <condition_variable>
#include <queue>
#include <atomic>

namespace OpenSim {

    class UKFThreadPool {
    public:
        //UKFThreadPool();
        UKFThreadPool(size_t num_threads);        

        template<class F>
        void enqueue(F f) {
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                tasks.emplace(std::function<void()>(f));
            }
            numTasksPending++;
            condition.notify_one();
        }

        void waitUntilCompleted();

        ~UKFThreadPool();

    private:
        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;
        std::mutex queue_mutex;
        std::condition_variable condition;
        std::atomic<int> numTasksPending;
        std::mutex main_mutex;
        std::condition_variable main_condition;
        bool stop;

        void workerThread();

    };  // END of class UKFThreadPool

}   // END of namespace OpenSim

#endif  // OPENSIM_UKFTHREADPOOL_H_