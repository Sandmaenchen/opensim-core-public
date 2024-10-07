#ifndef OPENSIM_UKF_IMU_INVERSE_KINEMATICS_TOOL_H_
#define OPENSIM_UKF_IMU_INVERSE_KINEMATICS_TOOL_H_
/* -------------------------------------------------------------------------- *
 *                    OpenSim:  UKFIMUInverseKinematicsTool.h                 *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Distribution Statement A – Approved for Public Release; Distribution is    *
 * Unlimited.                                                                 *
 *                                                                            *
 * UKFIMUInverseKinematicsTool is written by University of Eastern Finland    *
 * and bases on the IMUInverseKinematicsTool code written by Ajay Seth.       *
 *                                                                            *
 * The project was supported by the Research Council of Finland.              *
 *                                                                            *
 * Copyright (c) 2024 University of Eastern Finland                           *
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

#include "osimToolsDLL.h"
#include <OpenSim/Common/Object.h>
#include <OpenSim/Common/ModelDisplayHints.h>
#include <OpenSim/Common/Set.h>
#include <OpenSim/Common/TimeSeriesTable.h>
#include <OpenSim/Common/IO.h>
#include <OpenSim/Common/TableSource.h>
#include <OpenSim/Common/STOFileAdapter.h>
#include <OpenSim/Common/TRCFileAdapter.h>
#include <OpenSim/Common/Reporter.h>
#include <OpenSim/Simulation/Model/Point.h>
#include <OpenSim/Simulation/OrientationsReference.h>
#include <OpenSim/Simulation/OpenSense/OpenSenseUtilities.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PhysicalOffsetFrame.h>
#include <OpenSim/Simulation/InverseKinematicsSolver.h>
#include <OpenSim/Simulation/OrientationsReference.h>
#include <OpenSim/Tools/InverseKinematicsToolBase.h>
 //include to realize UKF
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/AssemblyCondition_OrientationSensors.h"
#include <OpenSim/Tools/IMUInverseKinematicsTool.h>
//#include <Eigen/Eigen>
//#include <Eigen/Eigenvalues>
//#include <Eigen/Cholesky>
#include "eigen/Eigen/Eigen"
#include "eigen/Eigen/Eigenvalues"
#include "eigen/Eigen/Cholesky"
#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>
#include <regex>
#include <functional>
#include <vector>
#include <condition_variable>
#include <queue>
#include <deque>



//#include <future>


namespace OpenSim {

//class Model;
// Only reason for this class to exist is to have nicer name in XML
//class OSIMTOOLS_API OrientationWeightSet : public Set<OrientationWeight> {
//    OpenSim_DECLARE_CONCRETE_OBJECT(
//            OrientationWeightSet, Set<OrientationWeight>);

//public:
    /** Use Super's constructors. */
//    using Super::Super;
    // default copy, assignment operator, and destructor
//    OrientationWeightSet() = default;
//    OrientationWeightSet(const OrientationWeightSet&) = default;
    //=============================================================================
//}; 
        //=============================================================================
//=============================================================================
/**
 * A Study that performs an Inverse Kinematics analysis with a given model.
 * Inverse kinematics is the solution of internal coordinates that poses
 * the model such that the body rotations (as measured by IMUs) affixed to the 
 * model minimize the weighted least-squares error with observations of IMU 
 * orientations in their spatial coordinates. 
 *
 * @author Ajay Seth
 */


class UKFThreadPool {
public:
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


class OSIMTOOLS_API UKFIMUInverseKinematicsTool
        : public InverseKinematicsToolBase {
    OpenSim_DECLARE_CONCRETE_OBJECT(
            UKFIMUInverseKinematicsTool, InverseKinematicsToolBase);

public:
    OpenSim_DECLARE_PROPERTY(orientations_file, std::string,
        "Name/path to a .sto file of sensor frame orientations as quaternions.");

    OpenSim_DECLARE_PROPERTY(sensor_to_opensim_rotations, SimTK::Vec3,
            "Space fixed Euler angles (XYZ order) from IMU Space to OpenSim."
            " Default to (0, 0, 0).");
    OpenSim_DECLARE_PROPERTY(orientation_weights, OrientationWeightSet,
            "Set of orientation weights identified by orientation name with "
            "weight being a positive scalar. If not provided, all IMU "
            "orientations are tracked with weight 1.0.");
    OpenSim_DECLARE_PROPERTY(write_UKF, bool, 
    "Write the means and covariances of Kalman smoother estimates up to 2nd order.");
    OpenSim_DECLARE_PROPERTY(num_threads, int, 
    "Number of threads to use in UKF forward filtering.");
    OpenSim_DECLARE_PROPERTY(alpha, double, 
    "UKF hyperparameter alpha.");
    OpenSim_DECLARE_PROPERTY(beta, double, 
    "UKF hyperparameter beta.");
    OpenSim_DECLARE_PROPERTY(kappa, double, 
    "UKF hyperparameter kappa. Special values: "
    "-1.337 => Set kappa equal to (3 - length of state vector)."
    "-4.337 => Set kappa equal to length of state vector.");
    OpenSim_DECLARE_PROPERTY(sgma2w, double, 
    "Scaling factor for process noise covariance matrix");
    OpenSim_DECLARE_PROPERTY(order, int, 
    "Largest order of position time derivatives to use in process model.");
    OpenSim_DECLARE_PROPERTY(lag_length, int, 
    "Number of future samples to use in backward pass of Kalman smoother.");
    OpenSim_DECLARE_PROPERTY(missing_data_scale, double, 
    "Scaling factor for observation noise covariance matrix elements in case of missing observations.");
    OpenSim_DECLARE_PROPERTY(imu_RMS_in_deg, SimTK::Vec3, 
    "RMS errors for each of 3 axes of orientation sensors.");
    OpenSim_DECLARE_PROPERTY(enable_resampling, bool, 
    "Resample sigma points after propagation through process model.");
    OpenSim_DECLARE_PROPERTY(process_covariance_method, int, 
    "Model for process noise covariance matrix. "
    "0 = classic white noise process (Fioretti and Jetto, 1989); "
    "1 = scaled remainders of Taylor series expansion.");


    //=============================================================================
// METHODS
//=============================================================================
    //--------------------------------------------------------------------------
    // CONSTRUCTION
    //--------------------------------------------------------------------------
public:
    virtual ~UKFIMUInverseKinematicsTool();
    UKFIMUInverseKinematicsTool();
    UKFIMUInverseKinematicsTool(const std::string &setupFile);
    //--------------------------------------------------------------------------
    // INTERFACE
    //--------------------------------------------------------------------------

    //The roll, heading, and pitch RMS default values as given for Xsens MTW2-3A7G6
    bool run(bool visualizeResults, SimTK::Vector_<double> processCovScales = SimTK::Vector_<double>()) SWIG_DECLARE_EXCEPTION;
    bool run() override SWIG_DECLARE_EXCEPTION { 
        return run(false);
    };

    static TimeSeriesTable_<SimTK::Vec3>
        loadMarkersFile(const std::string& markerFile);

    void runInverseKinematicsWithOrientationsFromFile(Model& model,
            const std::string& quaternionStoFileName, bool visualizeResults=false, SimTK::Vector_<double> processCovScales = SimTK::Vector_<double>());

    //template <class T>
    void UKFTool(Model& model, int nqf, int nuf, std::map<int, int> yMapFromSimbodyToEigen, std::map<int, int> yMapFromEigenToSimbody, 
            std::map<int, std::string> yMapFromSimbodyToOpenSim, std::map<int, int> oMapFromDataToModel, 
            std::queue<std::vector<Eigen::MatrixXd>>* priorStatsBuffer, std::mutex* fwdBwdMutex, 
            std::condition_variable* condVar, bool* fwdDone, SimTK::State& s, OpenSim::AnalysisSet& analysisSet,
            OpenSim::OrientationsReference oRefs, OpenSim::InverseKinematicsSolver& ikSolver,
            std::shared_ptr<OpenSim::TimeSeriesTable> modelOrientationErrors, bool visualizeResults, SimTK::Array_<double> orientationErrors,
            SimTK::Vector_<double> processCovScales = SimTK::Vector_<double>());

    std::tuple<std::map<std::string, int>, std::map<int, std::string>> CreateYMaps(Model model);

    double computeFactorial(int input);

    template <typename T>
    void deletePointers(std::vector<T*>& vec);

    //template <class T>
    void computeBackwardPass(std::queue<std::vector<Eigen::MatrixXd>>* priorStatsBuffer, std::mutex* fwdBwdMutex, 
    std::condition_variable* condVar, bool* fwdDone, std::map<int, int> yMapFromEigenToSimbody, 
    std::map<int, std::string> yMapFromSimbodyToOpenSim, int nqf, int nuf);

private:
    void constructProperties();
    

//=============================================================================
};  // END of class IMUInverseKinematicsTool
//=============================================================================
} // namespace

#endif // OPENSIM_IMU_INVERSE_KINEMATICS_TOOL_H_
