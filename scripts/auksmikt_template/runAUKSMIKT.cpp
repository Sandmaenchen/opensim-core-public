
#include "OpenSim/Common/STOFileAdapter.h"
#include "OpenSim/Simulation/OpenSense/OpenSenseUtilities.h"
#include "OpenSim/Simulation/OpenSense/IMUPlacer.h"
#include "OpenSim/Simulation/Model/Model.h"
#include "OpenSim/Tools/IMUInverseKinematicsTool.h"

#include <string>
#include <iostream>
#include <clocale>
#include <chrono> 
#include "OpenSim/Tools/UKFIMUInverseKinematicsTool.h"

int main() {
    
    // Inverse kinematics wit AUKSMIKT

    // path to MS model file
    auto modelFile = std::string("<path/modelfilename>.osim");
    
    // path to file containing IMU orientations as quaternions
    auto orientationsFile = std::string("<path/orientationfilename>.sto");

    // path to the results directory
    auto resultsDir = std::string("<path>/");

    // name for output filename; the estimated state mean and covariance are written to subdirectory called 'ukf'
    auto outputMotionFile = std::string("<filename>.sto");

    // create the tool
    auto ukfIK = OpenSim::UKFIMUInverseKinematicsTool();

    // configuration
    ukfIK.set_calibrate(true);                      // if true, runs IMUPlacer assuming the first time frame is in default model position
    ukfIK.set_abort_if_diverging(false);            // aborts the run if the AUKS solution starts to diverge based on difference to observations
    ukfIK.set_accuracy(0.0001);                     // accuracy used to solve the first time frame with Least Squares IK tool
    ukfIK.set_model_file(modelFile);
    ukfIK.set_orientations_file(orientationsFile);
    ukfIK.set_output_motion_file(outputMotionFile);
    ukfIK.set_results_directory(resultsDir);
    ukfIK.set_report_errors(false);
    ukfIK.set_sensor_to_opensim_rotations(
        SimTK::Vec3(-3.1416/2, 0, 0));              // set the relative orientation between IMUs and OpenSim ((-pi/2, 0, 0) for XSens Awinda MTw)
    ukfIK.set_base_imu_label("pelvis_imu");         // choose reference IMU
    ukfIK.set_base_heading_axis("-z");              // which axis of reference IMU is pointing towards 'front' in OpenSim (+x axis)

    bool visualizeResults = false;
    bool writeUKF = true;                           // if true, writes the estimated state mean and covariance files to the 'ukf' subdirectory
    bool enableClamping = false;                    // if true, enforces the solution to follow inequality constraints (clamping)
    bool zeroObsNoise = false;                      // if true, enforces the mean of observation noise is always zero
    bool zeroProcessNoise = true;                   // if true, enforces the mean of process noise is always zero
    bool whiteProcessNoise = true;                  // if true, enforces zero cross-covariances between different coordinates in the process noise covariance
    bool independentSensors = false;                // if true, enforces zero cross-covariances between different sensors in the observation noise covariance
    int order = 2;                                  // highest order of time derivatives used
    int lagLength = 6;                              // number of samples to use in backward smoothing
    int numCPUCores = 7;                            // number of threads used in the threadpool

    double alpha = 1.0;                             // UKF hyperparameter
    double beta = 2.0;                              // UKF hyperparameter
    double kappa = -1.337;                          // UKF hyperparameter; value -1.337 sets it equal to 3 - length_of_state_vector
    double sgma2w0 = std::pow(2.0, 10.0);           // initial value for process noise variances
    double sgma2wMin = std::pow(2.0, 10.0);         // minimum value for process noise variances
    double sgma2wMax = 3.0 * std::pow(10.0, 7.0);   // maximum value for process noise variances
    double processForgetFactor = 0.1;               // update rate for process noise
    double observationForgetFactor = 0.001;         // update rate for observation noise
    double missingDataScale = 1.0;                  // scaling factor for observation noise covariance components in case of missing data (IN PROGRESS; DO NOT CHANGE)
    double rollRMSinDeg = 0.75;                     // IMU roll RMS error in degrees
    double headingRMSinDeg = 1.50;                  // IMU heading/yaw RMS error in degrees
    double pitchRMSinDeg = 0.75;                    // IMU pitch RMS error in degrees
    auto imuRMSinDeg = SimTK::Vec3(rollRMSinDeg, headingRMSinDeg, pitchRMSinDeg);
    int processCovMethod = 0;                       // method to compute process noise covariance components (OBSOLETE; DO NOT CHANGE)
    auto processCovScales = SimTK::Vector(order+1, 
        1.0);                                       // optional scaling for different time derivative components

    ukfIK.set_alpha(alpha);
    ukfIK.set_beta(beta);
    ukfIK.set_kappa(kappa);
    ukfIK.set_order(order);
    ukfIK.set_lag_length(lagLength);
    ukfIK.set_num_threads(numCPUCores);
    ukfIK.set_processForgetFactor(processForgetFactor);
    ukfIK.set_observationForgetFactor(observationForgetFactor);
    ukfIK.set_write_UKF(writeUKF);
    ukfIK.set_sgma2w_0(sgma2w0);
    ukfIK.set_sgma2w_min(sgma2wMin);
    ukfIK.set_sgma2w_max(sgma2wMax);
    ukfIK.set_missing_data_scale(missingDataScale);
    ukfIK.set_imu_RMS_in_deg(imuRMSinDeg);
    ukfIK.set_enable_clamping(enableClamping);
    ukfIK.set_enforce_independent_sensors(independentSensors);
    ukfIK.set_enforce_white_process_noise(whiteProcessNoise);
    ukfIK.set_process_noise_zero(zeroProcessNoise);
    ukfIK.set_observation_noise_zero(zeroObsNoise);
    ukfIK.set_process_covariance_method(processCovMethod);

    ukfIK.run(visualizeResults, processCovScales);

    return 0;
}