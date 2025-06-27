
#include "AugAUKSMIKT.h"


// AugAUKSMIKT methods
OpenSim::AugAUKSMIKT::AugAUKSMIKT()
        : OpenSim::InverseKinematicsToolBase() {
    OpenSim::AugAUKSMIKT::constructProperties();
}

OpenSim::AugAUKSMIKT::AugAUKSMIKT(const std::string& setupFile)
        : OpenSim::InverseKinematicsToolBase(setupFile, true) {
    OpenSim::AugAUKSMIKT::constructProperties();
    updateFromXMLDocument();
}

OpenSim::AugAUKSMIKT::~AugAUKSMIKT()
{
}

void OpenSim::AugAUKSMIKT::constructProperties()
{
    OpenSim::AugAUKSMIKT::constructProperty_sensor_to_opensim_rotations(
            SimTK::Vec3(0));
    OpenSim::AugAUKSMIKT::constructProperty_base_imu_label("");
    OpenSim::AugAUKSMIKT::constructProperty_base_heading_axis("");
    OpenSim::AugAUKSMIKT::constructProperty_orientations_file("");
    OpenSim::OrientationWeightSet orientationWeights = OpenSim::OrientationWeightSet();
    OpenSim::AugAUKSMIKT::constructProperty_orientation_weights(orientationWeights);
    OpenSim::AugAUKSMIKT::constructProperty_calibrate(false);
    OpenSim::AugAUKSMIKT::constructProperty_alpha(1.0);
    OpenSim::AugAUKSMIKT::constructProperty_beta(2.0);
    OpenSim::AugAUKSMIKT::constructProperty_kappa(-1.337);
    OpenSim::AugAUKSMIKT::constructProperty_sgma2w_min(256.0);
    OpenSim::AugAUKSMIKT::constructProperty_sgma2w_max(1048576.0);
    OpenSim::AugAUKSMIKT::constructProperty_sgma2w_0(1024.0);
    OpenSim::AugAUKSMIKT::constructProperty_processForgetFactor(0.05);
    OpenSim::AugAUKSMIKT::constructProperty_observationForgetFactor(0.0005);
    OpenSim::AugAUKSMIKT::constructProperty_order(2);
    OpenSim::AugAUKSMIKT::constructProperty_lag_length(5);
    OpenSim::AugAUKSMIKT::constructProperty_missing_data_scale(1.0);
    OpenSim::AugAUKSMIKT::constructProperty_num_threads(3);
    OpenSim::AugAUKSMIKT::constructProperty_write_UKF(true);
    OpenSim::AugAUKSMIKT::constructProperty_enable_clamping(true);
    OpenSim::AugAUKSMIKT::constructProperty_enforce_independent_sensors(true);
    OpenSim::AugAUKSMIKT::constructProperty_enforce_white_process_noise(true);
    OpenSim::AugAUKSMIKT::constructProperty_process_noise_zero(true);
    OpenSim::AugAUKSMIKT::constructProperty_observation_noise_zero(false);
    OpenSim::AugAUKSMIKT::constructProperty_abort_if_diverging(false);
    OpenSim::AugAUKSMIKT::constructProperty_process_covariance_method(0);
    OpenSim::AugAUKSMIKT::constructProperty_imu_RMS_in_deg(SimTK::Vec3(0));
        
}

void OpenSim::AugAUKSMIKT::runInverseKinematicsWithOrientationsFromFile(
        OpenSim::Model& model, const std::string& orientationsFileName, bool visualizeResults, SimTK::Vector_<double> processCovScales) {

    // Ideally if we add a Reporter, we also remove it at the end for good hygiene but 
    // at the moment there's no interface to remove Reporter so we'll reuse one if exists
    const auto reporterExists = model.findComponent<OpenSim::TableReporter>("ik_reporter");

    bool reuse_reporter = true;
    OpenSim::TableReporter* ikReporter = nullptr;
    if (reporterExists == nullptr) {
        // Add a reporter to get IK computed coordinate values out
        ikReporter = new OpenSim::TableReporter();
        ikReporter->setName("ik_reporter");
        reuse_reporter = false;
    } 
	else {
		ikReporter = &model.updComponent<OpenSim::TableReporter>("ik_reporter");
	}
	
    auto coordinates = model.updComponentList<OpenSim::Coordinate>();

    // Hookup reporter inputs to the individual coordinate outputs
    // and lock coordinates that are translational since they cannot be
    for (auto& coord : coordinates) {
        ikReporter->updInput("inputs").connect(
                coord.getOutput("value"), coord.getName());
        if (coord.getMotionType() == OpenSim::Coordinate::Translational) {
            coord.setDefaultLocked(true);
        }
    }

    if (!reuse_reporter) {
        model.addComponent(ikReporter);
    }
    OpenSim::TimeSeriesTable_<SimTK::Quaternion> quatTable(orientationsFileName);
    log_info("Loading orientations as quaternions from '{}'...",
        orientationsFileName);
    // Will maintain only data in time range specified by the tool
    // If unspecified {-inf, inf} no trimming is done
    quatTable.trim(getStartTime(), getEndTime());
    // Convert to OpenSim Frame
    const SimTK::Vec3& rotations = OpenSim::AugAUKSMIKT::get_sensor_to_opensim_rotations();
    SimTK::Rotation sensorToOpenSim = SimTK::Rotation(
            SimTK::BodyOrSpaceType::SpaceRotationSequence, 
            rotations[0], SimTK::XAxis, rotations[1], SimTK::YAxis, 
            rotations[2], SimTK::ZAxis);

    // Rotate data so Y-Axis is up
    OpenSim::OpenSenseUtilities::rotateOrientationTable(quatTable, sensorToOpenSim);
    //Trim to time window required by Tool
    quatTable.trim(getStartTime(), getEndTime());

    OpenSim::TimeSeriesTable_<SimTK::Rotation> orientationsData =
        OpenSim::OpenSenseUtilities::convertQuaternionsToRotations(quatTable);

    OpenSim::OrientationsReference oRefs(orientationsData, &get_orientation_weights());

    SimTK::Array_<OpenSim::CoordinateReference> coordinateReferences;


    // visualize for debugging
    //if (visualizeResults)
    //    model.setUseVisualizer(true);
    SimTK::State& s0 = model.initSystem();

    OpenSim::AnalysisSet& analysisSet = model.updAnalysisSet();
    analysisSet.begin(s0);


    double t0 = s0.getTime();

    // create the solver given the input data
    const double accuracy = 1e-4;
    OpenSim::InverseKinematicsSolver ikSolver(model, nullptr,
            std::make_shared<OpenSim::OrientationsReference>(oRefs),
        coordinateReferences);
    ikSolver.setAccuracy(accuracy);

    auto& times = oRefs.getTimes();
    std::shared_ptr<OpenSim::TimeSeriesTable> modelOrientationErrors(
            get_report_errors() ? new OpenSim::TimeSeriesTable()
                                : nullptr);
    s0.updTime() = times[0];
    ikSolver.assemble(s0);
    ikSolver.track(s0); // solve the initial state; assemble() was called before
	log_info("Solved at time: {} s", times[0]);
    // Create place holder for orientation errors, populate based on user pref.
    // according to report_errors property
    int nos = ikSolver.getNumOrientationSensorsInUse();
    SimTK::Array_<double> orientationErrors(nos, 0.0);

    if (get_report_errors()) {
        SimTK::Array_<std::string> labels;
        for (int i = 0; i < nos; ++i) {
            labels.push_back(ikSolver.getOrientationSensorNameForIndex(i));
        }
        modelOrientationErrors->setColumnLabels(labels);
        modelOrientationErrors->updTableMetaData().setValueForKey<std::string>(
                "name", "OrientationErrors");
        ikSolver.computeCurrentOrientationErrors(orientationErrors);
    }
    //if (visualizeResults) {
    //    model.getVisualizer().show(s0);
    //    model.getVisualizer().getSimbodyVisualizer().setShowSimTime(true);
    //}

    // Solve the states with Unscented Kalman Filter
	//int step = 0;	
	
	//log_info("Managed to get to UKFTool.");

    // Eigen::MatrixXd meanVec;
    // Eigen::MatrixXd priorCovMatrix;
    // Eigen::MatrixXd stateCrossCovMatrix;
    // std::vector<Eigen::MatrixXd> backwardsPassElement;
    
    std::mutex* fwdBwdMutex = new std::mutex();
    std::condition_variable* condVarB = new std::condition_variable();
    bool* fwdDone = new bool;
    *(fwdDone) = false; // Signal when the producer is done
    std::queue<std::vector<Eigen::MatrixXd>>* priorStatsBuffer = new std::queue<std::vector<Eigen::MatrixXd>>();

    // Make mapping between the IMUs in data file and in the model
    int const ny = ikSolver.getNumOrientationSensorsInUse();    //number of sensors
    SimTK::Array_<std::string> oRefNames = oRefs.getNames();
    SimTK::Array_<std::string> oModelNames;
    for (int ii = 0; ii < ny; ii++) { 
        oModelNames.push_back(ikSolver.getOrientationSensorNameForIndex(ii));
    }
    std::map<int, int> oMapFromDataToModel;
    std::map<int, int> oMapFromModelToData;
    for (int ii = 0; ii < oRefs.getNumRefs(); ii++) {
        for (int jj = 0; jj < ny; jj++) {
            if (oRefNames[ii] == oModelNames[jj]) {
                oMapFromDataToModel.insert(std::pair<int, int>(ii, jj));
                oMapFromModelToData.insert(std::pair<int, int>(jj, ii));
                break;
            }
        }
    }

    // Create mappings between continuous state variables between OpenSim and Simbody
    std::tuple<std::map<std::string, int>, std::map<int, std::string>> mappings = OpenSim::AugAUKSMIKT::CreateYMaps(model);
    std::map<std::string, int> yMapFromOpenSimToSimbody = std::get<0>(mappings);
    std::map<int, std::string> yMapFromSimbodyToOpenSim = std::get<1>(mappings);

    // Find locked coordinates (we cannot change these in UKF)

    std::map<int, int> simbodyIndexTypeMap;
    std::map<int, int> yMapFromSimbodyToEigen;
    std::map<int, int> yMapFromEigenToSimbody;

    int iqx;
    int nqf = 0;    //number of free and clamped q's
    int nuf = 0;
    double deltaTime = 1.0 / oRefs.getSamplingFrequency();
    for (const auto& coord : model.getComponentList<OpenSim::Coordinate>()) {
        iqx = yMapFromOpenSimToSimbody[coord.getStateVariableNames()[0]];
        if (coord.getMotionType() == OpenSim::Coordinate::MotionType::Translational) {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 0));    //translational coordinate
            log_info("{} is translational coordinate, will be ignored in UKF.", coord.getName());
        }
        else if (coord.getLocked(s0)) {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 0));    //locked coordinate
            log_info("{} is locked, will be ignored in UKF.", coord.getName());
        }
        else if (coord.isDependent(s0)) {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 0));    //dependent coordinate
            log_info("{} is dependent coordinate, will be ignored in UKF.", coord.getName());
        }
        else if (coord.getClamped(s0)) {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 1));    //clamped coordinate
            nqf++;
        }
        else {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 2));    //free coordinate
            nqf++;
        }
    }
    nuf = nqf;
    log_info("{} coordinates will be used in UKF", nqf);

    log_info("Eigen world version is {}", EIGEN_WORLD_VERSION);
    log_info("Eigen major version is {}", EIGEN_MAJOR_VERSION);
    log_info("Eigen minor verison is {}", EIGEN_MINOR_VERSION);

    int ie = 0;
    for (std::map<int, int>::iterator it = simbodyIndexTypeMap.begin(); it != simbodyIndexTypeMap.end(); ++it) {
        if (it->second != 0) {
            yMapFromSimbodyToEigen.insert(std::pair<int, int>(it->first, ie));
            yMapFromEigenToSimbody.insert(std::pair<int, int>(ie, it->first));
            ie++;
        }
    }

    // Vector of structs holding value limits for clamped coordinates
    std::vector<OpenSim::UKFClampedCoordLimits> clampedCoordLimits;
    for (const auto& coord : model.getComponentList<OpenSim::Coordinate>()) {
        if (coord.get_clamped()) {
            clampedCoordLimits.emplace_back(coord.getStateVariableNames()[0], coord.getRangeMin(), coord.getRangeMax());            
        }
    }

    // Construct covariance and mean for process noise
    Eigen::MatrixXd Q((nqf + get_order()*nuf), (nqf + get_order()*nuf));
    Eigen::MatrixXd w((nqf + get_order()*nuf), 1);
    w.setZero();

    Eigen::MatrixXd F((nqf + get_order()*nuf), (nqf + get_order()*nuf));
    std::queue<std::vector<Eigen::MatrixXd>>* stateMeansBuffer = new std::queue<std::vector<Eigen::MatrixXd>>();
    //std::mutex* qMutex = new std::mutex();

    // If the process covariance scale factors are not provided, use simply ones
    if (processCovScales.size() == 0) {
        processCovScales = SimTK::Vector_<double>(get_order()+1, 1.0);
        log_info("Process covariance scales not provided, using ones instead.");
    }


    // RMS errors for observation noise covariance    
    SimTK::Vec3 imuRMSinDeg = get_imu_RMS_in_deg();    
    double xAxisRMS = imuRMSinDeg(0) * (SimTK::Pi / 180); // roll RMS error (x-axis)
    double yAxisRMS = imuRMSinDeg(1) * (SimTK::Pi / 180); // heading RMS error (y-axis)
    double zAxisRMS = imuRMSinDeg(2) * (SimTK::Pi / 180); // pitch RMS error (z-axis)

    // Construct mean vector for observation errors (gets updated)
    Eigen::MatrixXd v(3 * ny, 1);
    v.setZero();
    std::vector<Eigen::Quaternion<double>> vQuat;    
    for (int idx = 0; idx < ny; idx++) {
        vQuat.emplace_back(Eigen::Quaternion<double>(1, 0, 0, 0));
    }    

    // Construct covariance matrix for observation errors (gets updated)
    Eigen::MatrixXd R(3 * ny, 3 * ny);
    R.setZero();

    for (int ii = 0; ii < (3 * ny); ii++) {
        if (ii % 3 == 0) {      //x-axis
            R(ii, ii) = std::pow(xAxisRMS, 2);
        } 
		else if (ii % 3 == 1) { // y-axis
            R(ii, ii) = std::pow(yAxisRMS, 2);
        } 
		else {                //z-axis
            R(ii, ii) = std::pow(zAxisRMS, 2);        
        }
    }
    
    std::thread forwardThread([&] {OpenSim::AugAUKSMIKT::UKFTool(
        nqf, nuf, w, Q, v, vQuat, R, F, yMapFromSimbodyToEigen, yMapFromEigenToSimbody, yMapFromSimbodyToOpenSim, yMapFromOpenSimToSimbody, oMapFromDataToModel, clampedCoordLimits, 
        priorStatsBuffer, fwdBwdMutex, stateMeansBuffer, condVarB, fwdDone, s0, oRefs, ikSolver, 
        modelOrientationErrors, visualizeResults, orientationErrors, processCovScales);});
    
    
    std::thread backwardThread([&] {OpenSim::AugAUKSMIKT::computeBackwardPass(model, clampedCoordLimits, 
        priorStatsBuffer,fwdBwdMutex, condVarB, fwdDone, yMapFromEigenToSimbody, yMapFromSimbodyToEigen, yMapFromOpenSimToSimbody, analysisSet, yMapFromSimbodyToOpenSim, nqf, nuf);});
    

    forwardThread.join();
    backwardThread.join();
    
    log_info("Threads finished");

    delete fwdBwdMutex;
    delete condVarB;
    delete fwdDone;
    delete priorStatsBuffer;
    delete stateMeansBuffer;
    log_info("Deleted dynamically allocated stuff");

    /*
    for (auto time : times) {
        s0.updTime() = time;
        ikSolver.track(s0);
        if (get_report_errors()) {
            ikSolver.computeCurrentOrientationErrors(orientationErrors);
            modelOrientationErrors->appendRow(
                    s0.getTime(), orientationErrors);
        }
        if (visualizeResults)  
            model.getVisualizer().show(s0);
        else
            log_info("Solved at time: {} s", time);
        // realize to report to get reporter to pull values from model
        analysisSet.step(s0, step++);
        model.realizeReport(s0);
    }
    */

    auto report = ikReporter->getTable();
    // form resultsDir either from results_directory or output_motion_file
    auto resultsDir = get_results_directory();
    if (resultsDir.empty() && !get_output_motion_file().empty())
        resultsDir = OpenSim::IO::getParentDirectory(get_output_motion_file());
    if (!resultsDir.empty()) {
        OpenSim::IO::makeDir(resultsDir);
        // directory will be restored on block exit
        // by changing dir all other files are created in resultsDir
        auto cwd = OpenSim::IO::CwdChanger::changeTo(resultsDir);
        std::string outName = get_output_motion_file();
        outName = OpenSim::IO::GetFileNameFromURI(outName);
        if (outName.empty()) {
            bool isAbsolutePath;
            std::string directory, fileName, extension;
            SimTK::Pathname::deconstructPathname(orientationsFileName,
                    isAbsolutePath, directory, fileName, extension);
            outName = "ik_" + fileName;
        }
        std::string outputFile = outName;

        // Convert to degrees to compare with marker-based IK
        // but only for rotational coordinates
        model.getSimbodyEngine().convertRadiansToDegrees(report);
        report.updTableMetaData().setValueForKey<std::string>("name", outName);

        auto fullOutputFilename = outputFile;
        std::string::size_type extSep = fullOutputFilename.rfind(".");
        if (extSep == std::string::npos) { fullOutputFilename.append(".mot"); }
        OpenSim::STOFileAdapter_<double>::write(report, fullOutputFilename);

        log_info("Wrote IK with IMU tracking results to: '{}'.",
                fullOutputFilename);
        if (get_report_errors()) {
            OpenSim::STOFileAdapter_<double>::write(*modelOrientationErrors,
                    outName + "_orientationErrors.sto");
        }
    } 
    else
        log_info("AugAUKSMIKT: No output files were generated, "
            "set output_motion_file to generate output files.");
    // Results written to file, clear in case we run again
    ikReporter->clearTable();
}


// main driver
bool OpenSim::AugAUKSMIKT::run(bool visualizeResults, SimTK::Vector_<double> processCovScales)
{
    if (_model.empty()) {
        _model.reset(new OpenSim::Model(get_model_file()));
    }
    if (get_calibrate() == true) {
        OpenSim::IMUPlacer imuPlacer = OpenSim::IMUPlacer();
        _model->updForceSet().clearAndDestroy();
        _model->updControllerSet().clearAndDestroy();
        imuPlacer.setModel(*_model);
        imuPlacer.set_base_imu_label(get_base_imu_label());
        imuPlacer.set_base_heading_axis(get_base_heading_axis());
        imuPlacer.set_sensor_to_opensim_rotations(get_sensor_to_opensim_rotations());
        imuPlacer.set_orientation_file_for_calibration(get_orientations_file());
        bool success = imuPlacer.run();
        if (success) {
            log_info("managed to calibrate");
        }
        else {
            log_info("failed to calibrate");
        }
        log_info("trying to assign calibrated model next.");
        //(*_model) = imuPlacer.getCalibratedModel();
        _model.reset(new OpenSim::Model(imuPlacer.getCalibratedModel()));
        log_info("managed to assign calibrated model.");
        _model->finalizeFromProperties();
        log_info("model: finalized from properties.");
    }
    else {
        _model->updForceSet().clearAndDestroy();
        _model->updControllerSet().clearAndDestroy();
        try {
            _model->finalizeFromProperties();
        }
        catch(OpenSim::Exception &ex) {
            log_error("Could not finalize model from properties.");
            throw(ex);
        }        
    }

    OpenSim::AugAUKSMIKT::runInverseKinematicsWithOrientationsFromFile(*_model,
            get_orientations_file(), visualizeResults, processCovScales);

    return true;
}

OpenSim::TimeSeriesTable_<SimTK::Vec3> OpenSim::AugAUKSMIKT::loadMarkersFile(const std::string& markerFile)
{
    OpenSim::TimeSeriesTable_<SimTK::Vec3> markers(markerFile);
    log_info("'{}' loaded {} markers and {} rows of data.", markerFile,
        markers.getNumColumns(), markers.getNumRows());

    if (markers.hasTableMetaDataKey("Units")) {
        auto& value = markers.getTableMetaData().getValueForKey("Units");
        log_info("'{}' has Units meta data. Units are {}.", markerFile,
                value.getValue<std::string>());
        if (value.getValue<std::string>() == "mm") {
            log_info("Marker data in mm, converting to m.");
            for (size_t i = 0; i < markers.getNumRows(); ++i) {
                markers.updRowAtIndex(i) *= 0.001;
            }
            markers.updTableMetaData().removeValueForKey("Units");
            markers.updTableMetaData().setValueForKey<std::string>("Units", "m");
        }
    }
    auto& value = markers.getTableMetaData().getValueForKey("Units");
    log_info("'{}' Units are {}.", markerFile, value.getValue<std::string>());

    return markers;
}

// The actual workhorse of UKF-IK
//template <class T>
void OpenSim::AugAUKSMIKT::UKFTool(int nqf, int nuf, Eigen::MatrixXd& w, Eigen::MatrixXd& Q, Eigen::MatrixXd& v, 
        std::vector<Eigen::Quaternion<double>>& vQuat, Eigen::MatrixXd& R, Eigen::MatrixXd& F, std::map<int, int> yMapFromSimbodyToEigen, 
        std::map<int, int> yMapFromEigenToSimbody, std::map<int, std::string> yMapFromSimbodyToOpenSim, std::map<std::string, int> yMapFromOpenSimToSimbody,
        std::map<int, int> oMapFromDataToModel, std::vector<UKFClampedCoordLimits> clampedCoordLimits, std::queue<std::vector<Eigen::MatrixXd>>* priorStatsBuffer, 
        std::mutex* fwdBwdMutex, std::queue<std::vector<Eigen::MatrixXd>>* stateMeansBuffer, std::condition_variable* condVarB, bool* fwdDone, SimTK::State& s, 
        OpenSim::OrientationsReference oRefs, OpenSim::InverseKinematicsSolver& ikSolver,
        std::shared_ptr<OpenSim::TimeSeriesTable> modelOrientationErrors, bool visualizeResults,
        SimTK::Array_<double> orientationErrors, SimTK::Vector_<double> processCovScales) {

	log_info("Got inside UKFTOOL");

    SimTK::State ss;

    {
        std::unique_lock<std::mutex> lock(*fwdBwdMutex);
        ss = SimTK::State(s);
    }

    auto times = oRefs.getTimes();
    double alpha = get_alpha();
    double beta = get_beta();
    double kappa = get_kappa();
    int order = get_order();
    double sgma2w0 = get_sgma2w_0();
    double missingDataScale = get_missing_data_scale();

    
    int num_cores = get_num_threads();
    
   
    SimTK::Array_<SimTK::Rotation_<double>>
            osensorOrientations; // array for orientations computed by the
                                 // model; REMOVE _<double> IF CAUSES ERROR
    std::vector<SimTK::Array_<SimTK::Rotation_<double>>*> arr_osensorOrientations;
    SimTK::Array_<SimTK::Rotation_<double>>
            yArray; // array for observations (IMU orientations)
    //SimTK::SimbodyMatterSubsystem matterSubSys = model.getMatterSubsystem();    
    
    int const nq = ss.getNQ(); // number of generalized positions (joint angles)
    int const nu = ss.getNU(); // number of generalized velocities (joint angular
                        // velocities)
    int const nr = std::min(nq, nu); // probably not needed, usually nq >= nu

    // Coefficients for process model f()
    double deltaTime = 1.0 / oRefs.getSamplingFrequency();
    Eigen::MatrixXd fCoeffs(order+1, order+1);
    fCoeffs.setZero();
    for (int irow = 0; irow <= order; irow++) {
        for (int icol = irow; icol <= order; icol++) {
            double denum = OpenSim::AugAUKSMIKT::computeFactorial(icol-irow);
            fCoeffs(irow, icol) = std::pow(deltaTime, (icol-irow)) / denum;
        }
    }

    // Coefficients for process noise covariance matrix Q
    Eigen::MatrixXd QCoeffs(order+1, order+1);

    // Coefficients for process noise correction matrix Qw
    Eigen::MatrixXd QwCoeffs(order+1, order+1);
    Eigen::MatrixXd wCoeffs(order+1, 1);
    
    if (get_process_covariance_method() == 0) {
    // classic approach of Fioretti and Jetto, 1989 (no scaling tricks by default; should use ones)
        log_info("Using method of Fioretti and Jetto, 1989.");
        for (int irow = 0; irow <= order; irow++) {
            for (int icol = 0; icol <= order; icol++) {
                double denum1 = OpenSim::AugAUKSMIKT::computeFactorial(order-irow);
                double denum2 = OpenSim::AugAUKSMIKT::computeFactorial(order-icol);
                int deltaPower = (order-irow) + (order-icol);
                QCoeffs(irow, icol) = processCovScales(irow) * processCovScales(icol) * std::pow(deltaTime, (deltaPower+1)) / (denum1 * denum2 * (deltaPower+1));
            }
        }
        //for w
        for (int irow = 0; irow <= order; irow++) {
            double denum = OpenSim::AugAUKSMIKT::computeFactorial(order-irow);
            int deltaPower = order+1-irow;
            wCoeffs(irow, 0) = processCovScales(irow) * std::pow(deltaTime, deltaPower) / (denum * deltaPower);
        }
        QwCoeffs = wCoeffs * wCoeffs.transpose();
    }

    else if (get_process_covariance_method() == 1) {
    // approach using the Taylor remainders as white noise multipliers (plus additional scale factors)        
        log_info("Using Taylor remainders.");
        for (int irow = 0; irow <= order; irow++) {
            for (int icol = 0; icol <= order; icol++) {
                double denum1 = OpenSim::AugAUKSMIKT::computeFactorial(order+1-irow);
                double denum2 = OpenSim::AugAUKSMIKT::computeFactorial(order+1-icol);
                int deltaPower = (order+1-irow) + (order+1-icol);
                QCoeffs(irow, icol) = processCovScales(irow) * processCovScales(icol) *
                                    std::pow(deltaTime, (deltaPower)) / (denum1 * denum2);
            }
        }
    }

    // Construct covariance matrix Q for process noise
    //Eigen::MatrixXd Q((nqf + order*nuf), (nqf + order*nuf));

    Q.setZero();
    for (int irow = 0; irow <= order; irow++) {
        for (int icol = 0; icol <= order; icol++) {
            Q.block((irow*nuf), (icol*nuf), nuf, nuf) = sgma2w0 * QCoeffs(irow, icol) * Eigen::MatrixXd::Identity(nuf, nuf);
        }
    }
    Eigen::MatrixXd Qw((nqf + order*nuf), (nqf + order*nuf));
    
    // Store the original observation noise covariance
    int const ny = ikSolver.getNumOrientationSensorsInUse();
    Eigen::MatrixXd R0(3 * ny, 3 * ny);
    R0 = R;

    // Construct covariance matrix of state
    Eigen::MatrixXd P((nqf + (order*nuf)), (nqf + (order*nuf)));
    Eigen::MatrixXd P0((nqf + (order*nuf)), (nqf + (order*nuf)));
    P = Q;
    
    // Construct the state vector
    Eigen::MatrixXd x((nqf + (order*nuf)), 1);
    Eigen::MatrixXd x0((nqf + (order*nuf)), 1);
    Eigen::MatrixXd xx((nqf + (order*nuf)), 1);
    std::vector<Eigen::MatrixXd*> arr_xx;
    std::vector<std::vector<Eigen::Quaternion<double>>> arr_vQuat;
    //SimTK::Vector_<double> x(nq + nu);
    x.setZero();
    xx.setZero();
    //Eigen::MatrixXd qf(nqf, 1);
    //qf.setZero();
    SimTK::Vector_<double> q(nq);
    std::vector<SimTK::Vector_<double>*> arr_q;
    q = ss.getQ();
    SimTK::Vector_<double> u(nu);
    std::vector<SimTK::Vector_<double>*> arr_u;
    u = ss.getU();
    for (std::map<int, int>::iterator it = yMapFromSimbodyToEigen.begin(); it != yMapFromSimbodyToEigen.end(); ++it) {
        x(it->second) = q(it->first);
    }
    for (std::map<int, int>::iterator it = yMapFromSimbodyToEigen.begin();
            it != yMapFromSimbodyToEigen.end(); ++it) {
        x((it->second)+nqf) = u(it->first);
    }

    // Construct the process model matrix
    //Eigen::MatrixXd F((nqf + order*nuf), (nqf + order*nuf));
    F.setZero();
    for (int irow = 0; irow <= order; irow++) {
        for (int icol = 0; icol <= order; icol++) {
            F.block((irow*nuf), (icol*nuf), nuf, nuf) = fCoeffs(irow, icol) * Eigen::MatrixXd::Identity(nuf, nuf);
        }
    }

    // For data
    //SimTK::Vector_<double> y(3 * ny);   //probably not needed..
    Eigen::MatrixXd y(3 * ny, 1);
    std::vector<Eigen::MatrixXd*> arr_y;

    // If the user wants kappa to be equal to the augmented state vector length (or 3-length)
    if (std::abs(kappa+1.337) < 1e-5 ) {
        //kappa = -1.0 * ((double) nqf + (order*nuf) - 3);
        kappa = -1.0 * ((double) nqf + (order*nuf) + (3 * ny) - 3);
    }
    else if (std::abs(kappa+4.337) < 1e-5) {        
        //kappa = (double) nqf + (order*nuf);
        kappa = (double) nqf + (order*nuf) + (3 * ny);
    }
    log_info("kappa set to {}", kappa);

    double const lambda = (std::pow(alpha, 2) * (nqf + (order*nuf) + (3 * ny) + kappa)) - (nqf + (order*nuf) + (3 * ny));
    double const W0m = lambda / (nqf + (order*nuf) + (3 * ny) + lambda);
    double const W0c = W0m + (1 - std::pow(alpha, 2) + beta);
    double const Wi = 1 / (2 * (nqf + (order*nuf) + (3 * ny) + lambda));
    Eigen::VectorXd Wm(2 * (nqf + (order*nuf) + (3 * ny)) + 1);
    Eigen::VectorXd Wc(2 * (nqf + (order*nuf) + (3 * ny)) + 1);
    Wm.setConstant(Wi);
    Wc.setConstant(Wi);
    Wm(0) = W0m;
    Wc(0) = W0c;

    Eigen::MatrixXd SSx(nqf + (order*nuf), nqf + (order*nuf));
    Eigen::MatrixXd SSy((3 * ny), (3 * ny));
    Eigen::MatrixXd Sigmas(nqf + (order*nuf), 2 * (nqf + (order*nuf) + (3 * ny)) + 1);
    Eigen::MatrixXd Sigmas2(nqf + (order*nuf), 2 * (nqf + (order*nuf) + (3 * ny)) + 1);
    //Eigen::MatrixXd Sigmaprops(nqf + (order*nuf), 2 * (nqf + (order*nuf) + (3 * ny)) + 1);
    Eigen::MatrixXd Sigmas2props(3 * ny, 2 * (nqf + (order*nuf) + (3 * ny)) + 1);

    std::vector<SimTK::State*> arr_ss;
    SimTK::Vector_<double> qdot;
    Eigen::MatrixXd xsave(nqf + (order*nuf), 1);

    SimTK::Vector_<double> u_old(nu);
    u_old.setToZero();
    SimTK::Vector_<double> u_oldold(nu);
    u_oldold.setToZero();
    Eigen::MatrixXd ysave(3 * ny, 1);
    Eigen::MatrixXd ydiff(3 * ny, 1);
    Eigen::MatrixXd ydiff0(3 * ny, 1);
    Eigen::MatrixXd ypred(3 * ny, 1);

    int ycounter;
    Eigen::MatrixXd Py(3 * ny, 3 * ny);
    Eigen::MatrixXd Pxy(nqf + (order*nuf), 3 * ny);
    Eigen::MatrixXd K(nqf + (order*nuf), 3 * ny);
    Eigen::MatrixXd C((nqf + (order*nuf)), (nqf + (order*nuf)));
    Eigen::MatrixXd stateMean(nqf+(order*nuf), 1);
    SimTK::Matrix A;        //matrix for holonomic constraints
    SimTK::Vector_<double> qerr;
	SimTK::Vec3 angle_vector;
    SimTK::Vec3 angle_vector0;
    std::vector<SimTK::Vec3*> arr_angle_vector;
    Eigen::LLT<Eigen::MatrixXd> llt;    //construct LLT object

    // For concurrent computing
    
    if (num_cores == 0) {
        num_cores = std::thread::hardware_concurrency();
        if (num_cores == 0) {
            num_cores = 1;
            log_info("Could not compute hardware_concurrency(); using single core instead.");
        }
    }
    else if (std::thread::hardware_concurrency() > 0 && num_cores > (int)std::thread::hardware_concurrency()) {
        num_cores = std::thread::hardware_concurrency();
    }
    else if (num_cores < 0) {
        num_cores = 1;
    }
    log_info("Number of threads: {} ", num_cores);
    OpenSim::UKFThreadPool pool(num_cores);
    //std::vector<std::thread> threads;
    //threads.reserve(num_cores);
    //std::mutex modelMutex;    
    int elPerThread = 0;
    int num_cols = 0;
    std::vector<OpenSim::InverseKinematicsSolver*> solvers;
    //solvers.reserve(num_cores);
    std::vector<OpenSim::Model*> models;
    //models.reserve(num_cores);
    SimTK::Array_<OpenSim::CoordinateReference> coordRefArray;

    {
        std::unique_lock<std::mutex> lock(*fwdBwdMutex);
        for (int ii = 0; ii < num_cores; ii++) {
            //OpenSim::Model* model_clone = model.clone();
            //OpenSim::Model* model_clone = new OpenSim::Model(get_model_file());
            OpenSim::Model* model_clone = new OpenSim::Model(ikSolver.getModel());
            models.push_back(model_clone);  //Alternatively, emplace_back(), but that *should* be slower
            models[ii]->updForceSet().clearAndDestroy();
            models[ii]->updControllerSet().clearAndDestroy();
            models[ii]->initSystem();
            OpenSim::InverseKinematicsSolver* aSolver = new OpenSim::InverseKinematicsSolver(*(models[ii]), nullptr,
                    std::make_shared<OpenSim::OrientationsReference>(oRefs),
                    coordRefArray);
            aSolver->setAccuracy(1e-4);
            aSolver->assemble(ss);
            solvers.push_back(aSolver);     //Alternatively, emplace_back(), but that *should* be slower
            //delete aSolver;
            //delete model_clone;
        }
    }

    
    // Modifications to compute average orientations, comparisons between orientations, etc.
    SimTK::Quaternion_<double> simTKquat(1, 0, 0, 0);
    SimTK::Quaternion_<double> simTKquat0(1, 0, 0, 0);
    std::vector<SimTK::Quaternion_<double>*> arr_simTKquat;
    Eigen::Quaternion<double> dummyquat(1, 0, 0, 0);
    Eigen::Quaternion<double> dataMinusMeanQuat(1, 0, 0, 0);
    Eigen::Quaternion<double> dataMinusMeanQuat0(1, 0, 0, 0);
    std::vector<Eigen::Quaternion<double>> dataOVector(ny, dummyquat);
    std::vector<Eigen::Quaternion<double>> expectedOVector(ny, dummyquat);
    std::vector<Eigen::Quaternion<double>> expectedOVectorNoObsError(ny, dummyquat);
    std::vector<std::vector<Eigen::Quaternion<double>>> Sigmas2Orientations(2 * (nqf + (order*nuf) + (3 * ny)) + 1, dataOVector);
    std::vector<std::vector<Eigen::Quaternion<double>>> Sigmas2OrientationsNoObsError(2 * (nqf + (order*nuf) + (3 * ny)) + 1, dataOVector);
    std::vector<std::vector<Eigen::Quaternion<double>>> Sigmas2OminusMean(2 * (nqf + (order*nuf) + (3 * ny)) + 1, dataOVector);
    Eigen::MatrixXd M(4, 2 * (nqf + (order*nuf) + (3 * ny)) + 1);
    std::vector<Eigen::MatrixXd*> arr_M;
    Eigen::MatrixXd MM(4, 4);
    std::vector<Eigen::MatrixXd*> arr_MM;
    Eigen::EigenSolver<Eigen::MatrixXd> eigSolver;
    std::vector<Eigen::EigenSolver<Eigen::MatrixXd>*> arr_eigSolver;
    Eigen::VectorXd eigenVals(4);
    std::vector<Eigen::VectorXd*> arr_eigenVals;
    Eigen::MatrixXd eigenVecs(4, 4);
    std::vector<Eigen::MatrixXd*> arr_eigenVecs;
    Eigen::MatrixXd newSigmas(nqf, 1); newSigmas.setConstant(get_sgma2w_0());
    Eigen::MatrixXd oldSigmas(nqf, 1); oldSigmas.setConstant(get_sgma2w_0()); 
    int iMax = 0;
    double maxEigenVal = 0.0;
    SimTK::Vec4 dummy4vector = simTKquat.asVec4();
    std::vector<SimTK::Vec4*> arr_dummy4vector;
    std::vector<std::map<int, int>*> arr_yMapFromEigenToSimbody;
    int timeStep = 0;
    //Eigen::MatrixXd innovationValue(1,1); innovationValue.setZero();

    for (int ii = 0; ii < num_cores; ii++) {
        std::map<int, int>* aMap = new std::map<int,int>(yMapFromEigenToSimbody);
        arr_yMapFromEigenToSimbody.push_back(aMap);
        SimTK::State* aState = new SimTK::State(ss);
        arr_ss.push_back(aState);
        SimTK::Vector_<double>* aQ = new SimTK::Vector_<double>(ss.getQ());
        arr_q.push_back(aQ);
        SimTK::Vector_<double>* aU = new SimTK::Vector_<double>(ss.getU());
        arr_u.push_back(aU);
        Eigen::MatrixXd* aX = new Eigen::MatrixXd((nqf + (order*nuf)), 1);
        aX->setZero();
        arr_xx.push_back(aX);
        std::vector<Eigen::Quaternion<double>> obsErrQuats;    
        for (int idx = 0; idx < ny; idx++) {
            obsErrQuats.emplace_back(Eigen::Quaternion<double>(1, 0, 0, 0));
        }
        arr_vQuat.push_back(obsErrQuats);
        SimTK::Array_<SimTK::Rotation_<double>>* aOrientation = new SimTK::Array_<SimTK::Rotation_<double>>();
        arr_osensorOrientations.push_back(aOrientation);
        Eigen::MatrixXd* aY = new Eigen::MatrixXd(y);
        arr_y.push_back(aY);
        SimTK::Vec4* av4 = new SimTK::Vec4(dummy4vector);
        arr_dummy4vector.push_back(av4);
        Eigen::VectorXd* aVals = new Eigen::VectorXd(eigenVals);
        arr_eigenVals.push_back(aVals);
        Eigen::MatrixXd* aVecs = new Eigen::MatrixXd(eigenVecs);
        arr_eigenVecs.push_back(aVecs);
        Eigen::EigenSolver<Eigen::MatrixXd>* aSolver = new Eigen::EigenSolver<Eigen::MatrixXd>(eigSolver);
        arr_eigSolver.push_back(aSolver);
        Eigen::MatrixXd* aMM = new Eigen::MatrixXd(MM);
        arr_MM.push_back(aMM);
        Eigen::MatrixXd* aM = new Eigen::MatrixXd(M);
        arr_M.push_back(aM);
        SimTK::Quaternion_<double>* aQuat = new SimTK::Quaternion_<double>(simTKquat);
        arr_simTKquat.push_back(aQuat);
        SimTK::Vec3* aVec3 = new SimTK::Vec3(0,0,0);
        arr_angle_vector.push_back(aVec3);
    }

    // for noise updates
    Eigen::MatrixXd new_v(3 * ny, 1); 
    std::vector<Eigen::Quaterniond> new_vQuat(ny, dummyquat); 
    Eigen::MatrixXd new_w((nqf + (order*nuf)), 1); new_w.setZero();

    Eigen::MatrixXd z((order+1),1); z.setZero();
    Eigen::MatrixXd H((order+1),1); H.setZero(); 

    Eigen::MatrixXd newQ(nqf + (order * nuf), nqf + (order * nuf)); newQ.setZero();
    Eigen::MatrixXd tempQ(nqf + (order * nuf), nqf + (order * nuf)); newQ.setZero();
    Eigen::MatrixXd newR(3 * ny, 3 * ny); newR.setZero();
    Eigen::MatrixXd tempR(3 * ny, 3 * ny); tempR.setZero();
    Eigen::MatrixXd sigmas(nqf, 1);
    Eigen::MatrixXd vUpdateweightMat(2,2); vUpdateweightMat.setZero();
    vUpdateweightMat(0, 0) = (1.0 - get_observationForgetFactor());
    vUpdateweightMat(1, 1) = (get_observationForgetFactor());
    
	int step = 0;

    
	log_info("Got to the beginning of for loop.");

    //log_info("before for loop, x(0, 0) = {}", (double)x(0, 0));

    for (auto time : times) { 
        if (time == times[0]) { 
            std::vector<Eigen::MatrixXd> priorStatsVector;            
            Eigen::MatrixXd time_export(1,1);
            time_export(0,0) = time;
            priorStatsVector.emplace_back(time_export);

            // Apply inequality constraints (clamped coordinates)
            if (get_enable_clamping()) {                
                OpenSim::AugAUKSMIKT::clampCoordinates(x, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf);            
            }
            
            priorStatsVector.emplace_back(x);   // this won't be used in backward pass
            priorStatsVector.emplace_back(P);   // this won't be used in backward pass
            priorStatsVector.emplace_back(P);   // this won't be used in backward pass
            priorStatsVector.emplace_back(x);
            priorStatsVector.emplace_back(P);  
            {
                std::unique_lock<std::mutex> lock(*fwdBwdMutex);
                priorStatsBuffer->push(priorStatsVector);
            }
            (*condVarB).notify_one();
            //log_info("got through 1st frame");

        } 
        else {
            timeStep++;
            // Step 0. Cholesky factorization of state covariance            
            ss.updTime() = time;
            ss.updQ() = q;
            ss.updU() = u;

            // Steps 0-4 for linear model with Gaussian noise            
            C = P * F.transpose();
            //P0 = ((F) * P * (F.transpose())) - (2 * (w) * (x.transpose() * F.transpose()));
            P0 = (F * P * (F.transpose()));
            P = (F * P * (F.transpose())) + Q;
            x0 = F * x;          //a priori mean assuming zero noise
            x = F * x + w;    //a priori mean assuming non-zero noise
            //x = x0;
            //log_info("completed steps 0-4");

            // Apply inequality constraints (clamped coordinates)
            if (get_enable_clamping()) {
                OpenSim::AugAUKSMIKT::clampCoordinates(x, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf);     
                OpenSim::AugAUKSMIKT::clampCoordinates(x0, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf);         
            }
            //log_info("applied constraints");

            // Copy constrained a priori state mean to vector for buffer

            // Note! These are computed for time step 'k', not 'k+1'. Have to handle in the backward pass
            std::vector<Eigen::MatrixXd> priorStatsVector;              
            Eigen::MatrixXd time_export(1,1);
            time_export(0,0) = time;
            priorStatsVector.emplace_back(time_export);
            priorStatsVector.emplace_back(x);
            priorStatsVector.emplace_back(P);
            priorStatsVector.emplace_back(C);        
            //log_info("added stuff for bwdThread");

            // Sample the sigma points after propagating through process model (optional)
            // Due to diagonality, we can first sample for the state covariance
            llt.compute((nqf + (order*nuf) + (3 * ny) + lambda) * P);
            if (llt.info() == Eigen::Success) {
                SSx = llt.matrixL();
            }
            else {
                log_info("WARNING: state covariance matrix is NOT positive definite!");
                SSx = llt.matrixL();
                //log_info("Tried to Cholesky factorize the following matrix: \n{}", (nqf + (order*nuf) + (3 * ny) + lambda) * P);
                //break;
            }
            Sigmas2.col(0) = x;                
            for (int ii = 1; ii < (nqf + (order*nuf) + 1); ii++) {
                Sigmas2.col(ii) = x + SSx.col(ii - 1);
            }
            for (int ii = (nqf + (order*nuf) + 1); ii < (2 * (nqf + (order*nuf)) + 1); ii++) {
                Sigmas2.col(ii) = x - SSx.col(ii - (nqf + (order*nuf) + 1));
            }

            // Then, we sample for the observation noise covariance
            llt.compute((nqf + (order*nuf) + (3 * ny) + lambda) * R);
            if (llt.info() == Eigen::Success) {
                SSy = llt.matrixL();
            }
            else {
                log_info("WARNING: observation noise covariance matrix is NOT positive definite!");
                SSy = llt.matrixL();
                //log_info("Tried to Cholesky factorize the following matrix: \n{}", (nqf + (order*nuf) + (3 * ny) + lambda) * R);
                //break;
            }
            for (int ii = (2 * (nqf + (order*nuf)) + 1); ii < (2 * (nqf + (order*nuf) + (3 * ny)) + 1); ii++) {
                Sigmas2.col(ii) = x;
            }

            // Step 5. Propagate a priori sigma points of current step through
            // observation model
            xsave = x;
            num_cols = (2 * (nqf + (order*nuf) + (3 * ny)) + 1);
            elPerThread = num_cols / num_cores;
            for (int ithr = 0; ithr < num_cores; ++ithr) {
                int start = ithr * elPerThread;
                int end = (ithr == num_cores - 1) ? num_cols : (ithr + 1) * elPerThread;
                //threads.emplace_back([&, ithr, start, end, x, ss, q, u, yMapFromEigenToSimbody, osensorOrientations, dummy4vector]() mutable {
                pool.enqueue([&, ithr, start, end]() mutable {
                    for (int ii = start; ii < end; ii++) {
                        (*(arr_xx[ithr])) = Sigmas2.col(ii);
                        for (std::map<int, int>::iterator it = arr_yMapFromEigenToSimbody[ithr]->begin(); it != arr_yMapFromEigenToSimbody[ithr]->end(); ++it) {
                            (*(arr_q[ithr]))(it->second) = (*(arr_xx[ithr]))(it->first);
                        }
                        arr_ss[ithr]->updTime() = time;
                        arr_ss[ithr]->updQ() = (*(arr_q[ithr]));
                        if (order > 0) {
                            for (std::map<int, int>::iterator it =
                                    arr_yMapFromEigenToSimbody[ithr]->begin();
                                    it != arr_yMapFromEigenToSimbody[ithr]->end(); ++it) {
                                (*(arr_u[ithr]))(it->second) = (*(arr_xx[ithr]))((it->first)+nqf);
                            }
                            arr_ss[ithr]->updU() = (*(arr_u[ithr]));
                        }
                        solvers[ithr]->setState(*(arr_ss[ithr])); //method added to AssemblySolver (like in Kalman Smoother)
                        solvers[ithr]->computeCurrentSensorOrientations(*(arr_osensorOrientations[ithr]));   
                        auto noise = arr_vQuat[ithr];                        
                        if (ii >= ((2 * (nqf + (order*nuf)) + 1 + (3*ny)))) {
                            //log_info("trying to access SSy matrix column {}", (ii - ((2 * (nqf + (order*nuf)) + 1 + (3*ny)))));
                            auto sigma = SSy.col(ii - ((2 * (nqf + (order*nuf)) + 1 + (3*ny)))); //euler angles
                            for (int iy = 0; iy < ny; iy++) {
                                auto asd = SimTK::Quaternion_<double>(SimTK::Rotation(SimTK::BodyOrSpaceType::SpaceRotationSequence, 
                                (sigma)(iy*3+0), SimTK::XAxis, (sigma)(iy*3+1), SimTK::YAxis, (sigma)(iy*3+2), SimTK::ZAxis));
                                asd = asd.normalize(); 
                                auto sad = Eigen::Quaternion<double>(asd(0), asd(1), asd(2), asd(3));
                                noise[iy] = sad * noise[iy].conjugate();    // subtract the column
                            }
                        }
                        else if (ii >= ((2 * (nqf + (order*nuf)) + 1))) {
                            //log_info("trying to access SSy matrix column {}", (ii - (2 * (nqf + (order*nuf)) + 1)));
                            auto sigma = SSy.col(ii - (2 * (nqf + (order*nuf)) + 1)); //euler angles
                            for (int iy = 0; iy < ny; iy++) {
                                auto asd = SimTK::Quaternion_<double>(SimTK::Rotation(SimTK::BodyOrSpaceType::SpaceRotationSequence, 
                                (sigma)(iy*3+0), SimTK::XAxis, (sigma)(iy*3+1), SimTK::YAxis, (sigma)(iy*3+2), SimTK::ZAxis));
                                asd = asd.normalize(); 
                                auto sad = Eigen::Quaternion<double>(asd(0), asd(1), asd(2), asd(3));
                                noise[iy] = sad * noise[iy];    // add the column
                            }
                        }              
                        // Add the observation noise to the orientations
                        for (int yy = 0; yy < ny; yy++) {
                            (*(arr_dummy4vector[ithr])) = (*(arr_osensorOrientations[ithr]))[yy].convertRotationToQuaternion().asVec4();
                            Sigmas2Orientations[ii][yy] = Eigen::Quaternion<double>((*(arr_dummy4vector[ithr]))(0), (*(arr_dummy4vector[ithr]))(1), (*(arr_dummy4vector[ithr]))(2), (*(arr_dummy4vector[ithr]))(3));
                            Sigmas2OrientationsNoObsError[ii][yy] = Sigmas2Orientations[ii][yy];
                            Sigmas2Orientations[ii][yy] = noise[yy] * Sigmas2Orientations[ii][yy]; 
                        }
                        //Sigmas2props.col(ii) = y;
                    }
                });
            }
            pool.waitUntilCompleted();
            //log_info("completed step 5");

            // Step 6. Calculate expected observation of current step
            num_cols = ny;
            elPerThread = num_cols / num_cores;
            for (int ithr = 0; ithr < num_cores; ++ithr) {
                int start = ithr * elPerThread;
                int end = (ithr == num_cores - 1) ? num_cols : (ithr + 1) * elPerThread;
                //threads.emplace_back([&, start, end, M, MM, eigSolver, eigenVals, eigenVecs, iMax, maxEigenVal]() mutable {
                pool.enqueue([&, start, end, ithr, iMax, maxEigenVal]() mutable {
                    for (int iy = start; iy < end; iy++) {
                        //First, process the orientations with observation errors
                        for (int icol = 0; icol < (2 * (nqf + (order*nuf) + (3 * ny)) + 1); icol++) {
                            (*(arr_M[ithr]))(0, icol) = Sigmas2Orientations[icol][iy].w();
                            (*(arr_M[ithr]))(1, icol) = Sigmas2Orientations[icol][iy].x();
                            (*(arr_M[ithr]))(2, icol) = Sigmas2Orientations[icol][iy].y();
                            (*(arr_M[ithr]))(3, icol) = Sigmas2Orientations[icol][iy].z();
                        }
                        (*(arr_MM[ithr])) = (*(arr_M[ithr])) * Wm.asDiagonal() * (arr_M[ithr])->transpose();
                        arr_eigSolver[ithr]->compute((*(arr_MM[ithr])), true);
                        (*(arr_eigenVals[ithr])) = arr_eigSolver[ithr]->eigenvalues().real();
                        (*(arr_eigenVecs[ithr])) = arr_eigSolver[ithr]->eigenvectors().real();
                        iMax = 0;
                        maxEigenVal = (*(arr_eigenVals[ithr]))(iMax);
                        for (int ii = 0; ii < 4; ii++) {
                            if ((*(arr_eigenVals[ithr]))(ii) > maxEigenVal) {
                                iMax = ii;
                                maxEigenVal = (*(arr_eigenVals[ithr]))(iMax);
                            }
                        }
                        expectedOVector[iy] = Eigen::Quaternion<double>((*(arr_eigenVecs[ithr]))(0, iMax), 
                        (*(arr_eigenVecs[ithr]))(1, iMax), (*(arr_eigenVecs[ithr]))(2, iMax), 
                        (*(arr_eigenVecs[ithr]))(3, iMax));

                        //Second, process the orientations without observation errors
                        for (int icol = 0; icol < (2 * (nqf + (order*nuf) + (3 * ny)) + 1); icol++) {
                            (*(arr_M[ithr]))(0, icol) = Sigmas2OrientationsNoObsError[icol][iy].w();
                            (*(arr_M[ithr]))(1, icol) = Sigmas2OrientationsNoObsError[icol][iy].x();
                            (*(arr_M[ithr]))(2, icol) = Sigmas2OrientationsNoObsError[icol][iy].y();
                            (*(arr_M[ithr]))(3, icol) = Sigmas2OrientationsNoObsError[icol][iy].z();
                        }
                        (*(arr_MM[ithr])) = (*(arr_M[ithr])) * Wm.asDiagonal() * (arr_M[ithr])->transpose();
                        arr_eigSolver[ithr]->compute((*(arr_MM[ithr])), true);
                        (*(arr_eigenVals[ithr])) = arr_eigSolver[ithr]->eigenvalues().real();
                        (*(arr_eigenVecs[ithr])) = arr_eigSolver[ithr]->eigenvectors().real();
                        iMax = 0;
                        maxEigenVal = (*(arr_eigenVals[ithr]))(iMax);
                        for (int ii = 0; ii < 4; ii++) {
                            if ((*(arr_eigenVals[ithr]))(ii) > maxEigenVal) {
                                iMax = ii;
                                maxEigenVal = (*(arr_eigenVals[ithr]))(iMax);
                            }
                        }
                        expectedOVectorNoObsError[iy] = Eigen::Quaternion<double>((*(arr_eigenVecs[ithr]))(0, iMax), 
                        (*(arr_eigenVecs[ithr]))(1, iMax), (*(arr_eigenVecs[ithr]))(2, iMax), 
                        (*(arr_eigenVecs[ithr]))(3, iMax));
                    }
                });
            }
            pool.waitUntilCompleted();
            //log_info("completed step 6");

            //Step 7. Subtract expected mean orientation from propagated sigma points
            num_cols = (2 * (nqf + (order*nuf) + (3 * ny)) + 1);
            elPerThread = num_cols / num_cores;
            for (int ithr = 0; ithr < num_cores; ++ithr) {
                int start = ithr * elPerThread;
                int end = (ithr == num_cores - 1) ? num_cols : (ithr + 1) * elPerThread;
                //threads.emplace_back([&, start, end, simTKquat, angle_vector, y]() mutable {
                pool.enqueue([&, start, end, ithr]() mutable {
                    for (int icol = start; icol < end; icol++) {
                        for (int iy = 0; iy < ny; iy++) {
                            //Sigmas2OminusMean[icol][iy] = (vQuat[iy] *  expectedOVector[iy]) * Sigmas2Orientations[icol][iy].conjugate();
                            Sigmas2OminusMean[icol][iy] = expectedOVector[iy] * Sigmas2Orientations[icol][iy].conjugate();
                            //Sigmas2OminusMean[icol][iy] = Sigmas2Orientations[icol][iy].conjugate() * expectedOVector[iy];
                            (*(arr_simTKquat[ithr])) = SimTK::Quaternion_<double>(
                                    Sigmas2OminusMean[icol][iy].w(),
                                    Sigmas2OminusMean[icol][iy].x(),
                                    Sigmas2OminusMean[icol][iy].y(),
                                    Sigmas2OminusMean[icol][iy].z());
                            (*(arr_angle_vector[ithr])) = SimTK::Rotation_<double>((*(arr_simTKquat[ithr]))).convertThreeAxesRotationToThreeAngles(
                                    SimTK::BodyOrSpaceType::SpaceRotationSequence, SimTK::XAxis, SimTK::YAxis, SimTK::ZAxis);
                            for (int kk = 0; kk < 3; kk++) {
                                (*(arr_y[ithr]))(iy * 3 + kk) = (double) (*(arr_angle_vector[ithr]))(kk);
                            }
                        }
                        Sigmas2props.col(icol) = (*(arr_y[ithr]));
                    }
                });
            }
            pool.waitUntilCompleted();
            //log_info("completed step 7");

            // Step 8. Subtract expected mean orientation from data
            // NOTE: We need to update the observation noise; thus two different innovation terms ydiff and ydiff0

            if ((missingDataScale - 1.0) > std::pow(10.0, -4.0) && get_observationForgetFactor() > std::pow(10.0, -4.0)) {
                R = R0; //restore covariance of observations; gets modified if missing data
            }            

            oRefs.getValuesAtTime(time, yArray);
            for (std::map<int, int>::iterator it = oMapFromDataToModel.begin();
                    it != oMapFromDataToModel.end(); ++it) {
                dummy4vector = yArray[it->first].convertRotationToQuaternion().asVec4();
                if (SimTK::isNaN(dummy4vector(0))) {
                    dataOVector[it->second] = Eigen::Quaternion<double>(0, 0, 0, 0);
                } 
                else {
                    dataOVector[it->second] = Eigen::Quaternion<double>(
                            dummy4vector(0), dummy4vector(1), dummy4vector(2),
                            dummy4vector(3));
                }
            }
            for (int iy = 0; iy < ny; iy++) {
                if (dataOVector[iy].w() == 0 && dataOVector[iy].x() == 0 &&
                        dataOVector[iy].y() == 0 && dataOVector[iy].z() == 0) {
                    for (int kk = 0; kk < 3; kk++) {
                        ydiff(iy * 3 + kk) = 0.0;
                        ydiff0(iy * 3 + kk) = 0.0;
                        R(iy * 3 + kk, iy * 3 + kk) *= missingDataScale;
                    }
                } 
                else if (std::isnan(dataOVector[iy].w()) &&
                           std::isnan(dataOVector[iy].x()) &&
                           std::isnan(dataOVector[iy].y()) &&
                           std::isnan(dataOVector[iy].z())) {
                    for (int kk = 0; kk < 3; kk++) {
                        ydiff(iy * 3 + kk) = 0.0;
                        ydiff0(iy * 3 + kk) = 0.0;
                        R(iy * 3 + kk, iy * 3 + kk) *= missingDataScale;
                    }                    
                }
                else {
                    // dataMinusMeanQuat0 = (vQuat[iy] *
                    //         expectedOVector[iy]) * dataOVector[iy].conjugate(); // THIS AND NEXT TERM WERE SWAPPED (2025-06-02)
                    dataMinusMeanQuat0 = (vQuat[iy] *
                            expectedOVector[iy].conjugate()) * dataOVector[iy].conjugate();
                    //dataMinusMeanQuat0 = expectedOVectorNoObsError[iy] * dataOVector[iy].conjugate();
                    // dataMinusMeanQuat = (expectedOVector[iy] * vQuat[iy]) * dataOVector[iy].conjugate();
                    dataMinusMeanQuat = expectedOVector[iy] * dataOVector[iy].conjugate();
                    // dataMinusMeanQuat = dataOVector[iy].conjugate() *
                    // expectedOVector[iy];
                    simTKquat = SimTK::Quaternion_<double>(
                            dataMinusMeanQuat.w(), dataMinusMeanQuat.x(),
                            dataMinusMeanQuat.y(), dataMinusMeanQuat.z());
                    simTKquat0 = SimTK::Quaternion_<double>(
                            dataMinusMeanQuat0.w(), dataMinusMeanQuat0.x(),
                            dataMinusMeanQuat0.y(), dataMinusMeanQuat0.z());
                    angle_vector = SimTK::Rotation_<double>(simTKquat).convertThreeAxesRotationToThreeAngles(
                        SimTK::BodyOrSpaceType::SpaceRotationSequence,
                        SimTK::XAxis, SimTK::YAxis, SimTK::ZAxis);
                    angle_vector0 = SimTK::Rotation_<double>(simTKquat0).convertThreeAxesRotationToThreeAngles(
                        SimTK::BodyOrSpaceType::SpaceRotationSequence,
                        SimTK::XAxis, SimTK::YAxis, SimTK::ZAxis);
                    for (int kk = 0; kk < 3; kk++) {
                        ydiff(iy * 3 + kk) = (double)angle_vector(kk);
                        ydiff0(iy * 3 + kk) = (double)angle_vector0(kk);
                    }
                }
            }
            //log_info("completed step 8");

            // Step 9. Calculate covariance of expected orientations
            // Note: no addition of R here since this is augmented UKF

            Py = W0c * (Sigmas2props.col(0) * Sigmas2props.col(0).transpose());
            for (int ii = 1; ii < (2 * (nqf + (order*nuf) + (3 * ny)) + 1); ii++) {
                Py += Wi * (Sigmas2props.col(ii) * Sigmas2props.col(ii).transpose());
            }
            //log_info("completed step 9");
			
            // Step 10. Calculate cross-covariance of expected orientations and a
            // priori state mean
            for (int jj = 0; jj < (2 * (nqf + (order*nuf) + (3 * ny)) + 1); jj++) {
                Sigmas2.col(jj) -= xsave;
            }
            Pxy = W0c * (Sigmas2.col(0) * Sigmas2props.col(0).transpose());
            for (int ii = 1; ii < (2 * (nqf + (order*nuf) + (3 * ny)) + 1); ii++) {
                Pxy += Wi * (Sigmas2.col(ii) * Sigmas2props.col(ii).transpose());
            }
            //log_info("completed step 10");
			
            // Step 11. Calculate Kalman gain
            //K = Pxy * Py.inverse();
            K = Pxy * Py.completeOrthogonalDecomposition().pseudoInverse();            
			
            // Step 12. Calculate a posteriori state mean of current step
            x = xsave + K * ydiff;
            log_info("absolute difference in data, maxCoeff = {}", (double)ydiff.cwiseAbs().maxCoeff());
            //log_info("difference in data, minCoeff = {}", (double)ydiff.minCoeff());
            log_info("absolute difference in data, mean = {}", (double)ydiff.cwiseAbs().mean());
            // Step 12. Calculate a posteriori state covariance of current step
            P -= K * Py * K.transpose();
            //log_info("calculated P");

            // Apply inequality constraints (clamped coordinates)
            if (get_enable_clamping()) {
                OpenSim::AugAUKSMIKT::clampCoordinates(x, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf); 
            }


            // Step 13. Update noise means and covariances
            new_v = ydiff0;              
            new_w = x - x0; 
            z.setZero();
            H.setZero(); 

            //update noise covariances
            newQ.setZero();
            tempQ = (x - x0 - w) * (x - x0 - w).transpose(); // NOTE: we must use the old mean w here
            newR.setZero();
            
            //NOTE: we must use the old observation noise mean v here

            tempR = ydiff * ydiff.transpose();
            
            if (get_enforce_independent_sensors() == true) {
                for (int irow = 0; irow < ny; irow++) {
                    newR.block(irow*3, irow*3, 3, 3) = tempR.block(irow*3, irow*3, 3, 3);
                }
            }
            else {
                newR = tempR;
            }

            R = (1.0 - get_observationForgetFactor()) * R + get_observationForgetFactor() * (newR);
            
            //log_info("calculated newR");

            if (get_enforce_white_process_noise() == true) {
                for (int idx = 0; idx < (order+1); idx++) {
                    H(idx) = QCoeffs(idx, idx);
                }
                for (int cidx = 0; cidx < nqf; cidx++) {
                    for (int oidx = 0; oidx < (order+1); oidx++) {
                        z(oidx) = tempQ(cidx + oidx*nuf, cidx + oidx*nuf);
                    }
                    //sigmas(cidx, 0) = std::max(std::abs((H.transpose() * H).ldlt().solve(H.transpose() * z)(0,0)), get_sgma2w_min());
                    //sigmas(cidx, 0) = std::max((z(order, 0) / H(order, 0)), get_sgma2w_min());
                    //newSigmas(cidx, 0) = (1.0 - get_processForgetFactor()) * oldSigmas(cidx, 0) + get_processForgetFactor() * sigmas(cidx, 0);
                    sigmas(cidx, 0) = (z(order, 0) / H(order, 0));
                    newSigmas(cidx, 0) = std::min(std::max(((1.0 - get_processForgetFactor()) * oldSigmas(cidx, 0) + 
                    get_processForgetFactor() * sigmas(cidx, 0)), get_sgma2w_min()), get_sgma2w_max());
                }
                //log_info("calculated sigmas");

                //newSigmas = (1.0 - get_processForgetFactor()) * oldSigmas + get_processForgetFactor() * sigmas;

                log_info("new process noise variances, minCoeff = {}", (double)newSigmas.cwiseAbs().minCoeff());
                log_info("new process noise variances, maxCoeff = {}", (double)newSigmas.cwiseAbs().maxCoeff());
                log_info("new process noise variances, mean = {}", (double)newSigmas.cwiseAbs().mean());

                for (int cidx = 0; cidx < nqf; cidx++) {
                    for (int colidx = 0; colidx < (order+1); colidx++) {
                        for (int rowidx = 0; rowidx < (order+1); rowidx++) {
                            newQ(((rowidx*nuf) + cidx), ((colidx*nuf) + cidx)) = QCoeffs(rowidx, colidx) * newSigmas(cidx, 0);
                        } 
                    }                    
                }
                oldSigmas = newSigmas;
                Q = newQ;
            }
            else {
                Q = (1.0 - get_processForgetFactor()) * Q + get_processForgetFactor() * (tempQ);
            }

            if (get_observation_noise_zero() == false) {
                for (int iy = 0; iy < ny; iy++) {           
                    auto asd = SimTK::Quaternion_<double>(SimTK::Rotation(SimTK::BodyOrSpaceType::SpaceRotationSequence, 
                    (new_v)(iy*3+0), SimTK::XAxis, (new_v)(iy*3+1), SimTK::YAxis, (new_v)(iy*3+2), SimTK::ZAxis));    
                    asd = asd.normalize(); 
                    new_vQuat[iy] = Eigen::Quaternion<double>(asd(0), asd(1), asd(2), asd(3));
                    //new_vQuat.emplace_back(Eigen::Quaternion<double>(asd(0), asd(1), asd(2), asd(3)));
                }
                Eigen::MatrixXd MMM(4,2);
                for (int iy = 0; iy < ny; iy++) {
                    MMM(0,0) = vQuat[iy].w();
                    MMM(1,0) = vQuat[iy].x();
                    MMM(2,0) = vQuat[iy].y();
                    MMM(3,0) = vQuat[iy].z();
                    MMM(0,1) = new_vQuat[iy].w();
                    MMM(1,1) = new_vQuat[iy].x();
                    MMM(2,1) = new_vQuat[iy].y();
                    MMM(3,1) = new_vQuat[iy].z();                    
                    eigSolver.compute(MMM * vUpdateweightMat * MMM.transpose(), true);
                    auto eigVals = eigSolver.eigenvalues().real();
                    auto eigVecs = eigSolver.eigenvectors().real();
                    iMax = 0;
                    maxEigenVal = eigVals(iMax);
                    for (int ii = 0; ii < 4; ii++) {
                        if (eigVals(ii) > maxEigenVal) {
                            iMax = ii;
                            maxEigenVal = eigVals(iMax);
                        }
                    }
                    // the v quaternion form gets updated here
                    vQuat[iy] = Eigen::Quaternion<double>(eigVecs(0, iMax), 
                    eigVecs(1, iMax), eigVecs(2, iMax), eigVecs(3, iMax));
                    vQuat[iy].normalize(); //just ensure it's of unit length...
                }
                for (int ithr = 0; ithr < num_cores; ithr++) {
                    for (int iy = 0; iy < ny; iy++) {
                        arr_vQuat[ithr][iy] = vQuat[iy];
                    }                    
                }
            }   

            if (get_process_noise_zero() == false) {
                w = (1.0 - get_processForgetFactor()) * w + get_processForgetFactor() * (new_w);
                log_info("new process noise mean, minCoeff = {}", (double)(w).minCoeff());
                log_info("new process noise mean, maxCoeff = {}", (double)(w).maxCoeff());
                log_info("new process noise mean, mean = {}", (double)(w).mean());
            }       

            priorStatsVector.emplace_back(x);
            priorStatsVector.emplace_back(P);  
            {
                std::unique_lock<std::mutex> lock(*fwdBwdMutex);
                priorStatsBuffer->push(priorStatsVector);
            }
            condVarB->notify_one();

        }

        // Abort running inverse kinematics in case of significant numerical instabilities
        if (get_abort_if_diverging() == true) {
            if ((double)ydiff0.cwiseAbs().mean() > 1.0) {
                log_info("Significant numerical instabilities encountered, aborting...");
                break;
            }
            else if ((double)ydiff0.cwiseAbs().maxCoeff() > 3.0) {
                log_info("Significant numerical instabilities encountered, aborting...");
                break;
            }
        }
    }	//end of for

    {
        std::unique_lock<std::mutex> lock1(*fwdBwdMutex);
        *(fwdDone) = true;
    }
    condVarB->notify_all();  

    // Delete dynamically allocated stuff
    OpenSim::AugAUKSMIKT::deletePointers(solvers);
    OpenSim::AugAUKSMIKT::deletePointers(models);
    OpenSim::AugAUKSMIKT::deletePointers(arr_yMapFromEigenToSimbody);
    OpenSim::AugAUKSMIKT::deletePointers(arr_angle_vector);
    OpenSim::AugAUKSMIKT::deletePointers(arr_dummy4vector);
    OpenSim::AugAUKSMIKT::deletePointers(arr_eigenVals);
    OpenSim::AugAUKSMIKT::deletePointers(arr_eigenVecs);
    OpenSim::AugAUKSMIKT::deletePointers(arr_eigSolver);
    OpenSim::AugAUKSMIKT::deletePointers(arr_M);
    OpenSim::AugAUKSMIKT::deletePointers(arr_MM);
    OpenSim::AugAUKSMIKT::deletePointers(arr_osensorOrientations);
    OpenSim::AugAUKSMIKT::deletePointers(arr_q);
    OpenSim::AugAUKSMIKT::deletePointers(arr_simTKquat);
    OpenSim::AugAUKSMIKT::deletePointers(arr_ss);
    OpenSim::AugAUKSMIKT::deletePointers(arr_u);
    OpenSim::AugAUKSMIKT::deletePointers(arr_xx);
    OpenSim::AugAUKSMIKT::deletePointers(arr_y);

	log_info("Got out of for loop and end of UKFTool.");
	

}  //end of AugAUKSMIKT::UKFTool


std::tuple<std::map<std::string, int>, std::map<int, std::string>> OpenSim::AugAUKSMIKT::CreateYMaps(OpenSim::Model model) {
	model.initSystem();
    SimTK::State ss = model.getWorkingState();
    OpenSim::Array<std::string> modelStateVariableNames = model.getStateVariableNames();
    int numY = model.getNumStateVariables();
    ss.updY() = 0;
    std::map<std::string, int> yMapFromOpenSimToSimbody;
    std::map<int, std::string> yMapFromSimbodyToOpenSim;
    SimTK::Vector modelStateVariableValues;
    for (int iy = 0; iy < numY; iy++) { //this y-index runs for Simbody
        ss.updY()[iy] = SimTK::NaN;
        modelStateVariableValues = model.getStateVariableValues(ss);
        for (int ii = 0; ii < modelStateVariableNames.size(); ii++) {   //this index runs for OpenSim
            if (SimTK::isNaN(modelStateVariableValues[ii])) {
                yMapFromOpenSimToSimbody.insert(std::pair<std::string, int>(modelStateVariableNames[ii], iy));
                yMapFromSimbodyToOpenSim.insert(std::pair<int, std::string>(iy, modelStateVariableNames[ii]));
                ss.updY()[iy] = 0;
                break;
            }
        }
        if (SimTK::isNaN(ss.updY()[iy])) {
            // If we reach here, this is an unused slot for a quaternion (from Antoine Felisse code)
            ss.updY()[iy] = 0;
        }
    }
    std::tuple<std::map<std::string, int>, std::map<int, std::string>> mappings(yMapFromOpenSimToSimbody, yMapFromSimbodyToOpenSim);
    if (numY != (int)yMapFromOpenSimToSimbody.size()) {
        log_info("There were {} state variables, but got {} mappings from OpenSim to Simbody!", numY,
            (int)yMapFromOpenSimToSimbody.size());
    }
    return mappings;

}   //end of AugAUKSMIKT::CreateYMaps


double OpenSim::AugAUKSMIKT::computeFactorial(int input) {
    double fact = 1.0;
    for (int ii = 1; ii <= input; ii++) {
        fact *= ii;
    }
    return fact;
}

double OpenSim::AugAUKSMIKT::probWithinInterval(double x1, double x2) {
    return (std::erf(x2 / std::sqrt(2)) - std::erf(x1 / std::sqrt(2))) / 2;
}

template <typename T>
void OpenSim::AugAUKSMIKT::deletePointers(std::vector<T*>& vec) {
    for (T* ptr : vec) {
        delete ptr;
    }
    vec.clear();
}

//template <class T>
void OpenSim::AugAUKSMIKT::computeBackwardPass(OpenSim::Model& model, std::vector<OpenSim::UKFClampedCoordLimits> clampedCoordLimits, 
    std::queue<std::vector<Eigen::MatrixXd>>* priorStatsBuffer, std::mutex* fwdBwdMutex, std::condition_variable* condVar, bool* fwdDone, 
    std::map<int, int> yMapFromEigenToSimbody, std::map<int, int> yMapFromSimbodyToEigen, std::map<std::string, int> yMapFromOpenSimToSimbody, 
    OpenSim::AnalysisSet& analysisSet, std::map<int, std::string> yMapFromSimbodyToOpenSim, int nqf, int nuf) {

    SimTK::State ss;
    //log_info("bwd: managed to get in bwdThread");

    {
        std::unique_lock<std::mutex> lock(*fwdBwdMutex);
        ss = SimTK::State(model.getWorkingState());
    }
    //log_info("bwd: managed to get state");
    
    SimTK::Vector q = ss.getQ();
    SimTK::Vector u = ss.getU();
    int step = 0;
    std::vector<Eigen::MatrixXd> bwdPriorStats;
    std::deque<Eigen::MatrixXd> currentTimes;
    std::deque<Eigen::MatrixXd> currentPriorMeans;
    std::deque<Eigen::MatrixXd> currentPriorAutoCovs;
    std::deque<Eigen::MatrixXd> previousCrossCovs;
    std::deque<Eigen::MatrixXd> currentPosteriorMeans;
    std::deque<Eigen::MatrixXd> currentPosteriorAutoCovs;
    Eigen::MatrixXd D;
    Eigen::MatrixXd smoothPosteriorMean;
    Eigen::MatrixXd smoothPosteriorCov;
    double time;
    bool writeToFile = false;
    int lagLength = get_lag_length();
    // If lag_length is negative, we ignore backward pass and output the UKF (forward) estimate
    if (lagLength < 0) {
        lagLength = -1;
    }
    int order = get_order();

    // Files to store state means and covariances
    auto UKFresultsDir = get_results_directory();
    if (UKFresultsDir.empty() && !get_output_motion_file().empty()) {
        UKFresultsDir = OpenSim::IO::getParentDirectory(get_output_motion_file());
    }
    if (!UKFresultsDir.empty()) {
        OpenSim::IO::makeDir(UKFresultsDir);
    }
    UKFresultsDir.append("/ukf/");
    OpenSim::IO::makeDir(UKFresultsDir);
    auto UKFResultsDirX = UKFresultsDir;
    auto UKFResultsDirPx = UKFresultsDir;
    UKFResultsDirX.append("x.txt");
    UKFResultsDirPx.append("Px.txt");
    std::ofstream fileX(UKFResultsDirX, std::ios_base::trunc);
    std::ofstream filePx(UKFResultsDirPx, std::ios_base::trunc);
    fileX << "endheader" << std::endl << "time\t";
    filePx << "endheader" << std::endl << "time\t";


    std::regex rgx(".*/(\\w+)/value");
    std::smatch match;
    for (int iord = 0; iord <= order; iord++) {
        for (std::map<int, int>::iterator it = yMapFromEigenToSimbody.begin(); it != yMapFromEigenToSimbody.end(); ++it) {
            std::string s = yMapFromSimbodyToOpenSim.at(it->second);
            if (std::regex_search(s,match,rgx)) {
                fileX << match[1] << "_" << std::to_string(iord) << "\t";
            }
        }
    }

    for (int irow = 0; irow < ((order+1) * nqf); irow++) {
        for (int icol = irow; icol < ((order+1) * nqf); icol++) {
            filePx << "[" << irow << ", " << icol << "]" << "\t";
        }
    }

    fileX << std::endl;
    filePx << std::endl;

    while (true) {
        //log_info("bwd: start of loop");
        std::unique_lock<std::mutex> lock(*fwdBwdMutex);
        condVar->wait(lock, [&](){ return !(priorStatsBuffer->empty()) || (*fwdDone); });

        if (!priorStatsBuffer->empty()) {
            bwdPriorStats = priorStatsBuffer->front();
            priorStatsBuffer->pop();
            lock.unlock();
            currentTimes.push_back(bwdPriorStats[0]);
            currentPriorMeans.push_back(bwdPriorStats[1]);
            currentPriorAutoCovs.push_back(bwdPriorStats[2]);
            previousCrossCovs.push_back(bwdPriorStats[3]);
            currentPosteriorMeans.push_back(bwdPriorStats[4]);
            currentPosteriorAutoCovs.push_back(bwdPriorStats[5]);
        }

        if (lagLength < 0 && (int)currentPriorMeans.size() >= 1) {
            smoothPosteriorMean = currentPosteriorMeans.at(0);
            smoothPosteriorCov = currentPosteriorAutoCovs.at(0);
            time = currentTimes.at(0)(0,0);
            currentTimes.pop_front();
            currentPriorMeans.pop_front();
            currentPriorAutoCovs.pop_front();
            previousCrossCovs.pop_front();
            currentPosteriorMeans.pop_front();
            currentPosteriorAutoCovs.pop_front();
            writeToFile = true;
        }
        else if (lagLength >= 0 && (int)currentPriorMeans.size() >= (lagLength+2)) {
            if (lagLength < 0) {
                smoothPosteriorMean = currentPosteriorMeans.at(0);
                smoothPosteriorCov = currentPosteriorAutoCovs.at(0);
            }
            else {
                for (int ii = lagLength; ii >= 0; ii--) {
                    //D = previousCrossCovs.at(ii+1) * (currentPriorAutoCovs.at(ii+1).inverse());
                    D = previousCrossCovs.at(ii+1) * (currentPriorAutoCovs.at(ii+1).completeOrthogonalDecomposition().pseudoInverse());                    
                    if (ii == lagLength) {
                        smoothPosteriorMean = currentPosteriorMeans.at(ii) + D * (currentPosteriorMeans.at(ii+1) - currentPriorMeans.at(ii+1));
                        if (get_enable_clamping()) {
                            OpenSim::AugAUKSMIKT::clampCoordinates(smoothPosteriorMean, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf);
                        }
                        smoothPosteriorCov = currentPosteriorAutoCovs.at(ii) + D * (currentPosteriorAutoCovs.at(ii+1) - currentPriorAutoCovs.at(ii+1)) * D.transpose(); 
                    }
                    else if (ii < lagLength) {
                        smoothPosteriorMean = currentPosteriorMeans.at(ii) + D * (smoothPosteriorMean - currentPriorMeans.at(ii+1));
                        if (get_enable_clamping()) {
                            OpenSim::AugAUKSMIKT::clampCoordinates(smoothPosteriorMean, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf);
                        }
                        smoothPosteriorCov = currentPosteriorAutoCovs.at(ii) + D * (smoothPosteriorCov - currentPriorAutoCovs.at(ii+1)) * D.transpose(); 
                    }                    
                }
            }
            writeToFile = true;
            time = currentTimes.at(0)(0,0);
            currentTimes.pop_front();
            currentPriorMeans.pop_front();
            currentPriorAutoCovs.pop_front();
            previousCrossCovs.pop_front();
            currentPosteriorMeans.pop_front();
            currentPosteriorAutoCovs.pop_front();
        }
        else if ((*fwdDone)) {
            if (lagLength < 0 && (int)currentPriorMeans.size() >= 1) {
                time = currentTimes.at(0)(0,0);
                smoothPosteriorMean = currentPosteriorMeans.at(0);
                smoothPosteriorCov = currentPosteriorAutoCovs.at(0);
                writeToFile = true;
            }
            else {
                if ((int)currentPriorMeans.size() > 1) {
                    for (int ii = currentPriorMeans.size()-2; ii >= 0; ii--) {
                        //D = previousCrossCovs.at(ii+1) * (currentPriorAutoCovs.at(ii+1).inverse());
                        D = previousCrossCovs.at(ii+1) * (currentPriorAutoCovs.at(ii+1).completeOrthogonalDecomposition().pseudoInverse());                        
                        if (ii == (((int)currentPriorMeans.size())-2)) {
                            smoothPosteriorMean = currentPosteriorMeans.at(ii) + D * (currentPosteriorMeans.at(ii+1) - currentPriorMeans.at(ii+1));
                            if (get_enable_clamping()) {
                                OpenSim::AugAUKSMIKT::clampCoordinates(smoothPosteriorMean, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf);
                            }
                            smoothPosteriorCov = currentPosteriorAutoCovs.at(ii) + D * (currentPosteriorAutoCovs.at(ii+1) - currentPriorAutoCovs.at(ii+1)) * D.transpose(); 
                        }
                        else if (ii < (((int)currentPriorMeans.size())-2)) {
                            smoothPosteriorMean = currentPosteriorMeans.at(ii) + D * (smoothPosteriorMean - currentPriorMeans.at(ii+1));
                            if (get_enable_clamping()) {
                                OpenSim::AugAUKSMIKT::clampCoordinates(smoothPosteriorMean, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf);
                            }
                            smoothPosteriorCov = currentPosteriorAutoCovs.at(ii) + D * (smoothPosteriorCov - currentPriorAutoCovs.at(ii+1)) * D.transpose(); 
                        }
                    }
                    time = currentTimes.at(0)(0,0);
                    writeToFile = true;
                }
                else if (((int)currentPriorMeans.size() == 1) || ((int)currentPosteriorMeans.size() == 1)) {
                    smoothPosteriorMean = currentPosteriorMeans.at(0);
                    smoothPosteriorCov = currentPosteriorAutoCovs.at(0);
                    time = currentTimes.at(0)(0,0);
                    writeToFile = true;
                }
                else if ((int)currentPriorMeans.size() == 0) {
                    fileX.close();
                    filePx.close();
                    writeToFile = false;
                    break;
                }
                currentTimes.pop_front();
                currentPriorMeans.pop_front();
                currentPriorAutoCovs.pop_front();
                previousCrossCovs.pop_front();
                currentPosteriorMeans.pop_front();
                currentPosteriorAutoCovs.pop_front();
            }
        }
        
        // Send results to the reporter
        if (writeToFile) {
            // Constrain the smoothed mean
            if (get_enable_clamping()) {
                OpenSim::AugAUKSMIKT::clampCoordinates(smoothPosteriorMean, clampedCoordLimits, yMapFromOpenSimToSimbody, yMapFromSimbodyToEigen, order, nuf);
            }
            for (std::map<int, int>::iterator it = yMapFromEigenToSimbody.begin(); it != yMapFromEigenToSimbody.end(); ++it) {
                q(it->second) = smoothPosteriorMean(it->first);
            }
            for (std::map<int, int>::iterator it =
                            yMapFromEigenToSimbody.begin();
                    it != yMapFromEigenToSimbody.end(); ++it) {
                u(it->second) = smoothPosteriorMean((it->first)+nqf);
            }
            ss.updTime() = time;
            ss.updQ() = q;
            ss.updU() = u;
            log_info("Solved at time: {} s", time);
            // realize to report to get reporter to pull values from model
            analysisSet.step(ss, step++);
            model.realizeReport(ss);
        }        
            
        // Write the state to result files
        if (get_write_UKF() && writeToFile) {
            filePx << static_cast<float>(time) << "\t";
            for (int irow = 0; irow < ((order+1) * nqf); irow++) {
                for (int icol = irow; icol < ((order+1) * nqf); icol++) {
                    filePx << static_cast<float>(smoothPosteriorCov(irow, icol)) << "\t";
                }
            }
            fileX << static_cast<float>(time) << "\t";
            for (int irow = 0; irow < ((order+1) * nqf); irow++) {
                fileX << static_cast<float>(smoothPosteriorMean(irow)) << "\t";
            }
            fileX << std::endl;
            filePx << std::endl;
        }
    }   // end of while(true)
}   //end of AugAUKSMIKT::computeBackwardPass


void OpenSim::AugAUKSMIKT::clampCoordinates(Eigen::MatrixXd& stateMeans, std::vector<OpenSim::UKFClampedCoordLimits> clampedCoordLimits, 
std::map<std::string, int> yMapFromOpenSimToSimbody, std::map<int, int> yMapFromSimbodyToEigen, int order, int nuf) {
    for (auto coord : clampedCoordLimits) {
        if (stateMeans(yMapFromSimbodyToEigen[yMapFromOpenSimToSimbody[coord.stateVarName]]) > coord.rangeMax) {
            stateMeans(yMapFromSimbodyToEigen[yMapFromOpenSimToSimbody[coord.stateVarName]]) = coord.rangeMax;
            for (int ordidx = 1; ordidx <= order; ordidx++) {
                if (stateMeans((ordidx * nuf) + yMapFromSimbodyToEigen[yMapFromOpenSimToSimbody[coord.stateVarName]]) > 0.0) {
                    stateMeans((ordidx * nuf) + yMapFromSimbodyToEigen[yMapFromOpenSimToSimbody[coord.stateVarName]]) = 0.0;
                }
            }
        }
        else if (stateMeans(yMapFromSimbodyToEigen[yMapFromOpenSimToSimbody[coord.stateVarName]]) < coord.rangeMin) {
            stateMeans(yMapFromSimbodyToEigen[yMapFromOpenSimToSimbody[coord.stateVarName]]) = coord.rangeMin;
            for (int ordidx = 1; ordidx <= order; ordidx++) {
                if (stateMeans((ordidx * nuf) + yMapFromSimbodyToEigen[yMapFromOpenSimToSimbody[coord.stateVarName]]) < 0.0) {
                    stateMeans((ordidx * nuf) + yMapFromSimbodyToEigen[yMapFromOpenSimToSimbody[coord.stateVarName]]) = 0.0;
                }
            }
        }
    }
}
