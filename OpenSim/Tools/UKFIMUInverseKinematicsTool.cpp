
#include "UKFIMUInverseKinematicsTool.h"

//using namespace OpenSim;
//using namespace SimTK;
//using namespace std;


// UKFThreadPool methods


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




// UKFIMUInverseKinematicsTool methods
OpenSim::UKFIMUInverseKinematicsTool::UKFIMUInverseKinematicsTool()
        : OpenSim::InverseKinematicsToolBase() {
    OpenSim::UKFIMUInverseKinematicsTool::constructProperties();
}

OpenSim::UKFIMUInverseKinematicsTool::UKFIMUInverseKinematicsTool(const std::string& setupFile)
        : OpenSim::InverseKinematicsToolBase(setupFile, true) {
    OpenSim::UKFIMUInverseKinematicsTool::constructProperties();
    updateFromXMLDocument();
}

OpenSim::UKFIMUInverseKinematicsTool::~UKFIMUInverseKinematicsTool()
{
}

void OpenSim::UKFIMUInverseKinematicsTool::constructProperties()
{
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_sensor_to_opensim_rotations(
            SimTK::Vec3(0));
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_orientations_file("");
    OpenSim::OrientationWeightSet orientationWeights;
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_orientation_weights(
            orientationWeights);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_alpha(1.0);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_beta(2.0);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_kappa(-1.337);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_sgma2w(1000.0);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_order(2);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_lag_length(3);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_missing_data_scale(100.0);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_num_threads(3);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_write_UKF(true);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_enable_resampling(true);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_process_covariance_method(0);
    OpenSim::UKFIMUInverseKinematicsTool::constructProperty_imu_RMS_in_deg(SimTK::Vec3(0));
        
}

void OpenSim::UKFIMUInverseKinematicsTool::runInverseKinematicsWithOrientationsFromFile(
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
    const SimTK::Vec3& rotations = OpenSim::UKFIMUInverseKinematicsTool::get_sensor_to_opensim_rotations();
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

    std::mutex* fwdBwdMutex = new std::mutex();
    std::condition_variable* condVar = new std::condition_variable();
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
    std::tuple<std::map<std::string, int>, std::map<int, std::string>> mappings = OpenSim::UKFIMUInverseKinematicsTool::CreateYMaps(model);
    std::map<std::string, int> yMapFromOpenSimToSimbody = std::get<0>(mappings);
    std::map<int, std::string> yMapFromSimbodyToOpenSim = std::get<1>(mappings);

    // Find locked coordinates (we cannot change these in UKF)

    std::map<int, int> simbodyIndexTypeMap;
    std::map<int, int> yMapFromSimbodyToEigen;
    std::map<int, int> yMapFromEigenToSimbody;

    int iqx;
    int nqf = 0;    //number of free and clamped q's
    int nuf = 0;
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

    
    std::thread forwardThread([&] {OpenSim::UKFIMUInverseKinematicsTool::UKFTool(
        model, nqf, nuf, yMapFromSimbodyToEigen, yMapFromEigenToSimbody, yMapFromSimbodyToOpenSim, oMapFromDataToModel, 
        priorStatsBuffer, fwdBwdMutex, condVar, fwdDone, s0, analysisSet, oRefs, ikSolver, 
        modelOrientationErrors, visualizeResults, orientationErrors, processCovScales);});
    
    std::thread backwardThread([&] {OpenSim::UKFIMUInverseKinematicsTool::computeBackwardPass(
        priorStatsBuffer,fwdBwdMutex, condVar, fwdDone, yMapFromEigenToSimbody, yMapFromSimbodyToOpenSim, nqf, nuf);});
    

    forwardThread.join();
    backwardThread.join();

    delete fwdBwdMutex;
    delete condVar;
    delete fwdDone;
    delete priorStatsBuffer;

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
        log_info("UKFIMUInverseKinematicsTool: No output files were generated, "
            "set output_motion_file to generate output files.");
    // Results written to file, clear in case we run again
    ikReporter->clearTable();
}


// main driver
bool OpenSim::UKFIMUInverseKinematicsTool::run(bool visualizeResults, SimTK::Vector_<double> processCovScales)
{
    if (_model.empty()) {
        _model.reset(new Model(get_model_file()));
    }

    OpenSim::UKFIMUInverseKinematicsTool::runInverseKinematicsWithOrientationsFromFile(*_model,
            get_orientations_file(), visualizeResults, processCovScales);

    return true;
}

OpenSim::TimeSeriesTable_<SimTK::Vec3> OpenSim::UKFIMUInverseKinematicsTool::loadMarkersFile(const std::string& markerFile)
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
void OpenSim::UKFIMUInverseKinematicsTool::UKFTool(OpenSim::Model& model, int nqf, int nuf, std::map<int, int> yMapFromSimbodyToEigen, 
        std::map<int, int> yMapFromEigenToSimbody, std::map<int, std::string> yMapFromSimbodyToOpenSim, std::map<int, int> oMapFromDataToModel,
        std::queue<std::vector<Eigen::MatrixXd>>* priorStatsBuffer, 
        std::mutex* fwdBwdMutex, std::condition_variable* condVar, bool* fwdDone, SimTK::State& s, OpenSim::AnalysisSet& analysisSet,
        OpenSim::OrientationsReference oRefs, OpenSim::InverseKinematicsSolver& ikSolver,
        std::shared_ptr<OpenSim::TimeSeriesTable> modelOrientationErrors, bool visualizeResults,
        SimTK::Array_<double> orientationErrors, SimTK::Vector_<double> processCovScales) {

	log_info("Got inside UKFTOOL");


    // If the process covariance scale factors are not provided, use simply ones
    if (processCovScales.size() == 0) {
        processCovScales = SimTK::Vector_<double>(get_order(), 1.0);
        log_info("Process covariance scales not provided, using ones instead.");
    }

    auto times = oRefs.getTimes();
    double alpha = get_alpha();
    double beta = get_beta();
    double kappa = get_kappa();
    int order = get_order();
    double sgma2w = get_sgma2w();
    double missingDataScale = get_missing_data_scale();

    SimTK::Vec3 imuRMSinDeg = get_imu_RMS_in_deg();
    int num_cores = get_num_threads();
    
   
    SimTK::Array_<SimTK::Rotation_<double>>
            osensorOrientations; // array for orientations computed by the
                                 // model; REMOVE _<double> IF CAUSES ERROR
    std::vector<SimTK::Array_<SimTK::Rotation_<double>>*> arr_osensorOrientations;
    SimTK::Array_<SimTK::Rotation_<double>>
            yArray; // array for observations (IMU orientations)
    //SimTK::SimbodyMatterSubsystem matterSubSys = model.getMatterSubsystem();    
    int const ny = ikSolver.getNumOrientationSensorsInUse();    //number of sensors
    //int const ny = oRefs.getNumRefs(); // number of osensors, maybe same as above?


    int const nq = s.getNQ(); // number of generalized positions (joint angles)
    int const nu = s.getNU(); // number of generalized velocities (joint angular
                        // velocities)
    int const nr = std::min(nq, nu); // probably not needed, usually nq >= nu


    // RMS errors for observation noise covariance
    
    double xAxisRMS = imuRMSinDeg(0) * (SimTK::Pi / 180); // roll RMS error (x-axis)
    double yAxisRMS = imuRMSinDeg(1) * (SimTK::Pi / 180); // heading RMS error (y-axis)
    double zAxisRMS = imuRMSinDeg(2) * (SimTK::Pi / 180); // pitch RMS error (z-axis)

    // Coefficients for process model f()
    double deltaTime = 1.0 / oRefs.getSamplingFrequency();
    Eigen::MatrixXd fCoeffs(order+1, order+1);
    fCoeffs.setZero();
    for (int irow = 0; irow <= order; irow++) {
        for (int icol = irow; icol <= order; icol++) {
            double denum = OpenSim::UKFIMUInverseKinematicsTool::computeFactorial(icol-irow);
            fCoeffs(irow, icol) = std::pow(deltaTime, (icol-irow)) / denum;
        }
    }

    // Coefficients for process noise covariance matrix Q
    Eigen::MatrixXd QCoeffs(order+1, order+1);

    
    if (get_process_covariance_method() == 0) {
    // classic approach of Fioretti and Jetto, 1989 (no scaling tricks by default; should use ones)    
        log_info("Using method of Fioretti and Jetto, 1989.");
        for (int irow = 0; irow <= order; irow++) {
            for (int icol = 0; icol <= order; icol++) {
                double denum1 = OpenSim::UKFIMUInverseKinematicsTool::computeFactorial(order-irow);
                double denum2 = OpenSim::UKFIMUInverseKinematicsTool::computeFactorial(order-icol);
                int deltaPower = (order-irow) + (order-icol);
                QCoeffs(irow, icol) = processCovScales(irow) * processCovScales(icol) * std::pow(deltaTime, (deltaPower+1)) / (denum1 * denum2 * (deltaPower+1));
            }
        }
    }

    else if (get_process_covariance_method() == 1) {
    // approach using the Taylor remainders as white noise multipliers (plus additional scale factors)        
        log_info("Using Taylor remainders.");
        for (int irow = 0; irow <= order; irow++) {
            for (int icol = 0; icol <= order; icol++) {
                double denum1 = OpenSim::UKFIMUInverseKinematicsTool::computeFactorial(order+1-irow);
                double denum2 = OpenSim::UKFIMUInverseKinematicsTool::computeFactorial(order+1-icol);
                int deltaPower = (order+1-irow) + (order+1-icol);
                QCoeffs(irow, icol) = processCovScales(irow) * processCovScales(icol) *
                                    std::pow(deltaTime, (deltaPower)) / (denum1 * denum2);
            }
        }
    }

    // Construct covariance matrix Q for process noise
    Eigen::MatrixXd Q((nqf + order*nuf), (nqf + order*nuf));
    Q.setZero();
    for (int irow = 0; irow <= order; irow++) {
        for (int icol = 0; icol <= order; icol++) {
            Q.block((irow*nuf), (icol*nuf), nuf, nuf) = sgma2w * QCoeffs(irow, icol) * Eigen::MatrixXd::Identity(nuf, nuf);
        }
    }
        
    // Construct covariance matrix for observation errors
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
    Eigen::MatrixXd R0(3 * ny, 3 * ny);
    R0 = R;

    // Construct covariance matrix of state
    Eigen::MatrixXd P((nqf + (order*nuf)), (nqf + (order*nuf)));
    P = Q;
    
    // Construct the state vector
    Eigen::MatrixXd x((nqf + (order*nuf)), 1);
    Eigen::MatrixXd xx((nqf + (order*nuf)), 1);
    std::vector<Eigen::MatrixXd*> arr_xx;
    //SimTK::Vector_<double> x(nq + nu);
    x.setZero();
    xx.setZero();
    //Eigen::MatrixXd qf(nqf, 1);
    //qf.setZero();
    SimTK::Vector_<double> q(nq);
    std::vector<SimTK::Vector_<double>*> arr_q;
    q = s.getQ();
    SimTK::Vector_<double> u(nu);
    std::vector<SimTK::Vector_<double>*> arr_u;
    u = s.getU();
    for (std::map<int, int>::iterator it = yMapFromSimbodyToEigen.begin(); it != yMapFromSimbodyToEigen.end(); ++it) {
        x(it->second) = q(it->first);
    }
    for (std::map<int, int>::iterator it = yMapFromSimbodyToEigen.begin();
            it != yMapFromSimbodyToEigen.end(); ++it) {
        x((it->second)+nqf) = u(it->first);
    }

    // Construct the process model matrix
    Eigen::MatrixXd F((nqf + order*nuf), (nqf + order*nuf));
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

    // If the user wants kappa to be equal to state vector length (or 3-length)
    if (std::abs(kappa+1.337) < 1e-5 ) {
        kappa = -1.0 * ((double) nqf + (order*nuf) - 3);
    }
    else if (std::abs(kappa+4.337) < 1e-5) {        
        kappa = (double) nqf + (order*nuf);
    }
    log_info("kappa set to {}", kappa);

    double const lambda = (std::pow(alpha, 2) * (nqf + (order*nuf) + kappa)) - (nqf + (order*nuf));
    double const W0m = lambda / (nqf + (order*nuf) + lambda);
    double const W0c = W0m + (1 - std::pow(alpha, 2) + beta);
    double const Wi = 1 / (2 * (nqf + (order*nuf) + lambda));
    Eigen::VectorXd Wm(2 * (nqf + (order*nuf)) + 1);
    Eigen::VectorXd Wc(2 * (nqf + (order*nuf)) + 1);
    Wm.setConstant(Wi);
    Wc.setConstant(Wi);
    Wm(0) = W0m;
    Wc(0) = W0c;

    Eigen::MatrixXd SS(nqf + (order*nuf), nqf + (order*nuf));
    Eigen::MatrixXd Sigmas(nqf + (order*nuf), 2 * (nqf + (order*nuf)) + 1);
    Eigen::MatrixXd Sigmas2(nqf + (order*nuf), 2 * (nqf + (order*nuf)) + 1);
    Eigen::MatrixXd Sigmaprops(nqf + (order*nuf), 2 * (nqf + (order*nuf)) + 1);
    Eigen::MatrixXd Sigmas2props(3 * ny, 2 * (nqf + (order*nuf)) + 1);

    SimTK::State ss = s;
    std::vector<SimTK::State*> arr_ss;
    SimTK::Vector_<double> qdot;
    Eigen::MatrixXd xsave(nqf + (order*nuf), 1);

    SimTK::Vector_<double> u_old(nu);
    u_old.setToZero();
    SimTK::Vector_<double> u_oldold(nu);
    u_oldold.setToZero();
    Eigen::MatrixXd ysave(3 * ny, 1);
    Eigen::MatrixXd ydiff(3 * ny, 1);

    int ycounter;
    Eigen::MatrixXd Py(3 * ny, 3 * ny);
    Eigen::MatrixXd Pxy(nqf + (order*nuf), 3 * ny);
    Eigen::MatrixXd K(nqf + (order*nuf), 3 * ny);
    Eigen::MatrixXd C((nqf + (order*nuf)), (nqf + (order*nuf)));
    Eigen::MatrixXd stateMean(nqf+(order*nuf), 1);
    SimTK::Matrix A;        //matrix for holonomic constraints
    SimTK::Vector_<double> qerr;
	SimTK::Vec3 angle_vector;
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

    for (int ii = 0; ii < num_cores; ii++) {
        //OpenSim::Model* model_clone = model.clone();
        OpenSim::Model* model_clone = new Model(get_model_file());
        models.push_back(model_clone);  //Alternatively, emplace_back(), but that *should* be slower
        models[ii]->initSystem();
        OpenSim::InverseKinematicsSolver* aSolver = new OpenSim::InverseKinematicsSolver(*(models[ii]), nullptr,
                std::make_shared<OpenSim::OrientationsReference>(oRefs),
                coordRefArray);
        aSolver->setAccuracy(1e-4);
        aSolver->assemble(s);
        solvers.push_back(aSolver);     //Alternatively, emplace_back(), but that *should* be slower
        //delete aSolver;
        //delete model_clone;
    }

    // Modifications to compute average orientations, comparisons between orientations, etc.
    SimTK::Quaternion_<double> simTKquat(1, 0, 0, 0);
    std::vector<SimTK::Quaternion_<double>*> arr_simTKquat;
    Eigen::Quaternion<double> dummyquat(1, 0, 0, 0);
    Eigen::Quaternion<double> dataMinusMeanQuat(1, 0, 0, 0);
    std::vector<Eigen::Quaternion<double>> dataOVector(ny, dummyquat);
    std::vector<Eigen::Quaternion<double>> expectedOVector(ny, dummyquat);
    std::vector<std::vector<Eigen::Quaternion<double>>> Sigmas2Orientations(2 * (nqf + (order*nuf)) + 1, dataOVector);
    std::vector<std::vector<Eigen::Quaternion<double>>> Sigmas2OminusMean(2 * (nqf + (order*nuf)) + 1, dataOVector);
    Eigen::MatrixXd M(4, 2 * (nqf + (order*nuf)) + 1);
    std::vector<Eigen::MatrixXd*> arr_M;
    Eigen::MatrixXd MM(4, 4);
    std::vector<Eigen::MatrixXd*> arr_MM;
    Eigen::EigenSolver<Eigen::MatrixXd> eigSolver;
    std::vector<Eigen::EigenSolver<Eigen::MatrixXd>*> arr_eigSolver;
    Eigen::VectorXd eigenVals(4);
    std::vector<Eigen::VectorXd*> arr_eigenVals;
    Eigen::MatrixXd eigenVecs(4, 4);
    std::vector<Eigen::MatrixXd*> arr_eigenVecs;
    int iMax = 0;
    double maxEigenVal = 0;
    SimTK::Vec4 dummy4vector = simTKquat.asVec4();
    std::vector<SimTK::Vec4*> arr_dummy4vector;
    std::vector<std::map<int, int>*> arr_yMapFromEigenToSimbody;

    for (int ii = 0; ii < num_cores; ii++) {
        std::map<int, int>* aMap = new std::map<int,int>(yMapFromEigenToSimbody);
        arr_yMapFromEigenToSimbody.push_back(aMap);
        SimTK::State* aState = new SimTK::State(s);
        arr_ss.push_back(aState);
        SimTK::Vector_<double>* aQ = new SimTK::Vector_<double>(s.getQ());
        arr_q.push_back(aQ);
        SimTK::Vector_<double>* aU = new SimTK::Vector_<double>(s.getU());
        arr_u.push_back(aU);
        Eigen::MatrixXd* aX = new Eigen::MatrixXd((nqf + (order*nuf)), 1);
        aX->setZero();
        arr_xx.push_back(aX);
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
    
	int step = 0;

    
	log_info("Got to the beginning of for loop.");

    //log_info("before for loop, x(0, 0) = {}", (double)x(0, 0));

    for (auto time : times) { 
        if (time == times[0]) 
        { 
            std::vector<Eigen::MatrixXd> priorStatsVector;
            Eigen::MatrixXd time_export(1,1);
            time_export(0,0) = time;
            priorStatsVector.emplace_back(time_export);
            priorStatsVector.emplace_back(x);   // this won't be used in backward pass
            priorStatsVector.emplace_back(P);   // this won't be used in backward pass
            priorStatsVector.emplace_back(P);   // this won't be used in backward pass
            priorStatsVector.emplace_back(x);
            priorStatsVector.emplace_back(P);  
            {
                std::unique_lock<std::mutex> lock(*fwdBwdMutex);
                priorStatsBuffer->push(priorStatsVector);
            }
            (*condVar).notify_one();
        } 
        else {

            // Step 0. Cholesky factorization of state covariance            
            s.updTime() = time;
            s.updQ() = q;
            s.updU() = u;

            ss.updTime() = time;
            ss.updQ() = q;  //ss.setQ(q) failed to update some elements; is this supposed to happen?
            ss.updU() = u;
            /*
            llt.compute((nqf + (order*nuf) + lambda) * P);
            SS = llt.matrixL();

            // Step 1. Generate posterior sigma points of previous step
            Sigmas.col(0) = x;
            
            stateMean = x;
            for (int ii = 1; ii < (nqf + (order*nuf) + 1); ii++) {
                Sigmas.col(ii) = x + SS.col(ii - 1);
            }
            for (int ii = (nqf + (order*nuf) + 1); ii < (2 * (nqf + (order*nuf)) + 1); ii++) {
                Sigmas.col(ii) = x - SS.col(ii - (nqf + (order*nuf) + 1));
			}            

            // Step 2. Propagate posterior sigma points of previous step through
            // process model
            num_cols = (2 * (nqf + (order*nuf)) + 1);
            elPerThread = num_cols / num_cores;
            for (int ithr = 0; ithr < num_cores; ++ithr) {
                int start = ithr * elPerThread;
                int end = (ithr == num_cores - 1) ? num_cols : (ithr + 1) * elPerThread;
                //threads.emplace_back([&, start, end, x, q, qdot, u, yMapFromEigenToSimbody, yMapFromSimbodyToEigen]() mutable {
                //pool.enqueue([&, start, end, x, q, qdot, u, yMapFromEigenToSimbody, yMapFromSimbodyToEigen]() mutable {
                pool.enqueue([&, start, end, ithr]() mutable {
                    for (int ii = start; ii < end; ii++) {
                        //x = Sigmas.col(ii);
                        arr_xx[ithr]->setZero();
                        for (int irow = 0; irow <= order; irow++) {
                            for (int icol = irow; icol <= order; icol++) {
                                //xx((it->first)+(irow*nuf)) += fCoeffs(irow, icol) * x((it->first)+(icol*nuf));
                                arr_xx[ithr]->block((irow*nuf), 0, nuf, 1) += fCoeffs(irow, icol) * Sigmas.col(ii).block((icol*nuf), 0, nuf, 1);
                            }
                        }
                        
                        Sigmaprops.col(ii) = (*(arr_xx[ithr]));
                    }
                });
            }
            pool.waitUntilCompleted();

            if (get_enable_resampling() == false) {
                Sigmas2 = Sigmaprops;
            }           
            */ 

           /*

            // Step 3. Calculate a priori state mean of current step
            x = W0m * Sigmaprops.col(0);
            for (int jj = 1; jj < (2 * (nqf + (order*nuf)) + 1); jj++) {
                x += Wi * Sigmaprops.col(jj);
            }

            // Step 4. Calculate a priori state covariance of current step
            for (int jj = 0; jj < (2 * (nqf + (order*nuf)) + 1); jj++) {
                Sigmaprops.col(jj) -= x;
            }            
            P = W0c * (Sigmaprops.col(0) * Sigmaprops.col(0).transpose());
            for (int ii = 1; ii < (2 * (nqf + (order*nuf)) + 1); ii++) {
                P += Wi * (Sigmaprops.col(ii) * Sigmaprops.col(ii).transpose());
            }

            // Step 4b. Calculate cross-covariance between states of previous and current steps (for backwards pass)
            C = W0c * ((Sigmas.col(0) - stateMean) * Sigmaprops.col(0).transpose());
            for (int ii = 1; ii < (2 * (nqf + (order*nuf)) + 1); ii++) {
                C += Wi * ((Sigmas.col(ii) - stateMean) * Sigmaprops.col(ii).transpose());
            }

            */

            // Steps 0-4 for linear model with Gaussian noise
            x = F * x;
            C = P * F.transpose();
            P = (F * P * F.transpose()) + Q;
             

            // Note! These are computed for time step 'k', not 'k+1'. Have to handle in the backward pass
            std::vector<Eigen::MatrixXd> priorStatsVector;
            Eigen::MatrixXd time_export(1,1);
            time_export(0,0) = time;
            priorStatsVector.emplace_back(time_export);
            priorStatsVector.emplace_back(x);
            priorStatsVector.emplace_back(P);
            priorStatsVector.emplace_back(C);

            // Resample the sigma points after propagating through process model (optional)
            // if (get_enable_resampling() == true)
            if (true) {
                llt.compute((nqf + (order*nuf) + lambda) * P);
                SS = llt.matrixL();
                Sigmas2.col(0) = x;                
                for (int ii = 1; ii < (nqf + (order*nuf) + 1); ii++) {
                    Sigmas2.col(ii) = x + SS.col(ii - 1);
                }
                for (int ii = (nqf + (order*nuf) + 1); ii < (2 * (nqf + (order*nuf)) + 1); ii++) {
                    Sigmas2.col(ii) = x - SS.col(ii - (nqf + (order*nuf) + 1));
                }
            }

            // Step 5. Propagate a priori sigma points of current step through
            // observation model
            xsave = x;
            num_cols = (2 * (nqf + (order*nuf)) + 1);
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
                        models[ithr]->getMultibodySystem().realize(*(arr_ss[ithr]), SimTK::Stage::Velocity);
                        //model.getMultibodySystem().realize(ss, SimTK::Stage::Velocity);
                        solvers[ithr]->computeCurrentSensorOrientations(*(arr_osensorOrientations[ithr]));
                        //ikSolver.computeCurrentSensorOrientations(osensorOrientations);                        
                        for (int yy = 0; yy < ny; yy++) {
                            (*(arr_dummy4vector[ithr])) = (*(arr_osensorOrientations[ithr]))[yy].convertRotationToQuaternion().asVec4();
                            Sigmas2Orientations[ii][yy] = Eigen::Quaternion<double>((*(arr_dummy4vector[ithr]))(0), (*(arr_dummy4vector[ithr]))(1), (*(arr_dummy4vector[ithr]))(2), (*(arr_dummy4vector[ithr]))(3));
                        }
                        //Sigmas2props.col(ii) = y;
                    }
                });
            }
            pool.waitUntilCompleted();

            // Step 6. Calculate expected observation of current step
            num_cols = ny;
            elPerThread = num_cols / num_cores;
            for (int ithr = 0; ithr < num_cores; ++ithr) {
                int start = ithr * elPerThread;
                int end = (ithr == num_cores - 1) ? num_cols : (ithr + 1) * elPerThread;
                //threads.emplace_back([&, start, end, M, MM, eigSolver, eigenVals, eigenVecs, iMax, maxEigenVal]() mutable {
                pool.enqueue([&, start, end, ithr, iMax, maxEigenVal]() mutable {
                    for (int iy = start; iy < end; iy++) {
                        for (int icol = 0; icol < (2 * (nqf + (order*nuf)) + 1); icol++) {
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
                        expectedOVector[iy] = Eigen::Quaternion<double>((*(arr_eigenVecs[ithr]))(0, iMax), (*(arr_eigenVecs[ithr]))(1, iMax), (*(arr_eigenVecs[ithr]))(2, iMax), (*(arr_eigenVecs[ithr]))(3, iMax));
                    }
                });
            }
            pool.waitUntilCompleted();

            //Step 7. Subtract expected mean orientation from propagated sigma points
            num_cols = (2 * (nqf + (order*nuf)) + 1);
            elPerThread = num_cols / num_cores;
            for (int ithr = 0; ithr < num_cores; ++ithr) {
                int start = ithr * elPerThread;
                int end = (ithr == num_cores - 1) ? num_cols : (ithr + 1) * elPerThread;
                //threads.emplace_back([&, start, end, simTKquat, angle_vector, y]() mutable {
                pool.enqueue([&, start, end, ithr]() mutable {
                    for (int icol = start; icol < end; icol++) {
                        for (int iy = 0; iy < ny; iy++) {
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

            // Step 8. Subtract expected mean orientation from data
            R = R0; //restore covariance of observations; gets modified if missing data

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
                        R(iy * 3 + kk, iy * 3 + kk) *= missingDataScale;
                    }
                } 
                else if (std::isnan(dataOVector[iy].w()) &&
                           std::isnan(dataOVector[iy].x()) &&
                           std::isnan(dataOVector[iy].y()) &&
                           std::isnan(dataOVector[iy].z())) {
                    for (int kk = 0; kk < 3; kk++) {
                        ydiff(iy * 3 + kk) = 0.0;
                        R(iy * 3 + kk, iy * 3 + kk) *= missingDataScale;
                    }                    
                }
                else {
                    dataMinusMeanQuat =
                            expectedOVector[iy] * dataOVector[iy].conjugate();
                    // dataMinusMeanQuat = dataOVector[iy].conjugate() *
                    // expectedOVector[iy];
                    simTKquat = SimTK::Quaternion_<double>(
                            dataMinusMeanQuat.w(), dataMinusMeanQuat.x(),
                            dataMinusMeanQuat.y(), dataMinusMeanQuat.z());
                    angle_vector =
                            SimTK::Rotation_<double>(simTKquat)
                                    .convertThreeAxesRotationToThreeAngles(
                                            SimTK::BodyOrSpaceType::
                                                    SpaceRotationSequence,
                                            SimTK::XAxis, SimTK::YAxis,
                                            SimTK::ZAxis);
                    for (int kk = 0; kk < 3; kk++) {
                        ydiff(iy * 3 + kk) = (double)angle_vector(kk);
                    }
                }
            }

            // Step 9. Calculate covariance of expected orientations
            Py = R + W0c * (Sigmas2props.col(0) * Sigmas2props.col(0).transpose());
            for (int ii = 1; ii < (2 * (nqf + (order*nuf)) + 1); ii++) {
                Py += Wi * (Sigmas2props.col(ii) * Sigmas2props.col(ii).transpose());
            }
			
            // Step 10. Calculate cross-covariance of expected orientations and a
            // priori state mean
            for (int jj = 0; jj < (2 * (nqf + (order*nuf)) + 1); jj++) {
                Sigmas2.col(jj) -= xsave;
            }
            Pxy = W0c * (Sigmas2.col(0) * Sigmas2props.col(0).transpose());
            for (int ii = 1; ii < (2 * (nqf + (order*nuf)) + 1); ii++) {
                Pxy += Wi * (Sigmas2.col(ii) * Sigmas2props.col(ii).transpose());
            }
			
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

            priorStatsVector.emplace_back(x);
            priorStatsVector.emplace_back(P);  
            {
                std::unique_lock<std::mutex> lock(*fwdBwdMutex);
                priorStatsBuffer->push(priorStatsVector);
            }
            (*condVar).notify_one();

            // Step 12b. Update Kalman gain to satisfy position constraints (nope)
            for (std::map<int, int>::iterator it = yMapFromEigenToSimbody.begin(); it != yMapFromEigenToSimbody.end(); ++it) {
                q(it->second) = x(it->first);
            }
            for (std::map<int, int>::iterator it =
                            yMapFromEigenToSimbody.begin();
                    it != yMapFromEigenToSimbody.end(); ++it) {
                u(it->second) = x((it->first)+nqf);
            }
            ss.updQ() = q;
            ss.updU() = u;

            ikSolver.setState(s);
			model.getMultibodySystem().realize(s, SimTK::Stage::Velocity);
			//log_info("Managed to set q and u for state.");
        }


		if (get_report_errors()) {
            ikSolver.computeCurrentOrientationErrors(orientationErrors);
            modelOrientationErrors->appendRow(
                    s.getTime(), orientationErrors);
        }
        if (visualizeResults) {
            model.getVisualizer().show(s);
            log_info("Solved at time: {} s", time);
		}
        else {
            log_info("Solved at time: {} s", time);
		}
        // realize to report to get reporter to pull values from model
        analysisSet.step(s, step++);
        model.realizeReport(s);

        // Abort running inverse kinematics in case of significant numerical instabilities
        if ((double)ydiff.cwiseAbs().mean() > 0.75) {
            log_info("Significant numerical instabilities encountered, aborting...");
            break;
        }
    }	//end of for

    {
        std::unique_lock<std::mutex> lock(*fwdBwdMutex);
        *(fwdDone) = true;
    }
    condVar->notify_all();

    // Delete dynamically allocated stuff
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(solvers);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(models);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_yMapFromEigenToSimbody);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_angle_vector);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_dummy4vector);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_eigenVals);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_eigenVecs);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_eigSolver);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_M);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_MM);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_osensorOrientations);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_q);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_simTKquat);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_ss);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_u);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_xx);
    OpenSim::UKFIMUInverseKinematicsTool::deletePointers(arr_y);

	log_info("Got out of for loop and end of UKFTool.");
	

}  //end of IMUInverseKinematicsTool::UKFTool


std::tuple<std::map<std::string, int>, std::map<int, std::string>> OpenSim::UKFIMUInverseKinematicsTool::CreateYMaps(OpenSim::Model model) {
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
}


double OpenSim::UKFIMUInverseKinematicsTool::computeFactorial(int input) {
    double fact = 1.0;
    for (int ii = 1; ii <= input; ii++) {
        fact *= ii;
    }
    return fact;
}

template <typename T>
void OpenSim::UKFIMUInverseKinematicsTool::deletePointers(std::vector<T*>& vec) {
    for (T* ptr : vec) {
        delete ptr;
    }
    vec.clear();
}

//template <class T>
void OpenSim::UKFIMUInverseKinematicsTool::computeBackwardPass(std::queue<std::vector<Eigen::MatrixXd>>* priorStatsBuffer, 
    std::mutex* fwdBwdMutex, std::condition_variable* condVar, bool* fwdDone, std::map<int, int> yMapFromEigenToSimbody, 
    std::map<int, std::string> yMapFromSimbodyToOpenSim, int nqf, int nuf) {
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
    auto UKFResultsDirQ = UKFresultsDir;
    auto UKFResultsDirPq = UKFresultsDir;
    auto UKFResultsDirU = UKFresultsDir;
    auto UKFResultsDirPu = UKFresultsDir;
    auto UKFResultsDirUdot = UKFresultsDir;
    auto UKFResultsDirPudot = UKFresultsDir;
    UKFResultsDirQ.append("q.txt");
    UKFResultsDirPq.append("Pq.txt");
    std::ofstream fileQ(UKFResultsDirQ, std::ios_base::trunc);
    std::ofstream filePq(UKFResultsDirPq, std::ios_base::trunc);
    fileQ << "endheader" << std::endl << "time\t";
    filePq << "endheader" << std::endl << "time\t";
    UKFResultsDirU.append("u.txt");
    UKFResultsDirPu.append("Pu.txt");
    std::ofstream fileU(UKFResultsDirU, std::ios_base::trunc);
    std::ofstream filePu(UKFResultsDirPu, std::ios_base::trunc);
    fileU << "endheader" << std::endl << "time\t";
    filePu << "endheader" << std::endl << "time\t";
    UKFResultsDirUdot.append("udot.txt");
    UKFResultsDirPudot.append("Pudot.txt");
    std::ofstream fileUdot(UKFResultsDirUdot, std::ios_base::trunc);
    std::ofstream filePudot(UKFResultsDirPudot, std::ios_base::trunc);
    fileUdot << "endheader" << std::endl << "time\t";
    filePudot << "endheader" << std::endl << "time\t";

    //OpenSim::Array<std::string> modelStateVariableNames = model.getStateVariableNames();

    std::regex rgx(".*/(\\w+)/value");
    std::smatch match;
    for (std::map<int, int>::iterator it = yMapFromEigenToSimbody.begin(); it != yMapFromEigenToSimbody.end(); ++it) {
        std::string s = yMapFromSimbodyToOpenSim.at(it->second);
        if (std::regex_search(s,match,rgx)) {
            fileQ << match[1] << "\t";
            fileU << match[1] << "\t";
            fileUdot << match[1] << "\t";
        }
    }

    for (int irow = 0; irow < (nqf); irow++) {
        for (int icol = irow; icol < (nqf); icol++) {
            filePq << "[" << irow << ", " << icol << "]" << "\t";
        }
    }
    for (int irow = 0; irow < (nuf); irow++) {
        for (int icol = irow; icol < (nuf); icol++) {
            filePu << "[" << irow << ", " << icol << "]" << "\t";
        }
    }
    for (int irow = 0; irow < (nuf); irow++) {
        for (int icol = irow; icol < (nuf); icol++) {
            filePudot << "[" << irow << ", " << icol << "]" << "\t";
        }
    }

    fileQ << std::endl;
    filePq << std::endl;
    fileU << std::endl;
    filePu << std::endl;
    fileUdot << std::endl;
    filePudot << std::endl;

    while (true) {
        std::unique_lock<std::mutex> lock(*fwdBwdMutex);
        (*condVar).wait(lock, [&](){ return !(priorStatsBuffer->empty()) || (*fwdDone); });

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

        if ((int)currentPriorMeans.size() >= (lagLength+2)) {
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
                        smoothPosteriorCov = currentPosteriorAutoCovs.at(ii) + D * (currentPosteriorAutoCovs.at(ii+1) - currentPriorAutoCovs.at(ii+1)) * D.transpose(); 
                    }
                    else if (ii < lagLength) {
                        smoothPosteriorMean = currentPosteriorMeans.at(ii) + D * (smoothPosteriorMean - currentPriorMeans.at(ii+1));
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
            if (lagLength < 0) {
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
                            smoothPosteriorCov = currentPosteriorAutoCovs.at(ii) + D * (currentPosteriorAutoCovs.at(ii+1) - currentPriorAutoCovs.at(ii+1)) * D.transpose(); 
                        }
                        else if (ii < (((int)currentPriorMeans.size())-2)) {
                            smoothPosteriorMean = currentPosteriorMeans.at(ii) + D * (smoothPosteriorMean - currentPriorMeans.at(ii+1));
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
                    fileQ.close();
                    filePq.close();
                    fileU.close();
                    filePu.close();
                    fileUdot.close();
                    filePudot.close();
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
        // Write the state to result files
        if (get_write_UKF() && writeToFile) {
            if (order > 1) {
                filePq << static_cast<float>(time) << "\t";
                for (int irow = 0; irow < (nqf); irow++) {
                    for (int icol = irow; icol < (nqf); icol++) {
                        filePq << static_cast<float>(smoothPosteriorCov(irow, icol)) << "\t";
                    }
                }
                fileQ << static_cast<float>(time) << "\t";
                for (int irow = 0; irow < (nqf); irow++) {
                    fileQ << static_cast<float>(smoothPosteriorMean(irow)) << "\t";
                }
                filePu << static_cast<float>(time) << "\t";
                for (int irow = nqf; irow < (nqf+nuf); irow++) {
                    for (int icol = irow; icol < (nqf+nuf); icol++) {
                        filePu << static_cast<float>(smoothPosteriorCov(irow, icol)) << "\t";
                    }
                }
                fileU << static_cast<float>(time) << "\t";
                for (int irow = nqf; irow < (nqf+nuf); irow++) {
                    fileU << static_cast<float>(smoothPosteriorMean(irow)) << "\t";
                }
                filePudot << static_cast<float>(time) << "\t";
                for (int irow = (nqf+nuf); irow < (nqf+2*nuf); irow++) {
                    for (int icol = irow; icol < (nqf+2*nuf); icol++) {
                        filePudot << static_cast<float>(smoothPosteriorCov(irow, icol)) << "\t";
                    }
                }
                fileUdot << static_cast<float>(time) << "\t";
                for (int irow = (nqf+nuf); irow < (nqf+2*nuf); irow++) {
                    fileUdot << static_cast<float>(smoothPosteriorMean(irow)) << "\t";
                }
                fileQ << std::endl;
                filePq << std::endl;
                fileU << std::endl;
                filePu << std::endl;
                fileUdot << std::endl;
                filePudot << std::endl;
            }
            else if (order == 1) {
                filePq << static_cast<float>(time) << "\t";
                for (int irow = 0; irow < (nqf); irow++) {
                    for (int icol = irow; icol < (nqf); icol++) {
                        filePq << static_cast<float>(smoothPosteriorCov(irow, icol)) << "\t";
                    }
                }
                fileQ << static_cast<float>(time) << "\t";
                for (int irow = 0; irow < (nqf); irow++) {
                    fileQ << static_cast<float>(smoothPosteriorMean(irow)) << "\t";
                }
                filePu << static_cast<float>(time) << "\t";
                for (int irow = nqf; irow < (nqf+nuf); irow++) {
                    for (int icol = irow; icol < (nqf+nuf); icol++) {
                        filePu << static_cast<float>(smoothPosteriorCov(irow, icol)) << "\t";
                    }
                }
                fileU << static_cast<float>(time) << "\t";
                for (int irow = nqf; irow < (nqf+nuf); irow++) {
                    fileU << static_cast<float>(smoothPosteriorMean(irow)) << "\t";
                }
                fileQ << std::endl;
                filePq << std::endl;
                fileU << std::endl;
                filePu << std::endl;
            }
            else {
                filePq << static_cast<float>(time) << "\t";
                for (int irow = 0; irow < (nqf); irow++) {
                    for (int icol = irow; icol < (nqf); icol++) {
                        filePq << static_cast<float>(smoothPosteriorCov(irow, icol)) << "\t";
                    }
                }
                fileQ << static_cast<float>(time) << "\t";
                for (int irow = 0; irow < (nqf); irow++) {
                    fileQ << static_cast<float>(smoothPosteriorMean(irow)) << "\t";
                }
                fileQ << std::endl;
                filePq << std::endl;
            }
        }
    } 
}



