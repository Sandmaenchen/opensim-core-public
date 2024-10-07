
#include "IMUForwardKinematicsTool.h"

//using namespace OpenSim;
//using namespace SimTK;
//using namespace std;

OpenSim::IMUForwardKinematicsTool::IMUForwardKinematicsTool()
        : OpenSim::InverseKinematicsToolBase() {
    OpenSim::IMUForwardKinematicsTool::constructProperties();
}

OpenSim::IMUForwardKinematicsTool::IMUForwardKinematicsTool(const std::string& setupFile)
        : OpenSim::InverseKinematicsToolBase(setupFile, true) {
    OpenSim::IMUForwardKinematicsTool::constructProperties();
    updateFromXMLDocument();
}

OpenSim::IMUForwardKinematicsTool::~IMUForwardKinematicsTool()
{
}

void OpenSim::IMUForwardKinematicsTool::constructProperties()
{
    OpenSim::IMUForwardKinematicsTool::constructProperty_opensim_to_sensor_rotations(
            SimTK::Vec3(0));
    OpenSim::IMUForwardKinematicsTool::constructProperty_qStates_file("");
    OpenSim::OrientationWeightSet orientationWeights;
    OpenSim::IMUForwardKinematicsTool::constructProperty_orientation_weights(
            orientationWeights);
}

void OpenSim::IMUForwardKinematicsTool::runForwardKinematicsWithQStatesFromFile(
        OpenSim::Model& model, const std::string& qStatesFileName,
        double rollRMSinDeg, double headingRMSinDeg, double pitchRMSinDeg,
        double stateRMS, unsigned seed) {

    // Create noisy q vector for forward kinematics simu
    OpenSim::TimeSeriesTable_<SimTK::Real> qStatesTableDeg(qStatesFileName);
    OpenSim::TimeSeriesTable_<SimTK::Real> qStatesTableRad(qStatesFileName);
    log_info("Loading qStates as degrees from '{}'...",
            qStatesFileName);

    // Will maintain only data in time range specified by the tool
    // If unspecified {-inf, inf} no trimming is done
    qStatesTableDeg.trim(getStartTime(), getEndTime());
    qStatesTableRad.trim(getStartTime(), getEndTime());

    // Check if random seed has been set. if not, set by system clock
    if (seed == 0) {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    std::normal_distribution<double> normalDist;
    if (stateRMS <= 0.0) {
        normalDist = std::normal_distribution<double>(0.0, 0.00001);
    }
    else {
        normalDist = std::normal_distribution<double>(0.0, stateRMS);
    }

    // Create mappings between continuous state variables between OpenSim and Simbody
    std::tuple<std::map<std::string, int>, std::map<int, std::string>> mappings = OpenSim::IMUForwardKinematicsTool::CreateYMaps(model);
    std::map<std::string, int> yMapFromOpenSimToSimbody = std::get<0>(mappings);
    std::map<int, std::string> yMapFromSimbodyToOpenSim = std::get<1>(mappings);

    SimTK::State& s = model.initSystem();
    int nq = model.getNumCoordinates();

    // Create mapping from the data table to OpenSim model
    std::vector<std::string> qStateTableNames = qStatesTableDeg.getColumnLabels();
    OpenSim::Array<std::string> qStateModelNames = model.getStateVariableNames();
    std::map<int, std::string> yMapFromQStateTableToOpenSim;
    for (int icol = 0; icol < (int)qStateTableNames.size(); icol++) {
        for (int ii = 0; ii < qStateModelNames.size(); ii++) {
            if (qStateModelNames[ii].find(qStateTableNames[icol]) != std::string::npos) {
                yMapFromQStateTableToOpenSim.insert(std::pair<int, std::string>(icol, qStateModelNames[ii]));
            }
        }
    }

    // Find locked coordinates (we cannot change these in UKF)
    std::map<int, int> simbodyIndexTypeMap;
    //std::map<int, int> yMapFromSimbodyToEigen;
    //std::map<int, int> yMapFromEigenToSimbody;

    int iqx;
    int nqf = 0;    //number of free and clamped q's
    for (const auto& coord : model.getComponentList<OpenSim::Coordinate>()) {
        iqx = yMapFromOpenSimToSimbody[coord.getStateVariableNames()[0]];
        if (coord.getMotionType() == OpenSim::Coordinate::MotionType::Translational) {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 0));    //translational coordinate
        }
        else if (coord.getLocked(s)) {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 0));    //locked coordinate
        }
        else if (coord.getClamped(s)) {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 1));    //clamped coordinate
            nqf++;
        }
        else {
            simbodyIndexTypeMap.insert(std::pair<int, int>(iqx, 2));    //free coordinate
            nqf++;
        }
    }
    log_info("Found {} free coordinates from the model.", nqf);


    // Generate normally distributed noise to state q's
    SimTK::Matrix_<double> qStatesRand_matrix(qStatesTableDeg.getMatrix());
    std::default_random_engine generator(seed);

    // THIS COULD BE THE PROBLEMATIC PART! THE TIME COLUMN MAY NOT BELONG TO ALL
    if ((int)qStateTableNames.size() != qStatesRand_matrix.ncol()) {
        log_info("time is not included in the matrix; matrix ncol = {}.", qStatesRand_matrix.ncol());
        log_info("number of columns in data = {}", (int)qStateTableNames.size());
        log_info("number of columns in data tables = {}", qStatesTableRad.getNumColumns());
    }
    else {
        log_info("time is included in the matrix");
        log_info("number of columns in data tables = {}", qStatesTableRad.getNumColumns());
        log_info("number of columns in random matrix = {}", qStatesRand_matrix.ncol());
    }

    for (int irow=0; irow < (int)qStatesTableDeg.getNumRows(); irow++) {
        for (int icol=0; icol < (int)qStatesTableDeg.getNumColumns(); icol++) {
            if (simbodyIndexTypeMap.at(yMapFromOpenSimToSimbody[yMapFromQStateTableToOpenSim[icol]]) != 0) {
                qStatesRand_matrix(irow,icol) += normalDist(generator);
                qStatesRand_matrix(irow,icol) = (qStatesRand_matrix(irow,icol) * SimTK::Pi / 180);
            }
        }
        auto rowRef = qStatesTableRad.updNearestRow(qStatesTableDeg.getIndependentColumn()[irow], false);
        for (int icol=0; icol < (int)qStatesTableDeg.getNumColumns(); icol++) {
            rowRef(icol) = qStatesRand_matrix(irow,icol);
        }
    }

    // Find the names of IMUs from the model to create dummy data
    auto bodySetList = model.getComponentList<OpenSim::PhysicalOffsetFrame>();
    std::vector<std::string> imuNames;
    for (auto it = bodySetList.begin(); it != bodySetList.end(); it++) {
        if (it->getName().find("_imu_") != std::string::npos) {
            continue;
        }
        else if (it->getName().find("_imu") != std::string::npos) {
            imuNames.emplace_back(it->getName());
        }
        else {
            continue;
        }
    }
    log_info("Total of {} IMUs found in the model", imuNames.size());
    for (int ii = 0; ii < (int)imuNames.size(); ii++) {
        log_info(imuNames[ii]);
    }


    // Generate a dummy orientation reference so we can use the IK solver
    std::vector<double> tableTimes = qStatesTableRad.getIndependentColumn();
    //tableTimes.emplace_back(0.0);
    SimTK::Matrix_<SimTK::Rotation> dummyData((int)qStatesTableRad.getNumRows(), (int)imuNames.size());
    for (int irow = 0; irow < (int)qStatesTableRad.getNumRows(); irow++) {
        for (int ii = 0; ii < (int)imuNames.size(); ii++) {
            dummyData(irow, ii) = SimTK::Rotation(
                    SimTK::BodyOrSpaceType::SpaceRotationSequence,
                    0.0, SimTK::XAxis, 0.0, SimTK::YAxis, 0.0, SimTK::ZAxis);
        }
    }

    OpenSim::TimeSeriesTable_<SimTK::Rotation> dummyTable(tableTimes, dummyData,imuNames);
    OpenSim::OrientationsReference oRefs(dummyTable, &get_orientation_weights());
    std::shared_ptr<OpenSim::OrientationsReference> oRefPtr(&oRefs);
    SimTK::Array_<OpenSim::CoordinateReference> coordinateReferences;
    //log_info("Created dummy data");

    // Create IK solver to calculate orientations corresponding to state
    const double accuracy = 1e-4;
    //std::make_shared<OpenSim::OrientationsReference>(oRefs)
    OpenSim::InverseKinematicsSolver ikSolver(model, nullptr,
            oRefPtr,
            coordinateReferences);
    ikSolver.setAccuracy(accuracy);
    //log_info("Created ikSolver");

    int ny = (int)imuNames.size();
    std::vector<double> times = qStatesTableRad.getIndependentColumn();

    s.updTime() = times[0];
    ikSolver.assemble(s);

    // Get the IMU names from model (correct order)
    SimTK::Array_<std::string> oModelNames;
    for (int ii = 0; ii < ny; ii++) {
        oModelNames.push_back(ikSolver.getOrientationSensorNameForIndex(ii));
    }

    SimTK::Vector_<double> u(nq);
    u.setToZero();
    SimTK::Array_<SimTK::Rotation_<double>> osensorOrientations;
    SimTK::Vec3 dummy3vector;
    std::normal_distribution<double> XDist;
    std::normal_distribution<double> YDist;
    std::normal_distribution<double> ZDist;
    if (rollRMSinDeg > 0.0) {
        XDist = std::normal_distribution<double>(0.0, (rollRMSinDeg * SimTK::Pi / 180));
    }
    else {
        XDist = std::normal_distribution<double>(0.0, 0.0001);
    }
    if (headingRMSinDeg > 0.0) {
        YDist = std::normal_distribution<double>(0.0, (headingRMSinDeg * SimTK::Pi / 180));
    }
    else {
        YDist = std::normal_distribution<double>(0.0, 0.0001);
    }
    if (pitchRMSinDeg > 0.0) {
        ZDist = std::normal_distribution<double>(0.0, (pitchRMSinDeg * SimTK::Pi / 180));
    }
    else {
        ZDist = std::normal_distribution<double>(0.0, 0.0001);
    }
    //log_info("Created distributions");

    // FIX THE CORRECT DIRECTORIES ETC.
    auto resultsDir = get_results_directory();
    if (resultsDir.empty() && !get_output_motion_file().empty()) {
        resultsDir = OpenSim::IO::getParentDirectory(get_output_motion_file());
    }
    if (!resultsDir.empty()) {
        OpenSim::IO::makeDir(resultsDir);
    }
    //resultsDir.append("/simulated_orientations/");
    //OpenSim::IO::makeDir(resultsDir);

    resultsDir.append("/simulatedOrientations.sto");

    std::ofstream fileO(resultsDir, std::ios_base::out);

    fileO << "DataRate=" << static_cast<float>(std::round(1 / (times[1] - times[0]))) << ".0" << std::endl;
    fileO << "DataType=Quaternion" << std::endl;
    fileO << "random_seed=" << static_cast<float>(seed) << std::endl;
    fileO << "stateRMSinDeg=" << static_cast<float>(stateRMS) << std::endl;
    fileO << "rollRMSinDeg=" << static_cast<float>(rollRMSinDeg) << std::endl;
    fileO << "headingRMSinDeg=" << static_cast<float>(headingRMSinDeg) << std::endl;
    fileO << "pitchRMSinDeg=" << static_cast<float>(pitchRMSinDeg) << std::endl;
    fileO << "endheader" << std::endl;
    fileO << "time\t";
    for (int yy = 0; yy < ny; yy++) {
        if (yy < (ny-1)) {
            fileO << oModelNames[yy] << "\t";
        }
        else {
            fileO << oModelNames[yy];
        }
    }
    fileO << std::endl;

    SimTK::Vector_<double> q(nq);
    std::vector<SimTK::Quaternion> quaternions;
    quaternions.reserve(ny);
    for (int irow=0; irow < (int)qStatesTableDeg.getNumRows(); irow++) {
        s.updTime() = times[irow];
        for (int icol = 0; icol < (int)qStatesTableRad.getNumColumns(); icol++) {
            q(yMapFromOpenSimToSimbody[yMapFromQStateTableToOpenSim[icol]]) = qStatesTableRad.getNearestRow(times[irow], false)(icol);
        }
        s.updQ() = q;
        s.updU() = u;
        ikSolver.setState(s);
        model.getMultibodySystem().realize(s, SimTK::Stage::Velocity);
        ikSolver.computeCurrentSensorOrientations(osensorOrientations);
        for (int yy = 0; yy < ny; yy++) {
            dummy3vector = osensorOrientations[yy].convertThreeAxesRotationToThreeAngles(
                    SimTK::BodyOrSpaceType::SpaceRotationSequence, SimTK::XAxis, SimTK::YAxis, SimTK::ZAxis);
            if (rollRMSinDeg > 0.0) {
                dummy3vector(0) += XDist(generator);
            }
            if (headingRMSinDeg > 0.0) {
                dummy3vector(1) += YDist(generator);
            }
            if (pitchRMSinDeg > 0.0) {
                dummy3vector(2) += ZDist(generator);
            }
            quaternions[yy] = (SimTK::Rotation(SimTK::BodyOrSpaceType::SpaceRotationSequence,
                    dummy3vector(0), SimTK::XAxis, dummy3vector(1), SimTK::YAxis,
                    dummy3vector(2), SimTK::ZAxis).convertRotationToQuaternion());
        }
        fileO << static_cast<float>(times[irow]) << "\t";
        for (int yy = 0; yy < ny; yy++) {
            for (int icol = 0; icol < 4; icol++) {
                if (icol < 3) {
                    fileO << static_cast<float>(quaternions[yy](icol)) << ",";
                }
                else if (yy < (ny-1)) {
                    fileO << static_cast<float>(quaternions[yy](icol)) << "\t";
                }
                else {
                    fileO << static_cast<float>(quaternions[yy](icol));
                }
            }
        }
        fileO << std::endl;
    }
    log_info("Wrote all quaternions to the file.");
    fileO.close();
    log_info("Closed file.");
}


// main driver
bool OpenSim::IMUForwardKinematicsTool::run(double rollRMSinDeg, double headingRMSinDeg, double pitchRMSinDeg,
        double stateRMS, unsigned seed)
{
    if (_model.empty()) {
        _model.reset(new Model(get_model_file()));
    }

    OpenSim::IMUForwardKinematicsTool::runForwardKinematicsWithQStatesFromFile(*_model,
            get_qStates_file(), rollRMSinDeg, headingRMSinDeg, pitchRMSinDeg, stateRMS, seed);

    return true;
}

OpenSim::TimeSeriesTable_<SimTK::Vec3> OpenSim::IMUForwardKinematicsTool::loadMarkersFile(const std::string& markerFile)
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



std::tuple<std::map<std::string, int>, std::map<int, std::string>> OpenSim::IMUForwardKinematicsTool::CreateYMaps(OpenSim::Model model) {
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
