#include"CoordinateDerivativeTool.h"

OpenSim::CoordinateDerivativeTool::~CoordinateDerivativeTool() {
    if (_model != NULL) {
        delete _model;
        _model = NULL;
    }
    if (_statesValues != NULL) {
        delete _statesValues;
        _statesValues = NULL;
    }
}

OpenSim::CoordinateDerivativeTool::CoordinateDerivativeTool(std::string modelFile, std::string statesFile, std::string outputFile, std::string outputDirectory, double cutoffFreq, int splineDegree, int outputMaxDerOrder, bool outputDegrees) {
    _modelFileNameProp.setValue(modelFile);
    _modelFileName = _modelFileNameProp.getValueStr();
    _model = new OpenSim::Model(_modelFileName);
    _statesFileNameProp.setValue(statesFile);
    _statesFileName = _statesFileNameProp.getValueStr();
    _statesValues = new OpenSim::Storage(_statesFileName);
    _lowpassCutoffFrequencyProp.setValue(cutoffFreq);
    _lowpassCutoffFrequency = cutoffFreq;
    _outputFileName = outputFile;
    _outputDirectory = outputDirectory;
    _splineDegree = splineDegree;
    _outputMaxDerOrder = outputMaxDerOrder;
    _outputDegrees = outputDegrees;
}

/*
OpenSim::CoordinateDerivativeTool::CoordinateDerivativeTool() : CoordinateDerivativeTool(_modelFileName, _statesFileName, _outputFileName, _lowpassCutoffFrequency, _splineDegree, _outputMaxDerOrder) {
    //_modelFileName = _modelFileNameProp.getValueStr();
    if (_modelFileName != "" && _modelFileName != "Unassigned") {
        OPENSIM_THROW_IF_FRMOBJ(_modelFileName.empty(), OpenSim::Exception,
                "No model filename was provided.")
        _model = new OpenSim::Model(_modelFileName);
    }
    if (_statesFileName != "" && _statesFileName != "Unassigned") {
        OPENSIM_THROW_IF_FRMOBJ(_statesFileName.empty(), OpenSim::Exception,
                "No states filename was provided.")
        _statesValues = new OpenSim::Storage(_statesFileName);
    }
    if (_outputFileName == "" || _outputFileName == "Unassigned" || _outputFileName.empty()) {
        _outputFileName = "derivedKinematics";
    }
}
*/

void OpenSim::CoordinateDerivativeTool::setModelFile(std::string modelFileName) {
    _modelFileName = modelFileName;
    _modelFileNameProp.setValue(modelFileName);
}

void OpenSim::CoordinateDerivativeTool::setStatesFile(std::string statesFileName) {
    _statesFileName = statesFileName;
    _statesFileNameProp.setValue(statesFileName);
}

void OpenSim::CoordinateDerivativeTool::setOutputFileName(std::string outputFileName) {
    _outputFileName = outputFileName;
}

void OpenSim::CoordinateDerivativeTool::setOutputDirectory(std::string outputDirectory) {
    _outputDirectory = outputDirectory;
}

std::string OpenSim::CoordinateDerivativeTool::getOutputDirectory() {
    return _outputDirectory;
}

void OpenSim::CoordinateDerivativeTool::setCutoffFrequency(double cutoffFreq) {
    _lowpassCutoffFrequency = cutoffFreq;
}

void OpenSim::CoordinateDerivativeTool::setSplineDegree(int splineDegree) {
    _splineDegree = splineDegree;
}

void OpenSim::CoordinateDerivativeTool::setOutputMaxDerOrder(int outputMaxDerOrder) {
    _outputMaxDerOrder = outputMaxDerOrder;
}

void OpenSim::CoordinateDerivativeTool::setOutputDegrees(bool outputDegrees) {
    _outputDegrees = outputDegrees;
}


bool OpenSim::CoordinateDerivativeTool::run() {
    bool success = false;
    bool modelFromFile=true;
    try{
        //Load and create the indicated model
        if (_model == NULL) {
            OPENSIM_THROW_IF_FRMOBJ(_modelFileName.empty(), OpenSim::Exception,
                    "No model filename was provided.")
            _model = new OpenSim::Model(_modelFileName);
        }
        else {
            modelFromFile = false;
        }

        _model->finalizeFromProperties();
        _model->printBasicInfo();

        log_info("Running tool {}...", getName());
        // Do the maneuver to change then restore working directory
        // so that the parsing code behaves properly if called from a different directory.
        auto cwd = OpenSim::IO::CwdChanger::changeToParentOf(getDocumentFileName());

        SimTK::State& s = _model->initSystem();

               // OpenSim::Coordinates represent degrees of freedom for a model.
               // Each Coordinate's value and speed maps to an index
               // in the model's underlying SimTK::State (value to a slot in the
               // State's q, and speed to a slot in the State's u).
               // So we need to map each OpenSim::Coordinate value and speed to the
               // corresponding SimTK::State's q and u indices, respectively.
        auto coords = _model->getCoordinatesInMultibodyTreeOrder();
        int nq = s.getNQ();
        int nu = s.getNU();
        int nCoords = (int)coords.size();
        int intUnusedSlot = -1;

               // Create a vector mapCoordinateToQ whose i'th element provides
               // the index in vector 'coords' that corresponds with the i'th 'q' value.
               // Do the same for mapCoordinateToU, which tracks the order for
               // each 'u' value.
        auto coordMap = OpenSim::createSystemYIndexMap(*_model);
        std::vector<int> mapCoordinateToQ(nq, intUnusedSlot);
        std::vector<int> mapCoordinateToU(nu, intUnusedSlot);
        for (const auto& c : coordMap) {
            // SimTK::State layout is [q u z]
            // (i.e., all "q"s first then "u"s second).
            if (c.second < nq + nu) {
                std::string svName = c.first;

                // The state names corresponding to q's and u's will be of
                // the form:
                // /jointset/(joint)/(coordinate)/(value or speed).
                // So the corresponding coordinate name is second from the end.
                OpenSim::ComponentPath svPath(svName);
                std::string lastPathElement = svPath.getComponentName();
                std::string coordName = svPath.getSubcomponentNameAtLevel(
                        svPath.getNumPathLevels() - 2);

                for (int i = 0; i < nCoords; ++i) {
                    if (coordName == coords[i]->getName()) {
                        if (lastPathElement == "value") {
                            mapCoordinateToQ[c.second] = i;
                            break;
                        }

                        // Shift State/System indices by nq since u's follow q's
                        else if (lastPathElement == "speed") {
                            mapCoordinateToU[c.second - nq] = i;
                            break;
                        }

                        else {
                            throw Exception("Last element in state variable "
                                            " name " + svName + " is neither "
                                                     "'value' nor 'speed'");
                        }
                    }
                    if (i == nCoords - 1) {
                        throw Exception("Coordinate " + coordName +
                                        " not found in model.");
                    }
                }

            }
        }

               // Make sure that order of coordFunctions (which define splines for
               // State's q's) is in the same order as the State's q order.
               // Also make a new vector (coordinatesToSpeedsIndexMap) that, for each
               // u in the State, gives the corresponding index in q (which is same
               // order as  coordFunctions). This accounts for cases where qdot != u.
        OpenSim::FunctionSet coordFunctions;
        coordFunctions.ensureCapacity(nq);
        std::vector<int> coordinatesToSpeedsIndexMap(nu, intUnusedSlot);

        // Take care of the time (just in case filtering does the resampling)
        double first_time = _statesValues->getFirstTime();
        double last_time = _statesValues->getLastTime();

        // Determine the starting and final time for the Tool by comparing to what data is available
        double start_time = ( first_time > _timeRange[0]) ? first_time : _timeRange[0];
        double final_time = ( last_time < _timeRange[1]) ? last_time : _timeRange[1];
        int start_index = _statesValues->findIndex(start_time);
        int final_index = _statesValues->findIndex(final_time);

        int nt = final_index-start_index+1; //number of rows in data

        SimTK::Array_<double> times(nt, 0.0);
        for(int i=0; i<nt; i++){
            times[i]=_statesValues->getStateVector(start_index+i)->getTime();
        }

        if (loadCoordinatesFromFile()){
            if(_lowpassCutoffFrequency>=0) {
                log_info("Low-pass filtering coordinates data with a cutoff "
                         "frequency of {}...", _lowpassCutoffFrequency);
                _statesValues->pad(_statesValues->getSize()/2);
                _statesValues->lowpassIIR(_lowpassCutoffFrequency);
            }

            // Convert degrees to radian if indicated
            if(_statesValues->isInDegrees()){
                _model->getSimbodyEngine().convertDegreesToRadians(*_statesValues);
            }
            // Create differentiable splines of the coordinate data
            OpenSim::GCVSplineSet coordSplines(_splineDegree, _statesValues);

                   // Functions must correspond to model coordinates.
                   // Solver needs the order of Function's to be the same as order
                   // in State's q's.
            for (int i = 0; i < nq; i++) {
                int coordInd = mapCoordinateToQ[i];

                       // unused q slot
                if (coordInd == intUnusedSlot) {
                    coordFunctions.insert(i, new OpenSim::Constant(0));
                    continue;
                }

                const OpenSim::Coordinate& coord = *coords[coordInd];
                if (coordSplines.contains(coord.getName())) {
                    coordFunctions.insert(i, coordSplines.get(coord.getName()));
                }
                else{
                    coordFunctions.insert(i,new OpenSim::Constant(coord.getDefaultValue()));
                    log_info("InverseDynamicsTool: coordinate file does not "
                             "contain coordinate '{}'. Assuming default value.",
                            coord.getName());
                }

                       // Fill in coordinatesToSpeedsIndexMap as we go along to make
                       // sure we know which function corresponds to State's u's.
                for (int j = 0; j < nu; ++j) {
                    if (mapCoordinateToU[j] == coordInd) {
                        coordinatesToSpeedsIndexMap[j] = i;
                    }
                }
            }
            if(coordFunctions.getSize() > nq){
                coordFunctions.setSize(nq);
            }
        }
        else{
            throw OpenSim::Exception("InverseDynamicsTool: no coordinate file found, "
                            " or setCoordinateValues() was not called.");
        }        

        // We want to output the positons, velocities and accelerations

        OpenSim::Storage statesOut(nt);

        for (int it = 0; it < nt; it++) {

            //s.updTime() = times[it];
            SimTK::Vector_<double> stateCoords(nq + _outputMaxDerOrder*nu);
            stateCoords.setToZero();
            //SimTK::Vector &q = s.updQ();
            //SimTK::Vector &u = s.updU();
            //SimTK::Vector &udot = s.updUDot();

                   // Account for cases where qdot != u with coordinatesToSpeedsIndexMap
            for (int iord = 0; iord <= _outputMaxDerOrder; iord++) {
                if (iord == 0) {
                    for( int j = 0; j < nq; ++j) {
                        stateCoords(j) = coordFunctions.evaluate(j, 0, times[it]);
                    }
                }
                else {
                    for (int j = 0; j < nu; ++j) {
                        stateCoords(nq+((iord-1) * nu)+j) = coordFunctions.evaluate(
                            coordinatesToSpeedsIndexMap[j], iord, times[it]);
                    }
                }
            }

            OpenSim::StateVector state(times[it], stateCoords);
            statesOut.append(state);
        }

        OpenSim::Array<std::string> labels("time", nq + _outputMaxDerOrder*nu + 1);

        for (int iord = 0; iord <= _outputMaxDerOrder; iord++) {
            if (iord == 0) {
                for (int ii = 0; ii < nq; ii++) {
                    labels[ii+1] = coords[ii]->getName();
                    if (_outputMaxDerOrder > 0) {
                        labels[ii+1] += "_[";
                        labels[ii+1] += std::to_string(iord);
                        labels[ii+1] += "]";
                    }
                }
            }
            else {
                for (int ii = 0; ii < nu; ii++) {
                    labels[nq+((iord-1)*nu)+ii+1] = coords[coordinatesToSpeedsIndexMap[ii]]->getName();
                    labels[nq+((iord-1)*nu)+ii+1] += "_[";
                    labels[nq+((iord-1)*nu)+ii+1] += std::to_string(iord);
                    labels[nq+((iord-1)*nu)+ii+1] += "]";
                }
            }
        }        
        
        statesOut.setColumnLabels(labels);
        statesOut.setName("Positions, Velocities and Accelerations");

        // Convert results to degrees if required
        if (_outputDegrees) {
            if (!statesOut.isInDegrees()) {
                _model->getSimbodyEngine().convertRadiansToDegrees(statesOut);
            }
        }

        OpenSim::IO::makeDir(getOutputDirectory());
        OpenSim::Storage::printResult(&statesOut, _outputFileName, getOutputDirectory(), -1, ".sto");
        cwd.restore();
        success = true;
    }
    catch (const OpenSim::Exception& ex) {
        log_error("CoordinateDerivativeTool Failed: {}", ex.what());
        throw (Exception("CoordinateDerivativeTool Failed, please see messages window for details..."));
    }

    if (modelFromFile) {
        delete _model;
    }
    return success;
}

bool OpenSim::CoordinateDerivativeTool::loadCoordinatesFromFile() {
    if (_statesValues != NULL) {
        return true;
    }

    if (_statesFileName != "" && _statesFileName != "Unassigned") {
        _statesValues = new OpenSim::Storage(_statesFileName);
        _statesValues->setName(_statesFileName);
        return true;
    }
    return false;
}
