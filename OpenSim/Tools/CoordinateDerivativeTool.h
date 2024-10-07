#ifndef COORDINATEDERIVATIVETOOL_H
#define COORDINATEDERIVATIVETOOL_H

/*
 * CoordinateDerivativeTool is written by University of Eastern Finland  *
 * and bases on the InverseDynamicsTool code written by Ajay Seth.       *
 *                                                                       *
 * The project was supported by the Research Council of Finland.         *
 *                                                                       *
 * Copyright (c) 2024 University of Eastern Finland                      *
 * Author(s): Matti Kortelainen                                          */

#include <OpenSim/Common/Storage.h>
#include "DynamicsTool.h"
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/FunctionSet.h>
#include <OpenSim/Common/GCVSplineSet.h>
#include <OpenSim/Common/IO.h>
//#include <OpenSim/Common/Stopwatch.h>
#include <OpenSim/Common/XMLDocument.h>
#include <OpenSim/Simulation/InverseDynamicsSolver.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimulationUtilities.h>


#ifdef SWIG
    #ifdef OSIMTOOLS_API
        #undef OSIMTOOLS_API
        #define OSIMTOOLS_API
    #endif
#endif

namespace OpenSim {

class Model;
class JointSet;

class OSIMTOOLS_API CoordinateDerivativeTool: public DynamicsTool {
    OpenSim_DECLARE_CONCRETE_OBJECT(CoordinateDerivativeTool, DynamicsTool);

private:
    OpenSim::Storage* _statesValues;
    OpenSim::Model* _model;


public:
    virtual ~CoordinateDerivativeTool();
    CoordinateDerivativeTool(std::string modelFileName, std::string statesFileName, std::string outputFileName);
    CoordinateDerivativeTool();
    void setModelFile(std::string modelFileName);
    void setStatesFile(std::string statesFileName);
    void setOutputFileName(std::string outputFileName);
    bool run() override SWIG_DECLARE_EXCEPTION;

protected:
    std::string _statesFileName;
    OpenSim::PropertyStr _statesFileNameProp;
    std::string _modelFileName;
    OpenSim::PropertyStr _modelFileNameProp;
    std::string _outputFileName;
    //double _startTime;
    //double _endTime;
    bool loadCoordinatesFromFile();

};
}

#endif // COORDINATEDERIVATIVETOOL_H
