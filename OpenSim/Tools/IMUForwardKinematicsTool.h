#ifndef IMUFORWARDKINEMATICSTOOL_H
#define IMUFORWARDKINEMATICSTOOL_H

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
//
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/AssemblyCondition_OrientationSensors.h"
#include <OpenSim/Tools/IMUInverseKinematicsTool.h>
#include <OpenSim/Common/ComponentList.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>

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
 * This is a tool based on the IMU inverse kinematics tool source code originally 
 * written by Ajay Seth.
 * This tool computes (noisy) IMU observations for a time series of coordinate 
 * positions for a given MS model.
 *
 * @author Matti Kortelainen
 */
class OSIMTOOLS_API IMUForwardKinematicsTool
        : public InverseKinematicsToolBase {
    OpenSim_DECLARE_CONCRETE_OBJECT(
            IMUForwardKinematicsTool, InverseKinematicsToolBase);

public:
    OpenSim_DECLARE_PROPERTY(qStates_file, std::string,
            "Name/path to a .sto file of qStates as degrees.");

    OpenSim_DECLARE_PROPERTY(opensim_to_sensor_rotations, SimTK::Vec3,
            "Space fixed Euler angles (XYZ order) from OpenSim to IMU Space."
            " Default to (0, 0, 0).");

    OpenSim_DECLARE_PROPERTY(orientation_weights, OrientationWeightSet,
            "Set of orientation weights identified by orientation name with "
            "weight being a positive scalar. If not provided, all IMU "
            "orientations are tracked with weight 1.0.");

    //=============================================================================
    // METHODS
    //=============================================================================
    //--------------------------------------------------------------------------
    // CONSTRUCTION
    //--------------------------------------------------------------------------
public:
    virtual ~IMUForwardKinematicsTool();
    IMUForwardKinematicsTool();
    IMUForwardKinematicsTool(const std::string &setupFile);
    //--------------------------------------------------------------------------
    // INTERFACE
    //--------------------------------------------------------------------------

    bool run(double rollRMSinDeg, double headingRMSinDeg, double pitchRMSinDeg,
            double stateRMS = 0.0, unsigned seed = 0) SWIG_DECLARE_EXCEPTION;
    bool run() override SWIG_DECLARE_EXCEPTION {
        return run(0.75, 1.50, 0.75, 0.0, 0); //The roll, heading, and pitch RMS default values as given for Xsens MTW2-3A7G6
    };

    static TimeSeriesTable_<SimTK::Vec3>
    loadMarkersFile(const std::string& markerFile);

    void runForwardKinematicsWithQStatesFromFile(Model& model,
            const std::string& qStatesFileName, double rollRMSinDeg, double headingRMSinDeg, double pitchRMSinDeg,
            double stateRMS = 0.0, unsigned seed = 0);

    std::tuple<std::map<std::string, int>, std::map<int, std::string>> CreateYMaps(Model model);

private:
    void constructProperties();


    //=============================================================================
};  // END of class IMUInverseKinematicsTool
//=============================================================================
} // namespace






#endif // IMUFORWARDKINEMATICSTOOL_H
