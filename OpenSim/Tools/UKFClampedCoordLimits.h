#ifndef UKFCLAMPEDCOORDLIMITS_H
#define UKFCLAMPEDCOORDLIMITS_H

#include <string>
namespace OpenSim {
// Struct for holding the clamped coordinate limits
struct UKFClampedCoordLimits {
    public:
        std::string stateVarName;
        double rangeMin;
        double rangeMax;

        UKFClampedCoordLimits(std::string name, double min, double max) {
            stateVarName = name;
            rangeMin = min;
            rangeMax = max;
        }
}; // END of struct
}
#endif // UKFCLAMPEDCOORDLIMITS_H