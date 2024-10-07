
file(GLOB INCLUDES *.h)
file(GLOB SOURCES *.cpp)

OpenSimAddLibrary(
    KIT Tools
    AUTHORS "Frank_C_Anderson"
    LINKLIBS osimCommon osimSimulation osimActuators osimAnalyses Eigen3::Eigen
    INCLUDES ${INCLUDES}
    SOURCES ${SOURCES}
    TESTDIRS "Test"
    )
