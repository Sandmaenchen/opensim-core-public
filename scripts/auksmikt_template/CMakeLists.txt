cmake_minimum_required(VERSION 3.21)
project(Opensim_Analysis_Scripts)

# Settings.
# ---------
set(TARGET "runAUKSMIKT" CACHE STRING "runAUKSMIKT")

# OpenSim uses C++11 language features.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native")

# Find and hook up to OpenSim.
# ----------------------------
find_package(OpenSim REQUIRED PATHS "/home/<username>/opensim-core/")

# Configure this project.
# -----------------------
file(GLOB SOURCE_FILES *.h *.cpp)

add_executable(${TARGET} ${SOURCE_FILES})

target_link_libraries(${TARGET} ${OpenSim_LIBRARIES})

