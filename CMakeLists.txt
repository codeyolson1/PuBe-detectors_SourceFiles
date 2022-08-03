# Setup Project for CMAKE and find Components
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(Helium3Detector)

#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#find_package(ROOT REQUIRED COMPONENTS)
#include(${ROOT_USE_FILE})
#include_directories(${ROOT_INCLUDE_DIRS})
#set(CMAKE_CXX_COMPILER mpicxx)

# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
#if(ROOT_FOUND)
#add_definitions(-DG4ANALYSIS_USE_ROOT)
include(${Geant4_USE_FILE})
include_directories(${Geant4_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include ${G4mpi_INCLUDE_DIR})
#else()
#include(${Geant4_USE_FILE})
#include_directories(${PROJECT_SOURCE_DIR}/include)
#endif()
add_definitions(-DTOOLS_USE_NATIVE_MPI) 
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executable, and link it to the Geant4 libraries
#
add_executable(he3 main.cc ${sources} ${headers})
#if(ROOT_FOUND)
target_link_libraries(he3 ${Geant4_LIBRARIES} ${G4mpi_LIBRARIES}) #add ${ROOT_LIBRARIES} here if needed
#else()
#target_link_libraries(reactorBay ${Geant4_LIBRARIES})
#endif()

# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
#set(reactorBay_scripts
#  init.mac
#  vis.mac
#  run.mac
#)

#foreach(_script ${reactorBay_scripts})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/macros/${_script}
#    ${PROJECT_BINARY_DIR}/${_script}
#    COPYONLY
#    )
#endforeach()

# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(Helium3Detector DEPENDS he3)

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS he3 DESTINATION bin)