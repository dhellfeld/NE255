cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# setting a variable equal to the executable name
set(name PRISM)

# name the project
project(PRISM)

# find Geant4 (uses Geant4_DIR)
find_package(Geant4 REQUIRED ui_all vis_all)

# sets the Geant4_INCLUDE_DIR and Geant4_LIBRARIES
include(${Geant4_USE_FILE})

# include your headers and the Geant4 headers
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})

# get the implementation files
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

# get the header files
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# create the executable
add_executable(${name} ${name}.cc ${sources} ${headers})

# link the Geant4 libraries to the executable
target_link_libraries(${name} ${Geant4_LIBRARIES})
                    