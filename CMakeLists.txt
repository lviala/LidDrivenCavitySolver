cmake_minimum_required(VERSION 3.0.0)
project(LidDrivenCavitySolver VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)

include(CTest)
enable_testing()

find_package(Boost COMPONENTS system iostreams iostreams system program_options REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_executable(LidDrivenCavitySolver src/LidDrivenCavitySolver.cpp src/LidDrivenCavity.cpp src/LDCprogram_options)
TARGET_LINK_LIBRARIES( LidDrivenCavitySolver ${Boost_LIBRARIES} )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

