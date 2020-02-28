cmake_minimum_required(VERSION 3.0.0)
project(LidDrivenCavitySolver VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)

include(CTest)
enable_testing()

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(Boost COMPONENTS system iostreams iostreams system program_options REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

find_package(BLAS REQUIRED)
INCLUDE_DIRECTORIES(${BLAS_INCLUDE_DIR})

find_package(LAPACK REQUIRED)
INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIR})

add_executable(LidDrivenCavitySolver src/LidDrivenCavitySolver.cpp src/LidDrivenCavity.cpp src/LDCpoissonSolver.cpp src/LDCprogram_options.cpp src/LDCmngMPI.cpp)

TARGET_LINK_LIBRARIES( LidDrivenCavitySolver ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set_target_properties(LidDrivenCavitySolver PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

