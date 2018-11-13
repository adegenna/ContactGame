cmake_minimum_required(VERSION 2.6)
project(lagrangiansolver)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(GTest REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_definitions( -std=c++11 -c -pg -g -DBOOST_LOG_DYN_LINK -DNDEBUG -O3 )

include_directories( ${CMAKE_SOURCE_DIR}
                     ${Boost_INCLUDE_DIRS}
                     ${GTEST_INCLUDE_DIRS}
                     ${EIGEN3_INCLUDE_DIR}
)

add_executable( testall
  tests/tests.cpp
  src/Inputfile.hpp
  src/LagrangianState.cpp
  src/ParticlePhysics.cpp
  tests/inputs_test.cpp
  tests/lagrangian_test.cpp
  tests/physics_test.cpp
)

add_executable( solver
  src/driver.cpp
  src/Inputfile.hpp
  src/LagrangianState.cpp
  src/ParticlePhysics.cpp
  src/options.hpp
  src/options_parser.hpp
)


target_link_libraries( solver
  ${Boost_LIBRARIES}   
)

target_link_libraries(testall
  ${GTEST_LIBRARIES}
  ${Boost_LIBRARIES}
  /usr/lib/libgtest.a
  pthread
)
