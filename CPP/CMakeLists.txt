cmake_minimum_required(VERSION 2.6)
project(lagrangiansolver)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(GTest REQUIRED)
find_package(Threads)

add_definitions( -std=c++11 -c -pg -g -DBOOST_LOG_DYN_LINK -DNDEBUG -O3 )

include_directories( ${CMAKE_SOURCE_DIR}
                     ${Boost_INCLUDE_DIRS}
                     ${GTEST_INCLUDE_DIRS}
		     /usr/local/include/eigen/
)

add_executable( testall
  tests/tests.cpp
  src/Inputfile.hpp
  src/LagrangianSimulation.cpp
  src/ParticlePhysics.cpp
  tests/inputs_test.cpp
  tests/lagrangian_test.cpp
  tests/physics_test.cpp
)

add_executable( solver
  src/driver.cpp
  src/Inputfile.hpp
  src/LagrangianSimulation.cpp
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
