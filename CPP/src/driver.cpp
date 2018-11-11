#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <string>
#include <fstream>
#include <istream>
#include <random>
#include <cmath>
#include <eigen/Eigen/Dense>
#include <boost/program_options.hpp>
#include "options.hpp"
#include "options_parser.hpp"
#include "LagrangianSimulation.h"
#include "ParticlePhysics.h"

using namespace std;
using namespace Eigen;

// ***************************************************
// DRIVER PROGRAM FOR LAGRANGIAN SIMULATION
// ***************************************************

int main(int argc, char* argv[]) {
  printf("*********** DRIVER PROGRAM FOR LAGRANGIAN PARTICLE SOLVER ***********\n\n");

  // Parse input file options
  Options options;
  if (!parseOptions(argc,argv,options)){
   return 0;
  }
  cout << options << endl;

  // Pass parsed program options to simulation
  LagrangianSimulation simulation(options);
  ParticlePhysics physics(options, simulation);
  
  // Solve
  simulation.setupInitialConditions();
  physics.simulate();

  // Output
  simulation.writeXY("final");
  
  return 0;
}
