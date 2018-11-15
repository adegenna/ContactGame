#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <string>
#include <fstream>
#include <istream>
#include <random>
#include <cmath>
#include <Eigen/Dense>
#include <boost/program_options.hpp>
#include "options.hpp"
#include "options_parser.hpp"
#include "LagrangianState.h"
#include "ParticlePhysics.h"
#include "TimeIntegration.h"
#include "Inputfile.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
  printf("*********** DRIVER PROGRAM FOR LAGRANGIAN PARTICLE SOLVER ***********\n\n");

  // Parse input file options
  Options options;
  if (!parseOptions(argc,argv,options)){
   return 0;
  }
  cout << options << endl;

  // Load state
  MatrixXd input = load_csv<MatrixXd>(options.inputfile);
  
  // Pass parsed program options to simulation
  LagrangianState state(input);
  ParticlePhysics physics(options, state);
  TimeIntegration integrator(options, physics, state);
  
  // Solve
  integrator.euler();

  // Output
  state.writeXY(options.outputfile+"_final.csv");
  
  return 0;
}
